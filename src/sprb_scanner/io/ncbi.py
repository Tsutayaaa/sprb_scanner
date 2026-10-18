"""Single-assembly NCBI Datasets download with bounded, safe ZIP extraction."""

import json
import re
import shutil
import stat
import tempfile
import time
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zipfile import BadZipFile, ZipFile

from sprb_scanner.core.exceptions import DownloadError
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import write_json

ACCESSION = re.compile(r"^GC[AF]_[0-9]{9}\.[1-9][0-9]*$")


def extract_archive(archive, directory, max_bytes):
    with ZipFile(archive) as zf:
        total, seen = 0, set()
        members = zf.infolist()
        if len(members) > 10000:
            raise DownloadError("Archive contains too many files")
        for item in members:
            name = PurePosixPath(item.filename)
            total += item.file_size
            if (
                name.is_absolute()
                or ".." in name.parts
                or "\\" in item.filename
                or name.as_posix() in seen
                or stat.S_ISLNK(item.external_attr >> 16)
            ):
                raise DownloadError("Unsafe or duplicate archive member")
            seen.add(name.as_posix())
            if total > max_bytes:
                raise DownloadError("Archive exceeds configured extracted-byte limit")
        zf.extractall(directory)


def download(accession: str, output: Path, settings):
    if not ACCESSION.fullmatch(accession):
        raise DownloadError("Use one versioned assembly accession, e.g. GCF_000016645.1")
    if output.exists():
        raise DownloadError("Download destination exists; choose a new directory")
    query = urlencode(
        [
            ("include_annotation_type", x)
            for x in ("GENOME_FASTA", "PROT_FASTA", "GENOME_GFF", "GENOME_GBFF")
        ]
    )
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?{query}"
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".download-", dir=output.parent))
    try:
        archive = stage / "ncbi_dataset.zip"
        for attempt in range(settings.retries + 1):
            try:
                request = Request(
                    url, headers={"User-Agent": "sprb-scanner/2 genome-infrastructure"}
                )
                with (
                    urlopen(request, timeout=settings.timeout_seconds) as response,
                    archive.open("wb") as out,
                ):
                    size = 0
                    while chunk := response.read(1024 * 1024):
                        size += len(chunk)
                        if size > settings.max_download_bytes:
                            raise DownloadError("Download exceeds configured byte limit")
                        out.write(chunk)
                break
            except (URLError, TimeoutError) as exc:
                if attempt == settings.retries:
                    raise DownloadError(f"NCBI download failed: {exc}") from exc
                time.sleep(min(2**attempt, 8))
        extract_archive(archive, stage, settings.max_extracted_bytes)
        root = stage / "ncbi_dataset" / "data" / accession
        if not root.is_dir():
            raise DownloadError("NCBI package does not contain the requested assembly")
        required = [root / "genomic.gbff", root / "genomic.gff", root / "protein.faa"]
        required += sorted(root.glob("*_genomic.fna"))
        if len(required) != 4 or any(not p.is_file() or not p.stat().st_size for p in required):
            raise DownloadError("Requested annotation package is incomplete (GBFF/GFF/FNA/FAA)")
        reports = stage / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
        if not reports.is_file():
            raise DownloadError("Missing NCBI assembly metadata")
        assemblies = [json.loads(line) for line in reports.read_text().splitlines() if line.strip()]
        if len(assemblies) != 1 or assemblies[0].get("accession") != accession:
            raise DownloadError("Package contains unexpected assemblies")
        manifest = {
            "accession": accession,
            "url": url,
            "downloaded_at": datetime.now(UTC).isoformat(),
            "assembly_report": assemblies[0],
            "files": {
                str(p.relative_to(stage)): {"sha256": sha256(p), "bytes": p.stat().st_size}
                for p in sorted(stage.rglob("*"))
                if p.is_file()
            },
        }
        write_json(stage / "download_manifest.json", manifest)
        stage.rename(output)
        return manifest
    except (OSError, BadZipFile, ValueError) as exc:
        raise DownloadError(f"Invalid NCBI package: {exc}") from exc
    finally:
        if stage.exists():
            shutil.rmtree(stage)
