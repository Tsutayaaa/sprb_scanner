import io
import json
import zipfile
from urllib.error import URLError

import pytest

from sprb_scanner.core.config import DownloadConfig
from sprb_scanner.core.exceptions import DownloadError
from sprb_scanner.io import ncbi


def package():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        base = "ncbi_dataset/data/GCF_000016645.1/"
        for name in ("genomic.gbff", "genomic.gff", "protein.faa", "assembly_genomic.fna"):
            z.writestr(base + name, "fixture")
        z.writestr(
            "ncbi_dataset/data/assembly_data_report.jsonl",
            json.dumps({"accession": "GCF_000016645.1"}),
        )
    return buffer.getvalue()


def test_download_manifest_and_no_overwrite(monkeypatch, tmp_path):
    monkeypatch.setattr(ncbi, "urlopen", lambda *a, **k: io.BytesIO(package()))
    output = tmp_path / "download"
    result = ncbi.download("GCF_000016645.1", output, DownloadConfig())
    assert len(result["files"]) == 6
    assert (output / "download_manifest.json").is_file()
    with pytest.raises(DownloadError, match="exists"):
        ncbi.download("GCF_000016645.1", output, DownloadConfig())


@pytest.mark.parametrize("accession", ["GCF_000016645", "9606", "GCF_000016645.1,GCF_000001405.1"])
def test_accession_requires_one_version(accession, tmp_path):
    with pytest.raises(DownloadError):
        ncbi.download(accession, tmp_path / "out", DownloadConfig())


@pytest.mark.parametrize("name", ["../escape", "/absolute", "bad\\file"])
def test_unsafe_zip_rejected(name, tmp_path):
    p = tmp_path / "bad.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(name, "bad")
    with pytest.raises(DownloadError):
        ncbi.extract_archive(p, tmp_path / "out", 1000)


def test_size_limit_and_cleanup(monkeypatch, tmp_path):
    monkeypatch.setattr(ncbi, "urlopen", lambda *a, **k: io.BytesIO(package()))
    with pytest.raises(DownloadError, match="byte limit"):
        ncbi.download("GCF_000016645.1", tmp_path / "out", DownloadConfig(max_download_bytes=10))
    assert not list(tmp_path.iterdir())


def test_retry_failure_is_explicit(monkeypatch, tmp_path):
    calls = []

    def failure(*args, **kwargs):
        calls.append(1)
        raise URLError("offline")

    monkeypatch.setattr(ncbi, "urlopen", failure)
    monkeypatch.setattr(ncbi.time, "sleep", lambda _: None)
    with pytest.raises(DownloadError, match="failed"):
        ncbi.download("GCF_000016645.1", tmp_path / "out", DownloadConfig(retries=1))
    assert len(calls) == 2
    assert not list(tmp_path.iterdir())
