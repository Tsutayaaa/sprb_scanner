"""Generic MMseqs2 protein alignments and integrity-checked reusable databases."""

import csv
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
import tempfile
from collections import Counter
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator

from sprb_scanner import __version__
from sprb_scanner.core.exceptions import InputError, MMseqsError
from sprb_scanner.core.models import (
    DomainHit,
    Fraction,
    GenomeDataset,
    ModelProvenance,
    NonNegative,
    Positive,
    ProteinRecord,
    Record,
    internal_id,
)
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared, write_json, write_tsv
from sprb_scanner.io.fasta import parse_proteins

INT_MAX = 2147483647
FORMAT = (
    "query,target,qlen,tlen,qstart,qend,tstart,tend,evalue,bits,alnlen,nident,"
    "fident,qcov,tcov,qaln,taln"
)


class MMseqsConfig(Record):
    sensitivity: float = Field(default=5.7, ge=1, le=7.5, allow_inf_nan=False)
    evalue: float = Field(default=0.001, gt=0, allow_inf_nan=False)
    threads: Positive = 1
    # None removes the usual 300-target cap; alternative alignments are enabled by default.
    max_seqs: int | None = Field(default=None, gt=0, le=INT_MAX, strict=True)
    alternative_alignments: int | None = Field(default=None, ge=0, le=INT_MAX, strict=True)
    split_memory_limit: str = Field(default="2G", pattern=r"^[1-9][0-9]*[KMGT]$")
    timeout_seconds: Positive = 3600
    executable: str = Field(default="mmseqs", min_length=1)


class SequenceAlignment(DomainHit):
    """DomainHit-compatible sequence evidence; identity is a fraction, not a percent."""

    method: Literal["mmseqs2"] = "mmseqs2"
    evalue_kind: Literal["alignment"] = "alignment"
    sequence_identity: Fraction
    identical_residues: NonNegative
    alignment_length: Positive
    query_source_id: str
    target_source_id: str | None
    target_gene_id: str | None

    @model_validator(mode="after")
    def alignment_statistics(self):
        if self.hmm is not None:
            raise ValueError("Sequence alignments cannot carry HMM statistics")
        if self.identical_residues > self.alignment_length:
            raise ValueError("Identity count exceeds alignment length")
        if self.alignment_length < max(self.query_end - self.query_start, self.end - self.start):
            raise ValueError("Alignment length is shorter than the coordinate spans")
        for actual, expected in (
            (self.sequence_identity, self.identical_residues / self.alignment_length),
            (self.query_coverage, (self.query_end - self.query_start) / self.query_length),
            (self.target_coverage, (self.end - self.start) / self.protein_length),
        ):
            if not math.isclose(actual, expected, abs_tol=1e-12):
                raise ValueError("Alignment statistics disagree with counts/coordinates")
        return self


def fingerprint(proteins):
    return hashlib.sha256(
        json.dumps(
            [(p.protein_id, p.sequence) for p in proteins],
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def load_prepared(path):
    if validate_prepared(path)["errors"]:
        raise InputError("Prepared dataset has quality errors")
    dataset = GenomeDataset.model_validate_json((path / "dataset.json").read_text())
    if not dataset.proteins:
        raise InputError("Prepared proteome contains no proteins")
    return dataset.proteins, {
        "prepared": str(path.resolve()),
        "dataset_sha256": sha256(path / "dataset.json"),
        "prepared_metadata_sha256": sha256(path / "run_metadata.json"),
        "genome": dataset.genome.model_dump(mode="json"),
    }


def load_queries(path):
    checksum = sha256(path)
    dataset = parse_proteins(path, "mmseqs-query")
    if any(
        x.severity == "error" or x.code == "duplicate_source_protein_id"
        for x in dataset.quality_issues
    ):
        raise InputError("Query FASTA has invalid sequences or duplicate IDs")
    if not dataset.proteins:
        raise InputError("Query FASTA contains no valid proteins")
    if sha256(path) != checksum:
        raise InputError("Query FASTA changed while reading")
    return tuple(p.model_copy(update={"protein_id": p.source_id}) for p in dataset.proteins), {
        "path": str(path.resolve()),
        "sha256": checksum,
        "normalization_issues": [i.model_dump() for i in dataset.quality_issues],
    }


class Backend:
    def __init__(self, config, metadata, log):
        self.config, self.metadata, self.log = config, metadata, log
        executable = shutil.which(config.executable)
        if executable is None:
            raise MMseqsError(f"MMseqs2 executable not found: {config.executable}")
        self.executable = str(Path(executable).resolve())
        try:
            result = subprocess.run(
                [self.executable, "version"], check=True, capture_output=True, text=True, timeout=15
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise MMseqsError(f"Cannot read MMseqs2 version: {exc}") from exc
        self.version = result.stdout.strip()
        if not self.version:
            raise MMseqsError("MMseqs2 returned an empty version")
        metadata.update(backend_version=self.version, executable=self.executable, commands=[])

    def run(self, *args):
        command = [self.executable, *map(str, args)]
        self.metadata["commands"].append(command)
        with self.log.open("a", encoding="utf-8") as log:
            log.write(json.dumps(command) + "\n")
            log.flush()
            try:
                with subprocess.Popen(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=os.name == "posix",
                ) as process:
                    try:
                        code = process.wait(timeout=self.config.timeout_seconds)
                    except subprocess.TimeoutExpired:
                        # MMseqs workflows spawn children; stop the whole workflow on timeout.
                        if os.name == "posix":
                            os.killpg(process.pid, signal.SIGKILL)
                        else:
                            process.kill()
                        process.wait()
                        raise
                    if code:
                        raise subprocess.CalledProcessError(code, command)
            except (OSError, subprocess.SubprocessError) as exc:
                raise MMseqsError(f"MMseqs2 {args[0]} failed; see backend.log: {exc}") from exc

    def createdb(self, fasta, prefix):
        self.run(
            "createdb",
            fasta,
            prefix,
            "--dbtype",
            1,
            "--shuffle",
            0,
            "--threads",
            self.config.threads,
        )


@contextmanager
def run_directory(output, kind, config):
    output = Path(output)
    if output.exists():
        raise MMseqsError(f"Output already exists; choose a new directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".mmseqs-", dir=output.parent))
    metadata = {
        "schema_version": "1.0.0",
        "kind": kind,
        "software_version": __version__,
        "created_at": datetime.now(UTC).isoformat(),
        "configuration": config.model_dump(),
        "state": "running",
    }
    failure = None
    try:
        try:
            yield stage, metadata
            metadata["state"] = "succeeded"
        except (InputError, MMseqsError, OSError, ValueError) as exc:
            failure = MMseqsError(str(exc))
            metadata.update(state="failed", error=str(exc))
        write_json(
            stage / "status.json",
            {
                "state": metadata["state"],
                "message": metadata.get("error"),
            },
        )
        metadata["output_sha256"] = {
            str(p.relative_to(stage)): sha256(p) for p in sorted(stage.rglob("*")) if p.is_file()
        }
        write_json(stage / "run_metadata.json", metadata)
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    if failure:
        raise failure


def write_mapping(directory, proteins, prefix):
    lookup = {f"{prefix}{i}": p for i, p in enumerate(proteins)}
    write_json(
        directory / "records.json",
        {key: p.model_dump(mode="json", exclude_computed_fields=True) for key, p in lookup.items()},
    )
    (directory / "sequences.faa").write_text(
        "".join(f">{key}\n{p.sequence}\n" for key, p in lookup.items()), encoding="utf-8"
    )
    return lookup


def create_database(prepared: Path, output: Path, config=None):
    """Build a managed, portable database; never overwrite an existing directory."""
    config = config or MMseqsConfig()
    with run_directory(output, "mmseqs2-database", config) as (stage, metadata):
        proteins, source = load_prepared(prepared)
        metadata.update(source=source, proteins_sha256=fingerprint(proteins), count=len(proteins))
        backend = Backend(config, metadata, stage / "backend.log")
        write_mapping(stage, proteins, "t")
        backend.createdb(stage / "sequences.faa", stage / "target")
    return {"state": "succeeded", "proteins": len(proteins), "output": str(output.resolve())}


def load_database(directory, version):
    """Only managed databases can guarantee target identities and sequence provenance."""
    metadata = json.loads((directory / "run_metadata.json").read_text())
    if not isinstance(metadata, dict):
        raise InputError("Invalid MMseqs2 database metadata object")
    if (
        metadata.get("schema_version") != "1.0.0"
        or metadata.get("kind") != "mmseqs2-database"
        or metadata.get("state") != "succeeded"
    ):
        raise InputError("Not a successful managed MMseqs2 database")
    if metadata.get("backend_version") != version:
        raise InputError("MMseqs2 database version mismatch; rebuild with the current executable")
    files = metadata.get("output_sha256", {})
    required = {
        "records.json",
        "sequences.faa",
        "target",
        "target.index",
        "target.dbtype",
        "target_h",
        "target_h.index",
        "target_h.dbtype",
        "status.json",
    }
    if not isinstance(files, dict) or not required <= set(files):
        raise InputError("MMseqs2 database manifest is incomplete")
    for name, checksum in files.items():
        path = directory / name
        if (
            Path(name).name != name
            or path.is_symlink()
            or not path.is_file()
            or sha256(path) != checksum
        ):
            raise InputError(f"MMseqs2 database checksum mismatch or unsafe path: {name}")
    records = json.loads((directory / "records.json").read_text())
    if not isinstance(records, dict):
        raise InputError("Invalid MMseqs2 database records object")
    lookup = {k: ProteinRecord.model_validate(v) for k, v in records.items()}
    proteins = tuple(lookup.values())
    if (
        not proteins
        or len({p.protein_id for p in proteins}) != len(proteins)
        or len(proteins) != metadata.get("count")
        or fingerprint(proteins) != metadata.get("proteins_sha256")
        or list(lookup) != [f"t{i}" for i in range(len(proteins))]
    ):
        raise InputError("Invalid MMseqs2 database protein mapping")
    expected = "".join(f">{key}\n{p.sequence}\n" for key, p in lookup.items())
    if (directory / "sequences.faa").read_text() != expected:
        raise InputError("MMseqs2 database FASTA and mapping disagree")
    return lookup, metadata


def parse_alignments(path, queries, targets, version, query_source):
    """Convert every backend row; no best-hit collapsing, coverage filter or deduplication."""
    hits, counts = [], Counter()
    with path.open(newline="", encoding="utf-8") as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != len(FORMAT.split(",")):
                raise MMseqsError(f"Malformed MMseqs2 alignment row {number}")
            values = dict(zip(FORMAT.split(","), row, strict=True))
            if values["query"] not in queries or values["target"] not in targets:
                raise MMseqsError(f"Unknown query/target identity in alignment row {number}")
            q, t = queries[values["query"]], targets[values["target"]]
            qlen, tlen = int(values["qlen"]), int(values["tlen"])
            qs, qe, ts, te = (int(values[k]) for k in ("qstart", "qend", "tstart", "tend"))
            alnlen, nident = int(values["alnlen"]), int(values["nident"])
            qaln, taln = values["qaln"], values["taln"]
            if (
                qlen != q.length
                or tlen != t.length
                or alnlen <= 0
                or len(qaln) != alnlen
                or len(taln) != alnlen
                or len(qaln.replace("-", "")) != qe - qs + 1
                or len(taln.replace("-", "")) != te - ts + 1
            ):
                raise MMseqsError(f"Alignment length/identity mapping mismatch in row {number}")
            counts[q.protein_id, t.protein_id] += 1
            query_hash = hashlib.sha256(q.sequence.encode()).hexdigest()
            hits.append(
                SequenceAlignment(
                    hit_id=internal_id(
                        "mmseqs2",
                        q.protein_id,
                        t.protein_id,
                        str(counts[q.protein_id, t.protein_id]),
                    ),
                    protein_id=t.protein_id,
                    model=ModelProvenance(
                        model_id=q.protein_id,
                        version=f"sequence-sha256:{query_hash}",
                        source=query_source,
                        sha256=query_hash,
                    ),
                    query_id=q.protein_id,
                    query_source_id=q.source_id,
                    target_source_id=t.source_id,
                    target_gene_id=t.gene_id,
                    query_length=qlen,
                    protein_length=tlen,
                    query_start=qs - 1,
                    query_end=qe,
                    start=ts - 1,
                    end=te,
                    evalue=float(values["evalue"]),
                    bit_score=float(values["bits"]),
                    sequence_identity=nident / alnlen,
                    identical_residues=nident,
                    alignment_length=alnlen,
                    query_coverage=(qe - qs + 1) / qlen,
                    target_coverage=(te - ts + 1) / tlen,
                    backend_version=version,
                    raw_artifact="raw_alignments.tsv",
                )
            )
    return hits


def run_search(query: Path, output: Path, prepared=None, database=None, config=None):
    """Search one prepared proteome or reuse a managed database (optionally cross-check both)."""
    config = config or MMseqsConfig()
    with run_directory(output, "mmseqs2-search", config) as (stage, metadata):
        if prepared is None and database is None:
            raise InputError("Provide --prepared or --database")
        queries, source = load_queries(query)
        metadata["query_source"] = source
        backend = Backend(config, metadata, stage / "backend.log")
        with tempfile.TemporaryDirectory(prefix="mmseqs-work-", dir=stage) as scratch:
            scratch = Path(scratch)
            if database is None:
                database = scratch / "target_database"
                try:
                    create_database(prepared, database, config)
                except MMseqsError as exc:
                    raise MMseqsError(
                        f"Target database creation failed; see direct_database_backend.log "
                        f"and direct_database_metadata.json: {exc}"
                    ) from exc
                finally:
                    # Keep nested build failures/logs before the temporary database is removed.
                    for name, destination in (
                        ("backend.log", "direct_database_backend.log"),
                        ("run_metadata.json", "direct_database_metadata.json"),
                    ):
                        if (database / name).is_file():
                            shutil.copyfile(database / name, stage / destination)
                metadata["database_reused"] = False
            else:
                database = Path(database).resolve()
                metadata["database_reused"] = True
            targets, db_metadata = load_database(database, backend.version)
            metadata.update(
                target_database_metadata=db_metadata,
                database_metadata_sha256=sha256(database / "run_metadata.json"),
            )
            if prepared is not None:
                proteins, prepared_source = load_prepared(prepared)
                if fingerprint(proteins) != db_metadata["proteins_sha256"]:
                    raise InputError("Prepared proteome and reused MMseqs2 database disagree")
                metadata["prepared_source"] = prepared_source
            query_dir = scratch / "query"
            query_dir.mkdir()
            query_map = write_mapping(query_dir, queries, "q")
            shutil.copyfile(query_dir / "records.json", stage / "query_mapping.json")
            shutil.copyfile(database / "records.json", stage / "target_mapping.json")
            backend.createdb(query_dir / "sequences.faa", scratch / "query_db")
            max_seqs = config.max_seqs or len(targets)
            # MMseqs masks the previous target interval before the next local alignment.
            # One alternative per target residue is a finite, length-derived upper bound.
            alt = (
                config.alternative_alignments
                if config.alternative_alignments is not None
                else max(t.length for t in targets.values())
            )
            max_len = max(
                65535, *(p.length for p in queries), *(p.length for p in targets.values())
            )
            metadata["effective_parameters"] = {
                "max_seqs": max_seqs,
                "alternative_alignments": alt,
                "max_seq_len": max_len,
                "coverage_threshold": 0,
                "min_sequence_identity": 0,
                "max_accept": INT_MAX,
                "max_rejected": INT_MAX,
            }
            backend.run(
                "search",
                scratch / "query_db",
                database / "target",
                scratch / "result",
                scratch / "tmp",
                "--search-type",
                1,
                "-s",
                config.sensitivity,
                "-e",
                config.evalue,
                "--threads",
                config.threads,
                "-a",
                1,
                "--alignment-mode",
                3,
                "--alt-ali",
                alt,
                "--max-seqs",
                max_seqs,
                "--max-accept",
                INT_MAX,
                "--max-rejected",
                INT_MAX,
                "-c",
                0,
                "--min-seq-id",
                0,
                "--min-aln-len",
                0,
                "--seq-id-mode",
                0,
                "--max-seq-len",
                max_len,
                "--split-memory-limit",
                config.split_memory_limit,
            )
            backend.run(
                "convertalis",
                scratch / "query_db",
                database / "target",
                scratch / "result",
                stage / "raw_alignments.tsv",
                "--format-output",
                FORMAT,
                "--threads",
                config.threads,
            )
            hits = parse_alignments(
                stage / "raw_alignments.tsv", query_map, targets, backend.version, source["path"]
            )
            # Detect source/cache mutation rather than attach alignments to changed identities.
            if sha256(query) != source["sha256"]:
                raise InputError("Query FASTA changed during search")
            load_database(database, backend.version)
            if sha256(database / "run_metadata.json") != metadata["database_metadata_sha256"]:
                raise InputError("MMseqs2 database manifest changed during search")
            write_json(stage / "alignments.json", [h.model_dump(mode="json") for h in hits])
            write_tsv(stage / "alignments.tsv", hits, list(SequenceAlignment.model_fields))
            metadata.update(
                alignment_count=len(hits),
                query_count=len(queries),
                target_count=len(targets),
                query_proteins_sha256=fingerprint(queries),
                raw_columns=FORMAT.split(","),
                coordinates="0-based half-open; raw: 1-based inclusive",
                identity="identical_residues / alignment_length (including gap columns)",
                coverage="independent query/target coordinate span divided by its protein length",
                retention=(
                    "all convertalis rows, including alternative alignments; no downstream "
                    "filter/deduplication; MMseqs heuristic prefilter, E-value reporting and "
                    "recorded max_seqs/alternative limits still apply; not exhaustive HSPs"
                ),
                pairs_reaching_alternative_limit=sum(
                    count >= alt + 1
                    for count in Counter((h.query_id, h.protein_id) for h in hits).values()
                ),
            )
    return {"state": "succeeded", "alignment_count": len(hits), "output": str(output.resolve())}
