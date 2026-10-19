"""Real CPU MMseqs2 integration: full proteins, fragments, repeated local alignments."""

import csv
import json
import os
import random
import shutil

import pytest
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import MMseqsError
from sprb_scanner.io.common import sha256
from sprb_scanner.io.prepare import prepare
from sprb_scanner.scanners.mmseqs import MMseqsConfig, SequenceAlignment, run_search


@pytest.fixture
def real_data(tmp_path):
    if shutil.which("mmseqs") is None:
        if os.environ.get("SPRB_REQUIRE_MMSEQS") == "1":
            pytest.fail("CI requires the MMseqs2 binary")
        pytest.skip("install MMseqs2 to run real backend integration tests")
    rng = random.Random(42)
    seq = lambda n: "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=n))  # noqa: E731
    full, fragment, spacer, unrelated = seq(120), seq(90), seq(45), seq(150)
    targets = tmp_path / "target.faa"
    targets.write_text(
        f">sp|FULL:1 full self match\n{full}\n"
        f">ref|REPEAT:2 two copies\n{spacer}{fragment}{spacer}{fragment}{spacer}\n"
    )
    prepared = tmp_path / "prepared"
    prepare(targets, "protein_fasta", prepared, AppConfig(), "synthetic")
    query = tmp_path / "query.faa"
    query.write_text(f">full|query full length\n{full}\n>fragment|query\n{fragment}\n")
    empty = tmp_path / "unrelated.faa"
    empty.write_text(f">unrelated\n{unrelated}\n")
    return prepared, query, empty


def check_outputs(output):
    metadata = json.loads((output / "run_metadata.json").read_text())
    assert metadata["state"] == "succeeded" and metadata["backend_version"]
    assert metadata["configuration"]["threads"] == 1
    for name, checksum in metadata["output_sha256"].items():
        assert sha256(output / name) == checksum
    assert not list(output.glob("mmseqs-work-*"))
    records = json.loads((output / "alignments.json").read_text())
    for record in records:
        record.pop("ctd_distance")
        SequenceAlignment.model_validate(record)
    with (output / "alignments.tsv").open() as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    assert len(table) == len(records) == metadata["alignment_count"]
    assert len((output / "raw_alignments.tsv").read_text().splitlines()) == len(records)
    return records, metadata


def test_real_cli_full_fragment_multi_alignment_and_reuse(real_data, tmp_path):
    prepared, query, _ = real_data
    runner = CliRunner()
    database = tmp_path / "db"
    result = runner.invoke(
        app, ["mmseqs-createdb", "--prepared", str(prepared), "--output", str(database)]
    )
    assert result.exit_code == 0, result.output + (
        (database / "backend.log").read_text() if (database / "backend.log").exists() else ""
    )
    relocated = tmp_path / "relocated-db"
    database.rename(relocated)
    database = relocated
    db_checksums = {p.name: sha256(p) for p in database.iterdir()}
    output = tmp_path / "search"
    result = runner.invoke(
        app,
        [
            "mmseqs-search",
            "--query",
            str(query),
            "--database",
            str(database),
            "--prepared",
            str(prepared),
            "--output",
            str(output),
            "--evalue",
            "1e-8",
        ],
    )
    assert result.exit_code == 0, result.output
    hits, metadata = check_outputs(output)
    assert metadata["database_reused"]
    full = [
        h for h in hits if h["query_id"] == "full|query" and h["target_source_id"] == "sp|FULL:1"
    ]
    assert any(
        (h["start"], h["end"], h["query_start"], h["query_end"]) == (0, 120, 0, 120)
        and h["sequence_identity"] == h["query_coverage"] == h["target_coverage"] == 1
        for h in full
    )
    repeated = [
        h
        for h in hits
        if h["query_id"] == "fragment|query" and h["target_source_id"] == "ref|REPEAT:2"
    ]
    assert len(repeated) >= 2
    assert {(h["start"], h["end"]) for h in repeated} >= {(45, 135), (180, 270)}
    assert all(h["query_coverage"] == 1 and h["target_coverage"] < 0.3 for h in repeated)
    assert len({h["hit_id"] for h in hits}) == len(hits)
    rerun = tmp_path / "repeat-search"
    run_search(query, rerun, database=database, config=MMseqsConfig(evalue=1e-8))
    again, _ = check_outputs(rerun)
    assert hits == again
    assert db_checksums == {p.name: sha256(p) for p in database.iterdir()}
    # Broken caches and a different prepared proteome must fail, not remap hit identities.
    wrong = tmp_path / "wrong"
    prepare(query, "protein_fasta", wrong, AppConfig(), "different")
    with pytest.raises(MMseqsError, match="disagree"):
        run_search(query, tmp_path / "mismatch", prepared=wrong, database=database)
    manifest = database / "run_metadata.json"
    original = manifest.read_text()
    metadata = json.loads(original)
    metadata["backend_version"] = "incompatible"
    manifest.write_text(json.dumps(metadata))
    with pytest.raises(MMseqsError, match="version mismatch"):
        run_search(query, tmp_path / "version-mismatch", database=database)
    manifest.write_text(original)
    (database / "target.index").write_text("corrupted")
    with pytest.raises(MMseqsError, match="checksum mismatch"):
        run_search(query, tmp_path / "corrupt", database=database)


def test_real_direct_prepared_empty_results_and_failed_backend(real_data, tmp_path):
    prepared, query, unrelated = real_data
    empty = tmp_path / "empty"
    run_search(unrelated, empty, prepared=prepared, config=MMseqsConfig(evalue=1e-8))
    hits, metadata = check_outputs(empty)
    assert hits == [] and metadata["database_reused"] is False
    assert (empty / "alignments.tsv").read_text().startswith("start\tend\t")
    failed = tmp_path / "failed-backend"
    with pytest.raises(MMseqsError):
        run_search(query, failed, prepared=prepared, config=MMseqsConfig(executable="false"))
    assert json.loads((failed / "status.json").read_text())["state"] == "failed"


def test_createdb_failure_keeps_nested_diagnostics(real_data, tmp_path):
    prepared, query, _ = real_data
    fake = tmp_path / "mmseqs-fail"
    fake.write_text(
        '#!/bin/sh\nif [ "$1" = "version" ]; then echo 18-test; exit 0; fi\n'
        'if [ "$2" = "-h" ]; then echo " --threads INT"; exit 0; fi\n'
        'echo "deliberate database failure" >&2\nexit 23\n'
    )
    fake.chmod(0o755)
    failed = tmp_path / "failed-build"
    with pytest.raises(MMseqsError, match="Target database creation failed"):
        run_search(query, failed, prepared=prepared, config=MMseqsConfig(executable=str(fake)))
    assert "deliberate database failure" in (failed / "direct_database_backend.log").read_text()
    meta = json.loads((failed / "direct_database_metadata.json").read_text())
    assert meta["state"] == "failed"
    assert json.loads((failed / "status.json").read_text())["state"] == "failed"
