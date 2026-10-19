"""Normalization, identity integrity and failures without an external binary."""

import json
import sys

import pytest
from pydantic import ValidationError

from sprb_scanner.core.exceptions import InputError, MMseqsError
from sprb_scanner.core.models import DomainHit, ProteinRecord
from sprb_scanner.scanners.mmseqs import (
    Backend,
    MMseqsConfig,
    SequenceAlignment,
    load_database,
    load_queries,
    parse_alignments,
    run_search,
)


def protein(pid, sequence="ACDEFGHIKL"):
    return ProteinRecord(protein_id=pid, genome_id="test", source_id=pid, sequence=sequence)


@pytest.fixture
def raw_row():
    # Gapped alignment: query residues 2..6, target 4..9, six columns / four identities.
    return "q0\tt0\t10\t10\t2\t6\t4\t9\t1e-9\t42.5\t6\t4\t0.667\t0.5\t0.6\tCDE-FG\tCDAAFG\n"


def test_gapped_coordinates_independent_coverage_and_duplicate_rows(tmp_path, raw_row):
    raw = tmp_path / "raw.tsv"
    raw.write_text(raw_row * 2)
    q, t = protein("fragment|ref"), protein("genome::protein", "MKRCDAAFGW")
    hits = parse_alignments(raw, {"q0": q}, {"t0": t}, "18-test", "reference.faa")
    assert len(hits) == 2 and hits[0].hit_id != hits[1].hit_id
    h = hits[0]
    assert isinstance(h, DomainHit)
    assert (h.query_id, h.protein_id) == (q.protein_id, t.protein_id)
    assert (h.query_start, h.query_end, h.start, h.end) == (1, 6, 3, 9)
    assert (h.query_coverage, h.target_coverage) == (0.5, 0.6)
    assert h.sequence_identity == 4 / 6
    assert h.alignment_length == 6 and h.evalue_kind == "alignment"
    assert h.model.model_id == q.protein_id and h.model.source == "reference.faa"
    assert h.bit_score == 42.5 and h.evalue == 1e-9
    assert (
        SequenceAlignment.model_validate_json(h.model_dump_json(exclude_computed_fields=True)) == h
    )
    with pytest.raises(ValidationError):
        SequenceAlignment.model_validate(
            h.model_dump(exclude_computed_fields=True) | {"target_coverage": 1}
        )


@pytest.mark.parametrize(
    "change",
    [
        {0: "missing"},
        {1: "missing"},
        {2: "11"},
        {3: "9"},
        {4: "0"},
        {5: "11"},
        {6: "0"},
        {7: "11"},
        {8: "nan"},
        {9: "inf"},
        {10: "7"},
        {11: "7"},
        {15: "CDEF"},
    ],
)
def test_reject_corrupt_alignment_rows(tmp_path, raw_row, change):
    row = raw_row.rstrip().split("\t")
    for index, value in change.items():
        row[index] = value
    raw = tmp_path / "raw.tsv"
    raw.write_text("\t".join(row) + "\n")
    with pytest.raises((MMseqsError, ValueError)):
        parse_alignments(raw, {"q0": protein("q")}, {"t0": protein("t")}, "v", "query.faa")


def test_empty_and_malformed_alignments(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text("")
    assert parse_alignments(raw, {}, {}, "v", "query.faa") == []
    raw.write_text("incomplete\trow\n")
    with pytest.raises(MMseqsError, match="Malformed"):
        parse_alignments(raw, {}, {}, "v", "query.faa")


@pytest.mark.parametrize(
    "text", ["", "ACDE", ">q\n", ">q\nAC-DE\n", ">q\nACD*EF\n", ">q\nACDE\n>q\nACDE\n"]
)
def test_invalid_query_never_silently_dropped(tmp_path, text):
    path = tmp_path / "query.faa"
    path.write_text(text)
    with pytest.raises((InputError, ValueError)):
        load_queries(path)


def test_query_ids_and_explicit_normalization(tmp_path):
    path = tmp_path / "query.faa"
    path.write_text(">ref|fragment:1 description\nacdef*\n>full\nACDEFGHIKL\n")
    queries, meta = load_queries(path)
    assert [q.protein_id for q in queries] == ["ref|fragment:1", "full"]
    assert queries[0].sequence == "ACDEF"
    assert {x["code"] for x in meta["normalization_issues"]} == {
        "uppercase_normalization",
        "terminal_stop_removed",
    }
    assert len(meta["sha256"]) == 64


@pytest.mark.parametrize(
    "settings",
    [
        {"threads": 0},
        {"threads": True},
        {"sensitivity": 8},
        {"sensitivity": float("nan")},
        {"evalue": 0},
        {"evalue": float("inf")},
        {"max_seqs": 0},
        {"max_seqs": 2147483648},
        {"alternative_alignments": 2147483648},
        {"alternative_alignments": -1},
        {"split_memory_limit": "arbitrary"},
        {"timeout_seconds": 0},
    ],
)
def test_invalid_config(settings):
    with pytest.raises(ValidationError):
        MMseqsConfig(**settings)


def test_failed_runs_are_preserved_and_never_overwritten(tmp_path):
    output = tmp_path / "failed"
    with pytest.raises(MMseqsError, match="Provide"):
        run_search(tmp_path / "missing.faa", output)
    metadata = json.loads((output / "run_metadata.json").read_text())
    assert metadata["state"] == "failed" and metadata["output_sha256"]
    before = (output / "run_metadata.json").read_bytes()
    with pytest.raises(MMseqsError, match="already exists"):
        run_search(tmp_path / "missing.faa", output)
    assert (output / "run_metadata.json").read_bytes() == before
    q = tmp_path / "query.faa"
    q.write_text(">q\nACDEFGHIKL\n")
    with pytest.raises(MMseqsError, match="not found"):
        run_search(
            q,
            tmp_path / "no-executable",
            database=tmp_path,
            config=MMseqsConfig(executable="not-a-real-mmseqs-command"),
        )
    assert json.loads((tmp_path / "no-executable/status.json").read_text())["state"] == "failed"


def test_backend_failure_and_timeout_keep_log(tmp_path):
    # Invoke Python as a harmless external command; no mock of subprocess behavior.
    backend = object.__new__(Backend)
    backend.executable = sys.executable
    backend.config = MMseqsConfig(timeout_seconds=1)
    backend.metadata = {"commands": []}
    backend.log = tmp_path / "backend.log"
    with pytest.raises(MMseqsError, match="failed"):
        backend.run("-c", "raise SystemExit(23)")
    with pytest.raises(MMseqsError, match="timed out"):
        backend.run("-c", "import time; time.sleep(3)")
    assert len(backend.metadata["commands"]) == 2
    assert "SystemExit" in backend.log.read_text()


@pytest.mark.parametrize("manifest", [[], {}, {"kind": "mmseqs2-database", "state": "failed"}])
def test_invalid_database_metadata(tmp_path, manifest):
    (tmp_path / "run_metadata.json").write_text(json.dumps(manifest))
    with pytest.raises(InputError):
        load_database(tmp_path, "18-test")
