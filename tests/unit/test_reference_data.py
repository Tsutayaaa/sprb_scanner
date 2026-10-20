"""Offline checks for the curated RN28-first reference dataset and its source snapshot."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from sprb_scanner.io.common import sha256
from sprb_scanner.scanners.mmseqs import load_queries

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = ROOT / "reference_data"


def read_json(name):
    return json.loads((REFERENCES / name).read_text())


def read_tsv(name):
    with (REFERENCES / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_reference_file_integrity_and_source_snapshot():
    manifest = read_json("reference_manifest.json")
    actual = {str(p.relative_to(REFERENCES)) for p in REFERENCES.rglob("*") if p.is_file()}
    assert set(manifest["output_sha256"]) == actual - {"reference_manifest.json"}
    for name, checksum in manifest["output_sha256"].items():
        assert sha256(REFERENCES / name) == checksum, name
    for source in manifest["source_files"]:
        path = REFERENCES / source["path"]
        assert sha256(path) == source["sha256"] and path.stat().st_size == source["bytes"]
    fixture = REFERENCES / manifest["parent"]["validation_fixture"]
    assert sha256(fixture) == manifest["parent"]["validation_fixture_sha256"]


def test_primary_rn28_is_separate_from_reference_only_panel():
    primary, metadata = load_queries(REFERENCES / "RN28.fasta")
    other, _ = load_queries(REFERENCES / "reference_only.fasta")
    panel, _ = load_queries(REFERENCES / "archive/sprb_smf_panel.fasta")
    assert not metadata["normalization_issues"]
    assert len(primary) == 1 and len(other) == 73 and len(panel) == 74
    rn28 = primary[0]
    assert rn28.protein_id == "SMF_RN28__RN28" and rn28.length == 207
    expected = read_json("reference_manifest.json")["primary_query"]
    assert hashlib.sha256(rn28.sequence.encode()).hexdigest() == expected["sequence_sha256"]
    assert expected["path"] == "RN28.fasta" and expected["parent_interval"] == {
        "start": 681,
        "end": 888,
    }
    assert {p.protein_id for p in primary}.isdisjoint(p.protein_id for p in other)
    assert {p.protein_id: p.sequence for p in (*primary, *other)} == {
        p.protein_id: p.sequence for p in panel
    }


def test_all_catalog_intervals_and_legacy_group_mappings():
    manifest = read_json("reference_manifest.json")
    fixture = REFERENCES / manifest["parent"]["validation_fixture"]
    sprb = next(
        p for p in SeqIO.parse(fixture, "fasta") if p.id == manifest["parent"]["source_protein_id"]
    )
    sequence = str(sprb.seq)
    assert len(sequence) == manifest["parent"]["length"] == 6497
    assert hashlib.sha256(sequence.encode()).hexdigest() == manifest["parent"]["sequence_sha256"]
    panel, _ = load_queries(REFERENCES / "archive/sprb_smf_panel.fasta")
    by_id = {p.protein_id: p for p in panel}
    catalog = read_json("reference_catalog.json")
    table = read_tsv("reference_catalog.tsv")
    assert [r["query_id"] for r in catalog] == [p.protein_id for p in panel]
    assert len(catalog) == len(table) == 74
    module_rows = {r["output_label"]: r for r in read_tsv("archive/domain_table.tsv")}
    assignments = {r["sequence_id_full"]: r for r in read_tsv("archive/cluster_assignments.tsv")}
    groups = {r["cluster_id"]: r for r in read_tsv("archive/cluster_summary.tsv")}
    assert len(groups) == manifest["legacy_cluster_count"] == 46
    assert Counter(r["role"] for r in catalog) == {"primary_search_query": 1, "reference_only": 73}
    for row, tab in zip(catalog, table, strict=True):
        p = by_id[row["query_id"]]
        assert row["sequence"] == p.sequence == sequence[row["parent_start"] : row["parent_end"]]
        assert row["sequence_sha256"] == hashlib.sha256(p.sequence.encode()).hexdigest()
        assert row["length"] == p.length == row["parent_end"] - row["parent_start"]
        module = module_rows[row["member_label"]]
        assert row["source_module_row"] == module
        assert (row["parent_start"], row["parent_end"]) == (
            int(module["start"]) - 1,
            int(module["end"]),
        )
        assert int(module["length"]) == p.length
        assignment = assignments[row["member_label"]]
        assert row["source_assignment_row"] == assignment
        assert row["legacy_cluster_id"] == assignment["cluster"]
        assert row["legacy_group"] == groups[assignment["cluster"]]["group_name"]
        matches = [
            {"start": i, "end": i + p.length}
            for i in range(len(sequence))
            if sequence.startswith(p.sequence, i)
        ]
        assert matches == row["all_exact_parent_intervals"]
        assert len(matches) == row["exact_parent_match_count"]
        for key, value in tab.items():
            if key == "all_exact_parent_intervals":
                assert json.loads(value) == row[key]
            else:
                assert value == str(row[key])
        if row["role"] == "primary_search_query":
            assert matches == [{"start": 681, "end": 888}]
