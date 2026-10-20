"""Frozen migration snapshots: all legacy bytes, RN28 boundaries and model provenance."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from sprb_scanner.core.config import load_profile
from sprb_scanner.io.common import sha256
from sprb_scanner.references import ReferenceRegistry, RegistryDocument

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = ROOT / "references"


def test_migration_is_byte_identical_and_legacy_paths_still_work():
    migration = json.loads((REFERENCES / "migration.json").read_text())
    assert len(migration["files"]) == 16
    for row in migration["files"]:
        legacy, canonical = ROOT / row["legacy_path"], REFERENCES / row["artifact"]
        assert legacy.read_bytes() == canonical.read_bytes()
        assert sha256(canonical) == row["sha256"] and canonical.stat().st_size == row["bytes"]
    for canonical in (REFERENCES / "hmm").glob("*.profile.json"):
        legacy = ROOT / "profiles" / canonical.name
        # User HMM paths still resolve beside the supplied manifest, for either directory.
        old, new = load_profile(legacy), load_profile(canonical)
        assert old.provenance == new.provenance and old.method == new.method == "hmm"
        assert old.path.parent == legacy.parent and new.path.parent == canonical.parent
        assert old.path.name == new.path.name


def test_registry_completeness_schema_and_original_model_sources():
    registry = ReferenceRegistry.load(REFERENCES / "registry.json")
    assert len(registry.document.references) == 20
    assert (
        json.loads((REFERENCES / "registry.schema.json").read_text())
        == RegistryDocument.model_json_schema()
    )
    report = registry.validate()
    assert not [r for r in report["resources"] if r["state"] == "failed"]
    assert {r["reference_id"] for r in report["resources"] if r["state"] == "missing_external"} <= {
        "TIGR04131",
        "PF13585",
        "TIGR03519",
        "PF11751",
        "PF13573",
    }
    for model in ("TIGR04131", "PF13585", "TIGR03519", "PF11751", "PF13573"):
        config = json.loads((ROOT / "profiles" / f"{model}.profile.json").read_text())
        entry = registry.get(model)
        provenance = config["provenance"]
        assert entry.version == provenance["version"]
        assert entry.sha256 == provenance["sha256"]
        assert entry.source.locator == provenance["source"]
        assert entry.biological_target.model_id == provenance["model_id"]
        assert entry.biological_target.model_accession == provenance["version"]
        assert (
            entry.artifact.storage == "external" and entry.artifact.path == "hmm/" + config["path"]
        )
    # Registry coverage includes every canonical non-document payload, without confusing
    # the profile-configuration JSON checksum with the external HMM byte checksum.
    registered = {entry.artifact.path for entry in registry.document.references}
    data_files = {
        str(p.relative_to(REFERENCES))
        for folder in ("sprb", "hmm")
        for p in (REFERENCES / folder).rglob("*")
        if p.is_file() and p.suffix != ".md"
    }
    assert data_files <= registered


def test_canonical_sprb_panel_all_74_intervals_rn28_and_legacy_clusters():
    registry = ReferenceRegistry.load(REFERENCES / "registry.json")
    sequences = registry.sequences("SPRB_LEGACY_PANEL")
    parent = next(
        p
        for p in SeqIO.parse(ROOT / "tests/fixtures/fj_window.faa", "fasta")
        if p.id == "WP_012023065.1"
    )
    catalog = json.loads((REFERENCES / "sprb/reference_catalog.json").read_text())
    assert len(sequences) == len(catalog) == 74
    assert len(registry.sequences("SPRB_REFERENCE_ONLY")) == 73
    with (REFERENCES / "sprb/archive/cluster_assignments.tsv").open() as handle:
        assignments = {r["sequence_id_full"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with (REFERENCES / "sprb/archive/cluster_summary.tsv").open() as handle:
        clusters = {r["cluster_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert len(clusters) == 46
    counts = Counter()
    for protein, row in zip(sequences, catalog, strict=True):
        assert protein.protein_id == row["query_id"]
        assert (
            protein.sequence
            == row["sequence"]
            == str(parent.seq[row["parent_start"] : row["parent_end"]])
        )
        assert hashlib.sha256(protein.sequence.encode()).hexdigest() == row["sequence_sha256"]
        original = row["source_module_row"]
        assert (row["parent_start"], row["parent_end"]) == (
            int(original["start"]) - 1,
            int(original["end"]),
        )
        assignment = assignments[row["member_label"]]
        assert row["source_assignment_row"] == assignment
        assert row["legacy_group"] == clusters[assignment["cluster"]]["group_name"]
        counts[assignment["cluster"]] += 1
    assert all(counts[key] == int(row["n_members"]) for key, row in clusters.items())
    (rn28,) = registry.sequences("RN28")
    target = registry.get("RN28").biological_target
    assert rn28.protein_id == "SMF_RN28__RN28" and rn28.length == 207
    assert rn28.sequence == str(parent.seq[681:888])
    assert (target.protein_interval.start, target.protein_interval.end) == (681, 888)
    assert (
        target.sequence_sha256 == "0e92a4f69219d39a2c7a50b104c922fed707b998555378e4a95ea0be3dddb42c"
    )
    assert {p.reference_id for p in registry.lineage("RN28")} == {
        "SPRB_LEGACY_PANEL",
        "SPRB_DOMAIN_TABLE",
    }
