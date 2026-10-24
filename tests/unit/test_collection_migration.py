"""Version-aware collection migration never replaces previously successful evidence."""

import csv
import json
from pathlib import Path

import pytest

from sprb_scanner.batch.combine import SOURCE_FIELD, combine_collections
from sprb_scanner.batch.models import BatchEntry, BatchManifest, InputFile, SlurmResources
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256

A, B = "GCF_000000001.1", "GCF_000000002.1"


def table(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def snapshot(tmp_path, name, states):
    raw = tmp_path / "input.gff"
    if not raw.exists():
        raw.write_text("unchanged original input")
    batch = tmp_path / name
    batch.mkdir()
    resources = SlurmResources(
        partition="test",
        account="test",
        qos="test",
        cpus=1,
        memory="1G",
        time="00:01:00",
        concurrency=1,
        max_array_tasks=2,
        max_array_size=1001,
        python="/project/python",
    )
    entries = [
        BatchEntry(
            index=i,
            assembly_accession=gid,
            genome_id=gid,
            inventory_state="ready",
            files={
                k: InputFile(path=str(raw), size_bytes=raw.stat().st_size, sha256=sha256(raw))
                for k in ["gff", "fna", "faa"]
            },
            prepared_path=str(batch / "genomes" / gid / "prepared"),
        )
        for i, gid in enumerate(sorted(states))
    ]
    plan = BatchManifest(
        created_at="2026-10-24",
        output_root=str(batch),
        dataset_root=str(tmp_path),
        software={"version": name},
        pipeline_config=str(batch / "pipeline.yaml"),
        pinned_files={},
        resources=resources,
        entries=entries,
    )
    manifest = batch / "manifest.json"
    manifest.write_text(plan.model_dump_json())
    manifest.with_suffix(".sha256").write_text(sha256(manifest))
    root = tmp_path / (name + "-collection")
    root.mkdir()
    base = [
        {
            "genome_id": gid,
            "assembly_accession": gid,
            "array_index": str(i),
            "state": states[gid],
            "count": "1" if states[gid] == "succeeded" else "",
        }
        for i, gid in enumerate(sorted(states))
    ]
    table(root / "genome_summary.tsv", list(base[0]), base)
    records = [
        {
            "genome_id": gid,
            "hit_id": f"{gid}:{i}",
            "raw_evidence": json.dumps({"coordinate": [i, i + 10], "score": 17.6}),
        }
        for gid in sorted(states)
        if states[gid] == "succeeded"
        for i in range(2)
    ]
    table(root / "evidence_records.tsv", ["genome_id", "hit_id", "raw_evidence"], records)
    table(
        root / "source_status.tsv",
        ["genome_id", "state"],
        [
            {"genome_id": gid, "state": "succeeded" if state == "succeeded" else "unavailable"}
            for gid, state in states.items()
        ],
    )
    (root / "provenance.jsonl").write_text(
        "".join(
            json.dumps({"genome_id": gid, "raw": name}) + "\n"
            for gid, state in states.items()
            if state == "succeeded"
        )
    )
    (root / "collection_metadata.json").write_text(
        json.dumps(
            {
                "manifest": str(manifest),
                "manifest_sha256": sha256(manifest),
                "software": plan.software,
                "pinned_files": plan.pinned_files,
                "output_sha256": {p.name: sha256(p) for p in root.iterdir()},
            }
        )
    )
    (root / "status.json").write_text(json.dumps({"state": "succeeded"}))
    return root


def test_recovery_combines_versions_retaining_all_native_fields_and_original_snapshots(tmp_path):
    old = snapshot(tmp_path, "old", {A: "succeeded", B: "failed"})
    new = snapshot(tmp_path, "new", {B: "succeeded"})
    old_files = {p: sha256(p) for root in [old, new] for p in root.iterdir()}
    out = tmp_path / "combined"
    result = combine_collections([old, new], out)
    assert result["genome_states"] == {"succeeded": 2}
    assert result["genomes"] == 2
    combined = read_table(out / "evidence_records.tsv")
    originals = read_table(old / "evidence_records.tsv") + read_table(new / "evidence_records.tsv")
    assert len(combined) == len(originals) == 4
    assert [{k: v for k, v in r.items() if k != SOURCE_FIELD} for r in combined] == originals
    assert len({r[SOURCE_FIELD] for r in combined}) == 2
    assert all(sha256(p) == expected for p, expected in old_files.items())
    meta = json.loads((out / "collection_metadata.json").read_text())
    assert meta["sources"][0]["software"] != meta["sources"][1]["software"]
    assert meta["superseded_unavailable"][0]["previous"]["count"] == ""
    for name, checksum in meta["output_sha256"].items():
        assert sha256(out / name) == checksum
    with pytest.raises(InputError, match="exists"):
        combine_collections([old, new], out)


@pytest.mark.parametrize(
    "failure",
    [
        "tamper",
        "successful_duplicate",
        "changed_input",
        "provenance_mismatch",
        "duplicate_manifest",
    ],
)
def test_incompatible_or_ambiguous_recovery_is_rejected_before_output(tmp_path, failure):
    old = snapshot(tmp_path, "old", {A: "succeeded", B: "failed"})
    new = snapshot(
        tmp_path, "new", {A: "succeeded"} if failure == "successful_duplicate" else {B: "succeeded"}
    )
    if failure == "tamper":
        (new / "evidence_records.tsv").write_text("corruption")
    if failure in {"changed_input", "provenance_mismatch"}:
        manifest = tmp_path / "new/manifest.json"
        data = json.loads(manifest.read_text())
        if failure == "changed_input":
            data["entries"][0]["files"]["gff"]["sha256"] = "f" * 64
        else:
            data["software"] = {"fake": "mismatched"}
        manifest.write_text(json.dumps(data))
        manifest.with_suffix(".sha256").write_text(sha256(manifest))
        metadata = new / "collection_metadata.json"
        m = json.loads(metadata.read_text())
        m["manifest_sha256"] = sha256(manifest)
        metadata.write_text(json.dumps(m))
    if failure == "duplicate_manifest":
        new = old
    with pytest.raises(InputError):
        combine_collections([old, new], tmp_path / "combined")
    assert not (tmp_path / "combined").exists()


def test_concurrent_snapshot_metadata_change_fails_without_claiming_success(tmp_path, monkeypatch):
    import sprb_scanner.batch.combine as migration

    old = snapshot(tmp_path, "old", {A: "succeeded", B: "failed"})
    new = snapshot(tmp_path, "new", {B: "succeeded"})
    original = migration._source
    calls = 0

    def changed(path):
        nonlocal calls
        calls += 1
        if calls == 3:
            file = Path(path) / "collection_metadata.json"
            data = json.loads(file.read_text())
            data["created_at"] = "changed while copying"
            file.write_text(json.dumps(data))
        return original(path)

    monkeypatch.setattr(migration, "_source", changed)
    with pytest.raises(InputError, match="concurrently"):
        migration.combine_collections([old, new], tmp_path / "combined")
    assert json.loads((tmp_path / "combined/status.json").read_text())["state"] == "failed"
