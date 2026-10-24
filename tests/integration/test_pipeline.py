"""Small real-backend orchestration: no genomes/models/network are downloaded."""

import json
import os
import random
import shutil
from pathlib import Path

import pyhmmer
import pytest
import yaml
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import HMMError, InputError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import export_prepared, write_json
from sprb_scanner.pipeline import run_pipeline
from sprb_scanner.pipeline.identity import STAGE_MODELS
from sprb_scanner.pipeline.storage import inventory, workspace_lock


def read(path):
    return json.loads(Path(path).read_text())


@pytest.fixture
def small_pipeline(tmp_path):
    if shutil.which("mmseqs") is None:
        if os.environ.get("SPRB_REQUIRE_MMSEQS") == "1":
            pytest.fail("CI requires MMseqs2")
        pytest.skip("Install MMseqs2 for real pipeline integration tests")
    rng = random.Random(47)
    sequence = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=207))
    proteins = [sequence + "G" * 30 + sequence, "P" * 120]
    ds = GenomeDataset.model_validate(
        {
            "genome": {"genome_id": "synthetic", "input_format": "genbank"},
            "contigs": [
                {
                    "contig_id": "chr",
                    "genome_id": "synthetic",
                    "source_id": "chr",
                    "length": 5000,
                    "topology": "circular",
                }
            ],
            "genes": [
                {
                    "gene_id": f"g{i}",
                    "genome_id": "synthetic",
                    "contig_id": "chr",
                    "strand": -1,
                    "location_parts": [{"start": i * 2000, "end": i * 2000 + len(p) * 3}],
                }
                for i, p in enumerate(proteins)
            ],
            "proteins": [
                {
                    "protein_id": f"p{i}",
                    "gene_id": f"g{i}",
                    "genome_id": "synthetic",
                    "source_id": f"WP_{i}",
                    "sequence": p,
                }
                for i, p in enumerate(proteins)
            ],
        }
    )
    prep = tmp_path / "prepared"
    export_prepared(ds, prep, [], AppConfig())
    references = []
    abc = pyhmmer.easel.Alphabet.amino()
    for model in [m for ids in STAGE_MODELS.values() for m in ids]:
        if model == "RN28":
            path = tmp_path / "RN28.fasta"
            path.write_text(f">RN28_query\n{sequence}\n")
            kind, fmt = "sequence", "protein_fasta"
        else:
            seq = pyhmmer.easel.TextSequence(name=model, sequence=sequence[:90]).digitize(abc)
            hmm, _, _ = pyhmmer.plan7.Builder(abc).build(seq, pyhmmer.plan7.Background(abc))
            hmm.accession = model + ".1"
            hmm.cutoffs.gathering = (10.0, 10.0)
            path = tmp_path / (model + ".hmm")
            with path.open("wb") as handle:
                hmm.write(handle)
            kind, fmt = "profile_hmm", "hmm"
        references.append(
            {
                "reference_id": model,
                "version": "synthetic-v1",
                "resource_type": kind,
                "artifact": {"path": path.name, "format": fmt, "storage": "bundled"},
                "sha256": sha256(path),
                "source": {"provider": "test fixture", "locator": "synthetic"},
                "biological_target": {"description": "synthetic test resource"},
            }
        )
    registry = tmp_path / "registry.json"
    write_json(registry, {"schema_version": "1.0.0", "references": references})
    config = tmp_path / "pipeline.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0.0",
                "registry": registry.name,
                "mmseqs": {"sensitivity": 7.5, "threads": 1},
            }
        )
    )
    return {
        "prepared": prep,
        "dataset": ds,
        "config": config,
        "registry": registry,
        "directory": tmp_path,
    }


def actions(report):
    return {k: v["action"] for k, v in report["stages"].items()}


def edit_config(data, **updates):
    raw = yaml.safe_load(data["config"].read_text())
    raw.update(updates)
    data["config"].write_text(yaml.safe_dump(raw))


def evidence(report):
    return read(Path(report["stages"]["evidence"]["result_directory"]) / "evidence_profiles.json")


def semantic_evidence(report):
    result = evidence(report)
    for profile in result["profiles"]:
        for source in profile["sources"].values():
            for domain in source["domains"]:
                domain.pop("scan_directory")
    return result["profiles"]


def test_full_pipeline_native_exports_resume_and_semantic_equivalence(small_pipeline, tmp_path):
    data = small_pipeline
    output = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], output)
    assert first["state"] == "succeeded"
    assert set(actions(first).values()) == {"executed"}
    assert sum(len(p["sources"]["RN28"]["alignments"]) for p in evidence(first)["profiles"]) >= 2
    before = inventory(output)
    second = run_pipeline(data["prepared"], data["config"], output, resume=True)
    assert set(actions(second).values()) == {"reused"}
    for item in first["stages"].values():
        path = Path(item["execution_manifest"])
        assert sha256(path) == before[str(path.relative_to(output))]
        assert (path.parent / "stage.log").is_file()
    fresh = run_pipeline(data["prepared"], data["config"], tmp_path / "fresh")
    assert {k: v["fingerprint"] for k, v in first["stages"].items()} == {
        k: v["fingerprint"] for k, v in fresh["stages"].items()
    }
    assert semantic_evidence(first) == semantic_evidence(second) == semantic_evidence(fresh)
    with pytest.raises(InputError, match="--resume"):
        run_pipeline(data["prepared"], data["config"], output)


def test_mmseqs_parameter_change_invalidates_only_rn28_and_evidence(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    edit_config(data, mmseqs={"sensitivity": 6.5, "threads": 1})
    second = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert actions(second) == {
        "typeb_ctd": "reused",
        "porp_sprf": "reused",
        "sprb_repeat": "reused",
        "rn28": "executed",
        "evidence": "executed",
    }
    assert first["stages"]["rn28"]["fingerprint"] != second["stages"]["rn28"]["fingerprint"]


def test_rn28_reference_version_and_content_selective_invalidation(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    registry = read(data["registry"])
    entry = registry["references"][-1]
    original = entry.copy()
    path = data["directory"] / "RN28-v2.fasta"
    sequence = (data["directory"] / "RN28.fasta").read_text().splitlines()[1]
    sequence = ("A" if sequence[0] != "A" else "C") + sequence[1:]
    path.write_text(f">RN28_v2\n{sequence}\n")
    entry.update(
        version="synthetic-v2",
        artifact={**entry["artifact"], "path": path.name},
        sha256=sha256(path),
    )
    registry["references"].append(original)
    write_json(data["registry"], registry)
    edit_config(data, references={"RN28": {"version": "synthetic-v2"}})
    second = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert actions(second) == {
        "typeb_ctd": "reused",
        "porp_sprf": "reused",
        "sprb_repeat": "reused",
        "rn28": "executed",
        "evidence": "executed",
    }
    assert first["stages"]["rn28"]["fingerprint"] != second["stages"]["rn28"]["fingerprint"]
    assert (
        evidence(second)["runs"]["RN28"]["metadata"]["rn28_reference"]["version"] == "synthetic-v2"
    )


def test_prepared_change_invalidates_all_stages(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    ds = data["dataset"]
    p = ds.proteins[0].model_copy(update={"sequence": "A" + ds.proteins[0].sequence[1:]})
    changed = ds.model_copy(update={"proteins": (p, ds.proteins[1])})
    prep2 = tmp_path / "prepared-v2"
    export_prepared(changed, prep2, [], AppConfig())
    second = run_pipeline(prep2, data["config"], out, resume=True)
    assert set(actions(second).values()) == {"executed"}
    assert all(
        first["stages"][k]["fingerprint"] != v["fingerprint"] for k, v in second["stages"].items()
    )


@pytest.mark.parametrize("mode", ["strict", "partial"])
def test_failure_status_logs_safe_retry_and_strict_vs_partial(
    small_pipeline, tmp_path, monkeypatch, mode
):
    from sprb_scanner.pipeline import runner

    data = small_pipeline
    out = tmp_path / "run"
    edit_config(data, mode=mode)
    original = runner.run_search

    def broken(*args, **kwargs):
        raise HMMError("synthetic backend failure")

    monkeypatch.setattr(runner, "run_search", broken)
    failed = run_pipeline(data["prepared"], data["config"], out)
    assert failed["state"] == ("failed" if mode == "strict" else "partial")
    rn = failed["stages"]["rn28"]
    attempt = Path(rn["attempt_directory"])
    assert (
        rn["state"] == "failed"
        and "synthetic backend failure" in (attempt / "stage.log").read_text()
    )
    old_manifest = (attempt / "execution_manifest.json").read_bytes()
    if mode == "strict":
        assert failed["stages"]["evidence"]["state"] == "blocked"
    else:
        assert evidence(failed)["runs"]["RN28"]["state"] == "failed"
        assert evidence(failed)["runs"]["RN28"]["retained_hits"] is None
        assert all(p["sources"]["RN28"]["state"] == "failed" for p in evidence(failed)["profiles"])
    monkeypatch.setattr(runner, "run_search", original)
    resumed = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert resumed["state"] == "succeeded"
    assert resumed["stages"]["rn28"]["attempt_directory"] != rn["attempt_directory"]
    assert all(
        resumed["stages"][k]["action"] == "reused"
        for k in ["typeb_ctd", "porp_sprf", "sprb_repeat"]
    )
    assert (attempt / "execution_manifest.json").read_bytes() == old_manifest


def test_corrupted_success_is_not_reused_or_overwritten(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    old = first["stages"]["sprb_repeat"]
    result = Path(old["result_directory"])
    (result / "hits.json").write_text("[]")
    second = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert actions(second) == {
        "typeb_ctd": "reused",
        "porp_sprf": "reused",
        "sprb_repeat": "executed",
        "rn28": "reused",
        "evidence": "executed",
    }
    assert second["stages"]["sprb_repeat"]["fingerprint"] == old["fingerprint"]
    assert second["stages"]["sprb_repeat"]["attempt_directory"] != old["attempt_directory"]
    assert (result / "hits.json").read_text() == "[]"
    assert second["stages"]["sprb_repeat"]["cache_rejections"]


def test_backend_identity_change_invalidates_only_its_scanner(
    small_pipeline, tmp_path, monkeypatch
):
    from sprb_scanner.pipeline import identity

    data = small_pipeline
    out = tmp_path / "run"
    run_pipeline(data["prepared"], data["config"], out)
    original = identity.backend_identity

    def changed(stage, config):
        value = original(stage, config)
        if stage == "rn28":
            value["build_identity"] = "changed"
        return value

    monkeypatch.setattr(identity, "backend_identity", changed)
    result = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert actions(result) == {
        "typeb_ctd": "reused",
        "porp_sprf": "reused",
        "sprb_repeat": "reused",
        "rn28": "executed",
        "evidence": "executed",
    }


def test_interrupted_attempt_is_preserved_and_retried(small_pipeline, tmp_path, monkeypatch):
    from sprb_scanner.pipeline import runner

    data = small_pipeline
    out = tmp_path / "run"
    original = runner.run_search

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt("test interruption")

    monkeypatch.setattr(runner, "run_search", interrupted)
    with pytest.raises(KeyboardInterrupt):
        run_pipeline(data["prepared"], data["config"], out)
    crashed = read(read(out / "latest.json")["run_manifest"])
    assert crashed["stages"]["rn28"]["state"] == "interrupted"
    assert crashed["stages"]["evidence"]["state"] == "pending"
    files = list(out.rglob("execution_manifest.json"))
    interrupted_manifest = next(p for p in files if read(p)["state"] == "interrupted")
    old_bytes = interrupted_manifest.read_bytes()
    monkeypatch.setattr(runner, "run_search", original)
    result = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert result["state"] == "succeeded"
    assert result["stages"]["rn28"]["action"] == "executed"
    assert result["stages"]["typeb_ctd"]["action"] == "reused"
    assert interrupted_manifest.read_bytes() == old_bytes


def test_successful_zero_hits_distinct_from_disabled_sources(small_pipeline, tmp_path):
    data = small_pipeline
    ds = data["dataset"]
    # MMseqs2 15 can reject an entirely masked target database rather than
    # report an empty search. Use searchable unrelated sequences for the
    # successful-zero-hit case; backend failures must remain failures.
    rng = random.Random(999)
    unrelated = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=150))
    negative = ds.model_copy(
        update={
            "proteins": tuple(p.model_copy(update={"sequence": unrelated}) for p in ds.proteins)
        }
    )
    prep = tmp_path / "negative"
    export_prepared(negative, prep, [], AppConfig())
    report = run_pipeline(prep, data["config"], tmp_path / "run")
    logs = (
        ""
        if report["state"] == "succeeded"
        else "\n".join(p.read_text() for p in (tmp_path / "run").rglob("backend.log"))
    )
    assert report["state"] == "succeeded", {"report": report, "backend_logs": logs}
    merged = evidence(report)
    assert not merged["profiles"]
    assert all(
        r["state"] == "succeeded" and r["retained_hits"] == 0 for r in merged["runs"].values()
    )
    edit_config(data, stages={key: False for key in STAGE_MODELS}, neighborhood={"enabled": False})
    disabled = run_pipeline(data["prepared"], data["config"], tmp_path / "disabled")
    assert disabled["state"] == "succeeded"
    assert all(
        r["state"] == "not_run" and r["retained_hits"] is None
        for r in evidence(disabled)["runs"].values()
    )


def test_cli_and_concurrent_lock(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    cli = CliRunner()
    args = [
        "run",
        "--prepared",
        str(data["prepared"]),
        "--config",
        str(data["config"]),
        "--output",
        str(out),
    ]
    result = cli.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert cli.invoke(app, args).exit_code == 2
    assert cli.invoke(app, [*args, "--resume"]).exit_code == 0
    with workspace_lock(out), pytest.raises(InputError, match="Another pipeline"):
        run_pipeline(data["prepared"], data["config"], out, resume=True)


def test_missing_external_model_partial_is_failed_not_not_run(small_pipeline, tmp_path):
    data = small_pipeline
    edit_config(data, mode="partial", references={"PF13573": {"artifact": "absent.hmm"}})
    # Explicit overrides require external registry storage, rather than weakening Registry rules.
    raw = read(data["registry"])
    next(e for e in raw["references"] if e["reference_id"] == "PF13573")["artifact"]["storage"] = (
        "external"
    )
    write_json(data["registry"], raw)
    report = run_pipeline(data["prepared"], data["config"], tmp_path / "run")
    assert report["state"] == "partial"
    assert report["stages"]["sprb_repeat"]["state"] == "failed"
    assert evidence(report)["runs"]["PF13573"]["state"] == "failed"
    assert "missing external artifact" in report["stages"]["sprb_repeat"]["error"]


@pytest.mark.parametrize("all_failed", [False, True])
def test_partial_native_pair_preserves_independent_model_status(
    small_pipeline, tmp_path, monkeypatch, all_failed
):
    from sprb_scanner.pipeline import runner

    data = small_pipeline
    edit_config(data, mode="partial")
    original = runner.run_bctd

    def failing_pair(proteins, profiles, output, config):
        changed = tuple(
            p.model_copy(update={"path": tmp_path / "unavailable-model.hmm"})
            if all_failed or i == 1
            else p
            for i, p in enumerate(profiles)
        )
        return original(proteins, changed, output, config)

    monkeypatch.setattr(runner, "run_bctd", failing_pair)
    report = run_pipeline(data["prepared"], data["config"], tmp_path / "run")
    assert report["state"] == "partial"
    assert report["stages"]["typeb_ctd"]["state"] == ("failed" if all_failed else "partial")
    assert report["stages"]["typeb_ctd"]["result_usable"]
    sources = evidence(report)["runs"]
    assert sources["TIGR04131"]["state"] == ("failed" if all_failed else "succeeded")
    assert sources["PF13585"]["state"] == "failed"
    assert sources["PF13585"]["retained_hits"] is None


def test_incomplete_attempt_is_retried_without_overwriting(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    old = Path(first["stages"]["rn28"]["execution_manifest"])
    manifest = read(old)
    manifest["state"] = "running"
    write_json(old, manifest)
    before = old.read_bytes()
    second = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert second["stages"]["rn28"]["action"] == "executed"
    assert second["stages"]["rn28"]["execution_manifest"] != str(old)
    assert old.read_bytes() == before


def test_valid_checksums_do_not_hide_native_parameter_mismatch(small_pipeline, tmp_path):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    item = first["stages"]["rn28"]
    root = Path(item["result_directory"])
    metadata = read(root / "run_metadata.json")
    metadata["configuration"]["evalue"] = 1.0
    write_json(root / "run_metadata.json", metadata)
    manifest_path = Path(item["execution_manifest"])
    manifest = read(manifest_path)
    manifest["output_sha256"] = inventory(manifest_path.parent, ("execution_manifest.json",))
    write_json(manifest_path, manifest)
    second = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert second["stages"]["rn28"]["action"] == "executed"
    assert any(
        "different parameters" in x["reason"] for x in second["stages"]["rn28"]["cache_rejections"]
    )


@pytest.mark.parametrize(("mode", "exit_code"), [("strict", 2), ("partial", 3)])
def test_cli_exposes_incomplete_run_exit_codes(
    small_pipeline, tmp_path, monkeypatch, mode, exit_code
):
    from sprb_scanner.pipeline import runner

    data = small_pipeline
    edit_config(data, mode=mode)

    def failed(*args, **kwargs):
        raise HMMError("test failed stage")

    monkeypatch.setattr(runner, "run_search", failed)
    result = CliRunner().invoke(
        app,
        [
            "run",
            "--prepared",
            str(data["prepared"]),
            "--config",
            str(data["config"]),
            "--output",
            str(tmp_path / "run"),
        ],
    )
    assert result.exit_code == exit_code, result.output


def test_merger_failure_adapter_rejects_conflicting_success_output(small_pipeline, tmp_path):
    from sprb_scanner.evidence.models import SourceRun
    from sprb_scanner.evidence.service import run_merge

    data = small_pipeline
    report = run_pipeline(data["prepared"], data["config"], tmp_path / "run")
    with pytest.raises(InputError, match="conflicts"):
        run_merge(
            data["prepared"],
            tmp_path / "conflict",
            rn28_scan=Path(report["stages"]["rn28"]["result_directory"]),
            registry=data["registry"],
            source_states={
                "RN28": SourceRun(
                    source_id="RN28", state="failed", message="do not overwrite a valid source"
                )
            },
        )


def test_neighborhood_window_invalidates_only_final_evidence_and_keeps_config_snapshot(
    small_pipeline, tmp_path
):
    data = small_pipeline
    out = tmp_path / "run"
    first = run_pipeline(data["prepared"], data["config"], out)
    original = Path(first["configuration_snapshot"]).read_bytes()
    edit_config(data, neighborhood={"enabled": True, "window": 1})
    second = run_pipeline(data["prepared"], data["config"], out, resume=True)
    assert actions(second) == {
        "typeb_ctd": "reused",
        "porp_sprf": "reused",
        "sprb_repeat": "reused",
        "rn28": "reused",
        "evidence": "executed",
    }
    assert Path(first["configuration_snapshot"]).read_bytes() == original
    assert sha256(Path(first["configuration_snapshot"])) == first["configuration_sha256"]


def test_special_genome_id_cannot_escape_workspace(small_pipeline, tmp_path):
    data = small_pipeline
    ds = data["dataset"]
    ds = ds.model_copy(
        update={
            "genome": ds.genome.model_copy(update={"genome_id": ".."}),
            "contigs": tuple(c.model_copy(update={"genome_id": ".."}) for c in ds.contigs),
            "genes": tuple(g.model_copy(update={"genome_id": ".."}) for g in ds.genes),
            "proteins": tuple(p.model_copy(update={"genome_id": ".."}) for p in ds.proteins),
        }
    )
    prep = tmp_path / "dots"
    export_prepared(ds, prep, [], AppConfig())
    edit_config(data, stages={key: False for key in STAGE_MODELS}, neighborhood={"enabled": False})
    out = tmp_path / "run"
    report = run_pipeline(prep, data["config"], out)
    assert report["state"] == "succeeded"
    assert Path(report["stages"]["evidence"]["result_directory"]).is_relative_to(out / "genomes")


def test_source_changes_during_execution_fail_closed(small_pipeline, tmp_path, monkeypatch):
    from sprb_scanner.pipeline import runner

    data = small_pipeline
    out = tmp_path / "run"
    original = runner._execute
    changed = False

    def mutate(spec, *args, **kwargs):
        nonlocal changed
        result = original(spec, *args, **kwargs)
        if not changed:
            data["registry"].write_text(data["registry"].read_text() + "\n")
            changed = True
        return result

    monkeypatch.setattr(runner, "_execute", mutate)
    report = run_pipeline(data["prepared"], data["config"], out)
    assert report["state"] == "failed"
    assert not report["stages"]["typeb_ctd"]["result_usable"]
    assert report["stages"]["typeb_ctd"]["state"] == "failed"
    assert report["stages"]["evidence"]["state"] == "blocked"
