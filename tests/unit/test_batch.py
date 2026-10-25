"""Immutable inventories, scheduler mapping and raw data validation contracts."""

import json
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml
from typer.testing import CliRunner

from sprb_scanner.batch import batch_status, collect_batch, plan_batch, run_task, submit_batch
from sprb_scanner.batch.models import SlurmResources
from sprb_scanner.batch.plan import load_manifest
from sprb_scanner.batch.slurm import array_expression, parse_indexes, render_script
from sprb_scanner.cli import app
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.pipeline.storage import workspace_lock

ROOT = Path(__file__).parents[2]


@pytest.fixture
def batch(tmp_path, genome_files):
    data = tmp_path / "data"
    for accession in ["GCF_000000002.1", "GCF_000000001.1"]:
        folder = data / accession
        folder.mkdir(parents=True)
        for key in ["gff", "fna", "faa"]:
            shutil.copyfile(genome_files[key], folder / ("genomic." + key))
    config = tmp_path / "pipeline.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0.0",
                "registry": str(ROOT / "references/registry.json"),
                "stages": dict.fromkeys(["typeb_ctd", "porp_sprf", "sprb_repeat", "rn28"], False),
                "neighborhood": {"enabled": False},
            }
        )
    )
    resource = tmp_path / "resources.json"
    settings = {
        "partition": "test",
        "account": "test",
        "qos": "test",
        "cpus": 2,
        "memory": "4G",
        "time": "00:10:00",
        "concurrency": 1,
        "max_array_tasks": 2,
        "max_array_size": 1001,
        "python": "/test/python",
        "modules": ["mmseqs/14"],
    }
    resource.write_text(json.dumps(settings))
    output = tmp_path / "batch"

    def plan():
        return plan_batch(data, config, resource, output)

    return {
        "data": data,
        "config": config,
        "resources": resource,
        "settings": settings,
        "output": output,
        "manifest": output / "manifest.json",
        "plan": plan,
    }


def test_stable_inventory_and_config_snapshot(batch):
    batch["plan"]()
    plan = load_manifest(batch["manifest"])
    assert [e.assembly_accession for e in plan.entries] == ["GCF_000000001.1", "GCF_000000002.1"]
    assert [e.index for e in plan.entries] == [0, 1]
    assert all(len(e.files) == 3 and e.preparation_status == "not_run" for e in plan.entries)
    batch["config"].write_text("changed")
    assert load_manifest(batch["manifest"]).pipeline_config == str(
        batch["output"] / "pipeline.yaml"
    )
    with pytest.raises(InputError, match="already exists"):
        # Original source YAML changed above; use the frozen valid one.
        plan_batch(
            batch["data"], batch["output"] / "pipeline.yaml", batch["resources"], batch["output"]
        )


def test_manifest_tampering_is_rejected(batch):
    batch["plan"]()
    batch["manifest"].write_text(batch["manifest"].read_text() + " ")
    with pytest.raises(InputError, match="checksum"):
        load_manifest(batch["manifest"])


@pytest.mark.parametrize("index", [-1, 2, True])
def test_invalid_index_and_allocation_guard(batch, index):
    batch["plan"]()
    with pytest.raises(InputError):
        run_task(batch["manifest"], index, local=True)
    with patch.dict("os.environ", {}, clear=True), pytest.raises(InputError, match="allocation"):
        run_task(batch["manifest"], 0)


def test_missing_input_is_recorded_not_dropped(batch):
    (batch["data"] / "GCF_000000001.1/genomic.faa").unlink()
    batch["plan"]()
    assert load_manifest(batch["manifest"]).entries[0].inventory_state == "invalid"
    report = run_task(batch["manifest"], 0, local=True)
    assert report["state"] == "failed" and report["preparation_state"] == "failed"
    assert "Expected exactly one" in report["error"]
    assert run_task(batch["manifest"], 1, prepare_only=True, local=True)["state"] == "prepared"


def test_input_mutation_is_not_repaired_or_scanned(batch):
    batch["plan"]()
    raw = batch["data"] / "GCF_000000001.1/genomic.faa"
    raw.write_text(">p1\nINVALID*SEQUENCE\n")
    report = run_task(batch["manifest"], 0, local=True)
    assert report["state"] == "failed" and "input changed" in report["error"]
    assert not Path(load_manifest(batch["manifest"]).entries[0].prepared_path).exists()
    assert raw.read_text() == ">p1\nINVALID*SEQUENCE\n"


def test_annotation_error_isolated_and_strict(batch):
    path = batch["data"] / "GCF_000000001.1/genomic.gff"
    path.write_text(path.read_text().replace("protein_id=p1", "protein_id=missing"))
    batch["plan"]()
    failed = run_task(batch["manifest"], 0, local=True)
    assert failed["state"] == "failed"
    assert not Path(load_manifest(batch["manifest"]).entries[0].prepared_path).exists()
    assert run_task(batch["manifest"], 1, local=True)["state"] == "succeeded"


def test_prepare_resume_and_compact_zero_hit_collection(batch, tmp_path):
    batch["plan"]()
    first = run_task(batch["manifest"], 0, prepare_only=True, local=True)
    assert first["state"] == "prepared"
    prepared = Path(load_manifest(batch["manifest"]).entries[0].prepared_path)
    before = (prepared / "dataset.json").read_bytes()
    second = run_task(batch["manifest"], 0, local=True)
    assert second["state"] == "succeeded" and second["preparation_action"] == "reused"
    third = run_task(batch["manifest"], 0, local=True)
    from sprb_scanner.batch.status import pipeline_report

    plan = load_manifest(batch["manifest"])
    assert pipeline_report(plan, plan.entries[0], third)["stages"]["evidence"]["action"] == "reused"
    assert (prepared / "dataset.json").read_bytes() == before
    assert first["disk_usage"]["files"] > 0
    assert batch_status(batch["manifest"], scheduler=False, verify=True)["counts"] == {
        "succeeded": 1,
        "not_run": 1,
    }
    out = tmp_path / "collected"
    result = collect_batch(batch["manifest"], out)
    assert result["collected_genomes"] == 1 and result["excluded_genomes"] == 1
    assert result["profiles"] == 0
    assert (out / "evidence_records.tsv").exists()
    with pytest.raises(InputError, match="already exists"):
        collect_batch(batch["manifest"], out)


def test_corrupt_prepared_cannot_be_silently_overwritten(batch):
    batch["plan"]()
    run_task(batch["manifest"], 0, prepare_only=True, local=True)
    prepared = Path(load_manifest(batch["manifest"]).entries[0].prepared_path)
    (prepared / "proteome.faa").write_text("corrupt")
    result = run_task(batch["manifest"], 0, local=True)
    assert result["state"] == "failed" and "checksum" in result["error"]
    assert (prepared / "proteome.faa").read_text() == "corrupt"


def test_task_lock_prevents_duplicate_execution(batch):
    batch["plan"]()
    root = batch["output"] / "genomes/GCF_000000001.1"
    root.mkdir(parents=True)
    with workspace_lock(root), pytest.raises(InputError, match="Another pipeline"):
        run_task(batch["manifest"], 0, local=True)


@pytest.mark.parametrize("value", ["-1", "2", "1-0", "0;touch /tmp/x", "0\n1", "0-10000000000"])
def test_invalid_array_indexes(value):
    with pytest.raises(InputError):
        parse_indexes(value, 2)


def test_sparse_retry_mapping_and_explicit_resources(batch):
    batch["plan"]()
    manifest = load_manifest(batch["manifest"])
    assert parse_indexes("0,1,0", 2) == [0, 1]
    assert array_expression([0, 1, 2, 5, 9, 10]) == "0-2,5,9-10"
    script = render_script(manifest, batch["manifest"], [1])
    assert "#SBATCH --array=1%1" in script
    assert "--account=test" in script and "--cpus-per-task=2" in script
    assert '"$SLURM_ARRAY_TASK_ID"' in script and "TMPDIR=" in script
    run_task(batch["manifest"], 0, local=True)
    dry = submit_batch(batch["manifest"], dry_run=True)
    assert dry["indexes"] == [1]
    with pytest.raises(InputError, match="duplicate"):
        submit_batch(batch["manifest"], indexes="0", dry_run=True)
    assert submit_batch(batch["manifest"], retry=True, dry_run=True)["state"] == "nothing_to_submit"


def test_resource_contract_has_no_universal_defaults(batch):
    with pytest.raises(ValueError):
        SlurmResources.model_validate({})
    with pytest.raises(ValueError):
        SlurmResources.model_validate(batch["settings"] | {"partition": "x\n#SBATCH --x"})
    with pytest.raises(ValueError):
        SlurmResources.model_validate(batch["settings"] | {"max_array_size": 0})
    batch["resources"].write_text(json.dumps(batch["settings"] | {"max_array_size": 1}))
    with pytest.raises(ValueError, match="MaxArraySize"):
        batch["plan"]()


def test_cli_registration_and_local_guard(batch):
    runner = CliRunner()
    for name in [
        "batch-plan",
        "batch-task",
        "batch-status",
        "batch-submit",
        "batch-retry",
        "batch-collect",
        "batch-disk-usage",
    ]:
        assert runner.invoke(app, [name, "--help"]).exit_code == 0
    batch["plan"]()
    assert (
        runner.invoke(
            app, ["batch-task", "--manifest", str(batch["manifest"]), "--index", "0", "--local"]
        ).exit_code
        == 0
    )


def test_generated_script_executes_installed_cli(batch, tmp_path):
    import os
    import subprocess
    import sys

    settings = batch["settings"] | {"python": sys.executable, "modules": []}
    batch["resources"].write_text(json.dumps(settings))
    batch["plan"]()
    script = tmp_path / "job.sh"
    script.write_text(render_script(load_manifest(batch["manifest"]), batch["manifest"], [0]))
    result = subprocess.run(
        ["bash", str(script)],
        capture_output=True,
        text=True,
        env=os.environ
        | {
            "SLURM_JOB_ID": "123",
            "SLURM_ARRAY_JOB_ID": "123",
            "SLURM_ARRAY_TASK_ID": "0",
            "SLURM_CPUS_PER_TASK": "2",
        },
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["state"] == "succeeded"


def test_scheduler_uses_array_jobid_not_numeric_jobidraw(batch, monkeypatch):
    from types import SimpleNamespace

    from sprb_scanner.batch import status

    batch["plan"]()
    submit = batch["output"] / "submissions/submission-0001"
    submit.mkdir(parents=True)
    (submit / "submission.json").write_text(
        json.dumps({"indexes": [0, 1], "job_id": "100", "created_at": "2000-01-01"})
    )
    monkeypatch.setattr(status.shutil, "which", lambda _: "/usr/bin/sacct")

    def sacct(command, **kwargs):
        assert "--format=JobID%80,State%40,ExitCode" in command
        return SimpleNamespace(returncode=0, stdout="100_0|FAILED|1:0\n100_1|PENDING|0:0\n")

    monkeypatch.setattr(status.subprocess, "run", sacct)
    report = batch_status(batch["manifest"])
    assert report["counts"] == {"failed": 1, "pending": 1}
    retry = submit_batch(batch["manifest"], retry=True, dry_run=True)
    assert retry["indexes"] == [0]


def test_ambiguous_submission_never_duplicates(batch, monkeypatch):
    import subprocess

    from sprb_scanner.batch import slurm

    batch["plan"]()

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired("sbatch", 60)

    monkeypatch.setattr(slurm.subprocess, "run", timeout)
    with pytest.raises(subprocess.TimeoutExpired):
        submit_batch(batch["manifest"])
    receipt = json.loads(
        next((batch["output"] / "submissions").glob("*/submission.json")).read_text()
    )
    assert receipt["state"] == "submission_unknown"
    assert submit_batch(batch["manifest"], dry_run=True)["state"] == "active_wave"


def test_submission_success_and_terminal_scheduler_failure(batch, monkeypatch):
    from types import SimpleNamespace

    from sprb_scanner.batch import slurm
    from sprb_scanner.batch import status as service

    batch["plan"]()
    calls = []

    def dispatch(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=0, stdout="777;cluster\n", stderr="")

    monkeypatch.setattr(slurm.subprocess, "run", dispatch)
    receipt = submit_batch(batch["manifest"], indexes="1")
    assert receipt["job_id"] == "777" and receipt["indexes"] == [1]
    assert calls[0][0] == "sbatch"
    monkeypatch.setattr(service.shutil, "which", lambda _: "/sacct")
    monkeypatch.setattr(
        service.subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(
            returncode=0, stdout="777_1|OUT_OF_MEMORY|0:9\n", stderr=""
        ),
    )
    report = batch_status(batch["manifest"])
    assert report["entries"][1]["state"] == "failed"
    assert submit_batch(batch["manifest"], retry=True, dry_run=True)["indexes"] == [1]


def test_prepared_corruption_excluded_without_zero_filling(batch, tmp_path):
    import csv

    batch["plan"]()
    run_task(batch["manifest"], 0, local=True)
    prep = Path(load_manifest(batch["manifest"]).entries[0].prepared_path)
    (prep / "proteome.faa").write_text("corrupt")
    result = collect_batch(batch["manifest"], tmp_path / "excluded")
    assert result["collected_genomes"] == 0
    with (tmp_path / "excluded/genome_summary.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["state"] == "invalid_result"
    assert rows[0]["rn28_proteins"] == "" and rows[0]["TIGR04131_retained_records"] == ""


def test_zero_profiles_with_disabled_sources_keep_counts_unknown(batch, tmp_path):
    import csv

    batch["plan"]()
    run_task(batch["manifest"], 0, local=True)
    collect_batch(batch["manifest"], tmp_path / "disabled")
    with (tmp_path / "disabled/genome_summary.tsv").open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["evidence_profiles"] == "0"
    assert row["typeb_retained_proteins"] == ""
    assert row["local_hit_typeb_cargos"] == ""
    assert row["RN28_state"] == "not_run" and row["RN28_retained_records"] == ""


def test_scheduler_outage_remains_unknown(batch, monkeypatch):
    import subprocess

    from sprb_scanner.batch import status

    batch["plan"]()
    receipt = batch["output"] / "submissions/submission-0001"
    receipt.mkdir(parents=True)
    (receipt / "submission.json").write_text(
        json.dumps({"indexes": [0], "job_id": "123", "created_at": "2000-01-01"})
    )
    monkeypatch.setattr(status.shutil, "which", lambda _: "/sacct")

    def down(*args, **kwargs):
        raise subprocess.TimeoutExpired("sacct", 60)

    monkeypatch.setattr(status.subprocess, "run", down)
    assert batch_status(batch["manifest"])["entries"][0]["state"] == "submitted_unknown"


def test_later_single_task_supersedes_old_failed_array(batch, monkeypatch):
    from types import SimpleNamespace

    from sprb_scanner.batch import status

    batch["plan"]()
    receipt = batch["output"] / "submissions/submission-0001"
    receipt.mkdir(parents=True)
    (receipt / "submission.json").write_text(
        json.dumps({"indexes": [0], "job_id": "123", "created_at": "2000-01-01"})
    )
    run_task(batch["manifest"], 0, local=True)
    monkeypatch.setattr(status.shutil, "which", lambda _: "/sacct")
    monkeypatch.setattr(
        status.subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(returncode=0, stdout="123_0|FAILED|2:0\n"),
    )
    assert batch_status(batch["manifest"], verify=True)["entries"][0]["state"] == "succeeded"


def test_unknown_new_submission_does_not_reuse_old_failed_task(batch):
    from sprb_scanner.pipeline.storage import now

    batch["plan"]()
    raw = batch["data"] / "GCF_000000001.1/genomic.faa"
    raw.write_text("changed")
    assert run_task(batch["manifest"], 0, local=True)["state"] == "failed"
    receipt = batch["output"] / "submissions/submission-0001"
    receipt.mkdir(parents=True)
    (receipt / "submission.json").write_text(
        json.dumps({"indexes": [0], "state": "submission_unknown", "created_at": now()})
    )
    assert batch_status(batch["manifest"])["entries"][0]["state"] == "submitted_unknown"
    assert submit_batch(batch["manifest"], retry=True, dry_run=True)["state"] == "nothing_to_submit"
