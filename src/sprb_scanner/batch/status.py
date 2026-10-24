"""Filesystem and optional SLURM status; missing results never become zero hits."""

import shutil
import subprocess
from collections import Counter
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256
from sprb_scanner.pipeline.storage import verify_files

from .plan import load_manifest, read_json
from .task import disk_usage, validate_preparation

TERMINAL_FAILURES = {
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "OUT_OF_MEMORY",
    "NODE_FAIL",
    "PREEMPTED",
    "BOOT_FAIL",
    "DEADLINE",
    "REVOKED",
}
ACTIVE = {
    "PENDING",
    "RUNNING",
    "CONFIGURING",
    "COMPLETING",
    "SUSPENDED",
    "REQUEUED",
    "REQUEUE_FED",
    "RESIZING",
    "SIGNALING",
    "STAGE_OUT",
}


def task_record(manifest, entry):
    root = Path(manifest.output_root) / "genomes" / entry.genome_id
    pointer = root / "latest_task.json"
    if not pointer.exists():
        return None
    ref = read_json(pointer)
    path = Path(ref["task"])
    if not path.resolve().is_relative_to(root / "tasks"):
        raise InputError("Task record escapes its genome directory")
    record = read_json(path)
    if record["state"] != "running" and sha256(path) != ref["sha256"]:
        raise InputError("Task record checksum mismatch")
    if (
        record["index"] != entry.index
        or record["genome_id"] != entry.genome_id
        or record["manifest_sha256"] != sha256(Path(manifest.output_root) / "manifest.json")
    ):
        raise InputError("Task record identity mismatch")
    return record


def pipeline_report(manifest, entry, record, *, verify=True):
    root = Path(manifest.output_root) / "genomes" / entry.genome_id / "pipeline"
    path = Path(record["pipeline_latest"])
    if path != root / "latest.json" or sha256(path) != record["pipeline_latest_sha256"]:
        raise InputError("Pipeline pointer changed since task completed")
    ref = read_json(path)
    run = Path(ref["run_manifest"])
    if not run.resolve().is_relative_to(root) or sha256(run) != ref["run_manifest_sha256"]:
        raise InputError("Pipeline run manifest checksum/path mismatch")
    data = read_json(run)
    if data["genome_id"] != entry.genome_id or data["state"] != record["pipeline_state"]:
        raise InputError("Pipeline run identity/state mismatch")
    if verify:
        for item in data["stages"].values():
            if not item.get("execution_manifest"):
                continue
            p = Path(item["execution_manifest"])
            if (
                not p.resolve().is_relative_to(root)
                or sha256(p) != item["execution_manifest_sha256"]
            ):
                raise InputError("Pipeline stage identity/checksum mismatch")
            stage = read_json(p)
            verify_files(p.parent, stage["output_sha256"])
    return data


def scheduler_states(manifest):
    submissions = []
    for path in sorted((Path(manifest.output_root) / "submissions").glob("*/submission.json")):
        item = read_json(path)
        if item.get("job_id") or item.get("state") in {"submitting", "submission_unknown"}:
            submissions.append(item)
    if not any(s.get("job_id") for s in submissions) or shutil.which("sacct") is None:
        return submissions, {}
    try:
        result = subprocess.run(
            [
                "sacct",
                "-n",
                "-P",
                "-j",
                ",".join(s["job_id"] for s in submissions if s.get("job_id")),
                "--format=JobID%80,State%40,ExitCode",
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return submissions, {}
    if result.returncode:
        return submissions, {}
    states = {}
    for line in result.stdout.splitlines():
        columns = line.split("|")
        if len(columns) >= 3:
            states[columns[0]] = {
                "state": columns[1].split()[0].rstrip("+"),
                "exit_code": columns[2],
            }
    return submissions, states


def batch_status(path, *, scheduler=True, verify=False):
    manifest = load_manifest(path)
    submissions, scheduler_info = scheduler_states(manifest) if scheduler else ([], {})
    latest_submissions = {}
    for s in submissions:
        for index in s["indexes"]:
            latest_submissions[index] = s
    rows = []
    for entry in manifest.entries:
        row = {
            "index": entry.index,
            "assembly_accession": entry.assembly_accession,
            "genome_id": entry.genome_id,
            "inventory_state": entry.inventory_state,
            "preparation_state": "not_run",
            "pipeline_state": "not_run",
            "state": "not_run",
            "error": entry.error,
            "verified": False,
        }
        try:
            record = task_record(manifest, entry)
            if record:
                row.update(
                    {
                        k: record.get(k)
                        for k in (
                            "state",
                            "preparation_state",
                            "pipeline_state",
                            "error",
                            "disk_usage",
                        )
                    }
                )
            s = latest_submissions.get(entry.index)
            # A later explicit single-task execution supersedes an older array receipt.
            if (
                s
                and record
                and not record["slurm"].get("SLURM_ARRAY_JOB_ID")
                and record["started_at"] >= s["created_at"]
            ):
                s = None
            if s:
                jid = f"{s.get('job_id', 'unknown')}_{entry.index}"
                scheduled = scheduler_info.get(jid, {})
                row.update(
                    slurm_job=jid,
                    slurm_state=scheduled.get("state"),
                    slurm_exit_code=scheduled.get("exit_code"),
                )
                matches = (
                    s.get("job_id")
                    and record
                    and record["slurm"].get("SLURM_ARRAY_JOB_ID") == s["job_id"]
                )
                if (
                    not matches
                    or row["state"] == "running"
                    or scheduled.get("state") in TERMINAL_FAILURES
                ):
                    state = scheduled.get("state")
                    row["state"] = (
                        "failed"
                        if state in TERMINAL_FAILURES or state == "COMPLETED"
                        else "running"
                        if state == "RUNNING"
                        else "pending"
                        if state in ACTIVE
                        else "submitted_unknown"
                    )
                    if row["state"] == "failed" and not row.get("error"):
                        row["error"] = (
                            "SLURM failed or terminated without a completed matching task record"
                        )
            if verify and row["state"] in {"succeeded", "partial", "prepared"}:
                from .plan import verify_inputs

                verify_inputs(entry)
                validate_preparation(entry)
                if (
                    sha256(Path(entry.prepared_path) / "run_metadata.json")
                    != record["prepared_metadata_sha256"]
                ):
                    raise InputError("Prepared metadata changed since task completed")
                if row["state"] != "prepared":
                    pipeline_report(manifest, entry, record)
                row["verified"] = True
            if row["state"] == "running" and not s:
                row["state"] = "interrupted_or_running_unknown"
        except (OSError, ValueError, KeyError, TypeError, InputError) as exc:
            row.update(state="invalid_result", error=str(exc), verified=False)
        rows.append(row)
    return {
        "manifest": str(Path(path).resolve()),
        "genomes": len(rows),
        "counts": dict(Counter(r["state"] for r in rows)),
        "entries": rows,
    }


def measure_batch(path):
    manifest = load_manifest(path)
    return [
        {
            "genome_id": e.genome_id,
            **disk_usage(Path(manifest.output_root) / "genomes" / e.genome_id),
        }
        for e in manifest.entries
    ]
