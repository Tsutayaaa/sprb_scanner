"""Explicit resource arrays; original immutable indexes survive retry submissions."""

import re
import shlex
import subprocess
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256
from sprb_scanner.pipeline.storage import atomic_json, next_directory, now, workspace_lock

from .plan import load_manifest
from .status import batch_status


def parse_indexes(value, count):
    if value is None:
        return None
    indexes = set()
    if not re.fullmatch(r"[0-9]+(?:-[0-9]+)?(?:,[0-9]+(?:-[0-9]+)?)*", value):
        raise InputError("Indexes must be comma-separated integers or inclusive ranges")
    for group in value.split(","):
        parts = list(map(int, group.split("-")))
        lo, hi = parts[0], parts[-1]
        if hi < lo or hi >= count:
            raise InputError("Array index range outside manifest")
        indexes.update(range(lo, hi + 1))
    return sorted(indexes)


def array_expression(indexes):
    ranges = []
    start = end = indexes[0]
    for index in indexes[1:]:
        if index == end + 1:
            end = index
        else:
            ranges.append(str(start) if start == end else f"{start}-{end}")
            start = end = index
    ranges.append(str(start) if start == end else f"{start}-{end}")
    return ",".join(ranges)


def render_script(manifest, manifest_path, indexes, *, prepare_only=False):
    r, root = manifest.resources, Path(manifest.output_root)
    # Shell-quote commands; SBATCH paths use quoted directive arguments, with control chars banned.
    if any(c in str(root) for c in "\n\r\0%"):
        raise InputError("Unsafe batch root for SLURM output directives")
    directives = {
        "job-name": "sprb-prepare" if prepare_only else "sprb-pipeline",
        "partition": r.partition,
        "account": r.account,
        "qos": r.qos,
        "nodes": "1",
        "ntasks": "1",
        "cpus-per-task": str(r.cpus),
        "mem": r.memory,
        "time": r.time,
        "array": array_expression(indexes) + f"%{r.concurrency}",
        "output": str(root / "logs" / "%x_%A_%a.out"),
        "error": str(root / "logs" / "%x_%A_%a.err"),
    }
    lines = [
        "#!/bin/bash",
        *[f"#SBATCH --{k}={shlex.quote(v)}" for k, v in directives.items()],
        "set -euo pipefail",
        f"cd {shlex.quote(str(root))}",
        'export TMPDIR="$PWD/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"',
        'mkdir -p "$TMPDIR"',
        "export PYTHONDONTWRITEBYTECODE=1",
        f"export OMP_NUM_THREADS={r.cpus}",
        "export OPENBLAS_NUM_THREADS=1",
    ]
    if r.modules:
        lines.append("module load " + " ".join(shlex.quote(m) for m in r.modules))
    command = [
        r.python,
        "-m",
        "sprb_scanner",
        "batch-task",
        "--manifest",
        str(manifest_path),
        "--index",
    ]
    lines.append(
        "exec "
        + shlex.join(command)
        + ' "$SLURM_ARRAY_TASK_ID"'
        + (" --prepare-only" if prepare_only else "")
    )
    return "\n".join(lines) + "\n"


def submit_batch(path, *, indexes=None, retry=False, prepare_only=False, dry_run=False, wait=False):
    path = Path(path).resolve()
    manifest = load_manifest(path)
    root = Path(manifest.output_root)
    # Also serializes submissions, preventing two clients from dispatching duplicate genomes.
    with workspace_lock(root):
        status = batch_status(path)
        requested = parse_indexes(indexes, len(manifest.entries))
        eligible = {"failed", "interrupted", "partial", "invalid_result"} if retry else {"not_run"}
        if not prepare_only and not retry:
            eligible.add("prepared")
        rows = status["entries"]
        active = [
            r["index"]
            for r in rows
            if r["state"]
            in {"running", "pending", "submitted_unknown", "interrupted_or_running_unknown"}
        ]
        available = max(0, manifest.resources.max_array_tasks - len(active))
        if not available:
            return {"state": "active_wave", "indexes": active, "counts": status["counts"]}
        if requested is None:
            selected = [r["index"] for r in rows if r["state"] in eligible][:available]
        else:
            selected = requested
            if any(rows[i]["state"] not in eligible for i in selected):
                raise InputError(
                    "Selected tasks are successful, active or ineligible; "
                    "refusing duplicate submission"
                )
        if not selected:
            return {"state": "nothing_to_submit", "indexes": [], "counts": status["counts"]}
        if len(selected) > available:
            raise InputError("Submission exceeds site-configured max_array_tasks; submit in waves")
        script = render_script(manifest, path, selected, prepare_only=prepare_only)
        if dry_run:
            return {"state": "dry_run", "indexes": selected, "script": script}
        directory = next_directory(root / "submissions", "submission-")
        (directory / "job.sh").write_text(script)
        receipt = {
            "created_at": now(),
            "manifest_sha256": sha256(path),
            "indexes": selected,
            "prepare_only": prepare_only,
            "retry": retry,
            "state": "submitting",
            "script_sha256": sha256(directory / "job.sh"),
        }
        atomic_json(directory / "submission.json", receipt)
        # Do not hold the submission lock for a potentially long --wait.
        command = ["sbatch", "--parsable", str(directory / "job.sh")]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=60)
            receipt.update(stdout=result.stdout, stderr=result.stderr, returncode=result.returncode)
            if result.returncode:
                receipt["state"] = "submission_failed"
                raise InputError("sbatch rejected submission; see submission.json")
            if not re.fullmatch(r"[0-9]+(?:;[^\s]+)?\s*", result.stdout):
                raise InputError("Ambiguous sbatch response; inspect scheduler before resubmitting")
            receipt.update(state="submitted", job_id=result.stdout.strip().split(";")[0])
        except Exception as exc:
            # Timeout/connection interruption may occur after SLURM accepted the job.
            # Unknown dispatch is protected against duplicate retries until reconciled.
            receipt.update(
                state=receipt["state"]
                if receipt["state"] == "submission_failed"
                else "submission_unknown",
                error=str(exc),
            )
            raise
        finally:
            atomic_json(directory / "submission.json", receipt)
    if wait:
        import time

        while True:
            state = batch_status(path)
            pending = {
                r["index"]
                for r in state["entries"]
                if r["state"]
                in {"running", "pending", "submitted_unknown", "interrupted_or_running_unknown"}
            }
            if not pending.intersection(selected):
                break
            time.sleep(30)
        receipt["completed_counts"] = state["counts"]
    return receipt
