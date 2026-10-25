"""One isolated assembly: strict preparation then the unmodified Batch 6B runner."""

import os
import socket
import traceback
from pathlib import Path

from sprb_scanner.association.inputs import check_unchanged, read_manifest
from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared
from sprb_scanner.io.prepare import prepare
from sprb_scanner.pipeline import run_pipeline
from sprb_scanner.pipeline.storage import atomic_json, next_directory, now, workspace_lock

from .plan import load_manifest, verify_inputs, verify_pins


def disk_usage(root):
    """Current apparent/allocated bytes; no cleanup or archival, no hard-link double count."""
    total = {"apparent_bytes": 0, "allocated_bytes": 0, "files": 0}
    seen = set()
    for p in Path(root).rglob("*"):
        if p.is_symlink():
            raise InputError(f"Unexpected symlink in genome workspace: {p}")
        if p.is_file():
            s = p.stat()
            key = (s.st_dev, s.st_ino)
            if key not in seen:
                seen.add(key)
                total["files"] += 1
                total["apparent_bytes"] += s.st_size
                total["allocated_bytes"] += s.st_blocks * 512
    return total


def validate_preparation(entry):
    prepared = Path(entry.prepared_path)
    fingerprints = {}
    metadata = read_manifest(prepared, fingerprints)
    actual = {r["path"]: r["sha256"] for r in metadata["inputs"]}
    expected = {r.path: r.sha256 for r in entry.files.values()}
    if actual != expected:
        raise InputError("Prepared inputs differ from batch inventory")
    report = validate_prepared(prepared)
    if report["genome_id"] != entry.genome_id or report["errors"] or not report["proteins"]:
        raise InputError("Prepared identity, quality or protein count is invalid")
    check_unchanged(fingerprints)
    return report


def run_task(manifest_path, index, *, prepare_only=False, local=False):
    if type(index) is not int or index < 0:
        raise InputError("Task index must be a nonnegative integer")
    if not local and not os.environ.get("SLURM_JOB_ID"):
        raise InputError("batch-task requires a SLURM allocation; use --local only for local tests")
    manifest = load_manifest(manifest_path, pins=False)
    if index >= len(manifest.entries):
        raise InputError("Task index outside batch manifest")
    entry = manifest.entries[index]
    root = Path(manifest.output_root) / "genomes" / entry.genome_id
    root.mkdir(parents=True, exist_ok=True)
    if not root.resolve().is_relative_to(Path(manifest.output_root) / "genomes"):
        raise InputError("Genome workspace escapes batch directory")
    with workspace_lock(root):
        attempt = next_directory(root / "tasks", "attempt-")
        report = {
            "schema_version": "1.0.0",
            "index": index,
            "genome_id": entry.genome_id,
            "manifest_sha256": sha256(Path(manifest_path)),
            "state": "running",
            "preparation_state": "not_run",
            "pipeline_state": "not_run",
            "prepare_only": prepare_only,
            "started_at": now(),
            "host": socket.gethostname(),
            "slurm": {
                k: os.environ.get(k)
                for k in (
                    "SLURM_JOB_ID",
                    "SLURM_ARRAY_JOB_ID",
                    "SLURM_ARRAY_TASK_ID",
                    "SLURM_CPUS_PER_TASK",
                )
            },
        }
        path = attempt / "task.json"
        atomic_json(path, report)
        atomic_json(root / "latest_task.json", {"task": str(path), "sha256": sha256(path)})
        try:
            verify_pins(manifest)
            verify_inputs(entry)
            if (
                not local
                and int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) < manifest.resources.cpus
            ):
                raise InputError("Allocated CPUs are below the manifest resource request")
            expected_index = os.environ.get("SLURM_ARRAY_TASK_ID")
            if not local and expected_index is not None and int(expected_index) != index:
                raise InputError("SLURM array index differs from requested task")
            temporary = attempt / "tmp"
            temporary.mkdir()
            # Every child sees project-local temporary storage. Pipeline's own scratch is local too.
            previous_tmp = os.environ.get("TMPDIR")
            os.environ["TMPDIR"] = str(temporary)
            try:
                prepared = Path(entry.prepared_path)
                if not prepared.exists():
                    prepare(
                        Path(entry.files["gff"].path),
                        "gff3",
                        prepared,
                        AppConfig(),
                        entry.genome_id,
                        Path(entry.files["fna"].path),
                        Path(entry.files["faa"].path),
                    )
                    report["preparation_action"] = "executed"
                else:
                    report["preparation_action"] = "reused"
                report["preparation"] = validate_preparation(entry)
                report["preparation_state"] = "succeeded"
                report["prepared_metadata_sha256"] = sha256(prepared / "run_metadata.json")
                atomic_json(path, report)
                if prepare_only:
                    report["state"] = "prepared"
                else:
                    workspace = root / "pipeline"
                    result = run_pipeline(
                        prepared, manifest.pipeline_config, workspace, resume=workspace.exists()
                    )
                    report.update(
                        state=result["state"],
                        pipeline_state=result["state"],
                        pipeline_latest=str(workspace / "latest.json"),
                        pipeline_latest_sha256=sha256(workspace / "latest.json"),
                    )
            finally:
                if previous_tmp is None:
                    os.environ.pop("TMPDIR", None)
                else:
                    os.environ["TMPDIR"] = previous_tmp
            verify_inputs(entry)
            verify_pins(manifest)
        except BaseException as exc:
            report.update(
                state="interrupted" if not isinstance(exc, Exception) else "failed",
                error=str(exc),
                exception_type=type(exc).__name__,
            )
            if report["preparation_state"] != "succeeded":
                report["preparation_state"] = "failed"
            (attempt / "error.log").write_text(traceback.format_exc())
            if not isinstance(exc, Exception):
                raise
        finally:
            report.update(completed_at=now(), disk_usage=disk_usage(root))
            atomic_json(path, report)
            atomic_json(root / "latest_task.json", {"task": str(path), "sha256": sha256(path)})
        return report
