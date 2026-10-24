"""Local, locked, checkpointed orchestration of existing scientific service APIs."""

import json
import traceback
from pathlib import Path
from urllib.parse import quote

from sprb_scanner import __version__
from sprb_scanner.association.inputs import check_unchanged, read_manifest
from sprb_scanner.association.models import DomainEvidence
from sprb_scanner.core.exceptions import InputError, ScannerError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.evidence.inputs import load_hmm, load_pair, load_rn28
from sprb_scanner.evidence.models import EvidenceResult, SourceRun
from sprb_scanner.evidence.service import run_merge
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared
from sprb_scanner.pipeline.config import (
    PipelineConfig,
    ReferenceBinding,
    load_pipeline_config,
    resolve_paths,
)
from sprb_scanner.pipeline.identity import STAGE_MODELS, StageSpec, implementation, scanner_spec
from sprb_scanner.pipeline.storage import (
    PROTOCOL,
    atomic_json,
    digest,
    inventory,
    next_directory,
    now,
    verify_files,
    workspace_lock,
)
from sprb_scanner.references import ReferenceRegistry
from sprb_scanner.scanners.bctd import run_bctd
from sprb_scanner.scanners.mmseqs import run_search
from sprb_scanner.scanners.repeats import run_repeats
from sprb_scanner.scanners.sprf import run_sprf


def _read(path):
    return json.loads(path.read_text())


def _validate_result(spec, result, dataset, config):
    fingerprints = {}
    stage = spec.stage_id
    if stage in {"typeb_ctd", "porp_sprf"}:
        runs, _ = load_pair(result, STAGE_MODELS[stage], dataset, fingerprints)
        for key, run in runs.items():
            if run.metadata["configuration"] != spec.payload["parameters"]:
                raise InputError("Stage result used different HMM parameters")
            expected = spec.resources[key].as_profile().provenance.model_dump()
            if run.metadata["profile"]["provenance"] != expected:
                raise InputError("Stage result used a different reference")
            if (
                run.state == "succeeded"
                and run.metadata["backend"] != spec.payload["backend"]["version"]
            ):
                raise InputError("Stage result used a different HMM backend")
        successes = sum(r.state == "succeeded" for r in runs.values())
        state = "succeeded" if successes == 2 else "partial" if successes else "failed"
    elif stage == "sprb_repeat":
        run, _ = load_hmm(result, "PF13573", dataset, fingerprints)
        if run.metadata["configuration"] != spec.payload["parameters"]:
            raise InputError("Repeat stage used different parameters")
        if (
            run.metadata["profile"]["provenance"]
            != spec.resources["PF13573"].as_profile().provenance.model_dump()
        ):
            raise InputError("Repeat stage used a different reference")
        if (
            run.state == "succeeded"
            and run.metadata["backend"] != spec.payload["backend"]["version"]
        ):
            raise InputError("Repeat stage used a different backend")
        state = run.state
    elif stage == "rn28":
        binding = config.references.get("RN28", ReferenceBinding())
        run, _ = load_rn28(
            result, dataset, config.registry, fingerprints, binding.version, binding.artifact
        )
        if run.metadata["configuration"] != spec.payload["parameters"]:
            raise InputError("RN28 stage used different parameters")
        if (
            run.state == "succeeded"
            and run.metadata["backend_version"] != spec.payload["backend"]["version"]
        ):
            raise InputError("RN28 stage used a different backend")
        state = run.state
    else:
        read_manifest(result, fingerprints)
        data = _read(result / "evidence_profiles.json")
        # Existing exports include computed fields; remove only known computed properties
        # at their typed evidence locations, never from gene qualifiers or raw metadata.
        for profile in data["profiles"]:
            for evidence in profile["sources"].values():
                for domain in evidence["domains"]:
                    for key in DomainEvidence.model_computed_fields:
                        domain.pop(key, None)
                    domain["hit"].pop("ctd_distance", None)
                for alignment in evidence["alignments"]:
                    alignment.pop("ctd_distance", None)
        validated = EvidenceResult.model_validate(data)
        if validated.genome_id != dataset.genome.genome_id:
            raise InputError("Evidence result belongs to another genome")
        status = _read(result / "status.json")
        if status["state"] != "succeeded":
            raise InputError("Incomplete evidence merger output")
        state = (
            "partial"
            if config.neighborhood.enabled and validated.neighborhood_state != "succeeded"
            else "succeeded"
        )
    check_unchanged(fingerprints)
    return state


def _execute(spec, result, dataset, prepared, config, dependencies):
    stage = spec.stage_id
    if stage in {"typeb_ctd", "porp_sprf"}:
        service = run_bctd if stage == "typeb_ctd" else run_sprf
        profiles = tuple(spec.resources[m].as_profile() for m in STAGE_MODELS[stage])
        return service(dataset.proteins, profiles, result, config.hmm)
    if stage == "sprb_repeat":
        return run_repeats(
            dataset.proteins, spec.resources["PF13573"].as_profile(), result, config.hmm
        )
    if stage == "rn28":
        return run_search(
            spec.resources["RN28"].artifact_path, result, prepared=prepared, config=config.mmseqs
        )
    paths, source_states = {}, {}
    for dependency, models in STAGE_MODELS.items():
        item = dependencies[dependency]
        if item.get("result_usable"):
            paths[dependency] = Path(item["result_directory"])
        elif item["state"] != "not_run":
            for model in models:
                source_states[model] = SourceRun(
                    source_id=model,
                    state="failed",
                    message=item.get("error") or f"Pipeline stage {dependency} did not complete",
                    metadata={
                        "pipeline_stage": dependency,
                        "fingerprint": item.get("fingerprint"),
                        "execution_manifest": item.get("execution_manifest"),
                    },
                )
    binding = config.references.get("RN28", ReferenceBinding())
    return run_merge(
        prepared,
        result,
        bctd_scan=paths.get("typeb_ctd"),
        sprf_scan=paths.get("porp_sprf"),
        repeat_scan=paths.get("sprb_repeat"),
        rn28_scan=paths.get("rn28"),
        registry=config.registry,
        window=config.neighborhood.window,
        neighborhood=config.neighborhood.enabled,
        source_states=source_states,
        rn28_version=binding.version,
        rn28_artifact=binding.artifact,
    )


def _attempt_reference(directory, manifest, action, rejections=()):
    return {
        "stage_id": manifest["stage_id"],
        "state": manifest["state"],
        "action": action,
        "fingerprint": manifest["fingerprint"],
        "attempt_directory": str(directory),
        "execution_manifest": str(directory / "execution_manifest.json"),
        "execution_manifest_sha256": sha256(directory / "execution_manifest.json"),
        "result_directory": str(directory / "result") if manifest["result_usable"] else None,
        "result_usable": manifest["result_usable"],
        "error": manifest.get("error"),
        "cache_rejections": list(rejections),
    }


def _run_stage(spec, stages_root, dataset, prepared, config, dependencies):
    parent = stages_root / spec.stage_id / spec.fingerprint
    rejections = []
    check_unchanged(spec.inputs)
    if parent.exists():
        if not parent.resolve().is_relative_to(stages_root.resolve()):
            raise InputError("Attempt cache escapes pipeline workspace")
        for previous in sorted(parent.glob("attempt-*"), reverse=True):
            try:
                if previous.is_symlink() or not previous.resolve().is_relative_to(
                    stages_root.resolve()
                ):
                    raise InputError("Attempt path escapes pipeline workspace")
                manifest = _read(previous / "execution_manifest.json")
                if (
                    manifest["schema_version"] != PROTOCOL
                    or manifest["state"] != "succeeded"
                    or manifest["stage_id"] != spec.stage_id
                    or manifest["fingerprint"] != spec.fingerprint
                    or manifest["fingerprint_payload"] != spec.payload
                    or digest(manifest["fingerprint_payload"]) != spec.fingerprint
                    or manifest.get("dependency_bindings", {}) != spec.bindings
                ):
                    raise InputError("Stage identity/dependencies changed or attempt incomplete")
                verify_files(previous, manifest["output_sha256"])
                if _validate_result(spec, previous / "result", dataset, config) != "succeeded":
                    raise InputError("Previous stage is not fully successful")
                check_unchanged(spec.inputs)
                return _attempt_reference(previous, manifest, "reused", rejections)
            except (OSError, ValueError, KeyError, TypeError, AttributeError, ScannerError) as exc:
                rejections.append({"attempt": str(previous), "reason": str(exc)})
    directory = next_directory(parent, "attempt-")
    manifest = {
        "schema_version": PROTOCOL,
        "stage_id": spec.stage_id,
        "fingerprint": spec.fingerprint,
        "fingerprint_payload": spec.payload,
        "dependency_bindings": spec.bindings,
        "input_files_sha256": spec.inputs,
        "state": "running",
        "started_at": now(),
        "result_usable": False,
        "cache_rejections": rejections,
    }
    path = directory / "execution_manifest.json"
    atomic_json(path, manifest)
    interrupted = None
    with (directory / "stage.log").open("w", encoding="utf-8") as log:
        try:
            print(f"Starting {spec.stage_id}: {spec.fingerprint}", file=log, flush=True)
            if spec.errors:
                raise InputError("; ".join(spec.errors))
            _execute(spec, directory / "result", dataset, prepared, config, dependencies)
            state = _validate_result(spec, directory / "result", dataset, config)
            check_unchanged(spec.inputs)
            manifest.update(state=state, result_usable=True)
            if state != "succeeded":
                manifest["error"] = "Stage produced incomplete evidence; inspect native status"
        except Exception as exc:
            manifest.update(state="failed", error=str(exc), exception_type=type(exc).__name__)
            traceback.print_exc(file=log)
            # A native failed/partial output is usable only if its original loader
            # verifies it. Otherwise the merger receives an explicit failed source.
            if (directory / "result").is_dir() and not spec.errors:
                try:
                    _validate_result(spec, directory / "result", dataset, config)
                    check_unchanged(spec.inputs)
                    manifest["result_usable"] = True
                except Exception as validation_error:
                    manifest["result_validation_error"] = str(validation_error)
        except BaseException as exc:
            manifest.update(state="interrupted", error=str(exc), exception_type=type(exc).__name__)
            traceback.print_exc(file=log)
            interrupted = exc
        print(f"Completed {spec.stage_id}: {manifest['state']}", file=log, flush=True)
    manifest.update(
        completed_at=now(), output_sha256=inventory(directory, ("execution_manifest.json",))
    )
    atomic_json(path, manifest)
    if interrupted is not None:
        raise interrupted
    return _attempt_reference(directory, manifest, "executed", rejections)


def _merge_spec(config, registry, prepared_identity, prepared_files, dependencies):
    dependencies = {key: dependencies[key] for key in STAGE_MODELS}
    payload = {
        "protocol": PROTOCOL,
        "stage": "evidence",
        "prepared": prepared_identity,
        "implementation": implementation("evidence"),
        "parameters": config.neighborhood.model_dump(mode="json"),
        "dependencies": {
            k: {"fingerprint": d.get("fingerprint"), "state": d["state"]}
            for k, d in dependencies.items()
        },
    }
    inputs = dict(prepared_files)
    if config.stages.rn28:
        payload["registry_sha256"] = registry.registry_sha256
        inputs[str(registry.path)] = registry.registry_sha256
    bindings = {k: d.get("execution_manifest_sha256") for k, d in dependencies.items()}
    # Native dependency artifacts are inputs too: detect mutation between scanner
    # validation and merger execution/reuse, not just at the next invocation.
    for d in dependencies.values():
        if d.get("execution_manifest"):
            mp = Path(d["execution_manifest"])
            inputs[str(mp)] = d["execution_manifest_sha256"]
            native = _read(mp)
            for name, checksum in native["output_sha256"].items():
                inputs[str(mp.parent / name)] = checksum
    return StageSpec("evidence", payload, inputs=inputs, bindings=bindings)


def run_pipeline(prepared, config, output, *, resume=False):
    """Run enabled stages once; resume verifies caches and retries in new attempts."""
    if type(resume) is not bool:
        raise InputError("resume must be a boolean")
    config_path = None if isinstance(config, PipelineConfig) else Path(config).resolve()
    config_bytes = config_path.read_bytes() if config_path else None
    config_checksum = sha256(config_path) if config_path else None
    config = (
        resolve_paths(config, Path.cwd())
        if isinstance(config, PipelineConfig)
        else load_pipeline_config(config_path)
    )
    if config_path and (
        config_path.read_bytes() != config_bytes or sha256(config_path) != config_checksum
    ):
        raise InputError("Pipeline YAML changed during configuration loading")
    prepared, output = Path(prepared).resolve(), Path(output).resolve()
    files = {}
    read_manifest(prepared, files)
    if validate_prepared(prepared)["errors"]:
        raise InputError("Prepared dataset has quality errors")
    dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
    if not dataset.proteins:
        raise InputError("Pipeline requires a nonempty prepared proteome")
    if config.neighborhood.enabled and (not dataset.genes or not dataset.contigs):
        raise InputError("Neighborhood enabled but prepared data has no CDS coordinates")
    registry = ReferenceRegistry.load(config.registry)
    if output.exists() and not resume:
        raise InputError("Pipeline output already exists; use --resume or a new directory")
    output.mkdir(parents=True, exist_ok=True)
    with workspace_lock(output):
        workspace = output / "workspace.json"
        if workspace.exists():
            meta = _read(workspace)
            if (
                meta.get("schema_version") != PROTOCOL
                or meta.get("genome_id") != dataset.genome.genome_id
            ):
                raise InputError("Incompatible pipeline workspace or different genome ID")
        elif any(p.name != ".pipeline.lock" for p in output.iterdir()):
            raise InputError("Output is not a pipeline workspace")
        else:
            atomic_json(
                workspace,
                {
                    "schema_version": PROTOCOL,
                    "genome_id": dataset.genome.genome_id,
                    "created_at": now(),
                },
            )
        genome_key = quote(dataset.genome.genome_id, safe="._-")
        if genome_key in {".", ".."}:
            genome_key = "%2E" * len(genome_key)
        genome_root = output / "genomes" / genome_key
        if not genome_root.resolve().is_relative_to(output):
            raise InputError("Genome directory escapes pipeline workspace")
        for folder in (genome_root / "runs", genome_root / "stages"):
            if not folder.resolve().is_relative_to(output):
                raise InputError("Pipeline state directory escapes workspace")
        invocation = next_directory(genome_root / "runs", "run-")
        if config_bytes is not None:
            (invocation / "configuration.yaml").write_bytes(config_bytes)
        atomic_json(invocation / "effective_config.json", config.model_dump(mode="json"))
        journal = {
            "schema_version": PROTOCOL,
            "genome_id": dataset.genome.genome_id,
            "mode": config.mode,
            "state": "running",
            "started_at": now(),
            "configuration": config.model_dump(mode="json"),
            "configuration_file": str(config_path) if config_path else None,
            "configuration_sha256": config_checksum,
            "configuration_snapshot": str(invocation / "configuration.yaml")
            if config_bytes
            else None,
            "effective_config_sha256": sha256(invocation / "effective_config.json"),
            "prepared_files_sha256": files,
            "registry_sha256": registry.registry_sha256,
            "software_version": __version__,
            "stages": {
                stage: {"stage_id": stage, "state": "pending", "action": "not_run"}
                for stage in (*STAGE_MODELS, "evidence")
            },
        }
        prepared_identity = {Path(k).name: v for k, v in files.items()}
        journal_path = invocation / "run_manifest.json"
        atomic_json(journal_path, journal)
        try:
            blocked = False
            for stage in (*STAGE_MODELS, "evidence"):
                if blocked:
                    item = {
                        "stage_id": stage,
                        "state": "blocked",
                        "action": "not_run",
                        "error": "Strict mode stopped after an incomplete stage",
                    }
                elif stage != "evidence" and not getattr(config.stages, stage):
                    item = {"stage_id": stage, "state": "not_run", "action": "disabled"}
                else:
                    spec = (
                        _merge_spec(config, registry, prepared_identity, files, journal["stages"])
                        if stage == "evidence"
                        else scanner_spec(stage, config, registry, prepared_identity, files)
                    )
                    journal["active_stage"] = stage
                    atomic_json(journal_path, journal)
                    item = _run_stage(
                        spec, genome_root / "stages", dataset, prepared, config, journal["stages"]
                    )
                    if item["state"] != "succeeded" and config.mode == "strict":
                        blocked = True
                journal["stages"][stage] = item
                atomic_json(journal_path, journal)
            incomplete = any(
                d["state"] not in {"not_run", "succeeded"} for d in journal["stages"].values()
            )
            journal["state"] = "failed" if blocked else "partial" if incomplete else "succeeded"
        except BaseException as exc:
            active = journal.get("active_stage")
            if active and journal["stages"][active]["state"] == "pending":
                journal["stages"][active].update(
                    state="interrupted" if not isinstance(exc, Exception) else "failed",
                    error=str(exc),
                    exception_type=type(exc).__name__,
                )
            journal.update(
                state="interrupted" if not isinstance(exc, Exception) else "failed",
                error=str(exc),
                exception_type=type(exc).__name__,
            )
            raise
        finally:
            journal.pop("active_stage", None)
            journal["completed_at"] = now()
            atomic_json(journal_path, journal)
            atomic_json(
                output / "latest.json",
                {
                    "run_manifest": str(journal_path),
                    "run_manifest_sha256": sha256(journal_path),
                    "state": journal["state"],
                },
            )
    return journal
