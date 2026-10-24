"""Immutable inventory; annotation validation happens per compute task."""

import importlib
import json
import re
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256
from sprb_scanner.pipeline import load_pipeline_config
from sprb_scanner.pipeline.identity import STAGE_MODELS, implementation
from sprb_scanner.pipeline.storage import atomic_json, now
from sprb_scanner.references import ReferenceRegistry

from .models import BatchEntry, BatchManifest, InputFile, SlurmResources


def read_json(path):
    return json.loads(Path(path).read_text())


def software_identity():
    modules = (
        "io.prepare",
        "io.gff",
        "core.locations",
        "io.fasta",
        "pipeline.runner",
        "pipeline.config",
        "pipeline.identity",
        "pipeline.storage",
        "batch.task",
        "batch.plan",
        "batch.models",
    )
    return {
        "scanners": {stage: implementation(stage) for stage in (*STAGE_MODELS, "evidence")},
        "orchestration": {
            name: sha256(Path(importlib.import_module("sprb_scanner." + name).__file__))
            for name in modules
        },
    }


def verify_pins(manifest):
    if software_identity() != manifest.software:
        raise InputError("Batch implementation changed; create a new plan")
    for name, expected in manifest.pinned_files.items():
        if sha256(Path(name)) != expected:
            raise InputError(f"Pinned configuration/reference changed: {name}; create a new plan")


def load_manifest(path, *, pins=True):
    path = Path(path).resolve()
    expected = path.with_suffix(".sha256").read_text().strip()
    if sha256(path) != expected:
        raise InputError("Batch manifest checksum mismatch")
    manifest = BatchManifest.model_validate_json(path.read_text())
    if path != Path(manifest.output_root) / "manifest.json":
        raise InputError("Manifest moved; immutable absolute-path plans cannot be relocated")
    if pins:
        verify_pins(manifest)
    return manifest


def verify_inputs(entry):
    if entry.inventory_state != "ready":
        raise InputError(entry.error)
    for record in entry.files.values():
        path = Path(record.path)
        if path.stat().st_size != record.size_bytes or sha256(path) != record.sha256:
            raise InputError(f"Assembly input changed: {path}; create a new plan")


def plan_batch(dataset, config_path, resources_path, output, *, accessions=None):
    """Hash raw files and snapshot effective config/resources; no parsing/searching."""
    dataset, output = Path(dataset).resolve(), Path(output).resolve()
    config = load_pipeline_config(config_path)
    resources = SlurmResources.model_validate_json(Path(resources_path).read_text())
    if max(config.hmm.threads, config.mmseqs.threads) > resources.cpus:
        raise InputError("Scanner threads exceed requested SLURM CPUs")
    registry = ReferenceRegistry.load(config.registry)
    pins = {str(registry.path): registry.registry_sha256}
    for stage, models in STAGE_MODELS.items():
        if not getattr(config.stages, stage):
            continue
        for model in models:
            binding = config.references.get(model)
            resource = registry.resolve(
                model, binding.version if binding else None, binding.artifact if binding else None
            )
            pins[str(resource.artifact_path)] = resource.reference.sha256
    folders = sorted(
        p
        for p in dataset.iterdir()
        if p.is_dir() and re.fullmatch(r"GC[AF]_[0-9]+\.[0-9]+", p.name)
    )
    if accessions is not None:
        selected = set(accessions)
        if selected - {p.name for p in folders}:
            raise InputError("Requested assembly not present in dataset")
        folders = [p for p in folders if p.name in selected]
    reports = {}
    report = dataset / "assembly_data_report.jsonl"
    if report.exists():
        before = sha256(report)
        for line in report.read_text().splitlines():
            row = json.loads(line)
            key = row["accession"]
            if key in reports:
                raise InputError(f"Duplicate assembly report identity: {key}")
            reports[key] = row
        if before != sha256(report):
            raise InputError("Assembly report changed while planning")
    entries = []
    for i, folder in enumerate(folders):
        files, errors = {}, []
        for key in ("gff", "fna", "faa"):
            matches = sorted(folder.glob("*." + key))
            if len(matches) != 1:
                errors.append(f"Expected exactly one .{key} file; found {len(matches)}")
                continue
            path = matches[0].resolve()
            stat = path.stat()
            checksum = sha256(path)
            if path.stat().st_size != stat.st_size or path.stat().st_mtime_ns != stat.st_mtime_ns:
                raise InputError(f"Input changed during inventory: {path}")
            if not stat.st_size:
                errors.append(f"Empty {key} input")
            files[key] = InputFile(path=str(path), size_bytes=stat.st_size, sha256=checksum)
        metadata = reports.get(folder.name, {})
        entries.append(
            BatchEntry(
                index=i,
                assembly_accession=folder.name,
                genome_id=folder.name,
                files=files,
                inventory_state="invalid" if errors else "ready",
                error="; ".join(errors) or None,
                annotation_identity={
                    k: metadata.get(k)
                    for k in ("accession", "organism", "assemblyInfo", "annotationInfo")
                },
                prepared_path=str(output / "genomes" / folder.name / "prepared"),
            )
        )
    if output.exists():
        raise InputError("Batch output already exists; choose a new immutable plan directory")
    snapshot = output / "pipeline.yaml"
    manifest = BatchManifest(
        created_at=now(),
        output_root=str(output),
        dataset_root=str(dataset),
        software=software_identity(),
        pipeline_config=str(snapshot),
        pinned_files=pins,
        resources=resources,
        entries=entries,
    )
    output.mkdir(parents=True)
    import yaml

    snapshot.write_text(yaml.safe_dump(config.model_dump(mode="json"), sort_keys=True))
    manifest.pinned_files[str(snapshot)] = sha256(snapshot)
    atomic_json(output / "manifest.json", manifest.model_dump(mode="json"))
    (output / "manifest.sha256").write_text(sha256(output / "manifest.json") + "\n")
    (output / "logs").mkdir()
    return {
        "manifest": str(output / "manifest.json"),
        "assemblies": len(entries),
        "inventory_invalid": sum(e.inventory_state == "invalid" for e in entries),
    }
