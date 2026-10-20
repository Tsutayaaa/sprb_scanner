"""Relevant content/backend identities only; unrelated reference changes do not rescan HMMs."""

import importlib
import importlib.metadata
import platform
from dataclasses import dataclass, field
from pathlib import Path

import pyhmmer

from sprb_scanner import __version__
from sprb_scanner.core.exceptions import InputError, ScannerError
from sprb_scanner.io.common import sha256
from sprb_scanner.pipeline.config import ReferenceBinding
from sprb_scanner.pipeline.storage import PROTOCOL, digest
from sprb_scanner.scanners.mmseqs import Backend

STAGE_MODELS = {
    "typeb_ctd": ("TIGR04131", "PF13585"),
    "porp_sprf": ("TIGR03519", "PF11751"),
    "sprb_repeat": ("PF13573",),
    "rn28": ("RN28",),
}


@dataclass
class StageSpec:
    stage_id: str
    payload: dict
    inputs: dict = field(default_factory=dict)
    resources: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)
    bindings: dict = field(default_factory=dict)

    @property
    def fingerprint(self):
        return digest(self.payload)


def implementation(stage):
    modules = [
        "core.models",
        "io.common",
        "io.exporters",
        "references.registry",
        "references.models",
    ]
    modules += {
        "typeb_ctd": ["scanners.hmm", "scanners.pair", "scanners.bctd"],
        "porp_sprf": ["scanners.hmm", "scanners.pair", "scanners.sprf"],
        "sprb_repeat": ["scanners.hmm", "scanners.repeats"],
        "rn28": ["scanners.mmseqs", "io.fasta"],
        "evidence": [
            "evidence.inputs",
            "evidence.merger",
            "evidence.models",
            "evidence.service",
            "association.inputs",
            "association.analyzer",
            "association.models",
            "association.service",
            "neighborhood.models",
            "neighborhood.analyzer",
        ],
    }[stage]
    return {
        name: sha256(Path(importlib.import_module("sprb_scanner." + name).__file__))
        for name in modules
    }


def backend_identity(stage, config):
    common = {
        "python": platform.python_version(),
        "platform": platform.system(),
        "machine": platform.machine(),
        "pydantic": importlib.metadata.version("pydantic"),
        "biopython": importlib.metadata.version("biopython"),
    }
    if stage == "rn28":
        metadata = {}
        backend = Backend(config.mmseqs, metadata, Path("unused-probe.log"))
        return {
            **common,
            "version": backend.version,
            "executable": backend.executable,
            "executable_sha256": sha256(Path(backend.executable)),
            "capabilities": metadata["capabilities"],
        }
    return {
        **common,
        "version": "pyhmmer " + pyhmmer.__version__,
        "library_sha256": {
            str(Path(module.__file__).resolve()): sha256(Path(module.__file__))
            for module in (pyhmmer.easel, pyhmmer.plan7)
        },
    }


def scanner_spec(stage, config, registry, prepared_identity, prepared_files):
    parameters = config.mmseqs if stage == "rn28" else config.hmm
    spec = StageSpec(
        stage,
        {
            "protocol": PROTOCOL,
            "software_version": __version__,
            "stage": stage,
            "prepared": prepared_identity,
            "parameters": parameters.model_dump(mode="json"),
            "implementation": implementation(stage),
            "references": {},
        },
        inputs=dict(prepared_files),
    )
    # Registry bytes are checked for concurrent mutation but excluded from scan fingerprints:
    # only the selected entries and their pinned lineage belong to each scan's identity.
    spec.inputs[str(registry.path)] = registry.registry_sha256
    for model in STAGE_MODELS[stage]:
        binding = config.references.get(model, ReferenceBinding())
        info = {"binding": binding.model_dump(mode="json")}
        try:
            entry = registry.get(model, binding.version)
            info.update(
                reference=entry.model_dump(mode="json"),
                lineage=[e.model_dump(mode="json") for e in registry.lineage(*entry.key)],
            )
            resolved = registry.resolve(model, binding.version, binding.artifact)
            expected_type = "sequence" if model == "RN28" else "profile_hmm"
            if entry.resource_type != expected_type:
                raise InputError(f"{model} requires {expected_type}")
            if model != "RN28" and resolved.as_profile().provenance.model_id != model:
                raise InputError(f"Registry model does not match pipeline role: {model}")
            # Fingerprint the resolved identity, not redundant selection syntax:
            # explicit/default bindings of the same resource are equivalent.
            info.pop("binding")
            spec.resources[model] = resolved
            spec.inputs[str(resolved.artifact_path)] = entry.sha256
            info.update(artifact_sha256=entry.sha256, artifact_path=str(resolved.artifact_path))
        except (ScannerError, OSError, ValueError) as exc:
            info["error"] = str(exc)
            spec.errors.append(f"{model}: {exc}")
        spec.payload["references"][model] = info
    try:
        spec.payload["backend"] = backend_identity(stage, config)
        backend = spec.payload["backend"]
        if "executable" in backend:
            spec.inputs[backend["executable"]] = backend["executable_sha256"]
        spec.inputs.update(backend.get("library_sha256", {}))
    except (ScannerError, OSError, ValueError) as exc:
        spec.payload["backend"] = {"error": str(exc)}
        spec.errors.append(str(exc))
    return spec
