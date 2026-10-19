"""Fail-closed loading of complete paired HMM exports and exact proteome identity."""

import hashlib
import json
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import DomainHit
from sprb_scanner.io.common import sha256


def checked_path(root, name):
    path = root / name
    if Path(name).is_absolute() or not path.resolve().is_relative_to(root.resolve()):
        raise InputError(f"Scan manifest path escapes input directory: {name}")
    return path


def read_manifest(root, fingerprints):
    metadata_path = root / "run_metadata.json"
    before = sha256(metadata_path)
    meta = json.loads(metadata_path.read_text())
    key = str(metadata_path.resolve())
    if key in fingerprints and fingerprints[key] != before:
        raise InputError(f"Input changed during loading: {metadata_path}")
    fingerprints[key] = before
    checks = meta.get("output_sha256", {})
    if not checks:
        raise InputError(f"Missing scan output manifest: {root}")
    for name, checksum in checks.items():
        path = checked_path(root, name)
        if not path.is_file() or sha256(path) != checksum:
            raise InputError(f"Scan checksum mismatch: {path}")
        key = str(path.resolve())
        if key in fingerprints and fingerprints[key] != checksum:
            raise InputError(f"Conflicting input checksum: {path}")
        fingerprints[key] = checksum
    return meta


def require_manifest(meta, names):
    if not set(names) <= set(meta["output_sha256"]):
        raise InputError("Incomplete scan manifest for association inputs")


def load_pair(root, dataset, model_ids, fingerprints):
    """Require both successful models; failed searches cannot imply no_local_hit."""
    root = Path(root)
    meta = read_manifest(root, fingerprints)
    require_manifest(
        meta,
        [
            "summary.json",
            "hits.json",
            *[
                f"{m}/{n}"
                for m in model_ids
                for n in ("run_metadata.json", "hits.json", "status.json", "raw_domains.json")
            ],
        ],
    )
    report = json.loads((root / "summary.json").read_text())
    if report.get("state") != "succeeded" or set(report.get("models", {})) != set(model_ids):
        raise InputError("Association requires both models to have succeeded")
    if any(report["models"][m].get("state") != "succeeded" for m in model_ids):
        raise InputError("Association requires both models to have succeeded")
    proteins = {p.protein_id: p for p in dataset.proteins}
    merged = []
    for model in model_ids:
        directory = root / model
        native = read_manifest(directory, fingerprints)
        require_manifest(native, ["hits.json", "status.json", "raw_domains.json"])
        status = json.loads((directory / "status.json").read_text())
        if status.get("state") != "succeeded" or status.get("model_id") != model:
            raise InputError(f"Unsuccessful/mismatched model: {model}")
        ids = list(native.get("protein_ids", {}).values())
        if len(ids) != len(set(ids)) or set(ids) != set(proteins):
            raise InputError(f"Scan protein IDs do not match prepared genome: {model}")
        digest = hashlib.sha256(
            json.dumps(
                [(pid, proteins[pid].sequence) for pid in ids],
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        if native.get("input_proteins_sha256") != digest:
            raise InputError(f"Scan protein sequence SHA-256 mismatch or unavailable: {model}")
        if native.get("target_count") != len(proteins):
            raise InputError(f"Scan target count mismatch: {model}")
        provenance = meta["profiles"][model]["provenance"]
        if provenance != native["profile"]["provenance"] or provenance["model_id"] != model:
            raise InputError(f"Conflicting model provenance: {model}")
        for data in json.loads((directory / "hits.json").read_text()):
            data.pop("ctd_distance", None)
            hit = DomainHit.model_validate(data)
            if hit.model.model_dump() != provenance:
                raise InputError(f"Hit model provenance mismatch: {hit.hit_id}")
            if hit.raw_artifact != "raw_domains.json":
                raise InputError(f"Unsupported raw artifact reference: {hit.hit_id}")
            merged.append(hit.model_copy(update={"raw_artifact": f"{model}/raw_domains.json"}))
    expected = [h.model_dump(mode="json") for h in merged]
    if json.loads((root / "hits.json").read_text()) != expected:
        raise InputError("Paired and individual model hits disagree")
    return tuple(merged), meta


def check_unchanged(fingerprints):
    for path, checksum in fingerprints.items():
        if sha256(Path(path)) != checksum:
            raise InputError(f"Association input changed during loading: {path}")
