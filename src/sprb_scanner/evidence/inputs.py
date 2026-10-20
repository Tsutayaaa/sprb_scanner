"""Verify existing exports offline; no executable, model download or search is invoked."""

import hashlib
import json
import math
from pathlib import Path

from sprb_scanner.association.inputs import read_manifest, require_manifest
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import DomainHit, ProteinRecord
from sprb_scanner.evidence.models import SourceRun
from sprb_scanner.io.common import sha256
from sprb_scanner.references import ReferenceRegistry
from sprb_scanner.scanners.mmseqs import SequenceAlignment, fingerprint


def payload(path):
    return json.loads(path.read_text())


def protein_map(path):
    records = payload(path)
    return {
        key: ProteinRecord.model_validate({k: v for k, v in data.items() if k != "length"})
        for key, data in records.items()
    }


def verify_targets(records, dataset):
    proteins = {p.protein_id: p for p in dataset.proteins}
    ids = [p.protein_id for p in records]
    if len(ids) != len(set(ids)) or set(ids) != set(proteins):
        raise InputError("Evidence protein IDs do not match prepared proteome")
    if any(p != proteins[p.protein_id] for p in records):
        raise InputError("Evidence protein sequence/identity does not match prepared proteome")


def load_hmm(root, model, dataset, fingerprints):
    root = Path(root)
    meta = read_manifest(root, fingerprints)
    require_manifest(meta, ["status.json"])
    status = payload(root / "status.json")
    provenance = meta["profile"]["provenance"]
    if status.get("model_id") != model or provenance.get("model_id") != model:
        raise InputError(f"Unexpected model identity: {model}")
    state = status.get("state")
    if state not in {"succeeded", "failed"}:
        raise InputError(f"Incomplete HMM run: {root}")
    if state == "failed":
        if any(name in meta["output_sha256"] for name in ("hits.json", "raw_domains.json")):
            raise InputError(f"Failed HMM has ambiguous partial hit artifacts: {root}")
        return SourceRun(
            source_id=model,
            state="failed",
            directory=str(root.resolve()),
            message=status.get("message") or "HMM search failed",
            metadata=meta,
        ), ()
    require_manifest(meta, ["hits.json", "raw_domains.json"])
    proteins = {p.protein_id: p for p in dataset.proteins}
    ids = list(meta["protein_ids"].values())
    if len(ids) != len(set(ids)) or set(ids) != set(proteins):
        raise InputError(f"HMM protein IDs mismatch: {model}")
    if meta["input_proteins_sha256"] != fingerprint([proteins[pid] for pid in ids]):
        raise InputError(f"HMM protein sequence fingerprint mismatch: {model}")
    if meta["target_count"] != len(proteins):
        raise InputError(f"HMM target count mismatch: {model}")
    hits = []
    for data in payload(root / "hits.json"):
        hit = DomainHit.model_validate({k: v for k, v in data.items() if k != "ctd_distance"})
        p = proteins.get(hit.protein_id)
        if p is None or p.length != hit.protein_length:
            raise InputError(f"Hit target identity mismatch: {hit.hit_id}")
        if hit.model.model_dump() != provenance or hit.method != "hmm" or hit.hmm is None:
            raise InputError(f"Hit HMM provenance mismatch: {hit.hit_id}")
        if hit.hmm.ga is None or hit.raw_artifact != "raw_domains.json":
            raise InputError(f"Missing GA/raw evidence: {hit.hit_id}")
        if not math.isclose(
            hit.query_coverage, (hit.query_end - hit.query_start) / hit.query_length
        ):
            raise InputError(f"Invalid HMM coverage: {hit.hit_id}")
        if not math.isclose(hit.target_coverage, (hit.end - hit.start) / hit.protein_length):
            raise InputError(f"Invalid protein coverage: {hit.hit_id}")
        hits.append(hit)
    if len({h.hit_id for h in hits}) != len(hits):
        raise InputError(f"Duplicate domain hit IDs: {model}")
    if meta.get("domain_count", len(hits)) != len(hits):
        raise InputError(f"HMM domain count mismatch: {model}")
    return SourceRun(
        source_id=model,
        state="succeeded",
        directory=str(root.resolve()),
        retained_hits=len(hits),
        metadata=meta,
    ), tuple(hits)


def load_pair(root, models, dataset, fingerprints):
    root = Path(root)
    meta = read_manifest(root, fingerprints)
    require_manifest(
        meta,
        [
            "summary.json",
            "hits.json",
            *(
                f"{model}/{name}"
                for model in models
                for name in ("status.json", "run_metadata.json")
            ),
        ],
    )
    summary = payload(root / "summary.json")
    if set(summary["models"]) != set(models):
        raise InputError("Unexpected paired HMM models")
    runs, hits = {}, []
    for model in models:
        run, group = load_hmm(root / model, model, dataset, fingerprints)
        if summary["models"][model]["state"] != run.state:
            raise InputError("Conflicting paired model status")
        if meta["profiles"][model]["provenance"] != run.metadata["profile"]["provenance"]:
            raise InputError("Conflicting paired model provenance")
        runs[model] = run
        hits.extend(group)
    successes = sum(r.state == "succeeded" for r in runs.values())
    if summary["state"] != (
        "succeeded" if successes == 2 else "partial" if successes else "failed"
    ):
        raise InputError("Paired aggregate status disagrees with individual models")
    expected = [
        h.model_copy(update={"raw_artifact": f"{h.model.model_id}/raw_domains.json"}).model_dump(
            mode="json"
        )
        for h in hits
    ]
    if payload(root / "hits.json") != expected:
        raise InputError("Paired and individual model evidence disagree")
    return runs, tuple(hits)


def load_rn28(root, dataset, registry_path, fingerprints, version=None, artifact=None):
    root = Path(root)
    meta = read_manifest(root, fingerprints)
    require_manifest(meta, ["status.json"])
    state = payload(root / "status.json").get("state")
    if meta.get("kind") != "mmseqs2-search" or meta.get("state") != state:
        raise InputError("Not a complete MMseqs2 search export")
    if state == "failed":
        if "alignments.json" in meta["output_sha256"]:
            raise InputError("Failed MMseqs2 run has ambiguous partial alignment artifacts")
        return SourceRun(
            source_id="RN28",
            state="failed",
            directory=str(root.resolve()),
            message=meta.get("error") or "MMseqs2 search failed",
            metadata=meta,
        ), ()
    if state != "succeeded":
        raise InputError("Incomplete MMseqs2 search")
    require_manifest(
        meta, ["alignments.json", "target_mapping.json", "query_mapping.json", "raw_alignments.tsv"]
    )
    targets = protein_map(root / "target_mapping.json")
    verify_targets(tuple(targets.values()), dataset)
    database = meta["target_database_metadata"]
    if database["proteins_sha256"] != fingerprint(tuple(targets.values())):
        raise InputError("MMseqs2 database protein fingerprint mismatch")
    if meta["target_count"] != len(targets):
        raise InputError("MMseqs2 target count mismatch")
    queries = protein_map(root / "query_mapping.json")
    if meta["query_proteins_sha256"] != fingerprint(tuple(queries.values())):
        raise InputError("MMseqs2 query fingerprint mismatch")
    registry = ReferenceRegistry.load(registry_path)
    fingerprints[str(registry.path)] = registry.registry_sha256
    resolved = registry.resolve("RN28", version, artifact)
    fingerprints[str(resolved.artifact_path)] = sha256(resolved.artifact_path)
    reference = registry.sequences("RN28", version, artifact)[0]
    query_by_id = {q.protein_id: q for q in queries.values()}
    if len(query_by_id) != len(queries) or not queries:
        raise InputError("MMseqs2 queries have duplicate IDs or are empty")
    rn28_ids = {q.protein_id for q in queries.values() if q.sequence == reference.sequence}
    if not rn28_ids:
        raise InputError("MMseqs2 query does not contain the registered RN28 sequence")
    all_hits, hits = [], []
    by_id = {p.protein_id: p for p in targets.values()}
    for data in payload(root / "alignments.json"):
        hit = SequenceAlignment.model_validate(
            {k: v for k, v in data.items() if k != "ctd_distance"}
        )
        p = by_id.get(hit.protein_id)
        q = query_by_id.get(hit.query_id)
        if p is None or q is None or hit.protein_length != p.length or hit.query_length != q.length:
            raise InputError(f"MMseqs2 alignment identity mismatch: {hit.hit_id}")
        if (
            hit.target_source_id != p.source_id
            or hit.target_gene_id != p.gene_id
            or hit.query_source_id != q.source_id
            or hit.model.sha256 != hashlib.sha256(q.sequence.encode()).hexdigest()
            or hit.model.version != "sequence-sha256:" + hit.model.sha256
            or hit.raw_artifact != "raw_alignments.tsv"
            or hit.model.model_id != q.protein_id
            or hit.backend_version != meta["backend_version"]
        ):
            raise InputError(f"MMseqs2 alignment provenance mismatch: {hit.hit_id}")
        all_hits.append(hit)
        if hit.query_id in rn28_ids:
            hits.append(hit)
    raw_count = sum(bool(line.strip()) for line in (root / "raw_alignments.tsv").open())
    if len(all_hits) != meta["alignment_count"] or len(all_hits) != raw_count:
        raise InputError("MMseqs2 raw/exported alignment count mismatch")
    if len({h.hit_id for h in all_hits}) != len(all_hits):
        raise InputError("Duplicate alignment hit IDs")
    return SourceRun(
        source_id="RN28",
        state="succeeded",
        directory=str(root.resolve()),
        retained_hits=len(hits),
        metadata={
            **meta,
            "rn28_reference": resolved.reference.model_dump(mode="json"),
            "rn28_query_ids": sorted(rn28_ids),
        },
    ), tuple(hits)
