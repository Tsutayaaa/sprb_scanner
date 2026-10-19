"""Independent PorP/SprF family evidence via two existing HMM models."""

import csv
import json
import shutil
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from sprb_scanner import __version__
from sprb_scanner.core.exceptions import HMMError, InputError
from sprb_scanner.core.models import DomainHit
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import write_json
from sprb_scanner.scanners.hmm import HMMSearchConfig, run_hmm

MODEL_IDS = ("TIGR03519", "PF11751")


def collect_evidence(proteins, hits, statuses):
    """Group hits by protein without merging scores, coordinates or model identity."""
    grouped = {}
    by_id = {p.protein_id: p for p in proteins}
    for hit in hits:
        grouped.setdefault(hit.protein_id, []).append(hit)
    rows = []
    for protein_id, evidence in sorted(grouped.items()):
        hit_models = sorted({h.model.model_id for h in evidence})
        ga_models = sorted({h.model.model_id for h in evidence if h.hmm.ga.status == "passed_ga"})
        p = by_id[protein_id]
        rows.append(
            {
                "protein_id": protein_id,
                "source_id": p.source_id,
                "gene_id": p.gene_id,
                "protein_length": p.length,
                "hit_models": hit_models,
                "ga_models": ga_models,
                "support": (
                    "both_models"
                    if len(hit_models) == 2
                    else hit_models[0] + "_only"
                    if all(value["state"] == "succeeded" for value in statuses.values())
                    else hit_models[0] + "_observed_other_failed"
                ),
                "model_states": {key: value["state"] for key, value in statuses.items()},
                "ga_support": (
                    "both_models"
                    if len(ga_models) == 2
                    else ga_models[0] + "_only"
                    if ga_models
                    else "none"
                ),
                "domain_count": len(evidence),
                "hit_ids": [h.hit_id for h in evidence],
                "ga_by_model": {
                    model: (
                        "passed_ga"
                        if any(
                            h.model.model_id == model and h.hmm.ga.status == "passed_ga"
                            for h in evidence
                        )
                        else next(h.hmm.ga.status for h in evidence if h.model.model_id == model)
                    )
                    for model in hit_models
                },
            }
        )
    return rows


def run_sprf(proteins, profiles, output: Path, config=None):
    """Run both models independently; preserve partial success and all domain evidence."""
    config = config or HMMSearchConfig()
    proteins, profiles = tuple(proteins), tuple(profiles)
    if len(profiles) != 2 or {p.provenance.model_id for p in profiles} != set(MODEL_IDS):
        raise InputError("Provide exactly TIGR03519 and PF11751 profiles with canonical model IDs")
    if any(p.method != "hmm" for p in profiles):
        raise InputError("PorP/SprF profiles must use method hmm")
    if len({p.protein_id for p in proteins}) != len(proteins):
        raise InputError("Duplicate internal protein IDs")
    if output.exists():
        raise InputError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".sprf-", dir=output.parent))
    profiles = {p.provenance.model_id: p for p in profiles}
    hits, statuses, per_model = [], {}, {}
    try:
        for model_id in MODEL_IDS:
            try:
                result = run_hmm(proteins, profiles[model_id], stage / model_id, config)
                per_model[model_id] = {k: v for k, v in result.items() if k != "output"}
                payload = json.loads((stage / model_id / "hits.json").read_text())
                for data in payload:
                    data.pop("ctd_distance", None)
                    hit = DomainHit.model_validate(data)
                    hits.append(
                        hit.model_copy(update={"raw_artifact": f"{model_id}/raw_domains.json"})
                    )
                statuses[model_id] = {"state": "succeeded"}
            except HMMError as exc:
                statuses[model_id] = {"state": "failed", "message": str(exc)}
        successful = sum(s["state"] == "succeeded" for s in statuses.values())
        state = "succeeded" if successful == 2 else "partial" if successful else "failed"
        evidence = collect_evidence(proteins, hits, statuses)
        report = {
            "state": state,
            "input_proteins": len(proteins),
            "models": statuses,
            "per_model": per_model,
            "domain_count": len(hits),
            "hit_proteins": len(evidence),
            "both_models_hit_proteins": sum(x["support"] == "both_models" for x in evidence),
            "any_model_ga_proteins": sum(bool(x["ga_models"]) for x in evidence),
            "both_models_ga_proteins": sum(x["ga_support"] == "both_models" for x in evidence),
            "interpretation": (
                "family evidence only; no substrate pairing or final SprB classification"
            ),
            "missing_evidence": (
                "unlisted model hits are absent only when that model succeeded; "
                "consult models status"
            ),
        }
        write_json(stage / "summary.json", report)
        write_json(stage / "hits.json", [h.model_dump(mode="json") for h in hits])
        write_json(stage / "protein_evidence.json", evidence)
        fields = [
            "protein_id",
            "source_id",
            "gene_id",
            "protein_length",
            "hit_models",
            "ga_models",
            "support",
            "ga_support",
            "domain_count",
            "hit_ids",
            "ga_by_model",
            "model_states",
        ]
        with (stage / "protein_evidence.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in evidence:
                writer.writerow(
                    {k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()}
                )
        write_json(
            stage / "run_metadata.json",
            {
                "software_version": __version__,
                "created_at": datetime.now(UTC).isoformat(),
                "configuration": config.model_dump(),
                "profiles": {k: p.model_dump(mode="json") for k, p in profiles.items()},
                "output_sha256": {
                    str(p.relative_to(stage)): sha256(p)
                    for p in sorted(stage.rglob("*"))
                    if p.is_file()
                },
            },
        )
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return report
