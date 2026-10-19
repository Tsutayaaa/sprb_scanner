"""Preserve every PF13573 domain and summarize repeat evidence per protein."""

import csv
import json

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import write_json
from sprb_scanner.scanners.hmm import run_hmm


def run_repeats(proteins, profile, output, config=None):
    if profile.provenance.model_id != "PF13573":
        raise InputError("SprB repeat scanning requires the PF13573 model")
    proteins = tuple(proteins)
    result = run_hmm(proteins, profile, output, config)
    hits = json.loads((output / "hits.json").read_text())
    rows = []
    for protein in proteins:
        evidence = [h for h in hits if h["protein_id"] == protein.protein_id]
        if not evidence:
            continue
        rows.append(
            {
                "protein_id": protein.protein_id,
                "source_id": protein.source_id,
                "protein_length": protein.length,
                "all_repeat_domains": len(evidence),
                "passed_ga_domains": sum(h["hmm"]["ga"]["status"] == "passed_ga" for h in evidence),
                "below_ga_domains": sum(h["hmm"]["ga"]["status"] == "below_ga" for h in evidence),
                "ga_unavailable_domains": sum(
                    h["hmm"]["ga"]["status"] == "unavailable" for h in evidence
                ),
                "reported_domains": sum(h["hmm"]["domain_reported"] for h in evidence),
            }
        )
    write_json(output / "repeat_counts.json", rows)
    fields = [
        "protein_id",
        "source_id",
        "protein_length",
        "all_repeat_domains",
        "passed_ga_domains",
        "below_ga_domains",
        "ga_unavailable_domains",
        "reported_domains",
    ]
    with (output / "repeat_counts.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    meta = json.loads((output / "run_metadata.json").read_text())
    meta["repeat_count_unit"] = (
        "backend-returned domain alignments, not deduplicated physical repeats"
    )
    for name in ("repeat_counts.json", "repeat_counts.tsv"):
        meta["output_sha256"][name] = sha256(output / name)
    write_json(output / "run_metadata.json", meta)
    return result | {"repeat_hit_proteins": len(rows)}
