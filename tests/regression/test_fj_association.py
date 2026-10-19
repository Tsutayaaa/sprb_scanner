"""Opt-in regression using complete, previously generated FJ paired HMM exports."""

import json
import os
from pathlib import Path

import pytest

from sprb_scanner.association.service import run_association

REFERENCE = json.loads(
    (Path(__file__).parents[1] / "fixtures/fj_association_reference.json").read_text()
)


@pytest.mark.external_data
def test_fj_twelve_cargos_and_sprb_sprf_sprd(tmp_path):
    names = ("SPRB_FJ_PREPARED", "SPRB_FJ_BCTD_SCAN", "SPRB_FJ_SPRF_SCAN")
    if not all(os.environ.get(n) for n in names):
        pytest.skip("Set SPRB_FJ_PREPARED, SPRB_FJ_BCTD_SCAN and SPRB_FJ_SPRF_SCAN; offline only")
    out = tmp_path / "association"
    run_association(*(Path(os.environ[n]) for n in names), out)
    result = json.loads((out / "associations.json").read_text())
    by_old = {
        tag: c for c in result["cargos"] for tag in c["gene"]["qualifiers"].get("old_locus_tag", [])
    }
    domains = {d["hit"]["hit_id"]: d for d in result["domains"]}
    counts = {"adjacent": 0, "orphan": 0}
    for ref in REFERENCE["reference_cargos"]:
        cargo = by_old[ref["old_locus_tag"]]
        assert cargo["source_protein_id"] == ref["source_protein_id"]
        counts[ref["group"]] += 1
        if ref["group"] == "orphan":
            assert cargo["status"] == "no_local_hit" and not cargo["neighbors"]
        else:
            assert cargo["status"] == "local_hit"
            adjacent = [n for n in cargo["neighbors"] if n["relationship"]["directly_adjacent"]]
            observed = {
                tag
                for n in adjacent
                for tag in n["relationship"]["gene"]["qualifiers"]["old_locus_tag"]
            }
            assert observed == set(ref["adjacent_loci"])
            for n in adjacent:
                assert {
                    domains[h]["hit"]["model"]["model_id"]
                    for h in n["domain_hit_ids"]
                    if domains[h]["hit"]["hmm"]["ga"]["status"] == "passed_ga"
                } == {"TIGR03519", "PF11751"}
    assert counts == {"adjacent": 9, "orphan": 3}
    for tag, expected in [("Fjoh_0978", 43), ("Fjoh_3951", 61)]:
        cargo = by_old["Fjoh_0979" if tag == "Fjoh_0978" else "Fjoh_3952"]
        n = next(
            n
            for n in cargo["neighbors"]
            if tag in n["relationship"]["gene"]["qualifiers"]["old_locus_tag"]
        )
        assert n["relationship"]["distance_bp"] == expected
        assert n["relationship"]["relative_position"] == "downstream"
        assert n["relationship"]["same_strand"]
    sprb = by_old["Fjoh_0979"]
    refs = {}
    for n in sprb["neighbors"]:
        for tag in n["relationship"]["gene"]["qualifiers"]["old_locus_tag"]:
            refs[tag] = {
                domains[h]["hit"]["model"]["model_id"]: domains[h] for h in n["domain_hit_ids"]
            }
    for model in ("TIGR03519", "PF11751"):
        sprf, sprd = refs["Fjoh_0978"][model], refs["Fjoh_0980"][model]
        assert sprf["hmm_coverage"] > 0.99 and sprd["hmm_coverage"] > 0.96
        assert sprf["protein_coverage"] > 0.88 and sprd["protein_coverage"] < 0.20
        assert sprd["hit"]["end"] <= 301 and sprd["hit"]["protein_length"] == 1588
        assert sprd["hit"]["hmm"]["ga"]["status"] == "passed_ga"
    # Only a positional evidence schema is exported, no typical shuttle/false-positive label.
    assert set(refs["Fjoh_0980"]["TIGR03519"]) == {
        "hit",
        "scan_directory",
        "domain_order_by_protein_start",
        "domains_same_protein_model",
        "overlapping_same_model_hit_ids",
        "hmm_coverage",
        "protein_coverage",
        "domain_start_fraction",
        "domain_end_fraction",
        "n_terminal_flank_aa",
        "c_terminal_flank_aa",
    }
    for tag in ("Fjoh_1645", "Fjoh_3952"):
        hits = [domains[h]["hit"] for h in by_old[tag]["domain_hit_ids"]]
        assert any(
            h["model"]["model_id"] == "TIGR04131" and h["hmm"]["ga"]["status"] == "below_ga"
            for h in hits
        )
