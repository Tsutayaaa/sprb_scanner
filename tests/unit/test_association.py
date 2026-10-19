import copy
import hashlib
import json

import pytest
from typer.testing import CliRunner

from sprb_scanner.association import associate
from sprb_scanner.association.analyzer import CTD_MODELS, SPRF_MODELS
from sprb_scanner.association.inputs import load_pair
from sprb_scanner.association.service import run_association
from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import DomainHit, GenomeDataset
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import export_prepared, write_json


def hit(pid, model, start=0, end=90, length=100, number=1, below=False):
    score = 5 if below else 20
    return DomainHit.model_validate(
        {
            "start": start,
            "end": end,
            "hit_id": f"{model}::{pid}::{number}",
            "protein_id": pid,
            "method": "hmm",
            "model": {
                "model_id": model,
                "version": model + ".1",
                "source": "synthetic",
                "sha256": "a" * 64,
            },
            "query_id": model + ".1",
            "query_length": 100,
            "query_start": 0,
            "query_end": 90,
            "protein_length": length,
            "evalue": 0.002,
            "evalue_kind": "independent_domain",
            "bit_score": score,
            "query_coverage": 0.9,
            "target_coverage": (end - start) / length,
            "backend_version": "synthetic",
            "raw_artifact": "raw_domains.json",
            "hmm": {
                "sequence_evalue": 0.001,
                "sequence_bit_score": score,
                "independent_domain_evalue": 0.002,
                "conditional_domain_evalue": 0.0001,
                "envelope": {"start": start, "end": end},
                "sequence_reported": True,
                "sequence_included": True,
                "domain_reported": True,
                "domain_included": True,
                "ga": {
                    "status": "below_ga" if below else "passed_ga",
                    "sequence_cutoff": 10,
                    "domain_cutoff": 10,
                    "sequence_margin": score - 10,
                    "domain_margin": score - 10,
                    "sequence_passed": not below,
                    "domain_passed": not below,
                },
            },
        }
    )


@pytest.fixture
def data():
    return GenomeDataset.model_validate(
        {
            "genome": {"genome_id": "g", "input_format": "genbank"},
            "contigs": [{"contig_id": "c", "genome_id": "g", "source_id": "chr", "length": 1000}],
            "genes": [
                {
                    "gene_id": f"g{i}",
                    "genome_id": "g",
                    "contig_id": "c",
                    "location_parts": [{"start": i * 100, "end": i * 100 + 90}],
                    "strand": -1,
                }
                for i in range(3)
            ],
            "proteins": [
                {
                    "protein_id": f"p{i}",
                    "gene_id": f"g{i}",
                    "genome_id": "g",
                    "source_id": f"WP_{i}",
                    "sequence": "A" * (1000 if i == 2 else 100),
                }
                for i in range(3)
            ],
        }
    )


def scans(tmp_path, data, ctd=None, sprf=None):
    def pair(name, models, hits):
        root = tmp_path / name
        root.mkdir()
        profiles = {}
        for model in models:
            path = root / model
            path.mkdir()
            group = [h for h in hits if h.model.model_id == model]
            profile = {
                "provenance": {
                    "model_id": model,
                    "version": model + ".1",
                    "source": "synthetic",
                    "sha256": "a" * 64,
                }
            }
            profiles[model] = profile
            write_json(path / "hits.json", [h.model_dump(mode="json") for h in group])
            write_json(path / "raw_domains.json", [])
            write_json(path / "status.json", {"state": "succeeded", "model_id": model})
            write_json(
                path / "run_metadata.json",
                {
                    "profile": profile,
                    "target_count": len(data.proteins),
                    "protein_ids": {f"p{i}": p.protein_id for i, p in enumerate(data.proteins)},
                    "input_proteins_sha256": hashlib.sha256(
                        json.dumps(
                            [(p.protein_id, p.sequence) for p in data.proteins],
                            ensure_ascii=False,
                            separators=(",", ":"),
                        ).encode()
                    ).hexdigest(),
                    "output_sha256": {p.name: sha256(p) for p in path.iterdir()},
                },
            )
        write_json(
            root / "hits.json",
            [
                h.model_copy(
                    update={"raw_artifact": f"{h.model.model_id}/raw_domains.json"}
                ).model_dump(mode="json")
                for m in models
                for h in hits
                if h.model.model_id == m
            ],
        )
        write_json(
            root / "summary.json",
            {"state": "succeeded", "models": {m: {"state": "succeeded"} for m in models}},
        )
        write_json(
            root / "run_metadata.json",
            {
                "profiles": profiles,
                "output_sha256": {
                    str(p.relative_to(root)): sha256(p) for p in root.rglob("*") if p.is_file()
                },
            },
        )
        return root

    return (
        pair("bctd", CTD_MODELS, ctd if ctd is not None else [hit("p1", "TIGR04131", below=True)]),
        pair(
            "sprf",
            SPRF_MODELS,
            sprf
            if sprf is not None
            else [hit("p0", "TIGR03519"), hit("p2", "PF11751", length=1000)],
        ),
    )


def test_below_ga_one_to_many_positions_and_unmerged_domains(data):
    ctd = [hit("p1", "TIGR04131", below=True)]
    sprf = [
        hit("p0", "TIGR03519"),
        hit("p2", "PF11751", length=1000),
        hit("p2", "PF11751", start=50, end=140, length=1000, number=2, below=True),
    ]
    result = associate(data, ctd, sprf)
    (cargo,) = result.cargos
    assert cargo.ga_models == () and cargo.status == "local_hit"
    assert len(cargo.neighbors) == 2
    assert cargo.neighbors[0].relationship.relative_position == "downstream"
    assert cargo.neighbors[0].relationship.distance_bp == 10
    assert all(n.relationship.directly_adjacent for n in cargo.neighbors)
    domains = [d for d in result.domains if d.hit.protein_id == "p2"]
    assert len(domains) == 2 and domains[1].domain_order_by_protein_start == 2
    assert all(d.protein_coverage == 0.09 and d.hmm_coverage == 0.9 for d in domains)
    assert domains[0].overlapping_same_model_hit_ids == (domains[1].hit.hit_id,)
    assert all(d.domains_same_protein_model == 2 for d in domains)
    assert domains[1].hit.hmm.ga.status == "below_ga"
    assert not hasattr(cargo, "classification")
    empty = associate(data, ctd, sprf, window=0).cargos[0]
    assert empty.status == "no_local_hit" and not empty.neighbors
    assert not associate(data, [], sprf).cargos


@pytest.mark.parametrize("failure", ["unknown", "length", "coverage", "model", "duplicate"])
def test_invalid_hit_inputs(data, failure):
    h = hit("p1", "TIGR04131")
    if failure == "unknown":
        h = h.model_copy(update={"protein_id": "other"})
    if failure == "length":
        h = h.model_copy(update={"protein_length": 101})
    if failure == "coverage":
        h = h.model_copy(update={"target_coverage": 0.5})
    if failure == "model":
        h = hit("p1", "TIGR03519")
    with pytest.raises(InputError):
        associate(data, [h, h] if failure == "duplicate" else [h], [])


@pytest.mark.parametrize("window", [-1, True, 1.5])
def test_invalid_windows(data, window):
    with pytest.raises(InputError):
        associate(data, [], [], window)


def test_cli_exports_identity_integrity_and_no_overwrite(data, tmp_path):
    bctd, sprf = scans(tmp_path, data)
    prepared = tmp_path / "prepared"
    export_prepared(data, prepared, [], AppConfig())
    out = tmp_path / "out"
    args = [
        "typeb-sprf-associate",
        "--prepared",
        str(prepared),
        "--bctd-scan",
        str(bctd),
        "--sprf-scan",
        str(sprf),
        "--output",
        str(out),
    ]
    runner = CliRunner()
    response = runner.invoke(app, args)
    assert response.exit_code == 0, response.output
    payload = json.loads((out / "associations.json").read_text())
    assert payload["cargos"][0]["status"] == "local_hit"
    assert len(payload["domains"]) == 3
    import csv

    rows = list(csv.DictReader((out / "domain_evidence.tsv").open(), delimiter="\t"))
    assert len(rows) == 3 and all(r["hmm_coverage"] for r in rows)
    assert all(json.loads(r["overlapping_same_model_hit_ids"]) == [] for r in rows)
    for row in rows:
        h = json.loads(row["raw_hit"])
        assert float(row["sequence_evalue"]) == h["hmm"]["sequence_evalue"]
        assert float(row["independent_domain_evalue"]) == h["hmm"]["independent_domain_evalue"]
    meta = json.loads((out / "run_metadata.json").read_text())
    for name, digest in meta["output_sha256"].items():
        assert sha256(out / name) == digest
    assert meta["parameters"]["coverage_cutoff"] is None
    assert runner.invoke(app, args).exit_code == 2
    args[-1] = str(tmp_path / "zero")
    assert runner.invoke(app, args + ["--window", "0"]).exit_code == 0
    (bctd / "hits.json").write_text("[]")
    args[-1] = str(tmp_path / "tampered")
    assert runner.invoke(app, args).exit_code == 2
    assert not (tmp_path / "tampered").exists()


def test_successful_empty_scans(data, tmp_path):
    bctd, sprf = scans(tmp_path, data, ctd=[], sprf=[])
    prepared = tmp_path / "prepared"
    export_prepared(data, prepared, [], AppConfig())
    assert run_association(prepared, bctd, sprf, tmp_path / "empty")["cargo_count"] == 0


@pytest.mark.parametrize("issue", ["sequence", "id", "failed", "legacy", "path", "merged"])
def test_scan_input_rejection(data, tmp_path, issue):
    bctd, _ = scans(tmp_path, data)
    if issue in {"sequence", "id"}:
        p = data.proteins[0].model_copy(
            update={"sequence": "C" * 100} if issue == "sequence" else {"protein_id": "elsewhere"}
        )
        data = data.model_copy(update={"proteins": (p, *data.proteins[1:])})
    elif issue in {"failed", "path", "merged"}:
        meta = json.loads((bctd / "run_metadata.json").read_text())
        if issue == "failed":
            report = json.loads((bctd / "summary.json").read_text())
            report["state"] = "partial"
            write_json(bctd / "summary.json", report)
            meta["output_sha256"]["summary.json"] = sha256(bctd / "summary.json")
        elif issue == "path":
            meta["output_sha256"]["../outside"] = "a" * 64
        else:
            write_json(bctd / "hits.json", [])
            meta["output_sha256"]["hits.json"] = sha256(bctd / "hits.json")
        write_json(bctd / "run_metadata.json", meta)
    else:
        f = bctd / "TIGR04131/run_metadata.json"
        meta = json.loads(f.read_text())
        meta.pop("input_proteins_sha256")
        write_json(f, meta)
        top = json.loads((bctd / "run_metadata.json").read_text())
        top["output_sha256"]["TIGR04131/run_metadata.json"] = sha256(f)
        write_json(bctd / "run_metadata.json", top)
    with pytest.raises(InputError):
        load_pair(bctd, data, CTD_MODELS, {})


def test_circular_compound_overlap_and_contig_isolation(data):
    payload = copy.deepcopy(data.model_dump(mode="json", exclude_computed_fields=True))
    payload["contigs"][0]["topology"] = "circular"
    payload["genes"][0]["location_parts"] = [{"start": 0, "end": 20}, {"start": 950, "end": 1000}]
    payload["genes"][1]["location_parts"] = [{"start": 940, "end": 970}]
    result = associate(
        GenomeDataset.model_validate(payload), [hit("p0", "TIGR04131")], [hit("p1", "TIGR03519")], 1
    )
    n = result.cargos[0].neighbors[0].relationship
    assert n.overlap_bp == 20 and n.distance_bp == 0 and n.signed_gap_bp == -20
    assert n.crosses_origin
    payload["contigs"].append(
        {"contig_id": "other", "genome_id": "g", "source_id": "other", "length": 1000}
    )
    payload["genes"][1]["contig_id"] = "other"
    assert (
        associate(
            GenomeDataset.model_validate(payload),
            [hit("p0", "TIGR04131")],
            [hit("p1", "TIGR03519")],
        )
        .cargos[0]
        .status
        == "no_local_hit"
    )
