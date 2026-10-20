"""Offline merger regression: union, missing data, lossless records and provenance."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import DomainHit, GenomeDataset
from sprb_scanner.evidence import merge, run_merge
from sprb_scanner.evidence.models import SOURCES, EvidenceResult, NeighborhoodEvidence, SourceRun
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import export_prepared, write_json
from sprb_scanner.references import ReferenceRegistry
from sprb_scanner.scanners.mmseqs import SequenceAlignment, fingerprint

REGISTRY = Path(__file__).parents[2] / "references/registry.json"


@pytest.fixture
def dataset():
    return GenomeDataset.model_validate(
        {
            "genome": {"genome_id": "g", "input_format": "genbank"},
            "contigs": [
                {
                    "contig_id": "c",
                    "genome_id": "g",
                    "source_id": "chr",
                    "length": 20000,
                    "topology": "circular",
                }
            ],
            "genes": [
                {
                    "gene_id": f"g{i}",
                    "genome_id": "g",
                    "contig_id": "c",
                    "location_parts": [{"start": i * 3100, "end": i * 3100 + 3000}],
                    "strand": -1,
                }
                for i in range(4)
            ],
            "proteins": [
                {
                    "protein_id": f"p{i}",
                    "gene_id": f"g{i}",
                    "genome_id": "g",
                    "source_id": f"WP_{i}",
                    "sequence": "A" * 1000,
                }
                for i in range(4)
            ],
        }
    )


def domain(pid, model, start=10, number=1, below=False):
    score = 5 if below else 20
    return DomainHit.model_validate(
        {
            "hit_id": f"{model}::{pid}::{number}",
            "protein_id": pid,
            "method": "hmm",
            "model": {
                "model_id": model,
                "version": model + ".1",
                "source": "synthetic",
                "sha256": "a" * 64,
            },
            "query_id": model,
            "query_length": 100,
            "query_start": 0,
            "query_end": 90,
            "protein_length": 1000,
            "start": start,
            "end": start + 90,
            "evalue": 0.002,
            "evalue_kind": "independent_domain",
            "bit_score": score,
            "query_coverage": 0.9,
            "target_coverage": 0.09,
            "backend_version": "synthetic",
            "raw_artifact": "raw_domains.json",
            "hmm": {
                "sequence_evalue": 0.001,
                "sequence_bit_score": score,
                "independent_domain_evalue": 0.002,
                "conditional_domain_evalue": 0.0001,
                "envelope": {"start": start, "end": start + 90},
                "sequence_reported": True,
                "domain_reported": True,
                "sequence_included": True,
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


def alignment(pid="p3", number=1, start=300):
    query = ReferenceRegistry.load(REGISTRY).sequences("RN28")[0]
    digest = hashlib.sha256(query.sequence.encode()).hexdigest()
    return SequenceAlignment.model_validate(
        {
            "hit_id": f"rn28::{pid}::{number}",
            "protein_id": pid,
            "method": "mmseqs2",
            "query_id": query.protein_id,
            "query_source_id": query.source_id,
            "target_source_id": f"WP_{pid[1:]}",
            "target_gene_id": f"g{pid[1:]}",
            "model": {
                "model_id": query.protein_id,
                "version": "sequence-sha256:" + digest,
                "source": "synthetic",
                "sha256": digest,
            },
            "query_length": 207,
            "query_start": 0,
            "query_end": 207,
            "protein_length": 1000,
            "start": start,
            "end": start + 207,
            "query_coverage": 1,
            "target_coverage": 0.207,
            "sequence_identity": 1,
            "identical_residues": 207,
            "alignment_length": 207,
            "evalue": 1e-30,
            "evalue_kind": "alignment",
            "bit_score": 200,
            "backend_version": "synthetic",
            "raw_artifact": "raw_alignments.tsv",
        }
    )


def runs(hits=(), alignments=(), states=None):
    counts = Counter(h.model.model_id for h in hits)
    counts["RN28"] = len(alignments)
    states = states or {}
    return {
        s: SourceRun(
            source_id=s,
            state=states.get(s, "succeeded"),
            retained_hits=counts[s] if states.get(s, "succeeded") == "succeeded" else None,
            message="backend failed" if states.get(s) == "failed" else None,
        )
        for s in SOURCES
    }


def manifest(root, extra=None):
    write_json(
        root / "run_metadata.json",
        {
            **(extra or {}),
            "output_sha256": {
                str(p.relative_to(root)): sha256(p)
                for p in sorted(root.rglob("*"))
                if p.is_file() and p != root / "run_metadata.json"
            },
        },
    )


def native(root, model, dataset, hits=(), failed=False):
    root.mkdir(parents=True)
    profile = {"provenance": domain("p0", model).model.model_dump()}
    write_json(
        root / "status.json",
        {
            "model_id": model,
            "state": "failed" if failed else "succeeded",
            "message": "bad HMM" if failed else None,
        },
    )
    if not failed:
        write_json(root / "hits.json", [h.model_dump(mode="json") for h in hits])
        write_json(root / "raw_domains.json", [])
    manifest(
        root,
        {
            "profile": profile,
            "target_count": len(dataset.proteins),
            "protein_ids": {f"p{i}": p.protein_id for i, p in enumerate(dataset.proteins)},
            "input_proteins_sha256": fingerprint(dataset.proteins),
            "domain_count": len(hits),
        },
    )
    return profile


def pair(root, models, dataset, hits=(), failed=()):
    profiles = {
        m: native(
            root / m, m, dataset, [h for h in hits if h.model.model_id == m], failed=m in failed
        )
        for m in models
    }
    write_json(
        root / "hits.json",
        [
            h.model_copy(update={"raw_artifact": f"{m}/raw_domains.json"}).model_dump(mode="json")
            for m in models
            for h in hits
            if h.model.model_id == m
        ],
    )
    write_json(
        root / "summary.json",
        {
            "state": "partial" if len(failed) == 1 else "failed" if failed else "succeeded",
            "models": {m: {"state": "failed" if m in failed else "succeeded"} for m in models},
        },
    )
    manifest(root, {"profiles": profiles})
    return root


def sequence_run(root, dataset, hits=(), failed=False):
    root.mkdir()
    write_json(root / "status.json", {"state": "failed" if failed else "succeeded"})
    meta = {
        "kind": "mmseqs2-search",
        "state": "failed" if failed else "succeeded",
        "error": "backend failed" if failed else None,
    }
    if not failed:
        query = ReferenceRegistry.load(REGISTRY).sequences("RN28")[0]
        write_json(root / "query_mapping.json", {"q0": query.model_dump(mode="json")})
        write_json(
            root / "target_mapping.json",
            {f"t{i}": p.model_dump(mode="json") for i, p in enumerate(dataset.proteins)},
        )
        write_json(root / "alignments.json", [h.model_dump(mode="json") for h in hits])
        (root / "raw_alignments.tsv").write_text("".join("raw\n" for _ in hits))
        meta.update(
            backend_version="synthetic",
            query_proteins_sha256=fingerprint([query]),
            target_database_metadata={"proteins_sha256": fingerprint(dataset.proteins)},
            target_count=len(dataset.proteins),
            alignment_count=len(hits),
        )
    manifest(root, meta)
    return root


def prepared(tmp_path, dataset):
    root = tmp_path / "prepared"
    export_prepared(dataset, root, [], AppConfig())
    return root


def test_union_keeps_rn28_only_repeat_only_and_sprf_only(dataset):
    hits = [domain("p0", "TIGR04131", below=True), domain("p1", "PF13573"), domain("p2", "PF11751")]
    ali = [alignment(), alignment(number=2, start=600)]
    result = merge(dataset, runs(hits, ali), hits, ali)
    assert [p.protein_id for p in result.profiles] == ["p0", "p1", "p2", "p3"]
    p = result.profiles[-1]
    assert len(p.sources["RN28"].alignments) == 2
    assert p.sources["TIGR04131"].state == "succeeded_no_hits"
    assert p.sources["RN28"].alignments[1].start == 600
    assert result.profiles[0].sources["TIGR04131"].domains[0].hit.hmm.ga.status == "below_ga"
    assert p.neighborhood.outcome == "local_hit"
    assert p.neighborhood.neighbors[0].relationship.relative_position == "downstream"
    assert p.neighborhood.neighbors[0].relationship.directly_adjacent
    assert p.neighborhood.neighbors[0].relationship.distance_bp == 100
    EvidenceResult.model_validate(result.model_dump(exclude_computed_fields=True))


def test_domains_remain_independent_and_overlap_not_summed(dataset):
    hits = [domain("p0", "PF13573"), domain("p0", "PF13573", start=50, number=2, below=True)]
    p = merge(dataset, runs(hits), hits).profiles[0]
    ev = p.sources["PF13573"].domains
    assert len(ev) == 2 and [x.protein_coverage for x in ev] == [0.09, 0.09]
    assert ev[0].overlapping_same_model_hit_ids == (hits[1].hit_id,)
    assert ev[0].hit.hmm.sequence_evalue != ev[0].hit.evalue


@pytest.mark.parametrize("missing", ["not_run", "failed"])
def test_missing_not_negative_and_partial_neighborhood(dataset, missing):
    ali = [alignment()]
    states = {"TIGR03519": missing}
    p = merge(dataset, runs(alignments=ali, states=states), alignments=ali).profiles[0]
    assert p.sources["TIGR03519"].state == missing
    assert p.sources["PF11751"].state == "succeeded_no_hits"
    assert p.neighborhood.state == "partial" and p.neighborhood.outcome == "unknown"


def test_both_sprf_missing_complete_zero_and_disabled(dataset):
    ali = [alignment()]
    for s in ["not_run", "failed"]:
        p = merge(
            dataset,
            runs(alignments=ali, states={m: s for m in ["TIGR03519", "PF11751"]}),
            alignments=ali,
        ).profiles[0]
        assert p.neighborhood.state == s and p.neighborhood.outcome == "unknown"
    p = merge(dataset, runs(alignments=ali), alignments=ali).profiles[0]
    assert p.neighborhood.state == "succeeded" and p.neighborhood.outcome == "no_local_hit"
    assert (
        merge(dataset, runs(alignments=ali), alignments=ali, neighborhood=False)
        .profiles[0]
        .neighborhood.state
        == "not_run"
    )


def test_partial_positive_is_explicit_and_sprd_local_evidence_retained(dataset):
    hits = [domain("p2", "PF11751", below=True)]
    ali = [alignment()]
    p = merge(dataset, runs(hits, ali, {"TIGR03519": "failed"}), hits, ali).profiles[-1]
    assert p.neighborhood.state == "partial" and p.neighborhood.outcome == "local_hit"
    assert p.neighborhood.neighbors[0].domain_hit_ids == (hits[0].hit_id,)


def test_protein_only_neighborhood_unavailable(dataset):
    ds = GenomeDataset.model_validate(
        {
            "genome": {"genome_id": "g", "input_format": "protein_fasta"},
            "proteins": [
                p.model_copy(update={"gene_id": None}).model_dump(exclude_computed_fields=True)
                for p in dataset.proteins
            ],
        }
    )
    ali = [alignment()]
    p = merge(ds, runs(alignments=ali), alignments=ali).profiles[0]
    assert p.neighborhood.state == "unavailable" and p.neighborhood.outcome == "unknown"


@pytest.mark.parametrize("bad", [True, -1, 1.5])
def test_invalid_window(dataset, bad):
    with pytest.raises(InputError):
        merge(dataset, runs(), window=bad)


def test_duplicate_wrong_target_and_state_rejected(dataset):
    h = domain("p0", "PF13573")
    with pytest.raises(InputError, match="Duplicate"):
        merge(dataset, runs([h, h]), [h, h])
    with pytest.raises(InputError, match="protein mismatch"):
        merge(dataset, runs([h]), [h.model_copy(update={"protein_id": "missing"})])
    with pytest.raises(InputError, match="unsuccessful"):
        merge(dataset, runs(states={"PF13573": "failed"}), [h])
    with pytest.raises(InputError, match="count"):
        merge(dataset, runs([h]))


def test_circular_unknown_strand_and_geometry_error_state(dataset):
    hits = [domain("p0", "PF11751")]
    ali = [alignment()]
    p = merge(dataset, runs(hits, ali), hits, ali, window=1).profiles[-1]
    assert p.neighborhood.neighbors[0].relationship.index_crosses_origin
    bad = dataset.model_copy(
        update={
            "genes": (
                dataset.genes[0].model_copy(
                    update={"location_parts": (dataset.genes[0].location_parts[0],) * 2}
                ),
                *dataset.genes[1:],
            )
        }
    )
    p = merge(bad, runs(hits, ali), hits, ali).profiles[-1]
    assert p.neighborhood.state == "failed" and "slippage" in p.neighborhood.message
    assert p.neighborhood.outcome == "unknown"


def test_no_local_hit_requires_complete_and_source_status_contract():
    with pytest.raises(ValueError):
        NeighborhoodEvidence(state="partial", outcome="no_local_hit", window=5)
    with pytest.raises(ValueError):
        SourceRun(source_id="RN28", state="failed")
    with pytest.raises(ValueError):
        SourceRun(source_id="RN28", state="not_run", retained_hits=0)


def test_cli_partial_sources_complete_exports_and_checksums(tmp_path, dataset):
    prep = prepared(tmp_path, dataset)
    bctd = pair(
        tmp_path / "bctd",
        ("TIGR04131", "PF13585"),
        dataset,
        [domain("p0", "TIGR04131", below=True)],
        failed=("PF13585",),
    )
    sprf = pair(tmp_path / "sprf", ("TIGR03519", "PF11751"), dataset, [domain("p2", "PF11751")])
    repeat = tmp_path / "repeats"
    native(repeat, "PF13573", dataset, [domain("p1", "PF13573")])
    rn = sequence_run(tmp_path / "rn28", dataset, [alignment(), alignment(number=2, start=600)])
    out = tmp_path / "merged"
    result = CliRunner().invoke(
        app,
        [
            "merge-evidence",
            "--prepared",
            str(prep),
            "--output",
            str(out),
            "--bctd-scan",
            str(bctd),
            "--sprf-scan",
            str(sprf),
            "--repeat-scan",
            str(repeat),
            "--rn28-scan",
            str(rn),
            "--registry",
            str(REGISTRY),
        ],
    )
    assert result.exit_code == 0, result.output
    data = json.loads((out / "evidence_profiles.json").read_text())
    assert len(data["profiles"]) == 4
    assert data["profiles"][-1]["sources"]["TIGR04131"]["state"] == "succeeded_no_hits"
    assert data["profiles"][-1]["sources"]["PF13585"]["state"] == "failed"
    assert data["profiles"][-1]["sources"]["RN28"]["alignments"][0]["start"] == 300
    assert json.loads((out / "status.json").read_text()) == {
        "state": "succeeded",
        "evidence_complete": False,
    }
    meta = json.loads((out / "run_metadata.json").read_text())
    for name, digest in meta["output_sha256"].items():
        assert sha256(out / name) == digest
    assert meta["input_files_sha256"][str(REGISTRY)] == sha256(REGISTRY)
    for filename in [
        "evidence_profiles.tsv",
        "evidence_records.tsv",
        "neighborhood_associations.tsv",
    ]:
        assert (out / filename).read_text().strip()
    with pytest.raises(InputError, match="already exists"):
        run_merge(prep, out)


def test_empty_union_headers_and_not_run_states(tmp_path, dataset):
    out = tmp_path / "merged"
    report = run_merge(prepared(tmp_path, dataset), out)
    assert report["profiles"] == 0
    assert all(
        r["state"] == "not_run" and r["retained_hits"] is None
        for r in report["source_runs"].values()
    )
    assert len((out / "evidence_records.tsv").read_text().splitlines()) == 1


@pytest.mark.parametrize("failure_kind", ["hmm", "rn28"])
def test_failed_sources_do_not_become_zero_hits(tmp_path, dataset, failure_kind):
    kwargs = {}
    if failure_kind == "hmm":
        path = tmp_path / "repeat"
        native(path, "PF13573", dataset, failed=True)
        kwargs["repeat_scan"] = path
        source = "PF13573"
    else:
        path = sequence_run(tmp_path / "rn28", dataset, failed=True)
        kwargs["rn28_scan"] = path
        source = "RN28"
    report = run_merge(prepared(tmp_path, dataset), tmp_path / "merged", **kwargs)
    assert report["source_runs"][source]["state"] == "failed"
    assert report["source_runs"][source]["retained_hits"] is None


def test_success_zero_rn28_source_known_even_empty_union(tmp_path, dataset):
    rn = sequence_run(tmp_path / "rn28", dataset)
    report = run_merge(
        prepared(tmp_path, dataset), tmp_path / "merged", rn28_scan=rn, registry=REGISTRY
    )
    assert report["source_runs"]["RN28"]["state"] == "succeeded"
    assert report["source_runs"]["RN28"]["retained_hits"] == 0


@pytest.mark.parametrize(
    "damage",
    [
        "checksum",
        "protein_ids",
        "sequences",
        "wrong_model",
        "pair_hits",
        "status",
        "coverage",
        "manifest_path",
    ],
)
def test_corrupt_exports_fail_closed(tmp_path, dataset, damage):
    root = pair(tmp_path / "bctd", ("TIGR04131", "PF13585"), dataset, [domain("p0", "TIGR04131")])
    mp = root / "TIGR04131/run_metadata.json"
    meta = json.loads(mp.read_text())
    if damage == "checksum":
        (root / "TIGR04131/hits.json").write_text("[]")
    elif damage == "protein_ids":
        meta["protein_ids"]["p0"] = "wrong"
        write_json(mp, meta)
    elif damage == "sequences":
        meta["input_proteins_sha256"] = "b" * 64
        write_json(mp, meta)
    elif damage == "wrong_model":
        meta["profile"]["provenance"]["model_id"] = "wrong"
        write_json(mp, meta)
    elif damage == "pair_hits":
        write_json(root / "hits.json", [])
    elif damage == "status":
        write_json(root / "summary.json", {"models": {}, "state": "succeeded"})
    elif damage == "coverage":
        hp = root / "TIGR04131/hits.json"
        h = json.loads(hp.read_text())
        h[0]["query_coverage"] = 0.1
        write_json(hp, h)
        manifest(root / "TIGR04131", {k: v for k, v in meta.items() if k != "output_sha256"})
    else:
        meta["output_sha256"]["../../outside"] = "a" * 64
        write_json(mp, meta)
    if damage != "checksum":
        top = json.loads((root / "run_metadata.json").read_text())
        manifest(root, {k: v for k, v in top.items() if k != "output_sha256"})
    out = tmp_path / "merged"
    with pytest.raises(InputError):
        run_merge(prepared(tmp_path, dataset), out, bctd_scan=root)
    assert not out.exists()


@pytest.mark.parametrize(
    "damage",
    [
        "query_sequence",
        "target_sequence",
        "target_source",
        "hit_target",
        "raw_count",
        "query_fingerprint",
    ],
)
def test_mmseqs_identity_validation(tmp_path, dataset, damage):
    root = sequence_run(tmp_path / "rn28", dataset, [alignment()])
    meta = json.loads((root / "run_metadata.json").read_text())
    if damage in {"query_sequence", "target_sequence", "target_source"}:
        name = "query_mapping.json" if damage == "query_sequence" else "target_mapping.json"
        p = root / name
        d = json.loads(p.read_text())
        key = next(iter(d))
        d[key]["source_id" if damage == "target_source" else "sequence"] = "INVALID"
        write_json(p, d)
    elif damage == "hit_target":
        p = root / "alignments.json"
        d = json.loads(p.read_text())
        d[0]["target_source_id"] = "wrong"
        write_json(p, d)
    elif damage == "raw_count":
        (root / "raw_alignments.tsv").write_text("")
    else:
        meta["query_proteins_sha256"] = "b" * 64
    manifest(root, {k: v for k, v in meta.items() if k != "output_sha256"})
    with pytest.raises(InputError):
        run_merge(
            prepared(tmp_path, dataset), tmp_path / "merged", rn28_scan=root, registry=REGISTRY
        )


def test_provided_missing_path_is_error_not_not_run(tmp_path, dataset):
    with pytest.raises(InputError):
        run_merge(prepared(tmp_path, dataset), tmp_path / "merged", rn28_scan=tmp_path / "missing")


def test_three_rn28_positive_ctd_negative_proteins_are_not_lost(dataset):
    ali = [alignment(pid=f"p{i}") for i in range(3)]
    result = merge(dataset, runs(alignments=ali), alignments=ali)
    assert {p.protein_id for p in result.profiles} == {"p0", "p1", "p2"}
    assert all(p.sources["TIGR04131"].state == "succeeded_no_hits" for p in result.profiles)
    assert all(p.sources["PF13585"].state == "succeeded_no_hits" for p in result.profiles)


def test_all_failed_pair_preserves_both_failure_states(tmp_path, dataset):
    root = pair(
        tmp_path / "bctd", ("TIGR04131", "PF13585"), dataset, failed=("TIGR04131", "PF13585")
    )
    report = run_merge(prepared(tmp_path, dataset), tmp_path / "merged", bctd_scan=root)
    assert report["source_runs"]["TIGR04131"]["state"] == "failed"
    assert report["source_runs"]["PF13585"]["state"] == "failed"
    assert report["profiles"] == 0


def test_empty_union_preserves_neighborhood_failure(dataset):
    bad = dataset.model_copy(
        update={
            "genes": (
                dataset.genes[0].model_copy(
                    update={"location_parts": (dataset.genes[0].location_parts[0],) * 2}
                ),
                *dataset.genes[1:],
            )
        }
    )
    result = merge(bad, runs())
    assert not result.profiles
    assert result.neighborhood_state == "failed"
    assert "slippage" in result.neighborhood_message


def test_unreported_domains_retained_without_a_new_reporting_cutoff(dataset):
    h = domain("p0", "PF13573", below=True)
    h = h.model_copy(update={"hmm": h.hmm.model_copy(update={"domain_reported": False})})
    result = merge(dataset, runs([h]), [h])
    assert result.profiles[0].sources["PF13573"].domains[0].hit == h


@pytest.mark.parametrize(
    ("failed_genes", "sprf_failed"),
    [
        (("g0",), False),
        (("g1",), False),
        (("g3",), False),
        (("g0", "g1", "g2", "g3"), False),
        (("g1",), True),
    ],
)
def test_candidate_query_failures_are_isolated_and_exported(
    tmp_path, dataset, monkeypatch, failed_genes, sprf_failed
):
    import csv

    from sprb_scanner.neighborhood import NeighborhoodIndex

    original_query = NeighborhoodIndex.query
    queried = []

    def fail_selected(self, gene_id, window=5):
        # Use a real, successfully built index; fail only the requested candidate query.
        assert gene_id in self.by_gene
        queried.append(gene_id)
        if gene_id in failed_genes:
            raise InputError("synthetic candidate geometry error")
        return original_query(self, gene_id, window)

    monkeypatch.setattr(NeighborhoodIndex, "query", fail_selected)
    bctd = pair(
        tmp_path / "bctd",
        ("TIGR04131", "PF13585"),
        dataset,
        [domain("p0", "TIGR04131", below=True)],
    )
    sprf = pair(
        tmp_path / "sprf",
        ("TIGR03519", "PF11751"),
        dataset,
        [domain("p2", "PF11751")],
        failed=("TIGR03519",) if sprf_failed else (),
    )
    repeats = tmp_path / "repeats"
    native(repeats, "PF13573", dataset, [domain("p1", "PF13573")])
    rn28 = sequence_run(tmp_path / "rn28", dataset, [alignment()])
    output = tmp_path / "merged"
    report = run_merge(
        prepared(tmp_path, dataset),
        output,
        bctd_scan=bctd,
        sprf_scan=sprf,
        repeat_scan=repeats,
        rn28_scan=rn28,
        registry=REGISTRY,
    )
    result = json.loads((output / "evidence_profiles.json").read_text())
    assert queried == ["g0", "g1", "g2", "g3"]
    assert report["profiles"] == len(result["profiles"]) == 4
    expected_state = "failed" if len(failed_genes) == 4 else "partial"
    assert report["neighborhood_run_state"] == result["neighborhood_state"] == expected_state
    assert (
        f"{len(failed_genes)} candidate neighborhood queries failed"
        in (report["neighborhood_run_message"])
    )
    if sprf_failed:
        assert "Only one PorP/SprF model succeeded" in report["neighborhood_run_message"]
    expected_genes = {g.gene_id: g.model_dump(mode="json") for g in dataset.genes}
    for profile in result["profiles"]:
        context = profile["neighborhood"]
        gene_id = profile["gene"]["gene_id"]
        assert profile["gene"] == expected_genes[gene_id]
        assert any(e["domains"] or e["alignments"] for e in profile["sources"].values())
        if gene_id in failed_genes:
            assert context["state"] == "failed" and context["outcome"] == "unknown"
            assert not context["neighbors"]
            assert profile["protein_id"] in context["message"] and gene_id in context["message"]
            assert "synthetic candidate geometry error" in context["message"]
        else:
            assert context["state"] == ("partial" if sprf_failed else "succeeded")
            assert context["message"] is None or "Only one PorP/SprF" in context["message"]
    assert json.loads((output / "status.json").read_text()) == {
        "state": "succeeded",
        "evidence_complete": False,
    }
    with (output / "evidence_profiles.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row["gene_id"] in failed_genes:
            assert (
                row["neighborhood_state"] == "failed" and row["neighborhood_outcome"] == "unknown"
            )
            assert row["local_sprf_proteins"] == row["directly_adjacent_sprf_proteins"] == ""
    with (output / "neighborhood_associations.tsv").open() as handle:
        relations = list(csv.DictReader(handle, delimiter="\t"))
    assert all("g" + r["protein_id"][1:] not in failed_genes for r in relations)
    with (output / "evidence_records.tsv").open() as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 4
    metadata = json.loads((output / "run_metadata.json").read_text())
    for name, digest in metadata["output_sha256"].items():
        assert sha256(output / name) == digest


def test_unexpected_query_exception_is_not_silently_ignored(dataset, monkeypatch):
    from sprb_scanner.neighborhood import NeighborhoodIndex

    def unexpected(self, gene_id, window=5):
        raise RuntimeError("unexpected programming error")

    monkeypatch.setattr(NeighborhoodIndex, "query", unexpected)
    ali = [alignment()]
    with pytest.raises(RuntimeError, match="unexpected programming error"):
        merge(dataset, runs(alignments=ali), alignments=ali)
