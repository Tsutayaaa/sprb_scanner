"""Validated existing-scan association and atomic, lossless JSON/TSV exports."""

import csv
import json
import shutil
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from sprb_scanner import __version__
from sprb_scanner.association.analyzer import CTD_MODELS, SPRF_MODELS, associate
from sprb_scanner.association.inputs import check_unchanged, load_pair, read_manifest
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared, write_json
from sprb_scanner.neighborhood.analyzer import GEOMETRY_POLICY


def write_rows(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: json.dumps(v, ensure_ascii=False)
                    if isinstance(v, (list, dict, tuple))
                    else v
                    for k, v in row.items()
                }
            )


def gene_fields(gene, protein_id, source_id, length, prefix):
    return {
        prefix + "protein_id": protein_id,
        prefix + "source_protein_id": source_id,
        prefix + "protein_length": length,
        prefix + "gene_id": gene.gene_id,
        prefix + "locus_tag": gene.locus_tag,
        prefix + "old_locus_tags": list(gene.qualifiers.get("old_locus_tag", ())),
        prefix + "contig_id": gene.contig_id,
        prefix + "location_parts": [p.model_dump() for p in gene.location_parts],
        prefix + "strand": gene.strand,
        prefix + "product": gene.product,
    }


def export_tables(stage, result, dataset):
    domains = {e.hit.hit_id: e for e in result.domains}
    rows = []
    for cargo in result.cargos:
        base = gene_fields(
            cargo.gene, cargo.protein_id, cargo.source_protein_id, cargo.protein_length, "cargo_"
        ) | {
            "status": cargo.status,
            "cargo_hit_models": list(cargo.hit_models),
            "cargo_ga_models": list(cargo.ga_models),
            "cargo_domain_evidence": [
                domains[h].model_dump(mode="json") for h in cargo.domain_hit_ids
            ],
        }
        for neighbor in cargo.neighbors:
            relation = neighbor.relationship
            rows.append(
                base
                | gene_fields(
                    relation.gene,
                    neighbor.protein_id,
                    neighbor.source_protein_id,
                    neighbor.protein_length,
                    "neighbor_",
                )
                | relation.model_dump(mode="json", exclude={"gene"})
                | {
                    "neighbor_domain_evidence": [
                        domains[h].model_dump(mode="json") for h in neighbor.domain_hit_ids
                    ],
                }
            )
        if not cargo.neighbors:
            rows.append(base)
    identity = [
        "protein_id",
        "source_protein_id",
        "protein_length",
        "gene_id",
        "locus_tag",
        "old_locus_tags",
        "contig_id",
        "location_parts",
        "strand",
        "product",
    ]
    from sprb_scanner.neighborhood.models import NeighborRecord

    fields = [
        *("cargo_" + f for f in identity),
        "status",
        "cargo_hit_models",
        "cargo_ga_models",
        "cargo_domain_evidence",
        *("neighbor_" + f for f in identity),
        *(f for f in NeighborRecord.model_fields if f != "gene"),
        "neighbor_domain_evidence",
    ]
    write_rows(stage / "associations.tsv", rows, fields)
    proteins = {p.protein_id: p for p in dataset.proteins}
    genes = {g.gene_id: g for g in dataset.genes}
    evidence_rows = []
    extra_fields = [
        "hit_id",
        "model",
        "ga",
        "sequence_evalue",
        "sequence_bit_score",
        "independent_domain_evalue",
        "conditional_domain_evalue",
        "domain_bit_score",
        "hmm_length",
        "hmm_start",
        "hmm_end",
        "alignment_start",
        "alignment_end",
        "envelope_start",
        "envelope_end",
        "hmm_coverage",
        "protein_coverage",
        "domain_start_fraction",
        "domain_end_fraction",
        "n_terminal_flank_aa",
        "c_terminal_flank_aa",
        "domain_order_by_protein_start",
        "domains_same_protein_model",
        "overlapping_same_model_hit_ids",
        "scan_directory",
        "raw_artifact",
        "raw_hit",
    ]
    for e in result.domains:
        h = e.hit
        p = proteins[h.protein_id]
        evidence_rows.append(
            gene_fields(genes[p.gene_id], p.protein_id, p.source_id, p.length, "")
            | {
                "hit_id": h.hit_id,
                "model": h.model.model_dump(),
                "ga": h.hmm.ga.model_dump(),
                "sequence_evalue": h.hmm.sequence_evalue,
                "sequence_bit_score": h.hmm.sequence_bit_score,
                "independent_domain_evalue": h.hmm.independent_domain_evalue,
                "conditional_domain_evalue": h.hmm.conditional_domain_evalue,
                "domain_bit_score": h.bit_score,
                "hmm_length": h.query_length,
                "hmm_start": h.query_start,
                "hmm_end": h.query_end,
                "alignment_start": h.start,
                "alignment_end": h.end,
                "envelope_start": h.hmm.envelope.start,
                "envelope_end": h.hmm.envelope.end,
                **{f: getattr(e, f) for f in extra_fields[15:-2]},
                "raw_artifact": h.raw_artifact,
                "raw_hit": h.model_dump(mode="json"),
            }
        )
    write_rows(stage / "domain_evidence.tsv", evidence_rows, identity + extra_fields)


def run_association(prepared: Path, bctd_scan: Path, sprf_scan: Path, output: Path, window=5):
    if output.exists():
        raise InputError(f"Output already exists: {output}")
    fingerprints = {}
    try:
        prepared_meta = read_manifest(prepared, fingerprints)
        validate_prepared(prepared)
        dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
        ctd, ctd_meta = load_pair(bctd_scan, dataset, CTD_MODELS, fingerprints)
        sprf, sprf_meta = load_pair(sprf_scan, dataset, SPRF_MODELS, fingerprints)
        result = associate(dataset, ctd, sprf, window)
        result = result.model_copy(
            update={
                "domains": tuple(
                    e.model_copy(
                        update={
                            "scan_directory": str(
                                (
                                    bctd_scan if e.hit.model.model_id in CTD_MODELS else sprf_scan
                                ).resolve()
                            )
                        }
                    )
                    for e in result.domains
                )
            }
        )
        check_unchanged(fingerprints)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise InputError(f"Cannot associate scan inputs: {exc}") from exc
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".association-", dir=output.parent))
    try:
        write_json(stage / "associations.json", result.model_dump(mode="json"))
        export_tables(stage, result, dataset)
        summary = {
            "cargo_count": len(result.cargos),
            "local_hit": sum(c.status == "local_hit" for c in result.cargos),
            "no_local_hit": sum(c.status == "no_local_hit" for c in result.cargos),
            "relationships": sum(len(c.neighbors) for c in result.cargos),
            "retained_domains": len(result.domains),
        }
        write_json(stage / "summary.json", summary)
        check_unchanged(fingerprints)
        write_json(
            stage / "run_metadata.json",
            {
                "software_version": __version__,
                "created_at": datetime.now(UTC).isoformat(),
                "genome_id": result.genome_id,
                "parameters": {
                    "window": window,
                    "candidate_rule": result.candidate_rule,
                    "coverage_cutoff": None,
                    "classification": "none",
                },
                "coordinates": result.coordinates,
                "geometry_policy": GEOMETRY_POLICY,
                "distance_rule": (
                    "nonnegative minimum part-edge gap; overlap separate; shortest circle"
                ),
                "input_files_sha256": fingerprints,
                "original_inputs": prepared_meta.get("inputs", []),
                "scan_directories": {
                    "bctd": str(bctd_scan.resolve()),
                    "sprf": str(sprf_scan.resolve()),
                },
                "scan_metadata": {"bctd": ctd_meta, "sprf": sprf_meta},
                "evidence_semantics": (
                    "all retained domains, including below GA; "
                    "coverage is per-domain alignment span; "
                    "no domain/coverage merging; no_local_hit is not an orphan/function prediction"
                ),
                "output_sha256": {p.name: sha256(p) for p in stage.iterdir()},
            },
        )
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return summary
