"""Integrity-checked prepared input and atomic neighborhood export."""

import csv
import json
import shutil
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from sprb_scanner import __version__
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared, write_json
from sprb_scanner.neighborhood.analyzer import GEOMETRY_POLICY, NeighborhoodIndex
from sprb_scanner.neighborhood.models import NeighborRecord


def run_neighborhood(prepared: Path, gene_id: str, output: Path, window=5):
    if output.exists():
        raise InputError(f"Output already exists: {output}")
    before = {name: sha256(prepared / name) for name in ("dataset.json", "run_metadata.json")}
    validate_prepared(prepared)
    dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
    source_meta = json.loads((prepared / "run_metadata.json").read_text())
    if before != {name: sha256(prepared / name) for name in before}:
        raise InputError("Prepared input changed during loading")
    result = NeighborhoodIndex(dataset).query(gene_id, window)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".neighborhood-", dir=output.parent))
    try:
        write_json(stage / "neighborhood.json", result.model_dump(mode="json"))
        relation_fields = [name for name in NeighborRecord.model_fields if name != "gene"]
        fields = [
            "target_gene_id",
            "contig_id",
            "contig_source_id",
            "gene_id",
            "location_parts",
            "strand",
            "locus_tag",
            "gene_name",
            "product",
            "partial",
            "pseudogene",
            *relation_fields,
        ]
        with (stage / "neighborhood.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for neighbor in result.neighbors:
                gene = neighbor.gene
                row = {
                    key: getattr(gene, key)
                    for key in (
                        "gene_id",
                        "strand",
                        "locus_tag",
                        "gene_name",
                        "product",
                        "partial",
                        "pseudogene",
                    )
                }
                row.update(
                    target_gene_id=gene_id,
                    contig_id=result.contig.contig_id,
                    contig_source_id=result.contig.source_id,
                    location_parts=json.dumps([p.model_dump() for p in gene.location_parts]),
                    **{name: getattr(neighbor, name) for name in relation_fields},
                )
                writer.writerow(row)
        write_json(
            stage / "run_metadata.json",
            {
                "software_version": __version__,
                "created_at": datetime.now(UTC).isoformat(),
                "parameters": {"gene_id": gene_id, "window": window},
                "prepared_directory": str(prepared.resolve()),
                "prepared_sha256": before,
                "source_files": dataset.genome.source_files,
                "original_inputs": source_meta.get("inputs", []),
                "coordinate_convention": result.coordinates,
                "quality_issue_counts": {
                    severity: sum(q.severity == severity for q in dataset.quality_issues)
                    for severity in ("warning", "error")
                },
                "geometry_policy": GEOMETRY_POLICY,
                "distance_rule": (
                    "distance_bp=nonnegative minimum part-edge gap (0 for overlap); "
                    "overlap_bp=total coding intersection; "
                    "signed_gap_bp=distance_bp-overlap_bp; shortest circle gap"
                ),
                "output_sha256": {p.name: sha256(p) for p in stage.iterdir()},
            },
        )
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return result
