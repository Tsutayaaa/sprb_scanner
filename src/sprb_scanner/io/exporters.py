"""Atomic prepared-dataset export and integrity verification (JSON/TSV/FASTA)."""

import csv
import json
import shutil
import tempfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from sprb_scanner import __version__
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.io.common import sha256


def summary(dataset):
    return {
        "genome_id": dataset.genome.genome_id,
        "contigs": len(dataset.contigs),
        "genes": len(dataset.genes),
        "proteins": len(dataset.proteins),
        "quality_issue_counts": dict(Counter(x.code for x in dataset.quality_issues)),
        "errors": sum(x.severity == "error" for x in dataset.quality_issues),
        "coordinates": "0-based-half-open",
    }


def write_json(path, value):
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )


def write_tsv(path, records, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for record in records:
            data = record.model_dump(mode="json")
            writer.writerow(
                {
                    k: json.dumps(data[k], ensure_ascii=False)
                    if isinstance(data[k], (list, dict))
                    else data[k]
                    for k in fields
                }
            )


def export_prepared(dataset, output: Path, input_manifest, config):
    if output.exists():
        raise InputError(f"Output already exists; choose a new directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".prepare-", dir=output.parent))
    try:
        write_json(
            stage / "dataset.json", dataset.model_dump(mode="json", exclude_computed_fields=True)
        )
        write_json(stage / "genome_metadata.json", dataset.genome.model_dump(mode="json"))
        write_json(stage / "quality_report.json", [x.model_dump() for x in dataset.quality_issues])
        write_json(stage / "summary.json", summary(dataset))
        write_tsv(
            stage / "genes.tsv",
            dataset.genes,
            [
                "gene_id",
                "genome_id",
                "contig_id",
                "source_id",
                "locus_tag",
                "gene_name",
                "product",
                "location_parts",
                "strand",
                "partial",
                "pseudogene",
                "codon_start",
                "transl_table",
            ],
        )
        write_tsv(
            stage / "contigs.tsv",
            dataset.contigs,
            ["contig_id", "genome_id", "source_id", "length", "topology"],
        )
        write_tsv(
            stage / "proteins.tsv",
            dataset.proteins,
            ["protein_id", "genome_id", "gene_id", "source_id", "length", "description"],
        )
        with (stage / "proteome.faa").open("w", encoding="utf-8") as handle:
            for p in dataset.proteins:
                handle.write(f">{p.protein_id}\n")
                handle.write(
                    "\n".join(p.sequence[i : i + 80] for i in range(0, p.length, 80)) + "\n"
                )
        write_json(
            stage / "run_metadata.json",
            {
                "software_version": __version__,
                "schema_version": dataset.schema_version,
                "created_at": datetime.now(UTC).isoformat(),
                "coordinate_convention": "0-based-half-open",
                "inputs": input_manifest,
                "configuration": config.model_dump(mode="json"),
                "output_sha256": {p.name: sha256(p) for p in sorted(stage.iterdir())},
            },
        )
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return output


def validate_prepared(path: Path):
    meta = json.loads((path / "run_metadata.json").read_text())
    expected = {
        "dataset.json",
        "genome_metadata.json",
        "quality_report.json",
        "summary.json",
        "genes.tsv",
        "contigs.tsv",
        "proteins.tsv",
        "proteome.faa",
    }
    if set(meta.get("output_sha256", {})) != expected:
        raise InputError("Prepared output manifest is incomplete or unsupported")
    for name, checksum in meta["output_sha256"].items():
        if not (path / name).is_file() or sha256(path / name) != checksum:
            raise InputError(f"Output missing or checksum mismatch: {name}")
    dataset = GenomeDataset.model_validate_json((path / "dataset.json").read_text())
    return summary(dataset)
