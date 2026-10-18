"""Application boundary for offline input selection, normalization and export."""

import hashlib
from pathlib import Path

from pydantic import ValidationError

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import source_manifest
from sprb_scanner.io.exporters import export_prepared, summary
from sprb_scanner.io.fasta import parse_proteins
from sprb_scanner.io.genbank import parse_genbank
from sprb_scanner.io.gff import parse_gff


def prepare_dataset(
    input_path: Path, input_format: str, genome_id=None, fna=None, faa=None, default_table=11
):
    if input_format not in {"genbank", "gff3", "protein_fasta"}:
        raise InputError("Choose genbank, gff3 or protein_fasta explicitly")
    if input_format == "gff3" and (fna is None or faa is None):
        raise InputError("GFF3 requires both --fna and --faa from the same annotation release")
    if input_format != "gff3" and (fna is not None or faa is not None):
        raise InputError("--fna and --faa are only valid with GFF3 input")
    paths = [input_path] + ([fna, faa] if input_format == "gff3" else [])
    try:
        manifest = source_manifest(paths)
        if genome_id is None:
            digest = hashlib.sha256(
                (input_format + "|" + "|".join(x["sha256"] for x in manifest)).encode()
            ).hexdigest()
            genome_id = "local_" + digest
        if input_format == "genbank":
            dataset = parse_genbank(input_path, genome_id, default_table)
        elif input_format == "gff3":
            dataset = parse_gff(input_path, fna, faa, genome_id, default_table)
        else:
            dataset = parse_proteins(input_path, genome_id)
        if source_manifest(paths) != manifest:
            raise InputError("Inputs changed during parsing; retry from an immutable copy")
        return dataset, manifest
    except (OSError, EOFError, ValueError, ValidationError) as exc:
        raise InputError(f"Cannot prepare input: {exc}") from exc


def prepare(input_path, input_format, output, config, genome_id=None, fna=None, faa=None):
    dataset, manifest = prepare_dataset(
        input_path, input_format, genome_id, fna, faa, config.input.default_translation_table
    )
    report = summary(dataset)
    if report["errors"] and config.input.strict_quality:
        examples = "; ".join(
            f"{x.code}: {x.record_id}" for x in dataset.quality_issues if x.severity == "error"
        )
        raise InputError(
            f"{report['errors']} quality errors; no dataset exported. {examples[:1500]}. "
            "Use validate --input for the full report, or explicitly set "
            "input.strict_quality: false to export flagged evidence."
        )
    export_prepared(dataset, output, manifest, config)
    return report
