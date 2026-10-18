"""Shared input checks, provenance and CDS quality checks."""

import gzip
import hashlib
from pathlib import Path

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import QualityIssue


def open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8")
    )


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def source_manifest(paths):
    return [
        {"path": str(p.resolve()), "sha256": sha256(p), "bytes": p.stat().st_size} for p in paths
    ]


def first(attrs, key, default=None):
    values = attrs.get(key, ())
    return values[0] if values else default


def issue(issues, code, message, record_id=None, severity="warning"):
    issues.append(QualityIssue(code=code, message=message, record_id=record_id, severity=severity))


def normalize_protein(sequence, record_id, issues):
    normalized = sequence.upper()
    if sequence != normalized:
        issue(issues, "uppercase_normalization", "Protein converted to uppercase", record_id)
    if normalized.endswith("*"):
        normalized = normalized[:-1]
        issue(issues, "terminal_stop_removed", "One terminal stop symbol removed", record_id)
    if not normalized or set(normalized) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
        issue(
            issues,
            "invalid_protein",
            "Empty sequence or unsupported residue/gap/internal stop",
            record_id,
            "error",
        )
        return None
    if set(normalized) & set("BXZJ"):
        issue(issues, "ambiguous_residues", "Protein contains ambiguous residues", record_id)
    return normalized


def check_translation(gene, nucleotide, protein, issues, default_table=11):
    if (
        gene.pseudogene
        or gene.partial
        or gene.strand == 0
        or any(k in gene.qualifiers for k in ("transl_except", "exception", "ribosomal_slippage"))
    ):
        issue(
            issues,
            "translation_check_skipped",
            "Partial/pseudo/exception/unknown strand CDS",
            gene.gene_id,
        )
        return
    try:
        translated = str(
            Seq(str(nucleotide)[gene.codon_start - 1 :]).translate(
                table=gene.transl_table or default_table, cds=True
            )
        )
    except (TranslationError, ValueError) as exc:
        issue(issues, "invalid_cds_translation", str(exc), gene.gene_id, "error")
        return
    if translated != protein:
        issue(
            issues,
            "translation_mismatch",
            "Annotated protein differs from translated CDS",
            gene.gene_id,
            "error",
        )


def require_file(path: Path, suffixes):
    suffix = path.with_suffix("").suffix if path.suffix == ".gz" else path.suffix
    if suffix.lower() not in suffixes:
        raise InputError(
            f"Unsupported extension for {path.name}; expected {sorted(suffixes)} (+ .gz)"
        )
    if not path.is_file() or path.stat().st_size == 0:
        raise InputError(f"Missing or empty input: {path}")
