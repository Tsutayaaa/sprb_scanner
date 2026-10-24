"""Shared location structure and coordinate-bound provenance, without sequence inference."""

import hashlib
import json
import re

from sprb_scanner.core.exceptions import InputError

VERIFIABLE_GFF_EXCEPTIONS = {
    "unclassified translation discrepancy",
    "annotated by transcript or proteomic data",
}

GFF_LOCATION_POLICY = "ncbi-gff3-location-v1"


def location_digest(parts, strand):
    payload = {"strand": strand, "parts": [[p.start, p.end] for p in parts]}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def progressive_path(parts, strand, length, circular):
    """Validate a single forward genomic traversal; preserve biological part order."""
    forward = tuple(reversed(parts)) if strand == -1 else tuple(parts)
    offset, previous_start, previous_end, drops = 0, forward[0].start, forward[0].end, 0
    for previous, part in zip(forward, forward[1:], strict=False):
        if part.start < previous.start:
            drops += 1
            offset += length
        start, end = part.start + offset, part.end + offset
        if start <= previous_start or end <= previous_end:
            raise InputError("Non-progressive/duplicated CDS location parts")
        previous_start, previous_end = start, end
    if drops and not circular:
        raise InputError("Origin-spanning or unordered CDS on non-circular contig")
    if drops > 1:
        raise InputError("CDS makes multiple origin traversals")
    if drops and previous_end - forward[0].start > length:
        raise InputError("CDS span exceeds one circle")
    return forward, bool(drops)


def faa_verified_overlap(gene):
    """Recognize a parser proof bound to these exact parts; never infer slippage."""
    values = gene.qualifiers.get("gff_location_validation", ())
    if len(values) != 1:
        return False
    try:
        proof = json.loads(values[0])
        exceptions = {v.strip().casefold() for v in gene.qualifiers.get("exception", ())}
        expected = (
            "faa_verified_translation_discrepancy"
            if exceptions
            else "faa_verified_unannotated_overlap"
        )
        ordered = sorted(gene.location_parts, key=lambda p: (p.start, p.end))
        return (
            proof["policy"] == GFF_LOCATION_POLICY
            and proof["representation"] == expected
            and exceptions.issubset(VERIFIABLE_GFF_EXCEPTIONS)
            and not any(k in gene.qualifiers for k in ("ribosomal_slippage", "transl_except"))
            and proof["location_sha256"] == location_digest(gene.location_parts, gene.strand)
            and proof["phase_consistent"] is True
            and proof["joined_translation"] == "exact_faa_match"
            and proof["protein_mapping"] == "exact_identifier"
            and re.fullmatch(r"[0-9a-f]{64}", proof["protein_sha256"]) is not None
            and proof["protein_id"] in gene.qualifiers.get("protein_id", ())
            and proof["overlap_bp"]
            == [
                a.end - b.start
                for a, b in zip(
                    ordered,
                    ordered[1:],
                    strict=False,
                )
                if a.end > b.start
            ]
            and all(0 < x <= 2 for x in proof["overlap_bp"])
            and bool(proof["overlap_bp"])
            and not gene.partial
            and not gene.pseudogene
            and gene.strand in {-1, 1}
        )
    except (ValueError, KeyError, TypeError):
        return False
