"""GFF3 CDS parser with exact FASTA identifiers and explicit phase handling."""

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.locations import (
    GFF_LOCATION_POLICY,
    VERIFIABLE_GFF_EXCEPTIONS,
    location_digest,
    progressive_path,
)
from sprb_scanner.core.models import (
    ContigRecord,
    GeneRecord,
    GenomeDataset,
    GenomeRecord,
    Interval,
    ProteinRecord,
    internal_id,
)
from sprb_scanner.io.common import (
    check_translation,
    first,
    issue,
    normalize_protein,
    open_text,
    require_file,
)
from sprb_scanner.io.fasta import read_fasta


def attributes(text):
    result = {}
    if text == ".":
        return result
    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise InputError("GFF3 attributes must be key=value; GTF is unsupported")
        key, value = item.split("=", 1)
        key = unquote(key)
        if key in result:
            raise InputError(f"Repeated GFF attribute: {key}")
        result[key] = tuple(unquote(v) for v in value.split(","))
    return result


def _fasta_topology(record):
    # Only explicit sequence metadata; "complete genome" is not topology evidence.
    values = set(
        re.findall(r"(?:^|[\s\[])topology=(circular|linear)(?=[\s\]]|$)", record.description)
    )
    if len(values) > 1:
        raise InputError(f"Conflicting FASTA topology metadata: {record.id}")
    return next(iter(values), "unknown")


def _spans(start, end, length, circular, label):
    if start < 0 or start >= length or end <= start:
        raise InputError(f"Invalid GFF interval: {label}")
    if end <= length:
        return [(start, end)]
    if not circular or end - start > length:
        raise InputError(f"Out-of-bounds GFF feature without valid circular extent: {label}")
    return [(start, length), (0, end - length)]


def _joined_translation(gene, nucleotide, protein, default_table):
    if gene.partial or gene.strand == 0 or "transl_except" in gene.qualifiers:
        return (
            "not_verifiable",
            "Partial/unknown-strand/translation-exception CDS; FAA remains canonical",
        )
    try:
        translated = str(
            Seq(nucleotide[gene.codon_start - 1 :]).translate(
                table=gene.transl_table or default_table, cds=True
            )
        )
    except (TranslationError, ValueError) as exc:
        return "invalid_translation", str(exc)
    return (
        ("exact_faa_match", "Exact joined-CDS translation")
        if translated == protein
        else ("mismatch", "Joined CDS differs from the exact-identifier FAA record")
    )


def _location_validation(
    gene,
    rows,
    nucleotide,
    protein,
    length,
    circular,
    phase_ok,
    explicit_slippage,
    topology_source,
    default_table,
):
    """Validate accommodated joins; source intervals never become genomic unions here."""
    try:
        forward, crosses = progressive_path(gene.location_parts, gene.strand, length, circular)
    except InputError as exc:
        raise InputError(f"{exc}: {gene.source_id}") from exc
    ordered = sorted(gene.location_parts, key=lambda p: (p.start, p.end))
    overlaps = [
        a.end - b.start for a, b in zip(ordered, ordered[1:], strict=False) if a.end > b.start
    ]
    # Record all forward path gaps, unrolled once for an origin-spanning CDS.
    gaps = [
        (b.start - a.end) % length if b.start < a.start else b.start - a.end
        for a, b in zip(forward, forward[1:], strict=False)
    ]
    virtual = any(row[2] > length for row in rows)
    needs_proof = explicit_slippage or bool(overlaps) or virtual or crosses
    if not needs_proof:
        return None
    translation, reason = _joined_translation(gene, nucleotide, protein, default_table)
    representation = "circular_origin"
    if explicit_slippage:
        representation = "annotated_ribosomal_slippage"
        if translation not in {"exact_faa_match", "not_verifiable"}:
            raise InputError(
                "Annotated ribosomal slippage fails joined-CDS/FAA verification: "
                f"{gene.source_id}: {reason}"
            )
    elif overlaps:
        exceptions = {v.strip().casefold() for v in gene.qualifiers.get("exception", ())}
        representation = (
            "faa_verified_translation_discrepancy"
            if exceptions
            else "faa_verified_unannotated_overlap"
        )
        if (
            any(x > 2 for x in overlaps)
            or not phase_ok
            or translation != "exact_faa_match"
            or "protein_id" not in gene.qualifiers
            or not exceptions.issubset(VERIFIABLE_GFF_EXCEPTIONS)
            or any(k in gene.qualifiers for k in ("ribosomal_slippage", "transl_except"))
        ):
            raise InputError(
                "Unannotated overlapping CDS fails conservative structure/phase/FAA verification: "
                f"{gene.source_id}"
            )
    return {
        "policy": GFF_LOCATION_POLICY,
        "representation": representation,
        "biological_interpretation": "annotated_slippage" if explicit_slippage else "undetermined",
        "location_sha256": location_digest(gene.location_parts, gene.strand),
        "protein_mapping": "exact_identifier",
        "protein_id": gene.qualifiers.get("protein_id", (None,))[0],
        "protein_sha256": hashlib.sha256(protein.encode()).hexdigest(),
        "joined_translation": translation,
        "translation_detail": reason,
        "phase_consistent": phase_ok,
        "overlap_bp": overlaps,
        "join_gap_bp": gaps,
        "crosses_origin": crosses,
        "topology_source": topology_source,
        "source_coordinates": "1-based inclusive",
        "source_intervals": [
            {
                "seqid": r[0],
                "start": r[1] + 1,
                "end": r[2],
                "strand": r[3],
                "phase": r[4],
                "line": r[6],
            }
            for r in rows
        ],
    }


def parse_gff(path: Path, fna: Path, faa: Path, genome_id: str, default_table=11):
    require_file(path, {".gff", ".gff3"})
    dna = {rec.id: rec for rec in read_fasta(fna, nucleotide=True)}
    faa_records = read_fasta(faa)
    if len({x.id for x in faa_records}) != len(faa_records):
        raise InputError("Duplicate protein FASTA IDs make GFF mapping ambiguous")
    aa = {rec.id: rec for rec in faa_records}
    groups, features, all_rows = defaultdict(list), {}, []
    topology = {key: _fasta_topology(rec) for key, rec in dna.items()}
    topology_source = {
        key: "FASTA topology metadata" for key, value in topology.items() if value != "unknown"
    }
    issues, genes, proteins, used = [], [], [], set()
    with open_text(path) as handle:
        if handle.readline().strip() != "##gff-version 3":
            raise InputError("GFF input requires ##gff-version 3 header")
        for lineno, line in enumerate(handle, 2):
            if line.startswith("##FASTA"):
                raise InputError("Embedded GFF FASTA is unsupported; provide separate FNA and FAA")
            if line.startswith("##sequence-region"):
                _, seqid, start, end = line.split()
                seqid = unquote(seqid)
                if seqid not in dna or int(start) != 1 or int(end) != len(dna[seqid]):
                    raise InputError("GFF sequence-region does not match FNA")
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                raise InputError(f"GFF line {lineno}: expected 9 tab-separated columns")
            seqid, _, kind, a, b, _, strand, phase, raw = fields
            seqid = unquote(seqid)
            if seqid not in dna:
                raise InputError(f"GFF sequence ID missing from FNA: {seqid}")
            start, end = int(a) - 1, int(b)
            if start < 0 or end <= start or start >= len(dna[seqid]):
                raise InputError(f"Invalid GFF interval at line {lineno}")
            attrs = attributes(raw)
            fid = first(attrs, "ID")
            if "gff_location_validation" in attrs:
                raise InputError("Reserved parser provenance attribute: gff_location_validation")
            all_rows.append((seqid, start, end, lineno))
            if kind == "region":
                if start != 0 or end != len(dna[seqid]):
                    raise InputError(f"GFF region must span its sequence: {seqid}")
                flag = first(attrs, "Is_circular")
                if flag is not None:
                    if attrs["Is_circular"] not in {("true",), ("false",)}:
                        raise InputError(f"Invalid Is_circular metadata: {seqid}")
                    declared = "circular" if flag == "true" else "linear"
                    if topology[seqid] not in {"unknown", declared}:
                        raise InputError(f"Conflicting sequence circularity metadata: {seqid}")
                    topology[seqid] = declared
                    topology_source[seqid] = "GFF region Is_circular"
            if fid:
                if fid in features and (kind != "CDS" or features[fid][0] != "CDS"):
                    raise InputError(f"Duplicate non-CDS feature ID: {fid}")
                features[fid] = (kind, attrs, seqid, start, end)
            if kind != "CDS":
                continue
            if strand not in {"+", "-", ".", "?"} or phase not in {"0", "1", "2"}:
                raise InputError(f"CDS requires valid strand and phase 0/1/2 at line {lineno}")
            if not fid:
                raise InputError("Each GFF CDS requires an ID (repeated only for its parts)")
            groups[fid].append((seqid, start, end, strand, int(phase), attrs, lineno))
    # Region metadata may occur after a feature: validate extents after reading all rows.
    for seqid, start, end, lineno in all_rows:
        _spans(start, end, len(dna[seqid]), topology[seqid] == "circular", f"line {lineno}")
    if not groups:
        raise InputError("GFF has no CDS annotations")

    def inherited(attrs, seqid, visited=None):
        visited = set() if visited is None else visited
        parent = first(attrs, "Parent")
        if len(attrs.get("Parent", ())) > 1:
            raise InputError(
                "Multi-parent CDS/transcripts are unsupported in this prokaryotic contract"
            )
        if not parent:
            return dict(attrs)
        if parent in visited or parent not in features:
            raise InputError("Missing or cyclic GFF Parent relationship")
        if features[parent][2] != seqid:
            raise InputError("GFF Parent belongs to a different contig")
        merged = inherited(features[parent][1], seqid, visited | {parent})
        merged.pop("ID", None)
        merged.pop("Parent", None)
        return merged | attrs

    ordinal = defaultdict(int)
    sorted_groups = sorted(
        groups.items(),
        key=lambda x: (x[1][0][0], min(z[1] for z in x[1]), max(z[2] for z in x[1]), x[1][0][6]),
    )
    for fid, rows in sorted_groups:
        seqids, strands = {x[0] for x in rows}, {x[3] for x in rows}
        if len(seqids) != 1 or len(strands) != 1:
            raise InputError(f"CDS {fid} spans different contigs or strands")
        seqid, strand = rows[0][0], rows[0][3]
        if len(rows) > 1 and strand not in {"+", "-"}:
            raise InputError("Cannot order compound CDS with unknown strand")
        for key in ("protein_id", "Parent", "transl_table", "pseudo", "pseudogene"):
            if len({row[5].get(key) for row in rows}) > 1:
                raise InputError(f"Inconsistent {key} across CDS parts: {fid}")
        anchor = 0
        ancestor = first(rows[0][5], "Parent")
        visited = set()
        while ancestor and ancestor in features and ancestor not in visited:
            visited.add(ancestor)
            parent = features[ancestor]
            if parent[2] == seqid and parent[4] > len(dna[seqid]):
                anchor = parent[3]
            ancestor = first(parent[1], "Parent")
        rows.sort(key=lambda x: (x[1] - anchor) % len(dna[seqid]), reverse=strand == "-")
        # Explicit NCBI part=X/Y permits circular multipart ordering.
        if any("part" in row[5] for row in rows):
            try:
                parts = [tuple(map(int, first(row[5], "part", "").split("/"))) for row in rows]
                if sorted(parts) != [(i, len(rows)) for i in range(1, len(rows) + 1)]:
                    raise ValueError("Incomplete part numbering")
                rows = [row for _, row in sorted(zip(parts, rows, strict=True))]
            except (ValueError, TypeError) as exc:
                raise InputError(f"Invalid CDS part ordering for {fid}") from exc
        attrs = inherited(rows[0][5], seqid)
        for row in rows[1:]:
            for key in ("partial", "start_range", "end_range", "exception", "transl_except"):
                if key in row[5]:
                    attrs[key] = tuple(dict.fromkeys((*attrs.get(key, ()), *row[5][key])))
        ordinal[seqid] += 1
        gid = internal_id(genome_id, seqid, f"CDS_{ordinal[seqid]:06d}")
        locations, nucleotide, phases = [], [], []
        for row in rows:
            _, start, end, _, phase, _, _ = row
            length = len(dna[seqid])
            spans = _spans(start, end, length, topology[seqid] == "circular", fid)
            if strand == "-":
                spans.reverse()
            for a, b in spans:
                locations.append(Interval(start=a, end=b))
                piece = dna[seqid].seq[a:b]
                nucleotide.append(str(piece.reverse_complement() if strand == "-" else piece))
            phases.append(phase)
        # Structural validity is independent of exception labels and FAA identity.
        _, crosses_origin = progressive_path(
            tuple(locations),
            {"+": 1, "-": -1}.get(strand, 0),
            len(dna[seqid]),
            topology[seqid] == "circular",
        )
        explicit_slippage = "ribosomal slippage" in tuple(
            v.strip().casefold() for v in attrs.get("exception", ())
        )
        if "ribosomal_slippage" in attrs:
            if any(v.strip() for v in attrs["ribosomal_slippage"]):
                raise InputError(f"Invalid ribosomal_slippage flag: {fid}")
            explicit_slippage = True
        if explicit_slippage and any(
            v.strip().casefold() != "ribosomal slippage" for v in attrs.get("exception", ())
        ):
            raise InputError(f"Conflicting exceptions on ribosomal-slippage CDS: {fid}")
        if crosses_origin or any(r[2] > len(dna[seqid]) for r in rows):
            attrs["gff_source_intervals"] = tuple(
                json.dumps(
                    {
                        "seqid": r[0],
                        "start": r[1] + 1,
                        "end": r[2],
                        "strand": r[3],
                        "phase": r[4],
                        "line": r[6],
                    },
                    sort_keys=True,
                )
                for r in rows
            )
        attrs["gff_part_attributes"] = tuple(json.dumps(row[5], sort_keys=True) for row in rows)
        attrs["gff_phases"] = tuple(str(x) for x in phases)
        attrs["gff_source_rows"] = tuple(str(x[6]) for x in rows)
        partial = (
            first(attrs, "partial") == "true" or "start_range" in attrs or "end_range" in attrs
        )
        gene = GeneRecord(
            gene_id=gid,
            genome_id=genome_id,
            contig_id=internal_id(genome_id, seqid),
            source_id=fid,
            location_parts=tuple(locations),
            strand={"+": 1, "-": -1}.get(strand, 0),
            locus_tag=first(attrs, "locus_tag"),
            gene_name=first(attrs, "gene", first(attrs, "gene_name")),
            product=first(attrs, "product"),
            partial=partial,
            pseudogene=first(attrs, "pseudo") == "true" or "pseudogene" in attrs,
            transl_table=int(first(attrs, "transl_table")) if "transl_table" in attrs else None,
            codon_start=phases[0] + 1,
            qualifiers=attrs,
        )
        genes.append(gene)
        phase_ok = True
        coding_length = rows[0][2] - rows[0][1] - phases[0]
        for row in rows[1:]:
            if row[4] != (-coding_length) % 3:
                phase_ok = False
                issue(
                    issues,
                    "inconsistent_phase",
                    "CDS phases disagree with joined reading frame",
                    gid,
                    "warning" if explicit_slippage else "error",
                )
            coding_length += row[2] - row[1]
        if partial:
            issue(issues, "partial_cds", "Partial CDS retained", gid)
        # Matching is exact, with a documented NCBI cds- prefix fallback only.
        if "protein_id" in attrs and len(attrs["protein_id"]) != 1:
            raise InputError(f"Ambiguous protein_id attribute: {fid}")
        protein_key = first(attrs, "protein_id")
        if protein_key is None:
            aliases = [x for x in (fid, fid.removeprefix("cds-")) if x in aa]
            aliases = list(dict.fromkeys(aliases))
            if len(aliases) > 1:
                raise InputError(f"Ambiguous protein mapping for {fid}")
            protein_key = aliases[0] if aliases else None
        if gene.pseudogene:
            ordered = sorted(locations, key=lambda p: (p.start, p.end))
            if not explicit_slippage and any(
                a.end > b.start for a, b in zip(ordered, ordered[1:], strict=False)
            ):
                raise InputError(
                    f"Unannotated overlapping pseudogene cannot be FAA-verified: {fid}"
                )
            issue(issues, "pseudogene", "Pseudogene retained without scan protein", gid)
            if protein_key in aa:
                used.add(protein_key)
            continue
        if protein_key not in aa:
            issue(
                issues, "missing_translation", "No exact matching FAA record for CDS", gid, "error"
            )
            continue
        if protein_key in used:
            issue(issues, "duplicate_source_protein_id", protein_key, gid)
        used.add(protein_key)
        pid = gid + "::P1"
        seq = normalize_protein(str(aa[protein_key].seq), pid, issues)
        if seq is None:
            continue
        proof = _location_validation(
            gene,
            rows,
            "".join(nucleotide),
            seq,
            len(dna[seqid]),
            topology[seqid] == "circular",
            phase_ok,
            explicit_slippage,
            topology_source.get(seqid, "unknown"),
            default_table,
        )
        if proof is not None:
            proof["protein_id"] = protein_key
            # Preserve all parent intervals involved in virtual-space ordering as source evidence.
            proof["source_parent_intervals"] = [
                {"id": key, "start": features[key][3] + 1, "end": features[key][4]}
                for key in sorted(visited)
            ]
            gene = gene.model_copy(
                update={
                    "qualifiers": attrs
                    | {"gff_location_validation": (json.dumps(proof, sort_keys=True),)}
                }
            )
            genes[-1] = gene
            if proof["joined_translation"] == "not_verifiable":
                issue(issues, "gff_location_unverified", proof["translation_detail"], gid)
            if proof["representation"].startswith("faa_verified_"):
                issue(
                    issues,
                    "unannotated_join_uncertain",
                    "Exact FAA agreement supports representation; slippage remains unproven",
                    gid,
                )
        proteins.append(
            ProteinRecord(
                protein_id=pid,
                genome_id=genome_id,
                gene_id=gid,
                source_id=protein_key,
                sequence=seq,
                description=aa[protein_key].description,
            )
        )
        if not explicit_slippage:
            check_translation(gene, Seq("".join(nucleotide)), seq, issues, default_table)
    for key in sorted(aa.keys() - used):
        issue(issues, "unmapped_faa_record", f"FAA record has no CDS mapping: {key}", key, "error")
    return GenomeDataset(
        genome=GenomeRecord(
            genome_id=genome_id,
            accession=genome_id if genome_id.startswith(("GCF_", "GCA_")) else None,
            input_format="gff3",
            source_files=tuple(str(p.resolve()) for p in (path, fna, faa)),
        ),
        contigs=tuple(
            ContigRecord(
                contig_id=internal_id(genome_id, key),
                genome_id=genome_id,
                source_id=key,
                length=len(rec),
                topology=topology.get(key, "unknown"),
            )
            for key, rec in sorted(dna.items())
        ),
        genes=tuple(genes),
        proteins=tuple(proteins),
        quality_issues=tuple(issues),
    )
