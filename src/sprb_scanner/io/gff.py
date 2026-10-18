"""GFF3 CDS parser with exact FASTA identifiers and explicit phase handling."""

import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio.Seq import Seq

from sprb_scanner.core.exceptions import InputError
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


def parse_gff(path: Path, fna: Path, faa: Path, genome_id: str, default_table=11):
    require_file(path, {".gff", ".gff3"})
    dna = {rec.id: rec for rec in read_fasta(fna, nucleotide=True)}
    faa_records = read_fasta(faa)
    if len({x.id for x in faa_records}) != len(faa_records):
        raise InputError("Duplicate protein FASTA IDs make GFF mapping ambiguous")
    aa = {rec.id: rec for rec in faa_records}
    groups, features, topology, regions = defaultdict(list), {}, {}, {}
    issues, genes, proteins, used = [], [], [], set()
    with open_text(path) as handle:
        if handle.readline().strip() != "##gff-version 3":
            raise InputError("GFF input requires ##gff-version 3 header")
        for lineno, line in enumerate(handle, 2):
            if line.startswith("##FASTA"):
                raise InputError("Embedded GFF FASTA is unsupported; provide separate FNA and FAA")
            if line.startswith("##sequence-region"):
                _, seqid, start, end = line.split()
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
            if kind == "region":
                topology[seqid] = "circular" if first(attrs, "Is_circular") == "true" else "unknown"
                regions[seqid] = attrs
            if fid:
                if fid in features and (kind != "CDS" or features[fid][0] != "CDS"):
                    raise InputError(f"Duplicate non-CDS feature ID: {fid}")
                features[fid] = (kind, attrs, seqid)
            if kind != "CDS":
                if end > len(dna[seqid]) and kind not in {"gene", "mRNA"}:
                    raise InputError(f"Out-of-bounds GFF feature at line {lineno}")
                continue
            if strand not in {"+", "-", ".", "?"} or phase not in {"0", "1", "2"}:
                raise InputError(f"CDS requires valid strand and phase 0/1/2 at line {lineno}")
            if not fid:
                raise InputError("Each GFF CDS requires an ID (repeated only for its parts)")
            groups[fid].append((seqid, start, end, strand, int(phase), attrs, lineno))
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
        rows.sort(key=lambda x: x[1], reverse=strand == "-")
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
                    attrs[key] = row[5][key]
        ordinal[seqid] += 1
        gid = internal_id(genome_id, seqid, f"CDS_{ordinal[seqid]:06d}")
        locations, nucleotide, phases = [], [], []
        for row in rows:
            _, start, end, _, phase, _, _ = row
            length = len(dna[seqid])
            if end > length:
                if topology.get(seqid) != "circular" or end > start + length:
                    raise InputError(f"Out-of-bounds CDS on non-circular contig: {fid}")
                spans = [(start, length), (0, end - length)]
                if strand == "-":
                    spans.reverse()
            else:
                spans = [(start, end)]
            for a, b in spans:
                locations.append(Interval(start=a, end=b))
                piece = dna[seqid].seq[a:b]
                nucleotide.append(str(piece.reverse_complement() if strand == "-" else piece))
            phases.append(phase)
        flat = sorted((p.start, p.end) for p in locations)
        if any(a[1] > b[0] for a, b in zip(flat, flat[1:], strict=False)):
            raise InputError(f"Overlapping CDS parts require unsupported exception handling: {fid}")
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
        coding_length = rows[0][2] - rows[0][1] - phases[0]
        for row in rows[1:]:
            if row[4] != (-coding_length) % 3:
                issue(
                    issues,
                    "inconsistent_phase",
                    "CDS phases disagree with joined reading frame",
                    gid,
                    "error",
                )
            coding_length += row[2] - row[1]
        if partial:
            issue(issues, "partial_cds", "Partial CDS retained", gid)
        # Matching is exact, with a documented NCBI cds- prefix fallback only.
        protein_key = first(attrs, "protein_id")
        if protein_key is None:
            aliases = [x for x in (fid, fid.removeprefix("cds-")) if x in aa]
            aliases = list(dict.fromkeys(aliases))
            if len(aliases) > 1:
                raise InputError(f"Ambiguous protein mapping for {fid}")
            protein_key = aliases[0] if aliases else None
        if gene.pseudogene:
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
