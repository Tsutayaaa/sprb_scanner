"""Annotated GenBank CDS parser preserving joined and fuzzy locations."""

from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

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


def parse_genbank(path: Path, genome_id: str, default_table=11) -> GenomeDataset:
    require_file(path, {".gbff", ".gbk", ".gb"})
    issues, contigs, genes, proteins, seen_proteins = [], [], [], [], set()
    with open_text(path) as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records or len({x.id for x in records}) != len(records):
        raise InputError("GenBank is empty or contains duplicate contig IDs")
    organisms = {x.annotations.get("organism") for x in records} - {None}
    if len(organisms) > 1:
        raise InputError("GenBank contains multiple organisms; provide one assembly per input")
    taxids = set()
    for rec in sorted(records, key=lambda x: x.id):
        cid = internal_id(genome_id, rec.id)
        contigs.append(
            ContigRecord(
                contig_id=cid,
                genome_id=genome_id,
                source_id=rec.id,
                length=len(rec),
                topology=rec.annotations.get("topology", "unknown"),
            )
        )
        for source in (f for f in rec.features if f.type == "source"):
            taxids.update(
                int(x.split(":")[1])
                for x in source.qualifiers.get("db_xref", [])
                if x.startswith("taxon:")
            )
        cds = [(i, f) for i, f in enumerate(rec.features) if f.type == "CDS"]
        if any(f.location is None for _, f in cds):
            raise InputError(f"Unparseable CDS location on {rec.id}")
        cds.sort(key=lambda pair: (int(pair[1].location.start), int(pair[1].location.end), pair[0]))
        for ordinal, (_, feature) in enumerate(cds, 1):
            gid = internal_id(genome_id, rec.id, f"CDS_{ordinal:06d}")
            loc = feature.location
            if getattr(loc, "operator", "join") != "join":
                raise InputError(f"Unsupported location operator: {gid}")
            if any(p.ref or p.ref_db for p in loc.parts):
                raise InputError(f"Remote CDS locations are unsupported: {gid}")
            if len({p.strand for p in loc.parts}) != 1:
                raise InputError(f"Mixed-strand CDS is unsupported: {gid}")
            attrs = {k: tuple(v) for k, v in feature.qualifiers.items()}
            attrs["source_location"] = (str(loc),)
            partial = "partial" in attrs or any(
                not isinstance(p.start, ExactPosition) or not isinstance(p.end, ExactPosition)
                for p in loc.parts
            )
            gene = GeneRecord(
                gene_id=gid,
                genome_id=genome_id,
                contig_id=cid,
                location_parts=tuple(
                    Interval(start=int(p.start), end=int(p.end)) for p in loc.parts
                ),
                strand=loc.strand or 0,
                source_id=first(attrs, "protein_id"),
                locus_tag=first(attrs, "locus_tag"),
                gene_name=first(attrs, "gene"),
                product=first(attrs, "product"),
                partial=partial,
                pseudogene="pseudo" in attrs or "pseudogene" in attrs,
                transl_table=int(first(attrs, "transl_table")) if "transl_table" in attrs else None,
                codon_start=int(first(attrs, "codon_start", "1")),
                qualifiers=attrs,
            )
            genes.append(gene)
            if partial:
                issue(issues, "partial_cds", "Partial CDS retained", gid)
            if gene.pseudogene:
                issue(issues, "pseudogene", "Pseudogene retained without scan protein", gid)
                continue
            translation = first(attrs, "translation")
            if translation is None:
                issue(issues, "missing_translation", "CDS retained; no translation inferred", gid)
                continue
            pid = gid + "::P1"
            seq = normalize_protein(translation, pid, issues)
            if seq is None:
                continue
            source_id = first(attrs, "protein_id")
            if source_id is not None and source_id in seen_proteins:
                issue(issues, "duplicate_source_protein_id", source_id, pid)
            seen_proteins.add(source_id)
            proteins.append(
                ProteinRecord(
                    protein_id=pid,
                    genome_id=genome_id,
                    gene_id=gid,
                    source_id=source_id,
                    sequence=seq,
                    description=gene.product or "",
                )
            )
            check_translation(gene, feature.extract(rec.seq), seq, issues, default_table)
    if not genes:
        raise InputError("No CDS annotations found; bare genomes are unsupported")
    if len(taxids) > 1:
        raise InputError("Conflicting taxonomy IDs in one assembly")
    return GenomeDataset(
        genome=GenomeRecord(
            genome_id=genome_id,
            accession=genome_id if genome_id.startswith(("GCF_", "GCA_")) else None,
            organism=next(iter(organisms), None),
            taxon_id=next(iter(taxids), None),
            input_format="genbank",
            source_files=(str(path.resolve()),),
        ),
        contigs=tuple(contigs),
        genes=tuple(genes),
        proteins=tuple(proteins),
        quality_issues=tuple(issues),
    )
