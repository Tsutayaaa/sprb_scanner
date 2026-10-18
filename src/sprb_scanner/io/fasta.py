"""Strict FASTA reading and protein-only dataset normalization."""

from pathlib import Path

from Bio import SeqIO

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset, GenomeRecord, ProteinRecord, internal_id
from sprb_scanner.io.common import issue, normalize_protein, open_text, require_file

FASTA_EXTENSIONS = {".faa", ".fasta", ".fa", ".fas"}
NUCLEOTIDE_EXTENSIONS = {".fna", ".fasta", ".fa", ".fas"}


def read_fasta(path: Path, nucleotide=False):
    require_file(path, NUCLEOTIDE_EXTENSIONS if nucleotide else FASTA_EXTENSIONS)
    with open_text(path) as handle:
        if not handle.readline().startswith(">"):
            raise InputError(f"FASTA must start with a header: {path}")
        handle.seek(0)
        records = list(SeqIO.parse(handle, "fasta"))
    if not records or any(not x.id or not x.seq for x in records):
        raise InputError(f"FASTA contains no records or an empty record: {path}")
    if nucleotide:
        if len({x.id for x in records}) != len(records):
            raise InputError("Duplicate nucleotide FASTA IDs are ambiguous")
        if any(set(str(x.seq).upper()) - set("ACGTRYSWKMBDHVN") for x in records):
            raise InputError("Nucleotide FASTA must contain ungapped IUPAC DNA")
    return records


def parse_proteins(path: Path, genome_id: str) -> GenomeDataset:
    issues, proteins, seen = [], [], set()
    for index, rec in enumerate(read_fasta(path), 1):
        pid = internal_id(genome_id, f"PROTEIN_{index:06d}")
        if rec.id in seen:
            issue(issues, "duplicate_source_protein_id", rec.id, pid)
        seen.add(rec.id)
        seq = normalize_protein(str(rec.seq), pid, issues)
        if seq is not None:
            proteins.append(
                ProteinRecord(
                    protein_id=pid,
                    genome_id=genome_id,
                    source_id=rec.id,
                    sequence=seq,
                    description=rec.description,
                )
            )
    return GenomeDataset(
        genome=GenomeRecord(
            genome_id=genome_id, input_format="protein_fasta", source_files=(str(path.resolve()),)
        ),
        proteins=tuple(proteins),
        quality_issues=tuple(issues),
    )
