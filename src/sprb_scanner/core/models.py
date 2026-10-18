"""Versioned records shared by parsers, scanners and future evidence consumers."""

from typing import Annotated, Literal
from urllib.parse import quote

from pydantic import BaseModel, ConfigDict, Field, computed_field, field_validator, model_validator

Identifier = Annotated[str, Field(min_length=1, pattern=r"^\S+$")]
NonNegative = Annotated[int, Field(ge=0, strict=True)]
Positive = Annotated[int, Field(gt=0, strict=True)]
Fraction = Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)]
SCHEMA_VERSION = "1.0.0"


class Record(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class Interval(Record):
    start: NonNegative
    end: Positive

    @model_validator(mode="after")
    def ordered(self):
        if self.end <= self.start:
            raise ValueError("end must exceed start")
        return self


class QualityIssue(Record):
    code: Identifier
    message: str
    record_id: str | None = None
    severity: Literal["warning", "error"] = "warning"


class GenomeRecord(Record):
    genome_id: Identifier
    accession: str | None = None
    organism: str | None = None
    taxon_id: Positive | None = None
    input_format: Literal["genbank", "gff3", "protein_fasta"]
    source_files: tuple[str, ...] = ()


class ContigRecord(Record):
    contig_id: Identifier
    genome_id: Identifier
    source_id: str
    length: Positive
    topology: Literal["linear", "circular", "unknown"] = "unknown"


class GeneRecord(Record):
    gene_id: Identifier
    genome_id: Identifier
    contig_id: Identifier
    location_parts: tuple[Interval, ...] = Field(min_length=1)
    strand: Literal[-1, 0, 1]
    source_id: str | None = None
    locus_tag: str | None = None
    gene_name: str | None = None
    product: str | None = None
    partial: bool = False
    pseudogene: bool = False
    transl_table: Positive | None = None
    codon_start: Literal[1, 2, 3] = 1
    qualifiers: dict[str, tuple[str, ...]] = Field(default_factory=dict)


class ProteinRecord(Record):
    protein_id: Identifier
    genome_id: Identifier
    gene_id: Identifier | None = None
    source_id: str | None = None
    sequence: str = Field(min_length=1)
    description: str = ""

    @field_validator("sequence")
    @classmethod
    def valid_sequence(cls, value):
        # Ambiguity, selenocysteine and pyrrolysine are retained. Parsers handle stops.
        if set(value) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError("Expected uppercase ungapped protein residues without stop symbols")
        return value

    @computed_field
    @property
    def length(self) -> int:
        return len(self.sequence)


class ModelProvenance(Record):
    model_id: Identifier
    version: str = Field(min_length=1)
    source: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[a-f0-9]{64}$")


class DomainHit(Interval):
    hit_id: Identifier
    protein_id: Identifier
    method: Literal["hmm", "mmseqs2"]
    model: ModelProvenance
    query_id: Identifier
    query_length: Positive
    query_start: NonNegative
    query_end: Positive
    protein_length: Positive
    evalue: float = Field(ge=0, allow_inf_nan=False)
    evalue_kind: Literal["independent_domain", "conditional_domain", "sequence", "alignment"]
    bit_score: float = Field(allow_inf_nan=False)
    query_coverage: Fraction
    target_coverage: Fraction
    backend_version: str = Field(min_length=1)
    raw_artifact: str | None = None

    @model_validator(mode="after")
    def bounds(self):
        if self.end > self.protein_length:
            raise ValueError("Hit exceeds protein length")
        if not self.query_start < self.query_end <= self.query_length:
            raise ValueError("Invalid query coordinates")
        return self

    @computed_field
    @property
    def ctd_distance(self) -> int:
        return self.protein_length - self.end


class ScanStatus(Record):
    model_id: Identifier
    state: Literal["not_run", "succeeded", "failed", "unavailable"]
    message: str | None = None

    @model_validator(mode="after")
    def explain_failure(self):
        if self.state in {"failed", "unavailable"} and not self.message:
            raise ValueError("Failed/unavailable scanners require an explanation")
        return self


class CandidateRecord(Record):
    candidate_id: Identifier
    protein_id: Identifier
    domain_hits: tuple[DomainHit, ...] = ()
    support: Literal["dual", "general_bctd", "sprb_specific", "unassigned"] = "unassigned"
    evidence: dict[str, tuple[str, ...]] = Field(default_factory=dict)

    @model_validator(mode="after")
    def matching_protein(self):
        if any(h.protein_id != self.protein_id for h in self.domain_hits):
            raise ValueError("Candidate contains evidence for another protein")
        if len({h.hit_id for h in self.domain_hits}) != len(self.domain_hits):
            raise ValueError("Duplicate hit IDs")
        return self


class GenomeDataset(Record):
    schema_version: Literal["1.0.0"] = SCHEMA_VERSION
    genome: GenomeRecord
    contigs: tuple[ContigRecord, ...] = ()
    genes: tuple[GeneRecord, ...] = ()
    proteins: tuple[ProteinRecord, ...] = ()
    quality_issues: tuple[QualityIssue, ...] = ()

    @model_validator(mode="after")
    def relationships(self):
        for records, key in (
            (self.contigs, "contig_id"),
            (self.genes, "gene_id"),
            (self.proteins, "protein_id"),
        ):
            if len({getattr(x, key) for x in records}) != len(records):
                raise ValueError(f"Duplicate internal {key}")
            if any(x.genome_id != self.genome.genome_id for x in records):
                raise ValueError("Record belongs to another genome")
        contigs = {x.contig_id: x for x in self.contigs}
        genes = {x.gene_id: x for x in self.genes}
        for gene in self.genes:
            if gene.contig_id not in contigs:
                raise ValueError("Gene references missing contig")
            if any(p.end > contigs[gene.contig_id].length for p in gene.location_parts):
                raise ValueError("Gene exceeds contig bounds")
        for protein in self.proteins:
            if protein.gene_id is not None and protein.gene_id not in genes:
                raise ValueError("Protein references missing gene")
            if self.genome.input_format != "protein_fasta" and protein.gene_id is None:
                raise ValueError("Annotated proteins require a gene link")
        if self.genome.input_format == "protein_fasta" and (self.contigs or self.genes):
            raise ValueError("Protein-only input cannot imply genomic coordinates")
        return self


def internal_id(*components: str) -> str:
    """Escape delimiters and whitespace so original IDs cannot collide."""
    if not components or any(not c for c in components):
        raise ValueError("ID components cannot be empty")
    return "::".join(quote(c, safe="._-") for c in components)
