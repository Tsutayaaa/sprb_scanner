"""Structured coordinate relationships, without inferred biological function."""

from typing import Literal

from sprb_scanner.core.models import ContigRecord, GeneRecord, NonNegative, Record


class NeighborRecord(Record):
    gene: GeneRecord
    genomic_position: Literal["upstream", "downstream"]
    relative_position: Literal["upstream", "downstream"] | None
    genomic_cds_offset: int
    transcription_cds_offset: int | None
    intervening_cds: NonNegative
    distance_bp: int
    overlap_bp: NonNegative
    directly_adjacent: bool
    same_strand: bool | None
    overlaps: bool
    crosses_origin: bool
    index_crosses_origin: bool
    distance_crosses_origin: bool
    gene_crosses_origin: bool


class NeighborhoodResult(Record):
    schema_version: Literal["1.0.0"] = "1.0.0"
    genome_id: str
    contig: ContigRecord
    target: GeneRecord
    target_crosses_origin: bool
    window: NonNegative
    neighbors: tuple[NeighborRecord, ...] = ()
    coordinates: Literal["0-based-half-open"] = "0-based-half-open"
