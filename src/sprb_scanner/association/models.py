"""Objective local HMM evidence, without shuttle function predictions."""

from typing import Literal

from pydantic import computed_field

from sprb_scanner.core.models import DomainHit, GeneRecord, NonNegative, Positive, Record
from sprb_scanner.neighborhood.models import NeighborRecord


class DomainEvidence(Record):
    hit: DomainHit
    scan_directory: str | None = None
    domain_order_by_protein_start: Positive
    domains_same_protein_model: Positive
    overlapping_same_model_hit_ids: tuple[str, ...] = ()

    @computed_field
    @property
    def hmm_coverage(self) -> float:
        return self.hit.query_coverage

    @computed_field
    @property
    def protein_coverage(self) -> float:
        return self.hit.target_coverage

    @computed_field
    @property
    def domain_start_fraction(self) -> float:
        return self.hit.start / self.hit.protein_length

    @computed_field
    @property
    def domain_end_fraction(self) -> float:
        return self.hit.end / self.hit.protein_length

    @computed_field
    @property
    def n_terminal_flank_aa(self) -> int:
        return self.hit.start

    @computed_field
    @property
    def c_terminal_flank_aa(self) -> int:
        return self.hit.protein_length - self.hit.end


class LocalAssociation(Record):
    protein_id: str
    source_protein_id: str | None
    protein_length: Positive
    relationship: NeighborRecord
    domain_hit_ids: tuple[str, ...]


class CargoAssociation(Record):
    gene: GeneRecord
    protein_id: str
    source_protein_id: str | None
    protein_length: Positive
    domain_hit_ids: tuple[str, ...]
    hit_models: tuple[str, ...]
    ga_models: tuple[str, ...]
    status: Literal["local_hit", "no_local_hit"]
    neighbors: tuple[LocalAssociation, ...]


class AssociationResult(Record):
    schema_version: Literal["1.0.0"] = "1.0.0"
    genome_id: str
    window: NonNegative
    coordinates: Literal["0-based-half-open"] = "0-based-half-open"
    candidate_rule: Literal["any_retained_ctd_domain"] = "any_retained_ctd_domain"
    cargos: tuple[CargoAssociation, ...]
    domains: tuple[DomainEvidence, ...]
