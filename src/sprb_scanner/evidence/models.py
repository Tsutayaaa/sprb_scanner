"""Lossless per-protein evidence; data availability is separate from biological support."""

from typing import Literal

from pydantic import Field, model_validator

from sprb_scanner.association.models import DomainEvidence, LocalAssociation
from sprb_scanner.core.models import GeneRecord, NonNegative, Positive, Record
from sprb_scanner.scanners.mmseqs import SequenceAlignment

HMM_MODELS = ("TIGR04131", "PF13585", "PF13573", "TIGR03519", "PF11751")
SOURCES = (*HMM_MODELS, "RN28")


class SourceRun(Record):
    source_id: str
    state: Literal["not_run", "failed", "succeeded"]
    directory: str | None = None
    message: str | None = None
    retained_hits: NonNegative | None = None
    metadata: dict = Field(default_factory=dict)

    @model_validator(mode="after")
    def availability(self):
        if (self.state == "succeeded") != (self.retained_hits is not None):
            raise ValueError("Only successful sources have a known retained hit count")
        if self.state == "failed" and not self.message:
            raise ValueError("Failed sources require an explanation")
        return self


class ProteinEvidence(Record):
    source_id: str
    state: Literal["not_run", "failed", "succeeded_no_hits", "succeeded_with_hits"]
    domains: tuple[DomainEvidence, ...] = ()
    alignments: tuple[SequenceAlignment, ...] = ()

    @model_validator(mode="after")
    def consistent(self):
        if bool(self.domains or self.alignments) != (self.state == "succeeded_with_hits"):
            raise ValueError("Protein hit presence disagrees with evidence state")
        if self.source_id == "RN28":
            if self.domains:
                raise ValueError("RN28 evidence must contain sequence alignments")
        elif self.alignments or any(e.hit.model.model_id != self.source_id for e in self.domains):
            raise ValueError("Domain evidence belongs to another source")
        return self


class NeighborhoodEvidence(Record):
    state: Literal["not_run", "failed", "unavailable", "partial", "succeeded"]
    outcome: Literal["local_hit", "no_local_hit", "unknown"]
    window: NonNegative
    message: str | None = None
    sprf_sources: tuple[str, ...] = ()
    neighbors: tuple[LocalAssociation, ...] = ()

    @model_validator(mode="after")
    def missing_is_not_negative(self):
        if self.outcome == "no_local_hit" and self.state != "succeeded":
            raise ValueError("Incomplete neighborhood evidence cannot imply no_local_hit")
        if bool(self.neighbors) != (self.outcome == "local_hit"):
            raise ValueError("Local hit outcome disagrees with neighbor evidence")
        if self.state in {"not_run", "failed", "unavailable"} and self.outcome != "unknown":
            raise ValueError("Unavailable neighborhood evidence must remain unknown")
        return self


class EvidenceProfile(Record):
    protein_id: str
    source_protein_id: str | None
    protein_length: Positive
    protein_sequence_sha256: str
    description: str
    gene: GeneRecord | None
    sources: dict[str, ProteinEvidence]
    neighborhood: NeighborhoodEvidence

    @model_validator(mode="after")
    def consistent_identity(self):
        if set(self.sources) != set(SOURCES):
            raise ValueError("Evidence profiles require all six source states")
        for key, evidence in self.sources.items():
            if key != evidence.source_id:
                raise ValueError("Evidence source key disagrees with source ID")
            hits = tuple(e.hit for e in evidence.domains) + evidence.alignments
            if any(
                h.protein_id != self.protein_id or h.protein_length != self.protein_length
                for h in hits
            ):
                raise ValueError("Profile contains another protein's evidence")
        if not any(e.state == "succeeded_with_hits" for e in self.sources.values()):
            raise ValueError("Candidate union requires at least one retained hit")
        return self


class EvidenceResult(Record):
    schema_version: Literal["1.0.0"] = "1.0.0"
    genome_id: str
    coordinates: Literal["0-based-half-open"] = "0-based-half-open"
    candidate_rule: Literal["any_retained_scanner_hit"] = "any_retained_scanner_hit"
    runs: dict[str, SourceRun]
    neighborhood_state: Literal["not_run", "failed", "unavailable", "partial", "succeeded"]
    neighborhood_message: str | None = None
    profiles: tuple[EvidenceProfile, ...]
