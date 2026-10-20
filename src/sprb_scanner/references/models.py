"""Versioned, checksum-pinned resource identities and derivation DAGs."""

from collections import defaultdict, deque
from pathlib import PurePosixPath
from typing import Annotated, Literal

from pydantic import Field, field_validator, model_validator

from sprb_scanner.core.models import Identifier, Interval, Positive, Record

Checksum = Annotated[str, Field(pattern=r"^[a-f0-9]{64}$")]
Text = Annotated[str, Field(min_length=1, pattern=r"\S")]
ResourceType = Literal["sequence", "sequence_set", "msa", "profile_hmm", "annotation_catalog"]


class Artifact(Record):
    path: str
    format: Literal["protein_fasta", "aligned_fasta", "stockholm", "a3m", "hmm", "json", "tsv"]
    storage: Literal["bundled", "external"] = "bundled"

    @field_validator("path")
    @classmethod
    def portable_file_path(cls, value):
        path = PurePosixPath(value)
        if (
            not value
            or value.startswith("~")
            or path.is_absolute()
            or ".." in path.parts
            or ":" in value
            or "\\" in value
            or "\x00" in value
            or value != path.as_posix()
            or value == "."
        ):
            raise ValueError("artifact.path must be a canonical relative POSIX file path")
        return value


class Source(Record):
    provider: Text
    locator: Text
    release: Text | None = None
    notes: Text | None = None


class BiologicalTarget(Record):
    description: Text
    organism: Text | None = None
    protein_accession: Identifier | None = None
    protein_length: Positive | None = None
    protein_interval: Interval | None = None
    sequence_ids: tuple[Identifier, ...] = ()
    sequence_sha256: Checksum | None = None
    model_id: Identifier | None = None
    model_accession: Identifier | None = None

    @model_validator(mode="after")
    def bounds_and_ids(self):
        if len(set(self.sequence_ids)) != len(self.sequence_ids):
            raise ValueError("Duplicate biological_target sequence IDs")
        if (
            self.protein_interval is not None
            and self.protein_length is not None
            and self.protein_interval.end > self.protein_length
        ):
            raise ValueError("Biological target interval exceeds parent protein length")
        return self


class Dependency(Record):
    reference_id: Identifier
    version: Identifier
    sha256: Checksum

    @property
    def key(self):
        return self.reference_id, self.version


class Reference(Dependency):
    resource_type: ResourceType
    artifact: Artifact
    source: Source
    biological_target: BiologicalTarget
    derived_from: tuple[Dependency, ...] = ()

    @model_validator(mode="after")
    def type_and_format(self):
        formats = {
            "sequence": {"protein_fasta"},
            "sequence_set": {"protein_fasta"},
            "msa": {"aligned_fasta", "stockholm", "a3m"},
            "profile_hmm": {"hmm"},
            "annotation_catalog": {"json", "tsv"},
        }
        if self.artifact.format not in formats[self.resource_type]:
            raise ValueError(f"Invalid artifact format for {self.resource_type}")
        if len({p.key for p in self.derived_from}) != len(self.derived_from):
            raise ValueError("Duplicate derived_from dependency")
        if self.biological_target.sequence_sha256 and self.resource_type != "sequence":
            raise ValueError("sequence_sha256 describes a single sequence resource")
        return self


class RegistryDocument(Record):
    schema_version: Literal["1.0.0"] = "1.0.0"
    references: tuple[Reference, ...] = Field(min_length=1)

    @model_validator(mode="after")
    def dependency_graph(self):
        by_key = {entry.key: entry for entry in self.references}
        if len(by_key) != len(self.references):
            raise ValueError("Duplicate (reference_id, version) identity")
        indegree, children = {}, defaultdict(list)
        for entry in self.references:
            indegree[entry.key] = len(entry.derived_from)
            for parent in entry.derived_from:
                if parent.key not in by_key:
                    raise ValueError(f"Missing dependency {parent.reference_id}@{parent.version}")
                if parent.sha256 != by_key[parent.key].sha256:
                    raise ValueError(
                        f"Dependency checksum disagrees with registered parent: {parent.key}"
                    )
                children[parent.key].append(entry.key)
        ready = deque(key for key, count in indegree.items() if count == 0)
        visited = 0
        while ready:
            key = ready.popleft()
            visited += 1
            for child in children[key]:
                indegree[child] -= 1
                if indegree[child] == 0:
                    ready.append(child)
        if visited != len(by_key):
            raise ValueError("Cycle in reference derivation graph")
        return self
