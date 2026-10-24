"""Versioned, strict batch identities and explicit site resource requirements."""

from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

Token = Annotated[str, Field(pattern=r"^[A-Za-z0-9_.@/+:-]+$")]
SHA = Annotated[str, Field(pattern=r"^[0-9a-f]{64}$")]


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)


class SlurmResources(StrictModel):
    partition: Token
    account: Token
    qos: Token
    cpus: Annotated[int, Field(gt=0)]
    memory: Annotated[str, Field(pattern=r"^[1-9][0-9]*[MG]$")]
    time: Annotated[str, Field(pattern=r"^(?:[0-9]+-)?[0-9]{1,3}:[0-5][0-9]:[0-5][0-9]$")]
    concurrency: Annotated[int, Field(gt=0)]
    max_array_tasks: Annotated[int, Field(gt=0)]
    max_array_size: Annotated[int, Field(gt=0)]
    python: str = Field(min_length=1)
    modules: list[Token] = Field(default_factory=list)

    @model_validator(mode="after")
    def absolute_executable(self):
        if not Path(self.python).is_absolute() or any(c in self.python for c in "\n\r\0"):
            raise ValueError("SLURM python must be an absolute executable path")
        if self.concurrency > self.max_array_tasks:
            raise ValueError("concurrency exceeds max_array_tasks")
        return self


class InputFile(StrictModel):
    path: str
    size_bytes: Annotated[int, Field(ge=0)]
    sha256: SHA


class BatchEntry(StrictModel):
    index: Annotated[int, Field(ge=0)]
    assembly_accession: Annotated[str, Field(pattern=r"^GC[AF]_[0-9]+\.[0-9]+$")]
    genome_id: str
    files: dict[str, InputFile]
    annotation_identity: dict = Field(default_factory=dict)
    inventory_state: Literal["ready", "invalid"]
    error: str | None = None
    preparation_status: Literal["not_run"] = "not_run"
    prepared_path: str

    @model_validator(mode="after")
    def identity(self):
        if self.genome_id != self.assembly_accession:
            raise ValueError("Genome ID must equal versioned assembly accession")
        if self.inventory_state == "ready" and set(self.files) != {"gff", "fna", "faa"}:
            raise ValueError("Ready entry requires GFF/FNA/FAA")
        if self.inventory_state == "invalid" and not self.error:
            raise ValueError("Invalid inventory requires diagnostics")
        return self


class BatchManifest(StrictModel):
    schema_version: Literal["1.0.0"] = "1.0.0"
    created_at: str
    output_root: str
    dataset_root: str
    software: dict
    pipeline_config: str
    pinned_files: dict[str, SHA]
    resources: SlurmResources
    entries: list[BatchEntry]

    @model_validator(mode="after")
    def deterministic_mapping(self):
        ids = [e.assembly_accession for e in self.entries]
        if not ids or ids != sorted(set(ids)):
            raise ValueError("Assemblies must be unique, nonempty and sorted")
        if [e.index for e in self.entries] != list(range(len(ids))):
            raise ValueError("Array indexes must be contiguous and zero based")
        if len(ids) > self.resources.max_array_size:
            raise ValueError("Assembly indexes exceed configured SLURM MaxArraySize")
        root = Path(self.output_root)
        if not root.is_absolute():
            raise ValueError("Batch root must be absolute")
        for entry in self.entries:
            if Path(entry.prepared_path) != root / "genomes" / entry.genome_id / "prepared":
                raise ValueError("Prepared path does not match isolated genome directory")
        return self
