"""Versioned strict local pipeline settings; scanner defaults remain unchanged."""

from pathlib import Path
from typing import Literal

import yaml
from pydantic import ConfigDict, Field, StrictBool, field_validator, model_validator

from sprb_scanner.core.exceptions import ConfigurationError
from sprb_scanner.core.models import NonNegative, Record
from sprb_scanner.evidence.models import SOURCES
from sprb_scanner.scanners.hmm import HMMSearchConfig
from sprb_scanner.scanners.mmseqs import MMseqsConfig


class StrictHMMConfig(HMMSearchConfig):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)


class StrictMMseqsConfig(MMseqsConfig):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)


class ReferenceBinding(Record):
    version: str | None = Field(default=None, min_length=1)
    artifact: Path | None = None


class EnabledStages(Record):
    typeb_ctd: StrictBool = True
    porp_sprf: StrictBool = True
    sprb_repeat: StrictBool = True
    rn28: StrictBool = True


class NeighborhoodConfig(Record):
    enabled: StrictBool = True
    window: NonNegative = 5


class PipelineConfig(Record):
    schema_version: Literal["1.0.0"]
    registry: Path
    mode: Literal["strict", "partial"] = "strict"
    stages: EnabledStages = Field(default_factory=EnabledStages)
    references: dict[str, ReferenceBinding] = Field(default_factory=dict)
    hmm: StrictHMMConfig = Field(default_factory=StrictHMMConfig)
    mmseqs: StrictMMseqsConfig = Field(default_factory=StrictMMseqsConfig)
    neighborhood: NeighborhoodConfig = Field(default_factory=NeighborhoodConfig)

    @field_validator("references")
    @classmethod
    def known_references(cls, value):
        if set(value) - set(SOURCES):
            raise ValueError("Only the five canonical HMM models and RN28 may be bound")
        return value

    @model_validator(mode="after")
    def compatible_stages(self):
        if self.neighborhood.enabled and not self.stages.porp_sprf:
            raise ValueError("Neighborhood requires porp_sprf; disable neighborhood explicitly")
        return self


class UniqueKeyLoader(yaml.SafeLoader):
    pass


def unique_mapping(loader, node, deep=False):
    loader.flatten_mapping(node)
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str) or key in result:
            raise ValueError(f"Non-string or duplicate YAML key: {key}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)


def resolve_paths(config, base):
    def path(value):
        expanded = value.expanduser()
        return (expanded if expanded.is_absolute() else base / expanded).resolve()

    references = {
        key: binding.model_copy(update={"artifact": path(binding.artifact)})
        if binding.artifact
        else binding
        for key, binding in config.references.items()
    }
    executable = config.mmseqs.executable
    if "/" in executable or "\\" in executable or executable.startswith("~"):
        executable = str(path(Path(executable)))
    return config.model_copy(
        update={
            "registry": path(config.registry),
            "references": references,
            "mmseqs": config.mmseqs.model_copy(update={"executable": executable}),
        }
    )


def load_pipeline_config(path):
    path = Path(path)
    try:
        raw = yaml.load(path.read_text(encoding="utf-8"), Loader=UniqueKeyLoader)
        if not isinstance(raw, dict):
            raise ValueError("Pipeline configuration must be a mapping")
        config = PipelineConfig.model_validate(raw)
        return resolve_paths(config, path.resolve().parent)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"Cannot load pipeline configuration {path}: {exc}") from exc
