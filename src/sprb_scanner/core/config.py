"""Strict YAML configuration with paths resolved relative to its file."""

from pathlib import Path

import yaml
from pydantic import Field, ValidationError, model_validator

from sprb_scanner.core.exceptions import ConfigurationError
from sprb_scanner.core.models import ModelProvenance, Positive, Record


class ProfileConfig(Record):
    provenance: ModelProvenance
    method: str = Field(pattern=r"^(hmm|mmseqs2)$")
    path: Path


class BCTDConfig(Record):
    # Exploratory reporting thresholds, not biologically validated classification cutoffs.
    report_evalue: float = Field(default=10.0, gt=0, allow_inf_nan=False)
    c_terminal_distance: int = Field(default=100, ge=0, strict=True)
    min_overlap_residues: Positive = 1
    profiles: tuple[ProfileConfig, ...] = ()

    @model_validator(mode="after")
    def unique_profiles(self):
        ids = [p.provenance.model_id for p in self.profiles]
        if len(ids) != len(set(ids)):
            raise ValueError("Duplicate model IDs")
        return self


class AppConfig(Record):
    schema_version: str = Field(default="1.0.0", pattern=r"^1\.0\.0$")
    data_dir: Path = Path("~/.local/share/sprb-scanner/data")
    results_dir: Path = Path("~/.local/share/sprb-scanner/results")
    threads: Positive = 1
    bctd: BCTDConfig = Field(default_factory=BCTDConfig)


def load_config(path: Path | None = None) -> AppConfig:
    try:
        raw = {} if path is None else yaml.safe_load(path.read_text(encoding="utf-8"))
        if raw is None:
            raw = {}
        if not isinstance(raw, dict):
            raise ValueError("Configuration must be a mapping")
        config = AppConfig.model_validate(raw)
        base = Path.cwd() if path is None else path.resolve().parent

        def resolve(value: Path) -> Path:
            expanded = value.expanduser()
            return (expanded if expanded.is_absolute() else base / expanded).resolve()

        profiles = tuple(
            p.model_copy(update={"path": resolve(p.path)}) for p in config.bctd.profiles
        )
        return config.model_copy(
            update={
                "data_dir": resolve(config.data_dir),
                "results_dir": resolve(config.results_dir),
                "bctd": config.bctd.model_copy(update={"profiles": profiles}),
            }
        )
    except (OSError, ValueError, ValidationError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"Cannot load configuration: {exc}") from exc
