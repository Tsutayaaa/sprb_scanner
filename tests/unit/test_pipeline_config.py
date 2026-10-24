"""Strict versioned YAML, relative paths and source-state adapter contracts."""

from pathlib import Path

import pytest

from sprb_scanner.core.exceptions import ConfigurationError, InputError
from sprb_scanner.evidence.models import SourceRun
from sprb_scanner.evidence.service import run_merge
from sprb_scanner.pipeline import load_pipeline_config
from sprb_scanner.pipeline.storage import digest


def test_paths_are_relative_to_yaml_and_defaults_match_scanners(tmp_path):
    path = tmp_path / "pipeline.yaml"
    path.write_text(
        'schema_version: "1.0.0"\nregistry: refs/registry.json\n'
        "references:\n  TIGR04131:\n    artifact: models/test.hmm\n"
        "mmseqs:\n  executable: bin/mmseqs\n"
    )
    config = load_pipeline_config(path)
    assert config.registry == tmp_path / "refs/registry.json"
    assert config.references["TIGR04131"].artifact == tmp_path / "models/test.hmm"
    assert config.mmseqs.executable == str(tmp_path / "bin/mmseqs")
    assert config.mmseqs.evalue == 0.001 and config.mmseqs.sensitivity == 5.7
    assert config.hmm.sequence_report_evalue == config.hmm.domain_report_evalue == 10
    assert config.hmm.seed == 42


@pytest.mark.parametrize(
    "content",
    [
        "",
        "[]",
        "registry: refs.json",
        'schema_version: "2.0.0"\nregistry: refs.json',
        'schema_version: "1.0.0"\nregistry: refs.json\nunknown: true',
        'schema_version: "1.0.0"\nregistry: refs.json\nmode: other',
        'schema_version: "1.0.0"\nregistry: refs.json\nstages:\n  rn28: "false"',
        'schema_version: "1.0.0"\nregistry: refs.json\nhmm:\n  threads: true',
        'schema_version: "1.0.0"\nregistry: refs.json\nhmm:\n  max_sensitivity: 1',
        'schema_version: "1.0.0"\nregistry: refs.json\nmmseqs:\n  evalue: "0.001"',
        'schema_version: "1.0.0"\nregistry: refs.json\nmmseqs:\n  threshold: 1',
        'schema_version: "1.0.0"\nregistry: refs.json\nregistry: other.json',
        'schema_version: "1.0.0"\nregistry: refs.json\nreferences:\n  made_up: {}',
        'schema_version: "1.0.0"\nregistry: refs.json\nneighborhood:\n  window: -1',
        'schema_version: "1.0.0"\nregistry: refs.json\nstages:\n  porp_sprf: false',
    ],
)
def test_strict_configuration_rejects_ambiguous_or_unknown_settings(tmp_path, content):
    path = tmp_path / "pipeline.yaml"
    path.write_text(content)
    with pytest.raises(ConfigurationError):
        load_pipeline_config(path)


def test_fingerprint_is_canonical_and_sensitive_to_values():
    assert digest({"a": 1, "b": {"c": 2}}) == digest({"b": {"c": 2}, "a": 1})
    assert digest({"seed": 42}) != digest({"seed": 43})


@pytest.mark.parametrize(
    "override",
    [
        {"unexpected": SourceRun(source_id="unexpected", state="failed", message="failed")},
        {"RN28": SourceRun(source_id="RN28", state="succeeded", retained_hits=0)},
        {"RN28": SourceRun(source_id="PF13573", state="failed", message="failed")},
    ],
)
def test_merger_failure_adapter_cannot_inject_success_or_change_source_ids(tmp_path, override):
    with pytest.raises(InputError, match="Explicit source states"):
        run_merge(tmp_path / "unused", tmp_path / "out", source_states=override)


def test_committed_example_and_schema_are_current():
    import json

    from sprb_scanner.pipeline import PipelineConfig

    root = Path(__file__).parents[2]
    assert load_pipeline_config(root / "config/pipeline.yaml").schema_version == "1.0.0"
    assert (
        json.loads((root / "config/pipeline.schema.json").read_text())
        == PipelineConfig.model_json_schema()
    )


@pytest.mark.parametrize("resume", ["false", 1, None])
def test_python_api_rejects_nonboolean_resume(tmp_path, resume):
    from sprb_scanner.pipeline import run_pipeline

    with pytest.raises(InputError, match="boolean"):
        run_pipeline(tmp_path / "unused", tmp_path / "unused.yaml", tmp_path / "out", resume=resume)


def test_equivalent_reference_selection_has_same_effective_fingerprint(monkeypatch):
    from sprb_scanner.pipeline import PipelineConfig, identity
    from sprb_scanner.references import ReferenceRegistry

    root = Path(__file__).parents[2]
    registry = ReferenceRegistry.load(root / "references/registry.json")
    monkeypatch.setattr(identity, "backend_identity", lambda *args: {"version": "fixture"})
    first = PipelineConfig(schema_version="1.0.0", registry=registry.path)
    second = PipelineConfig.model_validate(
        {
            "schema_version": "1.0.0",
            "registry": registry.path,
            "references": {"RN28": {"version": registry.get("RN28").version}},
        }
    )
    assert identity.scanner_spec("rn28", first, registry, {}, {}).fingerprint == (
        identity.scanner_spec("rn28", second, registry, {}, {}).fingerprint
    )
