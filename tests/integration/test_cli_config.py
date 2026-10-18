import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import load_config
from sprb_scanner.core.exceptions import ConfigurationError

runner = CliRunner()


def test_help_version_and_default_config():
    assert runner.invoke(app, ["--help"]).exit_code == 0
    assert runner.invoke(app, ["version"]).stdout.strip() == "2.0.0a1"
    result = runner.invoke(app, ["validate"])
    assert result.exit_code == 0
    assert Path(json.loads(result.stdout)["data_dir"]).is_absolute()


def test_relative_paths(tmp_path):
    path = tmp_path / "settings.yaml"
    path.write_text("data_dir: input\nresults_dir: output\n")
    config = load_config(path)
    assert config.data_dir == tmp_path / "input"
    assert not config.data_dir.exists()


@pytest.mark.parametrize(
    "text",
    ["threads: 0", "typo: 1", "[]", "bctd: {report_evalue: .nan}", "threads: [", "threads: true"],
)
def test_invalid_yaml_and_config(tmp_path, text):
    path = tmp_path / "bad.yaml"
    path.write_text(text)
    with pytest.raises(ConfigurationError):
        load_config(path)
    assert runner.invoke(app, ["validate", "--config", str(path)]).exit_code == 2


def test_missing_config(tmp_path):
    with pytest.raises(ConfigurationError):
        load_config(tmp_path / "missing.yaml")


def test_example_configuration():
    assert load_config(Path("config/default.yaml")).threads == 1
