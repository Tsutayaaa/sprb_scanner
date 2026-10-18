import json

import pytest
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.exporters import validate_prepared
from sprb_scanner.io.prepare import prepare


def test_prepare_cli_and_integrity(genome_files, tmp_path):
    output = tmp_path / "prepared"
    runner = CliRunner()
    args = [
        "prepare",
        "--input",
        str(genome_files["gb"]),
        "--format",
        "genbank",
        "--output",
        str(output),
        "--genome-id",
        "g",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert validate_prepared(output)["proteins"] == 2
    assert (output / "proteome.faa").read_text().startswith(">g::chr1::CDS_000001::P1")
    assert json.loads((output / "run_metadata.json").read_text())["inputs"][0]["sha256"]
    assert runner.invoke(app, args).exit_code == 2  # never overwrite
    (output / "proteome.faa").write_text("corruption")
    assert runner.invoke(app, ["validate", "--prepared", str(output)]).exit_code == 2


def test_strict_quality_and_explicit_flagged_export(genome_files, tmp_path):
    f = genome_files
    f["faa"].write_text(">p1\nMM\n>p2\nMK\n")
    kwargs = dict(
        input_path=f["gff"], input_format="gff3", genome_id="g", fna=f["fna"], faa=f["faa"]
    )
    with pytest.raises(InputError, match="quality errors"):
        prepare(**kwargs, output=tmp_path / "strict", config=AppConfig())
    assert not (tmp_path / "strict").exists()
    config = AppConfig(input={"strict_quality": False})
    result = prepare(**kwargs, output=tmp_path / "flagged", config=config)
    assert result["errors"] == 1
    assert validate_prepared(tmp_path / "flagged")["errors"] == 1


def test_validate_reports_errors(genome_files):
    f = genome_files
    f["faa"].write_text(">p1\nMM\n>p2\nMK\n")
    result = CliRunner().invoke(
        app,
        [
            "validate",
            "--input",
            str(f["gff"]),
            "--format",
            "gff3",
            "--fna",
            str(f["fna"]),
            "--faa",
            str(f["faa"]),
        ],
    )
    assert result.exit_code == 2
    assert "translation_mismatch" in result.stdout
