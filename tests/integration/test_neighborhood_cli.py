import csv
import hashlib
import json
from pathlib import Path

from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig
from sprb_scanner.io.exporters import export_prepared
from sprb_scanner.io.genbank import parse_genbank

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/fj_window.gbff"


def prepare_window(tmp_path):
    dataset = parse_genbank(FIXTURE, "FJ")
    prepared = tmp_path / "prepared"
    export_prepared(dataset, prepared, [], AppConfig())
    sprb = next(g for g in dataset.genes if "Fjoh_0979" in g.qualifiers.get("old_locus_tag", ()))
    sprf = next(g for g in dataset.genes if "Fjoh_0978" in g.qualifiers.get("old_locus_tag", ()))
    return prepared, sprb, sprf


def test_cli_sprb_sprf_json_tsv_metadata_and_empty(tmp_path):
    prepared, sprb, sprf = prepare_window(tmp_path)
    output = tmp_path / "neighborhood"
    args = [
        "neighborhood",
        "--prepared",
        str(prepared),
        "--gene",
        sprb.gene_id,
        "--window",
        "5",
        "--output",
        str(output),
    ]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    assert "43 bp" in result.output and "downstream" in result.output
    data = json.loads((output / "neighborhood.json").read_text())
    assert len(data["neighbors"]) == 1
    neighbor = data["neighbors"][0]
    assert neighbor["gene"]["gene_id"] == sprf.gene_id
    assert neighbor["gene"]["strand"] == data["target"]["strand"] == -1
    assert neighbor["relative_position"] == "downstream"
    assert neighbor["genomic_position"] == "upstream"
    assert neighbor["distance_bp"] == neighbor["signed_gap_bp"] == 43
    assert neighbor["overlap_bp"] == neighbor["intervening_cds"] == 0
    assert neighbor["directly_adjacent"] and neighbor["same_strand"]
    with (output / "neighborhood.tsv").open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert json.loads(row["location_parts"]) == neighbor["gene"]["location_parts"]
    assert int(row["distance_bp"]) == int(row["signed_gap_bp"]) == 43
    assert int(row["overlap_bp"]) == 0
    meta = json.loads((output / "run_metadata.json").read_text())
    assert "distance_bp=nonnegative" in meta["distance_rule"]
    assert "signed_gap_bp=distance_bp-overlap_bp" in meta["distance_rule"]
    for name, digest in meta["output_sha256"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    args[args.index("5")] = "0"
    args[-1] = str(tmp_path / "empty")
    assert CliRunner().invoke(app, args).exit_code == 0
    with (tmp_path / "empty/neighborhood.tsv").open() as handle:
        assert not list(csv.DictReader(handle, delimiter="\t"))


def test_cli_invalid_target_window_input_and_no_overwrite(tmp_path):
    prepared, sprb, _ = prepare_window(tmp_path)
    output = tmp_path / "out"
    args = [
        "neighborhood",
        "--prepared",
        str(prepared),
        "--gene",
        "missing",
        "--output",
        str(output),
    ]
    runner = CliRunner()
    assert runner.invoke(app, args).exit_code == 2
    assert not output.exists()
    args[args.index("missing")] = sprb.gene_id
    assert runner.invoke(app, args + ["--window", "-1"]).exit_code == 2
    assert runner.invoke(app, args).exit_code == 0
    before = (output / "neighborhood.json").read_bytes()
    assert runner.invoke(app, args).exit_code == 2
    assert (output / "neighborhood.json").read_bytes() == before
    (prepared / "dataset.json").write_text("{}")
    args[-1] = str(tmp_path / "tampered")
    assert runner.invoke(app, args).exit_code == 2
    assert not (tmp_path / "tampered").exists()
