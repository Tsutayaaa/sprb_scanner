import csv
import hashlib
import json

import pyhmmer
import pytest
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import AppConfig, ProfileConfig, load_profile
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset, GenomeRecord, ModelProvenance, ProteinRecord
from sprb_scanner.io.exporters import export_prepared
from sprb_scanner.scanners.hmm import HMMSearchConfig
from sprb_scanner.scanners.sprf import run_sprf

SEQ = (
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEQHFKGLVLIAFSQYL"
    "QQCPFEDHVKLVNEVTEFAKTCVADESHAGCEKSLHTLFGDELCKVASL"
)


@pytest.fixture
def inputs(tmp_path):
    abc = pyhmmer.easel.Alphabet.amino()
    profiles = []
    for name in ("TIGR03519", "PF11751"):
        seq = pyhmmer.easel.TextSequence(name=name, sequence=SEQ).digitize(abc)
        hmm, _, _ = pyhmmer.plan7.Builder(abc).build(seq, pyhmmer.plan7.Background(abc))
        hmm.accession = name + ".1"
        hmm.cutoffs.gathering = (0, 0) if name == "TIGR03519" else (10000, 10000)
        path = tmp_path / (name + ".hmm")
        with path.open("wb") as handle:
            hmm.write(handle)
        profiles.append(
            ProfileConfig(
                method="hmm",
                path=path,
                provenance=ModelProvenance(
                    model_id=name,
                    version=name + ".1",
                    source="synthetic test fixture",
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                ),
            )
        )
    proteins = (
        ProteinRecord(
            protein_id="g::tandem",
            genome_id="g",
            source_id="original",
            sequence=SEQ + "G" * 60 + SEQ,
        ),
    )
    return proteins, profiles


def test_dual_models_preserve_independent_evidence(inputs, tmp_path):
    proteins, profiles = inputs
    out = tmp_path / "dual"
    result = run_sprf(proteins, profiles, out, HMMSearchConfig(max_sensitivity=True))
    assert result["state"] == "succeeded"
    assert result["both_models_hit_proteins"] == 1
    assert result["both_models_ga_proteins"] == 0
    assert result["any_model_ga_proteins"] == 1
    evidence = json.loads((out / "protein_evidence.json").read_text())
    assert evidence[0]["ga_by_model"] == {"TIGR03519": "passed_ga", "PF11751": "below_ga"}
    assert evidence[0]["source_id"] == "original"
    hits = json.loads((out / "hits.json").read_text())
    assert result["domain_count"] == len(hits) >= 4
    assert len({h["hit_id"] for h in hits}) == len(hits)
    for h in hits:
        assert (out / h["raw_artifact"]).is_file()
        assert h["evalue"] == h["hmm"]["independent_domain_evalue"]
        assert h["protein_id"] == proteins[0].protein_id
    for profile in profiles:
        model = profile.provenance.model_id
        native = json.loads((out / model / "hits.json").read_text())
        merged = [h for h in hits if h["model"]["model_id"] == model]
        for a, b in zip(native, merged, strict=True):
            a["raw_artifact"] = f"{model}/raw_domains.json"
            assert a == b
    metadata = json.loads((out / "run_metadata.json").read_text())
    for name, digest in metadata["output_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    with pytest.raises(InputError, match="already exists"):
        run_sprf(proteins, profiles, out)


@pytest.mark.parametrize("failed_models", [1, 2])
def test_model_failure_preserves_other_result(inputs, tmp_path, failed_models):
    proteins, profiles = inputs
    for profile in profiles[:failed_models]:
        profile.path.unlink()
    out = tmp_path / "failed"
    result = run_sprf(proteins, profiles, out)
    assert result["state"] == ("partial" if failed_models == 1 else "failed")
    assert result["models"]["TIGR03519"]["state"] == "failed"
    if failed_models == 1:
        assert result["domain_count"] >= 2
        row = json.loads((out / "protein_evidence.json").read_text())[0]
        assert row["support"] == "PF11751_observed_other_failed"
        assert row["model_states"]["TIGR03519"] == "failed"
    else:
        assert json.loads((out / "hits.json").read_text()) == []


def test_empty_success(inputs, tmp_path):
    _, profiles = inputs
    result = run_sprf((), profiles, tmp_path / "empty")
    assert result["state"] == "succeeded" and result["domain_count"] == 0
    with (tmp_path / "empty/protein_evidence.tsv").open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []


def test_profile_contract(inputs, tmp_path):
    proteins, profiles = inputs
    with pytest.raises(InputError, match="exactly"):
        run_sprf(proteins, [profiles[0], profiles[0]], tmp_path / "duplicate")
    with pytest.raises(InputError, match="Duplicate internal"):
        run_sprf(proteins * 2, profiles, tmp_path / "duplicate_proteins")


def test_cli_relative_manifest_and_partial_exit(inputs, tmp_path):
    proteins, profiles = inputs
    prepared = tmp_path / "prepared"
    export_prepared(
        GenomeDataset(
            genome=GenomeRecord(genome_id="g", input_format="protein_fasta"), proteins=proteins
        ),
        prepared,
        [],
        AppConfig(),
    )
    manifests = []
    for profile in profiles:
        path = tmp_path / (profile.provenance.model_id + ".json")
        data = profile.model_dump(mode="json")
        data["path"] = profile.path.name
        path.write_text(json.dumps(data))
        manifests.append(path)
        assert load_profile(path).path == profile.path
    args = [
        "sprf-search",
        "--prepared",
        str(prepared),
        "--tigr-profile",
        str(manifests[0]),
        "--pfam-profile",
        str(manifests[1]),
        "--output",
        str(tmp_path / "cli"),
    ]
    runner = CliRunner()
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    profiles[0].path.unlink()
    args[-1] = str(tmp_path / "partial_cli")
    result = runner.invoke(app, args)
    assert result.exit_code == 2, result.output
    assert json.loads(result.output)["state"] == "partial"
