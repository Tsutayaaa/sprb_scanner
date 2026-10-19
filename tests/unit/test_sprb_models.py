import hashlib
import json

import pyhmmer
import pytest

from sprb_scanner.core.config import ProfileConfig
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import ModelProvenance, ProteinRecord
from sprb_scanner.scanners.bctd import run_bctd
from sprb_scanner.scanners.repeats import run_repeats

SEQ = (
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEQHFKGLVLIAFSQYL"
    "QQCPFEDHVKLVNEVTEFAKTCVADESHAGCEKSLHTLFGDELCKVASL"
)


def model(tmp_path, model_id):
    abc = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(name=model_id, sequence=SEQ).digitize(abc)
    hmm, _, _ = pyhmmer.plan7.Builder(abc).build(seq, pyhmmer.plan7.Background(abc))
    hmm.accession = model_id + ".1"
    hmm.cutoffs.gathering = (0, 0)
    path = tmp_path / (model_id + ".hmm")
    with path.open("wb") as handle:
        hmm.write(handle)
    return ProfileConfig(
        method="hmm",
        path=path,
        provenance=ModelProvenance(
            model_id=model_id,
            version=model_id + ".1",
            source="synthetic fixture",
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        ),
    )


def proteins():
    return (ProteinRecord(protein_id="g::repeat", genome_id="g", sequence=SEQ + "G" * 60 + SEQ),)


def test_repeat_counts_preserve_all_domains_and_empty_result(tmp_path):
    profile = model(tmp_path, "PF13573")
    out = tmp_path / "repeats"
    result = run_repeats(proteins(), profile, out)
    hits = json.loads((out / "hits.json").read_text())
    counts = json.loads((out / "repeat_counts.json").read_text())
    assert len(hits) >= 2
    assert counts[0]["all_repeat_domains"] == result["domain_count"] == len(hits)
    assert counts[0]["passed_ga_domains"] == sum(
        h["hmm"]["ga"]["status"] == "passed_ga" for h in hits
    )
    assert len({h["hit_id"] for h in hits}) == len(hits)
    meta = json.loads((out / "run_metadata.json").read_text())
    for name, checksum in meta["output_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == checksum
    empty = run_repeats((), profile, tmp_path / "empty")
    assert empty["state"] == "succeeded" and empty["repeat_hit_proteins"] == 0
    with pytest.raises(InputError, match="PF13573"):
        run_repeats(proteins(), model(tmp_path, "PF13585"), tmp_path / "wrong")


def test_type_b_pair_uses_correct_models_and_raw_paths(tmp_path):
    profiles = [model(tmp_path, name) for name in ("TIGR04131", "PF13585")]
    out = tmp_path / "bctd"
    result = run_bctd(proteins(), profiles, out)
    assert result["state"] == "succeeded"
    assert result["both_models_ga_proteins"] == 1
    assert set(result["models"]) == {"TIGR04131", "PF13585"}
    for hit in json.loads((out / "hits.json").read_text()):
        assert (out / hit["raw_artifact"]).is_file()
        assert hit["protein_id"] == "g::repeat"
    with pytest.raises(InputError, match="exactly"):
        run_bctd(proteins(), [profiles[0], profiles[0]], tmp_path / "duplicate")
