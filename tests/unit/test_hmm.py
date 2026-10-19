import hashlib
import json

import pyhmmer
import pytest
from pydantic import ValidationError

from sprb_scanner.core.config import ProfileConfig
from sprb_scanner.core.exceptions import HMMError
from sprb_scanner.core.models import DomainHit, ModelProvenance, ProteinRecord
from sprb_scanner.scanners.hmm import HMMSearchConfig, load_hmm, run_hmm, scan

SEQ = (
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEQHFKGLVLIAFSQYL"
    "QQCPFEDHVKLVNEVTEFAKTCVADESHAGCEKSLHTLFGDELCKVASL"
)


@pytest.fixture
def profile(tmp_path):
    abc = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(name="synthetic", sequence=SEQ).digitize(abc)
    hmm, _, _ = pyhmmer.plan7.Builder(abc).build(seq, pyhmmer.plan7.Background(abc))
    path = tmp_path / "model.hmm"
    with path.open("wb") as handle:
        hmm.write(handle)
    return ProfileConfig(
        method="hmm",
        path=path,
        provenance=ModelProvenance(
            model_id="synthetic",
            version="test",
            source="generated synthetic fixture",
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        ),
    )


def proteins():
    return (
        ProteinRecord(protein_id="g::tandem", genome_id="g", sequence=SEQ + "G" * 60 + SEQ),
        ProteinRecord(protein_id="g::negative", genome_id="g", sequence="P" * 100),
    )


def test_all_domains_match_native_backend(profile):
    config = HMMSearchConfig(max_sensitivity=True)
    hits, raw, _, meta = scan(proteins(), profile, config)
    abc = pyhmmer.easel.Alphabet.amino()
    targets = pyhmmer.easel.DigitalSequenceBlock(
        abc,
        [
            pyhmmer.easel.TextSequence(name=f"p{i}", sequence=p.sequence).digitize(abc)
            for i, p in enumerate(proteins())
        ],
    )
    native = next(
        pyhmmer.hmmsearch(
            load_hmm(profile),
            targets,
            cpus=1,
            E=10,
            domE=10,
            seed=42,
            F1=1,
            F2=1,
            F3=1,
            bias_filter=False,
        )
    )
    domains = [(t, d) for t in native for d in t.domains]
    assert len(hits) == len(raw) == len(domains)
    assert sum(h.protein_id == "g::tandem" for h in hits) >= 2
    assert meta["Z"] == len(proteins())
    for h, (t, d) in zip(hits, domains, strict=True):
        a = d.alignment
        assert h.protein_id == proteins()[int(t.name[1:])].protein_id
        assert (h.start, h.end) == (a.target_from - 1, a.target_to)
        assert (h.query_start, h.query_end) == (a.hmm_from - 1, a.hmm_to)
        assert h.bit_score == d.score
        assert h.evalue == h.hmm.independent_domain_evalue == d.i_evalue
        assert h.hmm.sequence_evalue == t.evalue
        assert h.hmm.conditional_domain_evalue == d.c_evalue
        assert h.hmm.sequence_bit_score == t.score
        assert h.model == profile.provenance
        assert h.query_coverage == (a.hmm_to - a.hmm_from + 1) / h.query_length
        assert h.ctd_distance == h.protein_length - h.end
    assert any(h.hmm.sequence_evalue != h.evalue for h in hits)


def test_unreported_domains_are_retained(profile):
    hits, raw, _, _ = scan(
        proteins(),
        profile,
        HMMSearchConfig(domain_report_evalue=1e-300, max_sensitivity=True),
    )
    assert hits and len(hits) == len(raw)
    assert any(not h.hmm.domain_reported for h in hits)


@pytest.mark.parametrize(
    "kind",
    [
        "missing",
        "directory",
        "empty",
        "garbage",
        "truncated",
        "multiple",
        "checksum",
        "identity",
        "dna",
    ],
)
def test_invalid_models(profile, tmp_path, kind):
    path = profile.path
    if kind == "missing":
        path = tmp_path / "missing.hmm"
    elif kind == "directory":
        path = tmp_path
    elif kind == "empty":
        path.write_bytes(b"")
    elif kind == "garbage":
        path.write_text("not an HMM")
    elif kind == "truncated":
        path.write_bytes(path.read_bytes()[:300])
    elif kind == "multiple":
        path.write_bytes(path.read_bytes() * 2)
    elif kind == "dna":
        abc = pyhmmer.easel.Alphabet.dna()
        seq = pyhmmer.easel.TextSequence(name="synthetic", sequence="ACGT" * 20).digitize(abc)
        hmm, _, _ = pyhmmer.plan7.Builder(abc).build(seq, pyhmmer.plan7.Background(abc))
        with path.open("wb") as handle:
            hmm.write(handle)
    prov = profile.provenance.model_dump()
    if path.is_file() and kind != "checksum":
        prov["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    if kind == "checksum":
        prov["sha256"] = "0" * 64
    if kind == "identity":
        prov["model_id"] = "wrong"
    invalid = ProfileConfig(method="hmm", path=path, provenance=ModelProvenance(**prov))
    with pytest.raises(HMMError):
        load_hmm(invalid)


def test_empty_and_no_hit_success(profile, tmp_path):
    for index, seqs in enumerate(((), (proteins()[1],))):
        out = tmp_path / str(index)
        result = run_hmm(seqs, profile, out)
        assert result["state"] == "succeeded" and result["domain_count"] == 0
        assert json.loads((out / "hits.json").read_text()) == []
        assert len((out / "hits.tsv").read_text().splitlines()) == 1
        with pytest.raises(HMMError, match="already exists"):
            run_hmm(seqs, profile, out)


def test_failure_status_and_duplicate_ids(profile, tmp_path, monkeypatch):
    with pytest.raises(HMMError, match="Duplicate"):
        scan(proteins() * 2, profile, HMMSearchConfig())

    def broken(*args, **kwargs):
        raise RuntimeError("backend failure")

    monkeypatch.setattr(pyhmmer, "hmmsearch", broken)
    out = tmp_path / "failure"
    with pytest.raises(HMMError, match="backend failure"):
        run_hmm(proteins(), profile, out)
    assert json.loads((out / "status.json").read_text())["state"] == "failed"
    assert not (out / "hits.json").exists()


def test_export_and_statistics_validation(profile, tmp_path):
    out = tmp_path / "export"
    run_hmm(proteins(), profile, out)
    payload = json.loads((out / "hits.json").read_text())[0]
    payload.pop("ctd_distance")
    hit = DomainHit.model_validate(payload)
    assert hit.hmm is not None
    payload["evalue"] = 123.0
    with pytest.raises(ValidationError):
        DomainHit.model_validate(payload)
    meta = json.loads((out / "run_metadata.json").read_text())
    for name, checksum in meta["output_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == checksum
