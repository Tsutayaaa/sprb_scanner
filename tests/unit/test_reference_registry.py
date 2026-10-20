"""Registry format parsing, pinned lineage, precise failure states and legacy adapters."""

import hashlib
import json

import pyhmmer
import pytest
from typer.testing import CliRunner

from sprb_scanner.cli import app
from sprb_scanner.core.config import load_profile
from sprb_scanner.core.exceptions import ReferenceRegistryError, ReferenceUnavailableError
from sprb_scanner.references import ReferenceRegistry
from sprb_scanner.scanners.hmm import load_hmm

SEQ = (
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEQHFKGLVLIAFSQYL"
    "QQCPFEDHVKLVNEVTEFAKTCVADESHAGCEKSLHTLFGDELCKVASL"
)


def entry(
    tmp_path,
    rid="q",
    kind="sequence",
    fmt="protein_fasta",
    content=">original|id\nACDEFGHIKL\n",
    version="v1",
):
    path = tmp_path / f"{rid}-{version}.artifact"
    path.write_text(content)
    return dict(
        reference_id=rid,
        resource_type=kind,
        version=version,
        artifact=dict(path=path.name, format=fmt, storage="bundled"),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source=dict(provider="synthetic fixture", locator="local:generated"),
        biological_target=dict(description="test protein"),
        derived_from=[],
    )


def registry(tmp_path, entries):
    path = tmp_path / "registry.json"
    path.write_text(json.dumps(dict(schema_version="1.0.0", references=entries)))
    return ReferenceRegistry.load(path)


def dependency(rec):
    return {k: rec[k] for k in ("reference_id", "version", "sha256")}


def test_sequence_records_relative_paths_and_whole_file_sha(tmp_path, monkeypatch):
    rec = entry(tmp_path)
    resources = registry(tmp_path, [rec])
    monkeypatch.chdir(tmp_path.parent)
    assert resources.registry_sha256 == hashlib.sha256(resources.path.read_bytes()).hexdigest()
    resolved = resources.resolve("q")
    proteins = resources.sequences("q")
    assert resolved.artifact_path == tmp_path / rec["artifact"]["path"]
    assert proteins[0].protein_id == proteins[0].source_id == "original|id"
    assert proteins[0].sequence == "ACDEFGHIKL" and resolved.summary["lengths"] == [10]
    with pytest.raises(ReferenceRegistryError, match="Only profile_hmm"):
        resolved.as_profile()
    path = resolved.artifact_path
    path.write_text(path.read_text().replace("ACDE", "ACDF"))
    with pytest.raises(ReferenceRegistryError, match="SHA-256 mismatch"):
        resources.resolve("q")


@pytest.mark.parametrize(
    "kind,fmt,content,columns",
    [
        ("sequence_set", "protein_fasta", ">one\nACD\n>two\nEFGH\n", None),
        ("msa", "aligned_fasta", ">one\nAC-D\n>two\nACED\n", 4),
        ("msa", "a3m", ">one\nACqq-D\n>two\nACED\n", 4),
        ("msa", "stockholm", "# STOCKHOLM 1.0\none AC-D\ntwo ACED\n//\n", 4),
        ("annotation_catalog", "json", '{"raw": [1, 2]}', None),
        ("annotation_catalog", "tsv", "id\tlabel\na\traw\n", None),
    ],
)
def test_all_resource_formats(tmp_path, kind, fmt, content, columns):
    resources = registry(tmp_path, [entry(tmp_path, kind=kind, fmt=fmt, content=content)])
    summary = resources.resolve("q").summary
    if columns:
        assert summary["match_columns"] == columns
        assert summary["sequence_ids"] == ("one", "two")
    if kind == "sequence_set":
        assert len(resources.sequences("q")) == 2
    elif kind != "sequence":
        with pytest.raises(ReferenceRegistryError, match="ProteinRecord"):
            resources.sequences("q")


@pytest.mark.parametrize(
    "kind,fmt,content",
    [
        ("sequence", "protein_fasta", ">a\nACD\n>b\nACD\n"),
        ("sequence_set", "protein_fasta", ">a\nACD\n>a\nACD\n"),
        ("sequence", "protein_fasta", ">a\nAC-D\n"),
        ("sequence", "protein_fasta", ">a\nacde\n"),
        ("sequence", "protein_fasta", ">a\nACD*EF\n"),
        ("sequence", "protein_fasta", "not fasta"),
        ("msa", "aligned_fasta", ">a\nAC-D\n>b\nACD\n"),
        ("msa", "aligned_fasta", ">a\n----\n>b\nACDE\n"),
        ("msa", "a3m", ">a\nACqq-D\n>b\nACD\n"),
        ("msa", "stockholm", "broken"),
        ("annotation_catalog", "json", '{"a": 1, "a": 2}'),
        ("annotation_catalog", "json", "[NaN]"),
        ("annotation_catalog", "json", "true"),
        ("annotation_catalog", "tsv", "id\tid\na\tb\n"),
        ("annotation_catalog", "tsv", "id\tlabel\na\n"),
        ("profile_hmm", "hmm", "not a valid HMM"),
    ],
)
def test_type_failures_are_explicit(tmp_path, kind, fmt, content):
    resources = registry(tmp_path, [entry(tmp_path, kind=kind, fmt=fmt, content=content)])
    with pytest.raises(ReferenceRegistryError, match="invalid artifact"):
        resources.resolve("q")
    assert resources.validate()["state"] == "failed"


@pytest.mark.parametrize(
    "change",
    [
        {"resource_type": "unknown"},
        {"version": ""},
        {"sha256": "missing"},
        {"source": {"provider": " ", "locator": "source"}},
        {"resource_type": "msa"},
        {
            "biological_target": {
                "description": "t",
                "protein_length": 10,
                "protein_interval": {"start": 5, "end": 11},
            }
        },
    ],
)
def test_schema_failures(tmp_path, change):
    rec = entry(tmp_path) | change
    with pytest.raises(ReferenceRegistryError, match="Cannot load"):
        registry(tmp_path, [rec])


@pytest.mark.parametrize(
    "path",
    [
        "../outside",
        "/absolute",
        "https://example.org/model",
        "a/../b",
        "a//b",
        "./a",
        "~/.model",
        "C:\\model",
    ],
)
def test_artifact_paths_are_portable(tmp_path, path):
    rec = entry(tmp_path)
    rec["artifact"]["path"] = path
    with pytest.raises(ReferenceRegistryError, match="relative POSIX"):
        registry(tmp_path, [rec])


def test_missing_external_is_registered_but_never_downloaded(tmp_path):
    rec = entry(tmp_path)
    (tmp_path / rec["artifact"]["path"]).unlink()
    resources = registry(tmp_path, [rec])
    assert resources.validate()["state"] == "failed"  # Missing bundled artifact is never allowed.
    rec["artifact"]["storage"] = "external"
    resources = registry(tmp_path, [rec])
    with pytest.raises(ReferenceUnavailableError, match="no automatic download"):
        resources.resolve("q")
    assert resources.validate()["resources"][0]["state"] == "missing_external"
    assert resources.validate()["state"] == "partial"
    assert resources.validate(allow_missing_external=False)["state"] == "failed"
    external = tmp_path / "external.faa"
    external.write_text(">original|id\nACDEFGHIKL\n")
    assert resources.resolve("q", artifact=external).artifact_path == external
    assert resources.validate(artifact_overrides={("q", "v1"): external})["state"] == "validated"
    with pytest.raises(ReferenceRegistryError, match="unknown reference"):
        resources.validate(artifact_overrides={("absent", "v1"): external})


def test_bundled_symlink_escape_and_directory_artifact(tmp_path):
    inside = tmp_path / "registry"
    inside.mkdir()
    rec = entry(inside)
    external = tmp_path / "outside.faa"
    external.write_bytes((inside / rec["artifact"]["path"]).read_bytes())
    path = inside / rec["artifact"]["path"]
    path.unlink()
    path.symlink_to(external)
    resources = registry(inside, [rec])
    with pytest.raises(ReferenceRegistryError, match="escapes"):
        resources.resolve("q")
    with pytest.raises(ReferenceRegistryError, match="override requires external"):
        resources.resolve("q", artifact=external)
    path.unlink()
    path.mkdir()
    with pytest.raises(ReferenceRegistryError, match="not a file"):
        resources.resolve("q")


def test_lineage_version_pins_and_ambiguous_lookup(tmp_path):
    sequence = entry(tmp_path, rid="sequence", kind="sequence_set")
    msa = entry(tmp_path, rid="msa", kind="msa", fmt="aligned_fasta", content=">a\nACD\n>b\nACD\n")
    model = entry(tmp_path, rid="model", kind="profile_hmm", fmt="hmm", content="not built")
    msa["derived_from"] = [dependency(sequence)]
    model["derived_from"] = [dependency(msa), dependency(sequence)]
    resources = registry(tmp_path, [sequence, msa, model])
    assert {r.key for r in resources.lineage("model")} == {("sequence", "v1"), ("msa", "v1")}
    # Model artifact need not be present/valid to inspect the fixed derivation provenance.
    second = entry(tmp_path, rid="sequence", kind="sequence_set", version="v2")
    resources = registry(tmp_path, [sequence, msa, model, second])
    with pytest.raises(ReferenceRegistryError, match="Ambiguous"):
        resources.get("sequence")
    assert resources.get("sequence", "v1").sha256 == sequence["sha256"]
    with pytest.raises(ReferenceRegistryError, match="Unknown"):
        resources.get("no_such_id")


@pytest.mark.parametrize(
    "case", ["duplicate", "dangling", "checksum", "cycle", "self", "duplicate_parent"]
)
def test_invalid_derivation_graphs(tmp_path, case):
    a, b = entry(tmp_path, rid="a"), entry(tmp_path, rid="b")
    refs = [a, b]
    if case == "duplicate":
        refs.append(a)
    elif case == "dangling":
        a["derived_from"] = [dependency(b) | {"reference_id": "absent"}]
    elif case == "checksum":
        a["derived_from"] = [dependency(b) | {"sha256": "0" * 64}]
    elif case == "cycle":
        a["derived_from"], b["derived_from"] = [dependency(b)], [dependency(a)]
    elif case == "self":
        a["derived_from"] = [dependency(a)]
    else:
        a["derived_from"] = [dependency(b), dependency(b)]
    with pytest.raises(ReferenceRegistryError):
        registry(tmp_path, refs)


def test_biological_sequence_id_hash_and_interval_validation(tmp_path):
    rec = entry(tmp_path)
    for target in (
        {"sequence_ids": ["wrong"]},
        {"sequence_sha256": "0" * 64},
        {"protein_interval": {"start": 681, "end": 888}, "protein_length": 6497},
    ):
        rec["biological_target"] = {"description": "RN28", **target}
        with pytest.raises(ReferenceRegistryError):
            registry(tmp_path, [rec]).resolve("q")


@pytest.fixture
def hmm_registry(tmp_path):
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(name="demo", sequence=SEQ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(sequence, pyhmmer.plan7.Background(alphabet))
    hmm.accession = b"DEMO.1"
    hmm.cutoffs.gathering = (15.0, 7.0)
    path = tmp_path / "demo.hmm"
    with path.open("wb") as handle:
        hmm.write(handle)
    rec = entry(tmp_path, rid="demo", kind="profile_hmm", fmt="hmm", version="DEMO.1")
    rec["artifact"]["path"] = path.name
    rec["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    rec["biological_target"] = {
        "description": "demo",
        "model_id": "DEMO",
        "model_accession": "DEMO.1",
    }
    return registry(tmp_path, [rec]), path


def test_hmm_profile_adapter_and_ga_are_unchanged(hmm_registry):
    resources, path = hmm_registry
    before = path.read_bytes()
    result = resources.resolve("demo")
    assert result.summary["gathering_cutoffs"] == (15.0, 7.0)
    profile = resources.to_profile("demo")
    assert profile.provenance.source == "local:generated"
    assert load_hmm(profile).cutoffs.gathering == (15.0, 7.0)
    assert path.read_bytes() == before
    released = resources.get("demo").model_dump(mode="json")
    released["version"] = "curated-v2"
    adapter = registry(path.parent, [released]).to_profile("demo")
    assert adapter.provenance.version == "DEMO.1"
    assert load_hmm(adapter).cutoffs.gathering == (15.0, 7.0)
    rec = resources.get("demo").model_dump(mode="json")
    rec["biological_target"]["model_accession"] = "DEMO.2"
    with pytest.raises(ReferenceRegistryError, match="accession/version"):
        registry(path.parent, [rec]).resolve("demo")


def test_duplicate_json_keys_and_bad_registry_paths(tmp_path):
    with pytest.raises(ReferenceRegistryError, match="Cannot load"):
        ReferenceRegistry.load(tmp_path / "missing.json")
    path = tmp_path / "registry.json"
    path.write_text('{"references": [], "references": []}')
    with pytest.raises(ReferenceRegistryError, match="Duplicate JSON key"):
        ReferenceRegistry.load(path)


def test_cli_validation_resolution_profile_export_and_errors(hmm_registry, tmp_path):
    resources, _ = hmm_registry
    runner = CliRunner()
    base = ["--registry", str(resources.path)]
    result = runner.invoke(app, ["reference-list", *base])
    assert result.exit_code == 0 and len(json.loads(result.stdout)) == 1
    result = runner.invoke(app, ["reference-validate", *base])
    assert result.exit_code == 0 and json.loads(result.stdout)["state"] == "validated"
    output = tmp_path / "legacy.profile.json"
    args = ["reference-resolve", *base, "--reference-id", "demo", "--export-profile", str(output)]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert load_profile(output) == resources.to_profile("demo")
    assert runner.invoke(app, args).exit_code == 2  # No overwrite.
    result = runner.invoke(app, ["reference-resolve", *base, "--reference-id", "unknown"])
    assert result.exit_code == 2 and "Unknown reference" in result.output
    rec = resources.get("demo").model_dump(mode="json")
    rec["artifact"].update(storage="external", path="unavailable.hmm")
    registry(tmp_path, [rec])
    assert runner.invoke(app, ["reference-validate", *base]).exit_code == 0
    assert runner.invoke(app, ["reference-validate", *base, "--require-external"]).exit_code == 2


def test_symlink_loop_errors_include_registry_context(tmp_path):
    loop = tmp_path / "registry-loop.json"
    loop.symlink_to(loop.name)
    with pytest.raises(ReferenceRegistryError, match="Cannot load"):
        ReferenceRegistry.load(loop)
    rec = entry(tmp_path)
    resources = registry(tmp_path, [rec])
    path = tmp_path / rec["artifact"]["path"]
    path.unlink()
    path.symlink_to(path.name)
    with pytest.raises(ReferenceRegistryError, match="q@v1"):
        resources.resolve("q")
