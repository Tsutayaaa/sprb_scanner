import pytest
from pydantic import ValidationError

from sprb_scanner.core.models import (
    CandidateRecord,
    ContigRecord,
    DomainHit,
    GeneRecord,
    GenomeDataset,
    GenomeRecord,
    Interval,
    ModelProvenance,
    ProteinRecord,
    ScanStatus,
    internal_id,
)


def dataset_data():
    return dict(
        genome=GenomeRecord(genome_id="g", input_format="genbank"),
        contigs=(ContigRecord(contig_id="c", genome_id="g", source_id="chr", length=100),),
        genes=(
            GeneRecord(
                gene_id="x",
                genome_id="g",
                contig_id="c",
                strand=-1,
                location_parts=(Interval(start=90, end=100), Interval(start=0, end=20)),
            ),
        ),
        proteins=(ProteinRecord(protein_id="p", genome_id="g", gene_id="x", sequence="MAXUO"),),
    )


def test_compound_location_and_serialization():
    data = GenomeDataset(**dataset_data())
    restored = GenomeDataset.model_validate_json(data.model_dump_json(exclude_computed_fields=True))
    assert restored == data
    assert data.genes[0].strand == -1
    assert data.genes[0].location_parts[0].start == 90
    assert data.proteins[0].length == 5


@pytest.mark.parametrize("kind", ["duplicate", "orphan", "genome", "bounds", "missing_gene"])
def test_invalid_dataset_links(kind):
    data = dataset_data()
    if kind == "duplicate":
        data["proteins"] *= 2
    elif kind == "orphan":
        data["genes"] = ()
    elif kind == "genome":
        data["proteins"] = (data["proteins"][0].model_copy(update={"genome_id": "other"}),)
    elif kind == "bounds":
        data["contigs"] = (data["contigs"][0].model_copy(update={"length": 99}),)
    else:
        data["proteins"] = (data["proteins"][0].model_copy(update={"gene_id": None}),)
    with pytest.raises(ValidationError):
        GenomeDataset(**data)


def test_protein_only_and_missing_translation():
    GenomeDataset(
        genome=GenomeRecord(genome_id="g", input_format="protein_fasta"),
        proteins=(ProteinRecord(protein_id="p", genome_id="g", sequence="MXX"),),
    )
    data = dataset_data()
    data["proteins"] = ()
    GenomeDataset(**data)


@pytest.mark.parametrize("sequence", ["", "abc", "M*", "MA-A", "M A"])
def test_invalid_sequence(sequence):
    with pytest.raises(ValidationError):
        ProteinRecord(protein_id="p", genome_id="g", sequence=sequence)


def hit_data():
    return dict(
        hit_id="h",
        protein_id="p",
        start=90,
        end=100,
        method="hmm",
        model=ModelProvenance(
            model_id="TIGR04131", version="test", source="fixture", sha256="0" * 64
        ),
        query_id="q",
        query_length=10,
        query_start=0,
        query_end=10,
        protein_length=100,
        evalue=0,
        evalue_kind="independent_domain",
        bit_score=20,
        query_coverage=1,
        target_coverage=0.1,
        backend_version="test",
    )


def test_hit_terminal_and_candidate_links():
    hit = DomainHit(**hit_data())
    assert hit.ctd_distance == 0
    assert DomainHit(**(hit_data() | {"start": 80, "end": 90})).ctd_distance == 10
    with pytest.raises(ValidationError):
        CandidateRecord(candidate_id="c", protein_id="other", domain_hits=(hit,))
    with pytest.raises(ValidationError):
        CandidateRecord(candidate_id="c", protein_id="p", domain_hits=(hit, hit))


@pytest.mark.parametrize(
    "update",
    [
        {"end": 101},
        {"query_end": 11},
        {"evalue": float("nan")},
        {"target_coverage": 1.1},
        {"start": True},
    ],
)
def test_invalid_hit(update):
    with pytest.raises(ValidationError):
        DomainHit(**(hit_data() | update))


def test_status_and_identifiers():
    assert internal_id("a::b", "c") != internal_id("a", "b::c")
    assert " " not in internal_id("a b", "p")
    with pytest.raises(ValueError):
        internal_id("")
    with pytest.raises(ValidationError):
        ScanStatus(model_id="m", state="failed")
    assert ScanStatus(model_id="m", state="succeeded").message is None
