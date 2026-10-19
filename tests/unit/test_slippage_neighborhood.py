"""Original CDS provenance and exact union geometry for supported slippage joins."""

import json
from pathlib import Path

import pytest

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import ContigRecord, GeneRecord, GenomeDataset, GenomeRecord, Interval
from sprb_scanner.neighborhood import NeighborhoodIndex
from sprb_scanner.neighborhood.analyzer import genomic_union, geometry, index_gene

SLIP = {"ribosomal_slippage": ("",), "note": ("programmed frameshift",)}
REFERENCES = json.loads(
    (Path(__file__).parents[1] / "fixtures/slippage_reference.json").read_text()
)


def gene(name, parts, strand=1, qualifiers=None):
    return GeneRecord(
        gene_id=name,
        genome_id="g",
        contig_id="c",
        strand=strand,
        location_parts=tuple(Interval(start=a, end=b) for a, b in parts),
        qualifiers=qualifiers or {},
    )


def data(genes, circular=False):
    return GenomeDataset(
        genome=GenomeRecord(genome_id="g", input_format="genbank"),
        contigs=(
            ContigRecord(
                contig_id="c",
                genome_id="g",
                source_id="chr",
                length=100,
                topology="circular" if circular else "linear",
            ),
        ),
        genes=tuple(genes),
    )


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("overlap", [0, 1, 2, 5, 19])
def test_union_preserves_original_parts_direction_and_distances(strand, overlap):
    forward = [(0, 20), (20 - overlap, 40)]
    parts = forward if strand == 1 else list(reversed(forward))
    g = gene("slip", parts, strand, SLIP if overlap else {})
    before = g.model_dump(mode="json")
    ds = data([g, gene("right", [(50, 60)])])
    index = NeighborhoodIndex(ds)
    assert index.by_contig["c"][0].parts == (Interval(start=0, end=40),)
    n = index.query("slip").neighbors[0]
    assert n.distance_bp == 10 and n.overlap_bp == 0 and n.directly_adjacent
    assert n.relative_position == ("downstream" if strand == 1 else "upstream")
    assert index.query("right").neighbors[0].gene is g
    assert g.model_dump(mode="json") == before


@pytest.mark.parametrize(
    "qualifiers", [SLIP, {"exception": ("ribosomal slippage",)}, {"ribosomal_slippage": ()}]
)
def test_explicit_slippage_evidence_and_union_overlap_not_double_counted(qualifiers):
    a = gene("a", [(0, 20), (19, 30)], qualifiers=qualifiers)
    b = gene("b", [(18, 20), (19, 21)], qualifiers=qualifiers)
    index = NeighborhoodIndex(data([a, b]))
    for target in ["a", "b"]:
        n = index.query(target).neighbors[0]
        assert (n.distance_bp, n.overlap_bp, n.signed_gap_bp) == (0, 3, -3)
        assert n.overlaps and n.directly_adjacent
    x, y = index.by_contig["c"]
    assert geometry(x.parts, y.parts) == geometry(y.parts, x.parts) == (0, 3, False)


@pytest.mark.parametrize("strand", [1, -1])
def test_slippage_compound_gaps_are_not_filled(strand):
    parts = [(0, 10), (9, 20), (40, 50)]
    g = gene("slip", parts if strand == 1 else list(reversed(parts)), strand, SLIP)
    index = NeighborhoodIndex(data([g, gene("inside_gap", [(25, 30)])]))
    assert index.by_contig["c"][0].parts == (Interval(start=0, end=20), Interval(start=40, end=50))
    n = index.query("slip").neighbors[0]
    assert n.distance_bp == 5 and not n.overlaps


@pytest.mark.parametrize("strand", [1, -1])
def test_circular_origin_with_slippage_union_and_symmetric_overlap(strand):
    parts = [(90, 100), (0, 10), (9, 20)]
    g = gene("wrap", parts if strand == 1 else list(reversed(parts)), strand, SLIP)
    index = NeighborhoodIndex(data([g, gene("near", [(5, 15), (14, 25)], qualifiers=SLIP)], True))
    result = index.query("wrap")
    assert result.target_crosses_origin and result.target is g
    n = result.neighbors[0]
    assert n.overlap_bp == 15 and n.distance_bp == 0 and n.crosses_origin
    assert n.directly_adjacent
    x, y = index.by_contig["c"]
    assert geometry(x.parts, y.parts, 100) == geometry(y.parts, x.parts, 100) == (0, 15, False)
    assert g.location_parts == tuple(
        Interval(start=a, end=b) for a, b in (parts if strand == 1 else list(reversed(parts)))
    )


@pytest.mark.parametrize(
    "qualifiers",
    [
        {},
        {"note": ("programmed frameshift",)},
        {"exception": ("unclassified translation discrepancy",)},
        {"note": ("suspected assembly error",)},
        {"ribosomal_slippage": ("false",)},
        {**SLIP, "exception": ("unclassified translation discrepancy",)},
    ],
)
def test_unexplained_or_conflicting_overlap_still_errors(qualifiers):
    with pytest.raises(InputError):
        NeighborhoodIndex(data([gene("bad", [(0, 20), (19, 30)], qualifiers=qualifiers)]))


@pytest.mark.parametrize(
    "parts,circular",
    [
        ([(0, 20), (0, 20)], False),
        ([(0, 30), (5, 20)], False),
        ([(0, 20), (50, 70), (10, 30)], False),
        ([(90, 100), (0, 95)], True),  # more than one circle
        ([(90, 100), (0, 20), (80, 100), (0, 10)], True),  # second origin traversal
    ],
)
def test_flag_cannot_authorize_nonprogressive_or_multiturn_locations(parts, circular):
    with pytest.raises(InputError):
        NeighborhoodIndex(data([gene("bad", parts, qualifiers=SLIP)], circular))


@pytest.mark.parametrize("ref", REFERENCES["records"], ids=lambda r: r["gene_id"])
def test_all_48_real_slippage_locations(ref):
    parts = tuple(Interval.model_validate(p) for p in ref["location_parts"])
    contig = ContigRecord(
        contig_id="c", genome_id="g", source_id="original", length=max(p.end for p in parts) + 10
    )
    g = gene(ref["gene_id"], [(p.start, p.end) for p in parts], ref["strand"], ref["qualifiers"])
    before = g.model_dump(mode="json")
    indexed = index_gene(g, contig)
    expected = set().union(*(set(range(p.start, p.end)) for p in parts))
    observed = set().union(*(set(range(p.start, p.end)) for p in indexed.parts))
    assert observed == expected
    assert sum(p.end - p.start for p in indexed.parts) == len(expected)
    assert indexed.gene is g and g.model_dump(mode="json") == before
    # Reversing strand and biological order must give the same geometric union.
    flipped = g.model_copy(update={"strand": -g.strand, "location_parts": tuple(reversed(parts))})
    assert index_gene(flipped, contig).parts == indexed.parts
    assert genomic_union(parts) == indexed.parts
