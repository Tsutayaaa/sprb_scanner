import pytest

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import ContigRecord, GeneRecord, GenomeDataset, GenomeRecord, Interval
from sprb_scanner.neighborhood import NeighborhoodIndex


def gene(name, *parts, strand=1, contig="c", **kwargs):
    return GeneRecord(
        gene_id=name,
        genome_id="g",
        contig_id=contig,
        strand=strand,
        location_parts=tuple(Interval(start=a, end=b) for a, b in parts),
        **kwargs,
    )


def dataset(genes, circular=False):
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
            ContigRecord(contig_id="other", genome_id="g", source_id="chr2", length=100),
        ),
        genes=tuple(genes),
    )


@pytest.mark.parametrize("strand", [1, -1])
def test_transcription_direction_adjacency_and_contig_boundary(strand):
    data = dataset(
        [
            gene("a", (0, 10), pseudogene=True),
            gene("b", (20, 30), strand=strand),
            gene("c", (40, 50), partial=True),
            gene("other", (15, 25), contig="other"),
        ]
    )
    index = NeighborhoodIndex(data)
    result = index.query("b", 1)
    left, right = result.neighbors
    assert [x.gene.gene_id for x in result.neighbors] == ["a", "c"]
    assert left.genomic_position == "upstream" and right.genomic_position == "downstream"
    assert left.relative_position == ("upstream" if strand == 1 else "downstream")
    assert right.relative_position == ("downstream" if strand == 1 else "upstream")
    assert all(
        x.distance_bp == 10 and x.directly_adjacent and x.intervening_cds == 0
        for x in result.neighbors
    )
    assert left.gene.pseudogene and right.gene.partial
    assert len(index.query("a", 99).neighbors) == 2
    assert not index.query("other").neighbors
    assert index.query("b", 0).neighbors == ()


def test_circular_dedup_and_cross_origin_compound():
    # Negative-strand biological order: [0,10) followed by [90,100).
    data = dataset(
        [
            gene("wrap", (0, 10), (90, 100), strand=-1),
            gene("right", (20, 30)),
            gene("left", (70, 80)),
        ],
        True,
    )
    index = NeighborhoodIndex(data)
    result = index.query("wrap", 50)
    assert result.target_crosses_origin
    assert {n.gene.gene_id for n in result.neighbors} == {"left", "right"}
    assert len(result.neighbors) == 2
    assert all(n.distance_bp == 10 and n.crosses_origin for n in result.neighbors)
    assert (
        next(n for n in result.neighbors if n.gene.gene_id == "right").relative_position
        == "upstream"
    )
    assert result.target.location_parts == data.genes[0].location_parts
    # Positive-strand order for the same wrapped region gives the same coordinate anchor.
    plus = dataset([gene("wrap", (90, 100), (0, 10)), gene("right", (20, 30))], True)
    assert NeighborhoodIndex(plus).query("wrap").target_crosses_origin
    one = dataset([gene("only", (10, 20))], True)
    assert not NeighborhoodIndex(one).query("only", 500).neighbors


def test_overlap_nested_touching_and_compound_gap():
    data = dataset([gene("outer", (10, 60)), gene("nested", (20, 30)), gene("touch", (60, 70))])
    result = NeighborhoodIndex(data).query("outer", 2)
    assert result.neighbors[0].distance_bp == -10
    assert result.neighbors[0].overlap_bp == 10 and result.neighbors[0].overlaps
    assert result.neighbors[1].distance_bp == 0 and not result.neighbors[1].overlaps
    assert result.neighbors[1].intervening_cds == 1 and not result.neighbors[1].directly_adjacent
    split = dataset([gene("split", (0, 10), (90, 100)), gene("inside_gap", (40, 50))])
    neighbor = NeighborhoodIndex(split).query("split").neighbors[0]
    assert not neighbor.overlaps and neighbor.distance_bp == 30


def test_circular_minimum_gap_and_unknown_strand():
    data = dataset([gene("a", (0, 10), strand=0), gene("b", (90, 95))], True)
    result = NeighborhoodIndex(data).query("a", 1)
    neighbor = result.neighbors[0]
    assert neighbor.distance_bp == 5 and neighbor.distance_crosses_origin
    assert neighbor.relative_position is None and neighbor.same_strand is None
    assert neighbor.transcription_cds_offset is None


def test_invalid_inputs_and_index_reuse():
    index = NeighborhoodIndex(dataset([gene("a", (0, 10)), gene("b", (20, 30))]))
    for window in [-1, True, 1.5]:
        with pytest.raises(InputError, match="window"):
            index.query("a", window)
    with pytest.raises(InputError, match="Unknown"):
        index.query("missing")
    for bad in [gene("bad", (90, 100), (0, 10)), gene("bad", (0, 20), (10, 30))]:
        with pytest.raises(InputError):
            NeighborhoodIndex(dataset([bad]))

    class IndexedOnly(tuple):
        def __iter__(self):
            raise AssertionError("query must not iterate over the complete contig")

    index.by_contig["c"] = IndexedOnly(index.by_contig["c"])
    assert index.query("a", 1).neighbors[0].gene.gene_id == "b"
    assert index.query("b", 1).neighbors[0].gene.gene_id == "a"
