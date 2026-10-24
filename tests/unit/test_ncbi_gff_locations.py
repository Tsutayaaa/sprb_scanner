"""NCBI joins retain canonical FAA, source coordinates and conservative uncertainty."""

import json
from pathlib import Path

import pytest
from Bio.Seq import Seq

from sprb_scanner.core.config import AppConfig
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.io.gff import parse_gff
from sprb_scanner.io.prepare import prepare
from sprb_scanner.neighborhood import NeighborhoodIndex
from sprb_scanner.neighborhood.analyzer import geometry


def write_inputs(tmp_path, dna, rows, protein="MKK", header="c", preamble=""):
    fna, faa, gff = [tmp_path / name for name in ("test.fna", "test.faa", "test.gff3")]
    fna.write_text(f">{header}\n{dna}\n")
    faa.write_text(f">p\n{protein}\n")
    gff.write_text("##gff-version 3\n" + preamble + "".join(rows))
    return gff, fna, faa


def join_case(tmp_path, join=-1, strand="+", annotated=True, phase=0, protein="MKK"):
    first, last = "ATGAAA", "AAATAA"
    dna = first + ("C" * join + last if join >= 0 else last[-join:])
    intervals = [(0, 6), (6 + join, 12 + join)]
    if strand == "-":
        dna = str(Seq(dna).reverse_complement())
        intervals = [(len(dna) - b, len(dna) - a) for a, b in intervals]
    exception = ";exception=ribosomal slippage" if annotated else ""
    rows = [
        f"c\tRefSeq\tCDS\t{a + 1}\t{b}\t.\t{strand}\t{0 if i == 0 else phase}\t"
        f"ID=cds-p;protein_id=p{exception}\n"
        for i, (a, b) in enumerate(intervals)
    ]
    return write_inputs(tmp_path, dna, list(reversed(rows)), protein), intervals


def parse(paths):
    return parse_gff(*paths, "g")


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("join", [-2, -1, 0, 1, 2])
def test_annotated_slippage_preserves_parts_phases_and_canonical_faa(tmp_path, strand, join):
    paths, expected = join_case(tmp_path, join, strand)
    ds = parse(paths)
    g = ds.genes[0]
    assert len(ds.genes) == len(ds.proteins) == 1
    assert ds.proteins[0].sequence == "MKK" and ds.proteins[0].source_id == "p"
    assert [(p.start, p.end) for p in g.location_parts] == expected
    assert g.source_id == "cds-p" and g.qualifiers["exception"] == ("ribosomal slippage",)
    assert g.qualifiers["gff_phases"] == ("0", "0")
    proof = json.loads(g.qualifiers["gff_location_validation"][0])
    assert proof["representation"] == "annotated_ribosomal_slippage"
    assert proof["joined_translation"] == "exact_faa_match" and proof["phase_consistent"]
    assert [(p["start"] - 1, p["end"]) for p in proof["source_intervals"]] == expected
    assert not any(i.severity == "error" for i in ds.quality_issues)
    index = NeighborhoodIndex(ds)
    indexed = index.by_contig[g.contig_id][0]
    assert sum(p.end - p.start for p in indexed.parts) == 12 + min(join, 0)
    assert index.query(g.gene_id).neighbors == ()


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("overlap", [1, 2])
def test_unannotated_short_overlap_requires_exact_faa_and_keeps_uncertainty(
    tmp_path, strand, overlap
):
    paths, parts = join_case(tmp_path, -overlap, strand, False)
    ds = parse(paths)
    g = ds.genes[0]
    proof = json.loads(g.qualifiers["gff_location_validation"][0])
    assert proof["representation"] == "faa_verified_unannotated_overlap"
    assert proof["biological_interpretation"] == "undetermined"
    assert "exception" not in g.qualifiers and "ribosomal_slippage" not in g.qualifiers
    assert any(i.code == "unannotated_join_uncertain" for i in ds.quality_issues)
    index = NeighborhoodIndex(ds)
    x = index.by_contig[g.contig_id][0]
    assert sum(p.end - p.start for p in x.parts) == 12 - overlap
    # An old proof cannot authorize a different coordinate structure.
    altered = g.model_copy(update={"strand": -g.strand})
    with pytest.raises(InputError):
        NeighborhoodIndex(ds.model_copy(update={"genes": (altered,)}))


@pytest.mark.parametrize("annotated", [True, False])
def test_mismatch_is_not_accepted_merely_due_to_overlap_or_exception(tmp_path, annotated):
    paths, _ = join_case(tmp_path, annotated=annotated, protein="MKR")
    with pytest.raises(InputError, match="verification"):
        parse(paths)


def test_phase_discrepancy_requires_explicit_slippage_and_exact_faa(tmp_path):
    paths, _ = join_case(tmp_path, phase=1)
    ds = parse(paths)
    assert any(
        i.code == "inconsistent_phase" and i.severity == "warning" for i in ds.quality_issues
    )
    assert not json.loads(ds.genes[0].qualifiers["gff_location_validation"][0])["phase_consistent"]
    paths, _ = join_case(tmp_path, annotated=False, phase=1)
    with pytest.raises(InputError, match="verification"):
        parse(paths)


@pytest.mark.parametrize(
    "change",
    [
        "duplicated",
        "contained",
        "missing_faa",
        "wrong_id",
        "extra_faa",
        "conflicting_exception",
        "ambiguous_id",
        "large_overlap",
        "reserved_proof",
    ],
)
def test_invalid_joins_and_mapping_remain_explicit_errors(tmp_path, change):
    paths, _ = join_case(tmp_path, annotated=change != "large_overlap")
    gff, fna, faa = paths
    text = gff.read_text()
    if change == "duplicated":
        text += "\n".join(text.splitlines()[1:]) + "\n"
    elif change == "contained":
        text = text.replace("CDS\t6\t11", "CDS\t2\t5")
    elif change == "missing_faa":
        faa.write_text(">other\nMKK\n")
    elif change == "wrong_id":
        text = text.replace("protein_id=p", "protein_id=wrong")
    elif change == "extra_faa":
        faa.write_text(faa.read_text() + ">extra\nMK\n")
    elif change == "conflicting_exception":
        text = text.replace(
            "exception=ribosomal slippage",
            "exception=ribosomal slippage,unclassified translation discrepancy",
        )
    elif change == "ambiguous_id":
        text = text.replace("protein_id=p", "protein_id=p,q")
    elif change == "large_overlap":
        text = text.replace("CDS\t6\t11", "CDS\t3\t11")
    elif change == "reserved_proof":
        text = text.replace("protein_id=p", "protein_id=p;gff_location_validation=forged")
    gff.write_text(text)
    with pytest.raises(InputError):
        prepare(gff, "gff3", tmp_path / "out", AppConfig(), genome_id="g", fna=fna, faa=faa)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("metadata", ["region", "fasta", "late_region"])
def test_virtual_origin_with_noncoding_feature_and_source_coordinates(tmp_path, strand, metadata):
    # The biological CDS encodes MK across the origin, in either direction.
    forward = "AAATAACCCATG"
    dna = forward if strand == "+" else str(Seq(forward).reverse_complement())
    start, end = (10, 18) if strand == "+" else (7, 15)
    region = "c\tRefSeq\tregion\t1\t12\t.\t+\t.\tID=c;Is_circular=true\n"
    rows = [
        f"c\tRefSeq\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene\n",
        f"c\tRefSeq\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds-p;Parent=gene;protein_id=p\n",
        "c\tRefSeq\trRNA\t11\t15\t.\t+\t.\tID=rna\n",
    ]
    if metadata == "late_region":
        rows.append(region)
    paths = write_inputs(
        tmp_path,
        dna,
        rows,
        "MK",
        header="c [topology=circular]" if metadata == "fasta" else "c",
        preamble=region if metadata == "region" else "",
    )
    ds = parse(paths)
    gene = ds.genes[0]
    assert ds.contigs[0].topology == "circular"
    assert ds.proteins[0].sequence == "MK"
    assert sum(x.end - x.start for x in gene.location_parts) == 9
    proof = json.loads(gene.qualifiers["gff_location_validation"][0])
    assert proof["source_intervals"][0]["end"] == end
    assert proof["crosses_origin"] and proof["joined_translation"] == "exact_faa_match"
    assert NeighborhoodIndex(ds).query(gene.gene_id).target_crosses_origin


@pytest.mark.parametrize(
    "change",
    ["no_flag", "linear", "conflicting", "too_long", "start_outside", "bad_region", "bad_faa"],
)
def test_circularity_is_not_inferred_from_out_of_bounds(tmp_path, change):
    region = "c\tx\tregion\t1\t12\t.\t+\t.\tID=c;Is_circular=true\n"
    start, end = 10, 18
    if change == "no_flag":
        region = ""
    if change == "linear":
        region = region.replace("true", "false")
    if change == "too_long":
        end = 23
    if change == "start_outside":
        start, end = 13, 21
    if change == "bad_region":
        region = region.replace("1\t12", "2\t12")
    rows = [f"c\tx\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds-p;protein_id=p\n"]
    paths = write_inputs(
        tmp_path,
        "AAATAACCCATG",
        rows,
        "MM" if change == "bad_faa" else "MK",
        header="c [topology=linear]" if change == "conflicting" else "c",
        preamble=region,
    )
    with pytest.raises(InputError):
        prepare(
            paths[0],
            "gff3",
            tmp_path / "out",
            AppConfig(),
            genome_id="g",
            fna=paths[1],
            faa=paths[2],
        )


@pytest.mark.parametrize("strand", ["+", "-"])
def test_split_children_ordered_from_virtual_parent_and_adjacency(tmp_path, strand):
    # Raw child rows do not themselves extend beyond the contig.
    forward = "AAATAACCCATGCCCCCCCCCCCCC"
    length = len(forward)
    dna = forward if strand == "+" else str(Seq(forward).reverse_complement())
    parts = (
        [(9, 12), (0, 6)] if strand == "+" else [(length - 12, length - 9), (length - 6, length)]
    )
    parent_start = 10 if strand == "+" else length - 6 + 1
    parent_end = length + (6 if strand == "+" else length - 9)
    rows = [f"c\tx\tgene\t{parent_start}\t{parent_end}\t.\t{strand}\t.\tID=parent\n"]
    for a, b in sorted(parts):
        rows.append(
            f"c\tx\tCDS\t{a + 1}\t{b}\t.\t{strand}\t0\tID=cds-p;protein_id=p;Parent=parent\n"
        )
    paths = write_inputs(tmp_path, dna, rows, "MK", header="c [topology=circular]")
    ds = parse(paths)
    assert ds.proteins[0].sequence == "MK"
    assert not any(i.severity == "error" for i in ds.quality_issues)
    assert [(p.start, p.end) for p in ds.genes[0].location_parts] == parts
    idx = NeighborhoodIndex(ds)
    entry = idx.by_contig[ds.genes[0].contig_id][0]
    assert entry.crosses_origin
    assert geometry(entry.parts, entry.parts, length)[1] == 9


def test_prepared_roundtrip_retains_location_proof(tmp_path):
    paths, _ = join_case(tmp_path, annotated=False)
    output = tmp_path / "prepared"
    prepare(paths[0], "gff3", output, AppConfig(), genome_id="g", fna=paths[1], faa=paths[2])
    original = parse(paths)
    restored = GenomeDataset.model_validate_json((output / "dataset.json").read_text())
    assert restored == original
    assert NeighborhoodIndex(restored).query(restored.genes[0].gene_id).neighbors == ()


REAL_JOINS = json.loads((Path(__file__).parents[1] / "fixtures/ncbi_gff_joins.json").read_text())


@pytest.mark.parametrize("sample", REAL_JOINS, ids=lambda r: r["source_id"])
def test_real_ncbi_slippage_cds_examples(tmp_path, sample):
    gff, fna, faa = write_inputs(
        tmp_path, sample["dna"], [sample["gff"]], sample["protein"], header=sample["seqid"]
    )
    faa.write_text(f">{sample['protein_id']}\n{sample['protein']}\n")
    ds = parse((gff, fna, faa))
    assert len(ds.proteins) == 1 and ds.proteins[0].sequence == sample["protein"]
    assert ds.proteins[0].source_id == sample["protein_id"]
    proof = json.loads(ds.genes[0].qualifiers["gff_location_validation"][0])
    assert proof["joined_translation"] == "exact_faa_match"
    assert proof["phase_consistent"] and not any(i.severity == "error" for i in ds.quality_issues)
    reconstructed = [
        [r["start"] + sample["coordinate_offset"], r["end"] + sample["coordinate_offset"]]
        for r in proof["source_intervals"]
    ]
    assert reconstructed == sample["source_intervals"]
    NeighborhoodIndex(ds)


@pytest.mark.parametrize(
    "exception",
    ["unclassified translation discrepancy", "annotated by transcript or proteomic data"],
)
def test_translation_discrepancy_is_not_relabelled_as_slippage(tmp_path, exception):
    paths, _ = join_case(tmp_path, annotated=False)
    paths[0].write_text(
        paths[0].read_text().replace("protein_id=p", f"protein_id=p;exception={exception}")
    )
    ds = parse(paths)
    gene = ds.genes[0]
    proof = json.loads(gene.qualifiers["gff_location_validation"][0])
    assert proof["representation"] == "faa_verified_translation_discrepancy"
    assert proof["biological_interpretation"] == "undetermined"
    assert gene.qualifiers["exception"] == (exception,)
    assert not any(i.severity == "error" for i in ds.quality_issues)
    NeighborhoodIndex(ds)
