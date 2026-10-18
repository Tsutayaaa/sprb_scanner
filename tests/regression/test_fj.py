"""Offline FJ window regression and opt-in full assembly consistency check."""

import json
import os
from pathlib import Path

import pytest

from sprb_scanner.io.common import sha256
from sprb_scanner.io.fasta import parse_proteins
from sprb_scanner.io.genbank import parse_genbank
from sprb_scanner.io.gff import parse_gff

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
REFERENCE = json.loads((FIXTURES / "fj_reference.json").read_text())


def assert_reference_genes(dataset, offset):
    for ref in REFERENCE["reference_genes"].values():
        matches = [g for g in dataset.genes if g.locus_tag == ref["locus_tag"]]
        assert len(matches) == 1
        gene = matches[0]
        assert ref["old_locus_tag"] in gene.qualifiers["old_locus_tag"]
        assert [(p.start + offset, p.end + offset) for p in gene.location_parts] == [
            (ref["start"], ref["end"])
        ]
        assert gene.strand == ref["strand"]
        protein = next(p for p in dataset.proteins if p.gene_id == gene.gene_id)
        assert protein.source_id == ref["protein_id"]
        assert protein.length == ref["length"]


def assert_agree(a, b):
    assert [
        (g.locus_tag, g.location_parts, g.strand, g.partial, g.pseudogene) for g in a.genes
    ] == [(g.locus_tag, g.location_parts, g.strand, g.partial, g.pseudogene) for g in b.genes]
    assert [(p.protein_id, p.source_id, p.sequence) for p in a.proteins] == [
        (p.protein_id, p.source_id, p.sequence) for p in b.proteins
    ]
    assert not [i for d in (a, b) for i in d.quality_issues if i.severity == "error"]


def test_offline_sprb_sprf_window():
    a = parse_genbank(FIXTURES / "fj_window.gbff", "fj_window")
    b = parse_gff(
        FIXTURES / "fj_window.gff3",
        FIXTURES / "fj_window.fna",
        FIXTURES / "fj_window.faa",
        "fj_window",
    )
    assert len(a.genes) == len(a.proteins) == 2
    assert_reference_genes(a, 1094604)
    assert_reference_genes(b, 1094604)
    assert_agree(a, b)
    c = parse_proteins(FIXTURES / "fj_window.faa", "fj_window")
    assert {p.sequence for p in c.proteins} == {p.sequence for p in a.proteins}
    assert not c.genes and not c.contigs


@pytest.mark.external_data
def test_full_fj_assembly():
    location = os.environ.get("SPRB_FJ_DATA")
    if not location:
        pytest.skip(
            "Set SPRB_FJ_DATA to an already downloaded assembly directory; no auto-download"
        )
    root = Path(location)
    for name, record in REFERENCE["files"].items():
        assert sha256(root / name) == record["sha256"], "Annotation snapshot changed"
    a = parse_genbank(root / "genomic.gbff", REFERENCE["accession"])
    b = parse_gff(
        root / "genomic.gff",
        next(root.glob("*_genomic.fna")),
        root / "protein.faa",
        REFERENCE["accession"],
    )
    c = parse_proteins(root / "protein.faa", REFERENCE["accession"])
    assert len(a.genes) == REFERENCE["cds"]
    assert len(a.proteins) == REFERENCE["mapped_proteins"]
    assert len(c.proteins) == REFERENCE["faa_records"]
    assert {p.sequence for p in a.proteins} == {p.sequence for p in c.proteins}
    assert_agree(a, b)
    assert_reference_genes(a, 0)
    assert_reference_genes(b, 0)
