import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def genome_files(tmp_path):
    # + CDS 1..9 and - CDS 13..21 both encode MK, including a terminal stop.
    seq = "ATGAAATAA" + "CCC" + "TTATTTCAT" + "C" * 9
    rec = SeqRecord(Seq(seq), id="chr1", name="chr1", description="synthetic")
    rec.annotations = {
        "molecule_type": "DNA",
        "organism": "Synthetic bacterium",
        "topology": "linear",
    }
    rec.features = [
        SeqFeature(
            SimpleLocation(0, 9, strand=1),
            type="CDS",
            qualifiers={
                "protein_id": ["p1"],
                "locus_tag": ["L1"],
                "gene": ["sprB"],
                "translation": ["MK"],
                "transl_table": ["11"],
            },
        ),
        SeqFeature(
            SimpleLocation(12, 21, strand=-1),
            type="CDS",
            qualifiers={
                "protein_id": ["p2"],
                "locus_tag": ["L2"],
                "gene": ["sprF"],
                "translation": ["MK"],
                "transl_table": ["11"],
            },
        ),
    ]
    gb = tmp_path / "genomic.gbff"
    SeqIO.write(rec, gb, "genbank")
    fna = tmp_path / "genomic.fna"
    SeqIO.write(rec, fna, "fasta")
    faa = tmp_path / "protein.faa"
    faa.write_text(">p1 test\nMK\n>p2 test\nMK\n")
    gff = tmp_path / "genomic.gff"
    gff.write_text(
        "##gff-version 3\n##sequence-region chr1 1 30\n"
        "chr1\tx\tgene\t1\t9\t.\t+\t.\tID=g1;locus_tag=L1;gene=sprB\n"
        "chr1\tx\tCDS\t1\t9\t.\t+\t0\tID=cds-p1;Parent=g1;protein_id=p1\n"
        "chr1\tx\tCDS\t13\t21\t.\t-\t0\tID=cds-p2;protein_id=p2;locus_tag=L2;gene=sprF\n"
    )
    return {"gb": gb, "gff": gff, "fna": fna, "faa": faa, "record": rec}


@pytest.fixture
def complex_genbank(tmp_path):
    rec = SeqRecord(Seq("ATG" * 10), id="circ", name="circ")
    rec.annotations = {"molecule_type": "DNA", "topology": "circular"}
    rec.features = [
        SeqFeature(
            CompoundLocation([SimpleLocation(24, 30, strand=1), SimpleLocation(0, 6, strand=1)]),
            type="CDS",
            qualifiers={"partial": [""], "translation": ["MMMM"], "protein_id": ["same"]},
        ),
        SeqFeature(
            SimpleLocation(BeforePosition(6), 12, strand=-1),
            type="CDS",
            qualifiers={"translation": ["HH"], "protein_id": ["same"]},
        ),
        SeqFeature(SimpleLocation(12, 18, strand=1), type="CDS", qualifiers={"pseudo": [""]}),
        SeqFeature(SimpleLocation(18, 24, strand=1), type="CDS"),
    ]
    other = SeqRecord(Seq("ATGAAATAA"), id="other", name="other")
    other.annotations = {"molecule_type": "DNA"}
    path = tmp_path / "complex.gbk"
    SeqIO.write([rec, other], path, "genbank")
    return path
