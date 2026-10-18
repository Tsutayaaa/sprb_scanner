import gzip

import pytest
from Bio import SeqIO

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.fasta import parse_proteins
from sprb_scanner.io.genbank import parse_genbank
from sprb_scanner.io.gff import attributes, parse_gff
from sprb_scanner.io.prepare import prepare_dataset


def test_genbank_and_gff_agree(genome_files):
    f = genome_files
    a = parse_genbank(f["gb"], "g")
    b = parse_gff(f["gff"], f["fna"], f["faa"], "g")
    assert not a.quality_issues and not b.quality_issues
    assert [p.sequence for p in a.proteins] == [p.sequence for p in b.proteins] == ["MK", "MK"]
    assert [g.location_parts for g in a.genes] == [g.location_parts for g in b.genes]
    assert [g.strand for g in b.genes] == [1, -1]
    assert b.genes[0].gene_name == "sprB"  # inherited from parent
    assert a.proteins[1].protein_id == b.proteins[1].protein_id


def test_compound_partial_pseudo_and_multicontig(complex_genbank):
    d = parse_genbank(complex_genbank, "g")
    assert len(d.contigs) == 2 and len(d.genes) == 4 and len(d.proteins) == 2
    assert d.genes[0].location_parts[0].start == 24
    assert d.genes[0].location_parts[1].end == 6
    codes = {x.code for x in d.quality_issues}
    assert {
        "pseudogene",
        "partial_cds",
        "missing_translation",
        "duplicate_source_protein_id",
    } <= codes


def test_gzip_and_content_ids(genome_files, tmp_path):
    f = genome_files["gb"]
    compressed = tmp_path / "test.gbff.gz"
    with gzip.open(compressed, "wb") as out:
        out.write(f.read_bytes())
    assert len(parse_genbank(compressed, "g").proteins) == 2
    copy = tmp_path / "renamed.gbff"
    copy.write_bytes(f.read_bytes())
    a, _ = prepare_dataset(f, "genbank")
    b, _ = prepare_dataset(copy, "genbank")
    assert a.genome.genome_id == b.genome.genome_id


def test_protein_normalization_and_duplicates(tmp_path):
    p = tmp_path / "p.faa"
    p.write_text(">same\nmk*\n>same\nMX\n>bad\nM*X\n")
    d = parse_proteins(p, "g")
    assert len(d.proteins) == 2 and d.proteins[0].sequence == "MK"
    assert len({x.protein_id for x in d.proteins}) == 2
    assert not d.genes and not d.contigs
    assert {x.code for x in d.quality_issues} >= {
        "invalid_protein",
        "terminal_stop_removed",
        "uppercase_normalization",
        "ambiguous_residues",
    }


@pytest.mark.parametrize("text", ["", "MK", ">empty\n"])
def test_reject_empty_or_non_fasta(tmp_path, text):
    p = tmp_path / "bad.faa"
    p.write_text(text)
    with pytest.raises(InputError):
        parse_proteins(p, "g")


@pytest.mark.parametrize(
    "change",
    ["seqid", "phase", "missing_faa", "duplicate_faa", "gtf", "embedded", "parent", "bounds"],
)
def test_gff_contract_errors(genome_files, change):
    f = genome_files
    text = f["gff"].read_text()
    if change == "seqid":
        text = text.replace("chr1", "wrong")
    elif change == "phase":
        text = text.replace("+\t0", "+\t.")
    elif change == "missing_faa":
        f["faa"].write_text(">other\nMK\n")
        d = parse_gff(f["gff"], f["fna"], f["faa"], "g")
        assert any(x.code == "unmapped_faa_record" for x in d.quality_issues)
        assert any(x.code == "missing_translation" for x in d.quality_issues)
        return
    elif change == "duplicate_faa":
        f["faa"].write_text(">p1\nMK\n>p1\nMK\n")
    elif change == "gtf":
        text = text.replace("##gff-version 3", "##gff-version 2")
    elif change == "embedded":
        text += "##FASTA\n>chr1\nATG\n"
    elif change == "parent":
        text = text.replace("Parent=g1", "Parent=missing")
    elif change == "bounds":
        text = text.replace("CDS\t13\t21", "CDS\t13\t31")
    f["gff"].write_text(text)
    with pytest.raises(InputError):
        parse_gff(f["gff"], f["fna"], f["faa"], "g")


def test_translation_mismatch(genome_files):
    f = genome_files
    f["faa"].write_text(">p1\nMM\n>p2\nMK\n")
    d = parse_gff(f["gff"], f["fna"], f["faa"], "g")
    assert any(x.code == "translation_mismatch" and x.severity == "error" for x in d.quality_issues)


def test_gff_joined_phase_and_escape(tmp_path):
    fna, faa, gff = [tmp_path / name for name in ("x.fna", "x.faa", "x.gff3")]
    fna.write_text(">c\nATGACCAATAA\n")
    faa.write_text(">p\nMK\n")
    gff.write_text(
        "##gff-version 3\n"
        "c\tx\tCDS\t1\t4\t.\t+\t0\tID=p;product=a%3Bb%2Cc\n"
        "c\tx\tCDS\t7\t11\t.\t+\t2\tID=p;product=a%3Bb%2Cc\n"
    )
    d = parse_gff(gff, fna, faa, "g")
    assert d.genes[0].product == "a;b,c"
    assert not d.quality_issues
    gff.write_text(gff.read_text().replace("+\t2", "+\t0"))
    assert any(x.code == "inconsistent_phase" for x in parse_gff(gff, fna, faa, "g").quality_issues)


def test_circular_gff_virtual_end(tmp_path):
    fna, faa, gff = [tmp_path / name for name in ("x.fna", "x.faa", "x.gff3")]
    fna.write_text(">c\nAAATAACCCATG\n")
    faa.write_text(">p\nMK\n")
    gff.write_text(
        "##gff-version 3\n"
        "c\tx\tregion\t1\t12\t.\t+\t.\tID=c;Is_circular=true\n"
        "c\tx\tCDS\t10\t18\t.\t+\t0\tID=p\n"
    )
    d = parse_gff(gff, fna, faa, "g")
    assert [(p.start, p.end) for p in d.genes[0].location_parts] == [(9, 12), (0, 6)]
    assert not d.quality_issues


def test_attributes_and_unsupported_input(genome_files):
    assert attributes("product=some%3Btext%2Cmore")["product"] == ("some;text,more",)
    with pytest.raises(InputError):
        attributes("ID=a;ID=b")
    with pytest.raises(InputError):
        prepare_dataset(genome_files["gff"], "gff3")
    with pytest.raises(InputError):
        prepare_dataset(genome_files["fna"], "fna")


def test_genbank_missing_location_and_invalid_sequence(genome_files):
    rec = genome_files["record"]
    rec.features[0].qualifiers["translation"] = ["M*X"]
    SeqIO.write(rec, genome_files["gb"], "genbank")
    d = parse_genbank(genome_files["gb"], "g")
    assert any(x.code == "invalid_protein" for x in d.quality_issues)


def test_cross_contig_parent_is_rejected(genome_files):
    f = genome_files
    f["fna"].write_text(f["fna"].read_text() + ">chr2\nATGAAATAA\n")
    f["gff"].write_text(f["gff"].read_text().replace("chr1\tx\tgene", "chr2\tx\tgene"))
    with pytest.raises(InputError, match="different contig"):
        parse_gff(f["gff"], f["fna"], f["faa"], "g")


def test_partial_annotation_on_later_part(tmp_path):
    fna, faa, gff = [tmp_path / name for name in ("x.fna", "x.faa", "x.gff3")]
    fna.write_text(">c\nATGACCAATAA\n")
    faa.write_text(">p\nMK\n")
    gff.write_text(
        "##gff-version 3\n"
        "c\tx\tCDS\t1\t4\t.\t+\t0\tID=p\n"
        "c\tx\tCDS\t7\t11\t.\t+\t2\tID=p;partial=true;end_range=11,.\n"
    )
    dataset = parse_gff(gff, fna, faa, "g")
    assert dataset.genes[0].partial
    assert "end_range" in dataset.genes[0].qualifiers
    assert any(i.code == "translation_check_skipped" for i in dataset.quality_issues)
