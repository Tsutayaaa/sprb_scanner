# FJ fixture provenance

Source assembly: GCF_000016645.1, NC_009441.1, NCBI RefSeq annotation RS_2025_12_11.
Downloaded on 2026-10-18. Full source checksums and reference coordinates: fj_reference.json.
No whole genomes are committed. `fj_window.*` contains only the 20,512 nt interval
[1094604,1115116), enclosing sprF and sprB. Contig renamed NC_009441.1_window, coordinates
shifted by -1094604. GBFF/FNA were sliced with BioPython; GFF rows within the window were
selected and shifted; FAA contains the two annotated translations without modification.
The full package is external and is never downloaded automatically by tests.

Gene identity anchors: sprB = Fjoh_0979 (FJOH_RS26780); sprF = Fjoh_0978 (FJOH_RS05105).
Family product names alone are not identity evidence; Fjoh_0980 is SprD, not SprF.
Sources: [NCBI assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000016645.1/),
[SprB NCBI record](https://www.ncbi.nlm.nih.gov/protein/ABQ04012.1),
[De Smet et al. 2025, Fig. 6](https://www.nature.com/articles/s42003-025-07817-2),
[Rhodes et al., spr operon](https://pubmed.ncbi.nlm.nih.gov/21131497/).
This fixture checks parsing/identity/coordinates, not a new homology classification.
