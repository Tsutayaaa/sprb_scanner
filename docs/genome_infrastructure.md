# Genome infrastructure module contracts

| Module | Purpose | Input → Output | Dependencies | Parameters | Error handling | Tests |
|---|---|---|---|---|---|---|
| io.common | checksums and sequence QC | paths/sequences → manifest/issues | BioPython, stdlib | translation table | InputError or recorded quality issue | parser tests |
| io.fasta | protein-only normalization | FAA → GenomeDataset | BioPython/core | genome ID | empty/malformed input fails; invalid residues recorded | protein-only and gzip cases |
| io.genbank | annotated CDS parsing | GBFF → GenomeDataset | BioPython/core | genome ID, default table | unsupported locations fail; missing translation retained as warning | compound, reverse, pseudo, FJ |
| io.gff | exact annotation/protein association | GFF3 + FNA + FAA → GenomeDataset | BioPython/core | genome ID, default table | ambiguous IDs fail; translation/mapping problems recorded | phases, missing files, circular, FJ |
| io.prepare | preparation service | input specification/config → prepared dataset directory | parsers/exporter | format, output, strict quality | strict quality stops export; source changes detected | CLI integration |
| io.exporters | atomic output and integrity | dataset → JSON/TSV/FASTA | core/stdlib | new output path | existing output rejected; output hash verification | round trip/corruption |
| io.ncbi | independent acquisition | one versioned accession → package/manifest | stdlib, NCBI Datasets REST | timeout/retries/byte limits | bounded retries, unsafe/incomplete ZIP rejection, cleanup | mocked offline failures + recorded real download |

The downloader uses NCBI Datasets v2 REST rather than requiring the standalone `datasets` CLI.
It downloads one assembly per invocation, with GBFF/GFF3/FNA/FAA, rejects an existing destination,
and saves source URL, NCBI report, UTC download time, file sizes and SHA-256. Parsing never uses network.
Defaults: timeout 120 seconds, two retries, 100 MB compressed and 500 MB extracted limits.
These are configurable byte limits, not a claim to support arbitrary large/eukaryotic assemblies.

`validate --prepared` checks checksums and dataset schema/relationships. It detects accidental
corruption; the manifest is not a digital signature or a replacement for biological validation.
Dataset JSON carries all CDS qualifiers; metadata TSV is a convenient subset.
The NCBI assembly report preserves organism/annotation details even when GFF alone lacks a name.
