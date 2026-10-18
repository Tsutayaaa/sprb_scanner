# Phase 1 roadmap

Each batch is independently reviewed and accepted by Tsutayaaa before merging.

1. **Project Foundation (merged PR #1):** src package, uv lock, CLI foundation, six core records,
   dataset validation, coordinates, strict configuration, automated unit/integration checks and CI.
2. **Genome Infrastructure (this batch):** BioPython GBFF/GFF3+FNA+FAA/protein-only parsers, deterministic
   IDs, annotation/translation QC, NCBI Datasets downloader and normalized metadata/proteome.
3. **BCTD Discovery:** TIGR04131 HMM provenance verification, PyHMMER, reference CTD MMseqs2
   and custom HMM, merger, full CLI, TSV/Parquet/JSON/FASTA, real FJ regression validation.

Batch 1 acceptance: uv sync, pytest, Ruff, CLI help and package build. GitHub CI must pass too.
Batch 2 adds offline FJ sprB/sprF parsing regression plus opt-in full assembly validation.
No stub scanner commands report success; BCTD biological validation remains batch 3.

## Decisions before biological acceptance

- Parser reference established: F. johnsoniae UW101, GCF_000016645.1, RS_2025_12_11.
- Obtain the owner's SprB CTD reference FASTA/HMM and verified boundary/source annotations.
- Establish literature-supported sprB/sprF locations and known Type B cargo reference set.
- Confirm experimental reporting/terminal-distance thresholds against that reference set.
- Record model download URLs, versions, checksums and redistribution terms before bundling.

Actual downloaded genomes, models and outputs remain outside Git. Small deterministic fixtures
will cover multi-contig, reverse-strand, partial/pseudo CDS, compound locations, repeated IDs,
mismatched annotations, terminal/internal/multiple hits, model failure, and empty outputs.
