# Phase 1 roadmap

Each batch is independently reviewed and accepted by Tsutayaaa before merging.

1. **Project Foundation (merged PR #1):** src package, uv lock, CLI foundation, six core records,
   dataset validation, coordinates, strict configuration, automated unit/integration checks and CI.
2. **Genome Infrastructure (merged PR #2):** BioPython GBFF/GFF3+FNA+FAA/protein-only parsers, deterministic
   IDs, annotation/translation QC, NCBI Datasets downloader and normalized metadata/proteome.
3. **HMM evidence (merged PR #3):** reusable PyHMMER engine, GA tiers, TIGR04131/PF13585 CTD,
   TIGR03519/PF11751 PorP/SprF family and PF13573 repeat evidence.
4. **Batch 4A (this batch):** independent indexed CDS neighborhoods, exact part geometry,
   circular handling, prepared-data CLI and FJ SprB–SprF validation.
5. **Batch 4B (next):** associate Type B CTD and PorP/SprF evidence with objective neighborhoods.

MMseqs2/custom HMM orchestration, evidence classification and Parquet remain future work.

Batch 1 acceptance: uv sync, pytest, Ruff, CLI help and package build. GitHub CI must pass too.
Batch 2 adds offline FJ sprB/sprF parsing regression plus opt-in full assembly validation.
No stub commands report success; final biological classification remains future work.

## Decisions before biological acceptance

- Parser reference established: F. johnsoniae UW101, GCF_000016645.1, RS_2025_12_11.
- Obtain the owner's SprB CTD reference FASTA/HMM and verified boundary/source annotations.
- Establish literature-supported sprB/sprF locations and known Type B cargo reference set.
- Confirm experimental reporting/terminal-distance thresholds against that reference set.
- Record model download URLs, versions, checksums and redistribution terms before bundling.

Actual downloaded genomes, models and outputs remain outside Git. Small deterministic fixtures
will cover multi-contig, reverse-strand, partial/pseudo CDS, compound locations, repeated IDs,
mismatched annotations, terminal/internal/multiple hits, model failure, and empty outputs.
