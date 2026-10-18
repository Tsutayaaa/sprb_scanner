# Architecture

## Implemented foundation

`core` contains records, coordinate conversion and configuration. It imports no scanners.
The CLI only invokes application functions and reports errors. There is no scanning logic yet.
Record validation rejects structural corruption; biological annotation problems will be recorded
as `QualityIssue` by parsers, retaining original annotations and source files.

| Module | Purpose | Input | Output | Dependencies | Parameters | Error handling | Test cases |
|---|---|---|---|---|---|---|---|
| core.models | shared schema and referential integrity | parsed fields | validated records/dataset | Pydantic | fields in schema | ValidationError | orphan links, bounds, duplicate IDs, sequence-only mode |
| core.coordinates | boundary conversions | integer intervals | normalized intervals/overlap | standard library | start/end | ValueError | single residue, touching regions, invalid coordinates |
| core.config | load reproducible settings | YAML/path or defaults | resolved AppConfig | PyYAML/Pydantic | paths, threads, BCTD | ConfigurationError | unknown keys, invalid values, relative paths |
| cli | user entry point | arguments | version/config JSON | Typer/core | --config | exit 2 for config errors | help, version, valid/invalid configs |

Records use frozen Pydantic models, but nested dictionaries are not deeply immutable: consumers
must not mutate shared qualifiers/evidence. JSON serialization is supported. Derived fields
(length and CTD distance) are emitted for users; use `exclude_computed_fields=True` for validated
JSON round trips. Batch 3 exporters will explicitly define storage schemas rather than infer them
from empty tables or dump arbitrary nested objects into Parquet.

## Planned interfaces (not implemented)

- Parsers: input paths + explicit genome ID → GenomeDataset + quality issues.
- Downloader: assembly accession + output directory → local file manifest and checksums.
  Parsers remain fully usable offline, without the downloader.
- Scanner: protein records + validated model/config → raw DomainHits + ScanStatus + artifacts.
- Merger: hits + per-model status → CandidateRecords and overlap relationships; preserve all hits,
  do not average or rank HMM and MMseqs2 scores together. Report status even for zero hits.
- Exporter: dataset + hits + candidates + run metadata → separate typed tables and FASTA/JSON.
- SprF/RN28/CBD: proteins → model-specific hits, depending only on core.
- Locus: genes/contigs/candidates → neighborhoods respecting contig topology and strand.
- Motility: genome/metadata → species evidence. Classifier: all evidence → ranked candidates.
- Structure: candidate proteins → structural evidence. Viewers: records/evidence → visual output.

Future modules are documented rather than represented by empty packages. BioPython, PyHMMER,
MMseqs2 and Parquet dependencies are introduced with their implementations, not in this batch.
