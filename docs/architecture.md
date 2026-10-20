# Architecture

## Implemented foundation

`core` contains records, coordinate conversion and configuration. It imports no scanners.
The CLI invokes preparation/download/validation services and reports errors. The single-profile HMM engine is implemented; see hmm_engine.md.
See genome_infrastructure.md for implemented parser, preparation, downloader and export contracts.
Record validation rejects structural corruption; biological annotation problems are recorded
as `QualityIssue` by parsers, retaining original annotations and source files.

| Module | Purpose | Input | Output | Dependencies | Parameters | Error handling | Test cases |
|---|---|---|---|---|---|---|---|
| core.models | shared schema and referential integrity | parsed fields | validated records/dataset | Pydantic | fields in schema | ValidationError | orphan links, bounds, duplicate IDs, sequence-only mode |
| core.coordinates | boundary conversions | integer intervals | normalized intervals/overlap | standard library | start/end | ValueError | single residue, touching regions, invalid coordinates |
| core.config | load reproducible settings | YAML/path or defaults | resolved AppConfig | PyYAML/Pydantic | paths, threads, BCTD | ConfigurationError | unknown keys, invalid values, relative paths |
| cli | user entry point | arguments/files | version/config/quality JSON and prepared files | Typer/core/io | --config, input format and paths | exit 2 for invalid input/config/quality errors | help, configuration, prepare/validate integrity |

Records use frozen Pydantic models, but nested dictionaries are not deeply immutable: consumers
must not mutate shared qualifiers/evidence. JSON serialization is supported. Derived fields
(length and CTD distance) are emitted for users; use `exclude_computed_fields=True` for validated
JSON round trips. Batch 3 exporters will explicitly define storage schemas rather than infer them
from empty tables or dump arbitrary nested objects into Parquet.

## Planned interfaces (not implemented)

- Scanner: protein records + validated model/config → raw DomainHits + ScanStatus + artifacts.
- Future classifier: Evidence Profiles → separately defined biological classifications; no such scoring exists yet.
- Exporter: dataset + hits + candidates + run metadata → separate typed tables and FASTA/JSON.
- RN28/CBD: proteins → model-specific hits, depending only on core.
- Candidate-to-neighborhood association: candidate evidence + objective CDS neighborhoods → context evidence (future Batch 4B).
- Motility: genome/metadata → species evidence. Classifier: all evidence → ranked candidates.
- Structure: candidate proteins → structural evidence. Viewers: records/evidence → visual output.

Future modules are documented rather than represented by empty packages. BioPython is used by batch 2. PyHMMER powers the HMM engine. MMseqs2 is implemented in scanners.mmseqs; Parquet remains future work.

PorP/SprF dual-model orchestration is implemented in scanners.sprf and calls the existing HMM
engine independently. It preserves model-specific results, partial failures and per-protein evidence;
see sprf_scanner.md. No substrate/shuttle pairing or score fusion is performed.

Shared pair orchestration now lives in scanners.pair; scanners.sprf and scanners.bctd
select their canonical model pairs and interpretation. scanners.repeats adds per-protein PF13573
domain counts on top of run_hmm without replacing the search algorithm or filtering raw evidence.

Batch 4A implements an independent neighborhood package over core records, without scanner imports.
See neighborhood.md for index reuse, compound/circular geometry and CLI contracts.

Batch 5A adds scanners.mmseqs: external CPU MMseqs2, managed database manifests, and
SequenceAlignment as a DomainHit subclass. No shared HMM/neighborhood behavior is changed.
See mmseqs_engine.md for the database/CLI contracts and retention limits.


## Unified Reference Registry

`sprb_scanner.references` imports core/io and reads existing artifact bytes; it imports no scanners.
It owns typed reference identities, version/checksum-pinned derivation DAGs, offline resolution and
format inspection (protein sequences/sets, aligned FASTA/Stockholm/A3M, HMM, JSON/TSV catalogs).
`ResolvedReference.as_profile()` adapts verified HMMs to the unchanged core ProfileConfig;
`sequences()` exposes existing ProteinRecord IDs without changing the search algorithms.
Core does not depend on the Registry. The CLI adds only listing/validation/resolution/export adapters.

`references/registry.json` is authoritative; `references/migration.json` records exact byte migration.
Legacy `profiles/` and `reference_data/` are fixed compatibility snapshots, checked in CI against
canonical resources. External model absence is explicit and never triggers download or a negative
biological conclusion. Future MSA-to-HMM construction can publish new versioned resources with
pinned `derived_from`; no such building workflow exists in this batch.
See [resource contracts](../references/README.md) for schema, errors, supported formats and tests.


## Unified evidence integration (Batch 6A)

`evidence.inputs` verifies existing paired/native HMM and MMseqs2 exports against prepared protein
identities, sequence fingerprints and manifests. Failed sources are explicit data states;
corrupt/incompatible exports are errors. `evidence.merger` unions all retained hits into typed
Evidence Profiles, reuses per-domain architecture evidence and queries NeighborhoodIndex for the union.
`evidence.service` orchestrates offline loading and atomic JSON/TSV/provenance exports.
No executable search, threshold change, score fusion, coverage aggregation or functional inference occurs.
Core, scanners and neighborhood remain unchanged; CLI adds a thin `merge-evidence` adapter.
See [evidence merger contracts](evidence_merger.md) for Purpose, Input, Output, Dependencies,
Parameters, Error Handling, Test Cases and real-data acceptance.


## Single-genome orchestration (Batch 6B)

`pipeline.config` owns the strict versioned YAML contract. `pipeline.identity` builds per-stage
content/reference/parameter/backend identities. `pipeline.storage` provides canonical hashes,
atomic checkpoints, immutable attempt directories and an exclusive workspace lock.
`pipeline.runner` invokes the existing four scanner service APIs and Batch 6A Merger, verifies
native outputs on both execution and reuse, and records dependency bindings and run history.
Pipeline policy belongs here; core/scanner algorithms and neighborhood logic remain unchanged.

The Merger gains backward-compatible API arguments for explicit failed/not-run source states
when a pipeline stage fails before producing native output, plus RN28 resource version/external
artifact selection. These adapters preserve the existing typed evidence/status semantics;
they cannot inject success or override a supplied successful scanner directory.
See [pipeline runner contracts](pipeline_runner.md) for configuration, storage, resume and validation.
