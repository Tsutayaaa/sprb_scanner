# Batch 6C MVP — HPC validation

Validation date: 2026-10-24. Deployed scientific/execution code: `c3bdd4a20fdfecedc87ae5abe55e1d8e42895dbc`.
See [SLURM usage and contracts](slurm_batch.md) for purpose, inputs, outputs, parameters, errors and tests.

## Environment and reproducibility

- XJTLU Slurm 23.02.2, Linux x86_64; searches ran only in compute allocations.
- Project-local environment: Python 3.12.9, PyHMMER 0.12.3, Biopython 1.88, MMseqs2 14-7e284.
- Five HMM artifacts matched the existing Reference Registry SHA-256 values. RN28 and GA thresholds unchanged.
- Explicit resources: partition `cpu6348`, account `chunchan_igem`, QOS `52cores`, 4 CPU / 8G / 30 min per task;
  concurrency 6, at most 40 array entries per wave. These are site-specific settings, not universal defaults.
- Immutable 893-entry manifest SHA-256:
  `93c99477eed2e5a4b2d6a8bf1a6344845c07ac037303a69b233a06b52a154c7a`.
- All 2,679 source files (7,404,933,695 bytes) were rehashed unchanged after execution.
- Jobs: preflight `3079546`, inventory `3079604`, controller `3079607`, annotation diagnostics `3080964`,
  independent final acceptance `3080683`. Acceptance state: `passed`.

## Outcome and scope

| Measure | Result |
|---|---:|
| Assemblies inventoried and preparation attempted | 893 |
| Prepared and completed with all evidence sources | 698 |
| Explicit preparation rejections | 195 |
| Protein Evidence Profiles | 16,125 |
| Proteins with any retained Type B CTD domain evidence | 10,509 |
| Local neighborhood association rows | 4,721 |
| RN28 alignments | 372 |

The 195 preparation rejections comprise 165 explicitly annotated ribosomal-slippage overlapping CDS and
30 circular-origin-spanning GFF representations beyond the sequence-region boundary. The first fatal feature
and its qualifiers were inspected for each rejected assembly; this is not an exhaustive audit of every CDS.
These are limitations of the current GFF adapter, not evidence of assembly damage. The owner explicitly chose
to defer parser adaptation to the next batch. No coordinates, proteins or annotations were repaired or deleted.
Rejected assemblies remain in the summary with unavailable evidence and blank unknown counts; they are
excluded from the denominator of successful analyses. All 698 prepared assemblies completed their pipelines.

Native domain/alignment rows and fields were checked against collected evidence and subtables, without
collapsing hits or changing GA/coverage semantics. Raw retained and below-GA records are not confirmed family
members. `no_local_hit` does not establish orphan status; PorP/SprF evidence does not prove shuttle function.

## Reference and resume checks

The HPC FJ assembly is `GCF_034479105.1` (5,134 prepared proteins), distinct from the earlier `NC_009441.1`
regression assembly. All 12 literature-reference protein source IDs were recovered: nine have directly adjacent
PorP/SprF-like evidence and three have `no_local_hit`. SprB `WP_012023065.1` (6,497 aa) retains the RN28 match
at `[681,888)`, 207/207 identical residues, query coverage 1.0 and full-protein coverage 0.0318608589.
SprB has two directly adjacent retained neighbor candidates; this is evidence, not two established shuttles.

A separate FJ resume reused all five verified stages: Type B CTD, PorP/SprF, repeat, RN28 and evidence.
The pilot also exercised failed-task isolation and retry with original array indexes. Local integration tests
cover checksum corruption, unchanged native evidence, scanner results and verified reuse.

## Retained disk measurement

Decimal GB (10^9 bytes), measured after completion; these are current retained sizes, not peak usage.

| Scope | Apparent GB | Allocated GB |
|---|---:|---:|
| Full batch workspace | 14.870 | 16.066 |
| Scientific collection | 1.518 | 1.518 |
| Entire project, including environment/pilot/logs | 16.634 | 17.969 |

No cleanup or archival mechanism was added. Existing scanner temporary-file lifecycles are unchanged.

## Artifacts and checks

On HPC, under `/gpfs/work/bio/shuleihe23/sprb_scanner/`:

- `batches/bacteroidota-893-v1/`: immutable manifest/config, prepared data, attempts, independent scanner results,
  final status, submissions/logs and per-genome disk measurement.
- `collections/bacteroidota-893-v1/`: genome/protein/cargo/association/architecture/RN28/source tables,
  all domain records, provenance, output SHA-256 and collection status.
- `reports/batch6c-893-v1/`: `acceptance_summary.json` and `failures.tsv`.
- `reports/preparation-diagnostics-v1/`: original first-fatal-feature annotations and category summary.

Large results, genomes, models and databases are not committed. Seven key scientific tables were copied locally
and checked against the original collection SHA-256; full evidence/provenance remains on HPC.

Local full pytest: **435 passed**. Ruff check/format, CLI and wheel/sdist build passed. GitHub Actions run
`40077723463` passed Python 3.11, 3.12 and 3.13 on the deployed code commit. Documentation-only follow-up
commits do not change the deployed code identity; their CI is checked separately on the PR.
