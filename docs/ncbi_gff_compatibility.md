# Batch 6D — NCBI GFF compatibility

## Purpose / inputs / dependencies

Support NCBI biological CDS locations while retaining strict annotation validation. Inputs remain the same
annotation-release GFF3 + FNA + FAA; no source files, CDS coordinates, FAA proteins, scanner thresholds or
reference sequences are repaired or rewritten. Dependencies are the existing Biopython/Pydantic stack.

NCBI documents overlapping CDS rows and micro-introns for both programmed slippage and putative assembly
errors; an overlap alone is not a diagnosis. The parser follows the
[NCBI GFF3 conventions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).

## Join validation / parameters / errors

- Multiple rows with one CDS ID become one GeneRecord and at most one scan ProteinRecord. Distinct CDS
  copies referencing the same exact FAA ID remain separate, as before; duplicate FAA IDs are ambiguous errors.
- Preserve original intervals, biological part order, strand, phase, row attributes, ID and exceptions.
  Genomic interval union is used **only** for neighborhood geometry; translation uses biological extraction
  parts, retaining repeated bases at an overlap and omitting annotated gaps.
- Explicit `exception=ribosomal slippage` (or an empty INSDC-style flag) permits progressive joins. Full,
  known-strand CDS without `transl_except` must reconstruct the exact supplied FAA protein. Phase
  discrepancies remain visible warnings only for explicit slippage with independently checked FAA support.
- Partial/unknown-strand/translation-exception slippage remains FAA-canonical with explicit
  `not_verifiable` provenance and a warning; this is not asserted to be a verified genomic reconstruction.
- Without a slippage annotation, accept only progressive **1–2 bp overlaps**, consistent phase, an explicit
  exact `protein_id` mapping, and complete joined translation matching FAA, including valid start/stop.
  Partial, ambiguous or longer overlaps are rejected. This conservative structural bound is not a scientific
  search cutoff and is not applied to HMM/MMseqs2 hits.
- Known `unclassified translation discrepancy` / `annotated by transcript or proteomic data` exceptions
  can use the same strict FAA verification, but are recorded separately. Do not relabel them as slippage
  or assert that an assembly error occurred. Unsupported/conflicting exceptions remain errors.
- Ordinary non-overlapping compound CDS/gaps retain their existing phase/translation checks. No global
  strict-quality bypass was introduced. Duplicate/contained/non-progressive segments or multiple circular
  traversals fail even when an exception is present.

FAA is canonical after exact identifier mapping; the existing documented `cds-` fallback is retained for
ordinary/explicitly annotated CDS. No fuzzy protein-ID matching, regenerated FAA or sequence substitution.

## Circular locations and neighborhood review

Circularity must be explicit: GFF `region;Is_circular=true` or FNA header `[topology=circular]` (also the
unbracketed `topology=circular` token). `linear` declarations are honored; conflicting metadata is an error.
`complete genome`, sequence length or out-of-bounds coordinates alone do not establish circularity.

Landmark region coordinates must match FNA length. An origin-spanning feature may end in virtual space
if its start lies within the contig and its span does not exceed one circle. This applies to RNA and other
feature rows as well as CDS/parents. Normalize CDS into in-bounds 0-based half-open circular segments;
never clip the virtual end or duplicate the protein. Negative-strand segment extraction reverses biological
order correctly. Region metadata may appear later in the file.

For multipart children, explicit `part=X/Y` has precedence; otherwise a virtual-space parent supplies the
circular anchor before strand-aware sorting. Without explicit part/parent information, no circular ordering
is inferred from FAA alone. Source virtual intervals, including parent intervals used for ordering, remain
in provenance. Pseudogene CDS retains source intervals even without an exported scan protein.

`core.locations.progressive_path` shares the former NeighborhoodIndex structural checks. `geometry()`,
union, coordinate/transcription directions, tie policy, nonnegative `distance_bp`, `overlap_bp`,
`signed_gap_bp`, adjacency by CDS rank and ±window semantics are unchanged. NeighborhoodIndex additionally
recognizes only a parser FAA-verification proof bound to the exact coordinates/strand for unannotated joins.
It never guesses slippage from geometry.

## Outputs / provenance

Prepared dataset schema and existing CLI remain compatible. Per-part attributes/phases/rows remain in
qualifiers. New accommodated CDS carry `gff_location_validation` JSON containing policy version,
representation, uncertainty, source intervals, normalized-location SHA-256, canonical protein SHA-256,
FAA ID, phase/translation check and topology source. Origin-spanning CDS also carry `gff_source_intervals`,
including those without proteins. Reserved proof attributes cannot be supplied in the source GFF.

Ordinary errors, missing FAA, unmapped FAA, incompatible coordinates or malformed IDs still fail strict
preparation. Partial/unverifiable and conservative unannotated joins retain explicit warnings. Software
fingerprints include the shared location helper; old execution plans are not silently resumed under new code.

## Versioned HPC recovery / scientific collection

Create a **new** plan selecting only Batch 6C's 195 failed assembly accessions. Preserve the old manifest,
prepared datasets, successful attempts and collection checksums. Use an isolated code deployment and pinned
resource/config files; do not modify the old installation while it validates legacy plans. Recheck all original
inputs and validate old successful outputs read-only. Reparse the 698 successful inputs for parser regression
without rerunning their scientific searches. Record all newly surfaced failures after resolving the first
blocker, rather than dropping another feature or weakening strict preparation.

The recovery plan has its own stable sorted indexes. Save a mapping to original array indexes and both
manifest SHA-256 values. After only recovered tasks run, collect that plan using its original collector.
Combine immutable collection snapshots with:

```bash
python scripts/hpc/combine_collections.py \
  --collection /project/collections/bacteroidota-893-v1 \
  --collection /project/collections/bacteroidota-recovery-195-v1 \
  --output /project/collections/bacteroidota-combined-v2
```

The Python service is `sprb_scanner.batch.combine.combine_collections`. Every combined TSV/JSONL row
retains its native fields and gains `collection_manifest_sha256`, resolving its original array/software
identity. Metadata records each source snapshot's checksums, immutable manifest, software fingerprints,
references and superseded unavailable rows. Only prior unavailable genomes may be superseded; previously
successful or partial scientific results cannot be replaced. Different source files, corrupt artifacts,
provenance/schema mismatches and concurrent changes fail explicitly. Failed recovery rows stay visible
with unknown counts blank. Combining snapshots performs no scientific searches and never rewrites inputs.

## Test cases

`tests/unit/test_ncbi_gff_locations.py` covers both strands, 1/2 bp overlaps/gaps, exact FAA mapping,
uncertainty, phase, real small NCBI CDS examples, parent-based circular order, virtual RNA rows, topology
conflicts, failed strict preparation, prepared roundtrip and neighborhood geometry.
`tests/unit/test_collection_migration.py` checks raw multi-hit retention, software provenance, old-output
immutability and rejection of duplicate success, incompatible inputs or corrupt tables.
Run existing parser/neighborhood/FJ tests and the full pytest/Ruff/CLI/package checks. HPC acceptance records
195 recovery outcomes, protein/CDS mapping, further blockers, old checksums and the unchanged FJ reference.
