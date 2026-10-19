# Data schema 1.0.0

| Record | Identity and key fields |
|---|---|
| GenomeRecord | genome_id, accession, organism, taxon_id, input_format, source_files |
| ContigRecord | contig_id, genome_id, source_id, length, topology |
| GeneRecord | gene_id, genome_id, contig_id, ordered location_parts, strand, source_id, locus_tag, gene_name, product, partial, pseudogene, transl_table, codon_start, qualifiers |
| ProteinRecord | protein_id, genome_id, nullable gene_id, source_id, sequence, derived length |
| DomainHit | hit_id, protein_id, method, model provenance, query coordinates/length, target coordinates/length, typed E-value, bit score, coverage, CTD distance, backend version, raw artifact |
| CandidateRecord | candidate_id, protein_id, all domain_hits, support category, extensible evidence references |
| ScanStatus | model_id, not_run/succeeded/failed/unavailable, explanation |
| GenomeDataset | schema_version, genome, contigs, genes, proteins, quality_issues |

Internal ID components are percent-escaped before joining with `::`; source IDs remain intact.
Examples: `GCF_x.1::NC_x::CDS_000182`, then `::P1` for its protein.
Parsers will number CDS deterministically in a documented contig/coordinate order, tie-breaking
with source order. Local genome IDs must be supplied or content-derived; file basenames alone
are not globally unique. Source protein IDs may repeat; internal IDs may not.

Protein-only records have no fabricated genes/contigs. Missing translations and pseudogenes may
have a GeneRecord without any ProteinRecord. Protein alphabet accepts canonical/ambiguous and
rare residues. Parsers must report normalization (e.g. terminal-stop removal), invalid residues
and missing sequences as quality issues instead of silently repairing biological content.

Structural validation covers ID uniqueness, genome ownership, contig bounds, gene/protein links,
query and target hit bounds, finite statistics and candidate/protein identity.
Schema validation does not establish biological correctness or classify SprB-like proteins.
`support` is evidence grouping only; its calculation belongs to the future merger.

Model provenance requires source, version and SHA-256. HMM domain independent/conditional
E-values, sequence E-values and MMseqs2 alignment E-values remain explicitly distinguishable.
Run-level statuses are separate from candidates, so failure/unavailability never implies a negative.
Future run metadata must record schema/software/backend versions, complete resolved configuration,
input/model checksums, timestamps, commands and raw artifact locations.

HMM Engine adds optional `DomainHit.hmm: HMMStatistics` for backward-compatible records.
New engine hits always populate separate sequence/independent-domain/conditional-domain E-values,
sequence score, envelope and reporting/inclusion flags. See [hmm_engine.md](hmm_engine.md).
