# Coordinate conventions

All internal genomic and protein intervals are **0-based, half-open** `[start, end)`.
Intervals are nonempty; length = end - start. Start is never negative.
GFF3 and HMMER 1-based inclusive `[a,b]` become `[a-1,b)` at input boundaries.
Output conversion is explicit; never infer a coordinate system from a filename.
For example `[1,1]` becomes `[0,1)`, and `[101,200]` becomes `[100,200)`.
Adjacent `[0,10)` and `[10,20)` have zero overlap.

Gene `location_parts` are in biological extraction order, preserving joined CDS and
origin-spanning circular CDS. Every part uses ascending coordinates, even on the minus strand.
Strand is 1, -1 or 0 (unknown). Do not collapse compound locations to a min/max envelope.
Parsers must preserve fuzzy-location details in qualifiers and set partial appropriately.
GFF phase/codon_start and translation table affect translation, not genomic interval conversion.
Future parser tests must cover reverse strands and circular origin crossings explicitly.

DomainHit start/end are target protein residues; query_start/query_end are HMM match positions
or reference sequence residues. Both use this same convention. CTD distance = protein_length - end.
Coverage values are backend-reported aligned coverage fractions and must document their precise
backend definition; they must not become full-protein-length exclusion filters.
