"""Preindexed CDS neighborhoods; never scan the genome again for a query."""

from dataclasses import dataclass

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.locations import faa_verified_overlap, progressive_path
from sprb_scanner.core.models import GeneRecord, GenomeDataset, Interval
from sprb_scanner.neighborhood.models import NeighborhoodResult, NeighborRecord

GEOMETRY_POLICY = (
    "cached genomic interval union for annotated slippage or exactly FAA-verified CDS; "
    "original biological location_parts and qualifiers preserved"
)


@dataclass(frozen=True)
class IndexedCDS:
    gene: GeneRecord
    parts: tuple
    anchor: int
    crosses_origin: bool


def genomic_union(parts):
    """Return disjoint genomic intervals; never change biological extraction parts."""
    merged = []
    for part in sorted(parts, key=lambda p: (p.start, p.end)):
        if merged and part.start <= merged[-1].end:
            previous = merged.pop()
            merged.append(Interval(start=previous.start, end=max(previous.end, part.end)))
        else:
            merged.append(part)
    return tuple(merged)


def index_gene(gene, contig):
    # Keep biological order for provenance; reverse only the coordinate-forward view.
    ordered = tuple(sorted(gene.location_parts, key=lambda p: (p.start, p.end)))
    overlapping = any(a.end > b.start for a, b in zip(ordered, ordered[1:], strict=False))
    if overlapping:
        # INSDC flag or the exact NCBI GFF3 exception, never free-text notes.
        flag = "ribosomal_slippage" in gene.qualifiers
        flag_valid = flag and all(not v.strip() for v in gene.qualifiers["ribosomal_slippage"])
        exceptions = tuple(v.strip().casefold() for v in gene.qualifiers.get("exception", ()))
        if flag and not flag_valid:
            raise InputError(f"Invalid value for /ribosomal_slippage: {gene.gene_id}")
        verified = faa_verified_overlap(gene)
        if not (flag_valid or "ribosomal slippage" in exceptions or verified):
            raise InputError(
                "Overlapping location parts require ribosomal slippage annotation "
                f"or exact FAA verification: {gene.gene_id}"
            )
        if not verified and any(value and value != "ribosomal slippage" for value in exceptions):
            raise InputError(f"Unexplained exception on overlapping CDS: {gene.gene_id}")
    try:
        forward, crosses = progressive_path(
            gene.location_parts, gene.strand, contig.length, contig.topology == "circular"
        )
    except InputError as exc:
        raise InputError(f"{exc}: {gene.gene_id}") from exc
    return IndexedCDS(gene, genomic_union(ordered), forward[0].start, crosses)


def geometry(a, b, circle_length=None):
    """Return (nonnegative gap, coding overlap bp, shortest-gap origin crossing).

    Inputs must be nonempty, coordinate-sorted, internally disjoint Interval tuples
    (NeighborhoodIndex caches each validated CDS genomic union at construction).
    Swapping a and b preserves all three values. Overlap and touching both have
    distance zero; overlap is distinguished by its positive intersection length.
    """
    i = j = overlap = 0
    best = None
    while i < len(a) and j < len(b):
        x, y = a[i], b[j]
        overlap += max(0, min(x.end, y.end) - max(x.start, y.start))
        gap = max(0, x.start - y.end, y.start - x.end)
        best = gap if best is None else min(best, gap)
        if x.end <= y.end:
            i += 1
        else:
            j += 1
    if overlap:
        return 0, overlap, False
    crossed = False
    if circle_length is not None:
        wrapped = min(
            circle_length - a[-1].end + b[0].start, circle_length - b[-1].end + a[0].start
        )
        if wrapped < best:
            best, crossed = wrapped, True
    return best, 0, crossed


class NeighborhoodIndex:
    """O(N log N) CDS indexing; O(k) selection plus returned location-part processing."""

    def __init__(self, dataset: GenomeDataset):
        self.genome_id = dataset.genome.genome_id
        self.contigs = {c.contig_id: c for c in dataset.contigs}
        grouped = {key: [] for key in self.contigs}
        for gene in dataset.genes:
            grouped[gene.contig_id].append(index_gene(gene, self.contigs[gene.contig_id]))
        self.by_contig = {}
        self.by_gene = {}
        for contig_id, entries in grouped.items():
            entries.sort(key=lambda x: (x.anchor, x.parts[-1].end, x.gene.gene_id))
            self.by_contig[contig_id] = tuple(entries)
            for position, entry in enumerate(entries):
                self.by_gene[entry.gene.gene_id] = (contig_id, position)

    def query(self, gene_id: str, window: int = 5):
        if isinstance(window, bool) or not isinstance(window, int) or window < 0:
            raise InputError("window must be a non-negative integer")
        if gene_id not in self.by_gene:
            raise InputError(f"Unknown CDS gene_id: {gene_id}")
        contig_id, position = self.by_gene[gene_id]
        contig, entries = self.contigs[contig_id], self.by_contig[contig_id]
        target = entries[position]
        n = len(entries)
        seen, neighbors = {gene_id}, []
        circular = contig.topology == "circular"
        # Visit nearest CDS ranks first; coordinate-upstream wins equal-rank circle ties.
        for steps in range(1, min(window, n - 1) + 1):
            for direction in (-1, 1):
                raw = position + direction * steps
                if not circular and not 0 <= raw < n:
                    continue
                entry = entries[raw % n]
                if entry.gene.gene_id in seen:
                    continue
                seen.add(entry.gene.gene_id)
                distance, overlap, distance_crosses = geometry(
                    target.parts, entry.parts, contig.length if circular else None
                )
                trans_offset = (
                    direction * steps * target.gene.strand if target.gene.strand else None
                )
                index_crosses = circular and not 0 <= raw < n
                neighbors.append(
                    NeighborRecord(
                        gene=entry.gene,
                        genomic_position="upstream" if direction < 0 else "downstream",
                        relative_position=(
                            None
                            if trans_offset is None
                            else "upstream"
                            if trans_offset < 0
                            else "downstream"
                        ),
                        genomic_cds_offset=direction * steps,
                        transcription_cds_offset=trans_offset,
                        intervening_cds=steps - 1,
                        distance_bp=distance,
                        signed_gap_bp=distance - overlap,
                        overlap_bp=overlap,
                        directly_adjacent=steps == 1,
                        same_strand=(
                            target.gene.strand == entry.gene.strand
                            if target.gene.strand and entry.gene.strand
                            else None
                        ),
                        overlaps=overlap > 0,
                        crosses_origin=(
                            index_crosses
                            or distance_crosses
                            or target.crosses_origin
                            or entry.crosses_origin
                        ),
                        index_crosses_origin=index_crosses,
                        distance_crosses_origin=distance_crosses,
                        gene_crosses_origin=entry.crosses_origin,
                    )
                )
            if len(seen) == n:
                break
        return NeighborhoodResult(
            genome_id=self.genome_id,
            contig=contig,
            target=target.gene,
            target_crosses_origin=target.crosses_origin,
            window=window,
            neighbors=tuple(neighbors),
        )
