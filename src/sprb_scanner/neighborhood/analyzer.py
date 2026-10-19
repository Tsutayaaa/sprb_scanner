"""Preindexed CDS neighborhoods; never scan the genome again for a query."""

from dataclasses import dataclass

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GeneRecord, GenomeDataset
from sprb_scanner.neighborhood.models import NeighborhoodResult, NeighborRecord


@dataclass(frozen=True)
class IndexedCDS:
    gene: GeneRecord
    parts: tuple
    anchor: int
    crosses_origin: bool


def index_gene(gene, contig):
    # location_parts are in biological extraction order, including on the negative strand.
    forward = gene.location_parts if gene.strand != -1 else tuple(reversed(gene.location_parts))
    parts = tuple(sorted(gene.location_parts, key=lambda p: (p.start, p.end)))
    if any(a.end > b.start for a, b in zip(parts, parts[1:], strict=False)):
        raise InputError(f"Overlapping location parts within CDS: {gene.gene_id}")
    drops = sum(b.start < a.start for a, b in zip(forward, forward[1:], strict=False))
    circular = contig.topology == "circular"
    if drops and not circular:
        raise InputError(f"Origin-spanning or unordered CDS on non-circular contig: {gene.gene_id}")
    if drops > 1:
        raise InputError(f"CDS makes multiple origin traversals: {gene.gene_id}")
    crosses = bool(drops)
    if crosses and forward[-1].end > forward[0].start:
        raise InputError(f"CDS span exceeds one circle: {gene.gene_id}")
    return IndexedCDS(gene, parts, forward[0].start, crosses)


def geometry(a, b, circle_length=None):
    """Exact union-of-parts intersection and minimum edge gap; touching means zero."""
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
        return -overlap, overlap, False
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
