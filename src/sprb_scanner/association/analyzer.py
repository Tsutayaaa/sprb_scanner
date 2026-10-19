"""Join independent model evidence to an existing CDS neighborhood index."""

import math
from collections import defaultdict

from sprb_scanner.association.models import (
    AssociationResult,
    CargoAssociation,
    DomainEvidence,
    LocalAssociation,
)
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.neighborhood import NeighborhoodIndex

CTD_MODELS = ("TIGR04131", "PF13585")
SPRF_MODELS = ("TIGR03519", "PF11751")


def domain_evidence(hits):
    """Never merge domains or sum overlapping/cross-model coverage fractions."""
    grouped = defaultdict(list)
    for hit in hits:
        grouped[hit.protein_id, hit.model.model_id].append(hit)
    result = []
    for key in sorted(grouped):
        group = sorted(grouped[key], key=lambda h: (h.start, h.end, h.hit_id))
        for ordinal, hit in enumerate(group, 1):
            result.append(
                DomainEvidence(
                    hit=hit,
                    domain_order_by_protein_start=ordinal,
                    domains_same_protein_model=len(group),
                    overlapping_same_model_hit_ids=tuple(
                        h.hit_id
                        for h in group
                        if h.hit_id != hit.hit_id and max(h.start, hit.start) < min(h.end, hit.end)
                    ),
                )
            )
    return tuple(result)


def associate(dataset, ctd_hits, sprf_hits, window=5):
    """Any retained CTD domain creates a cargo, including below-GA-only evidence."""
    if type(window) is not int or window < 0:
        raise InputError("window must be a nonnegative integer")
    if not dataset.genes or not dataset.contigs:
        raise InputError("Association requires an annotated genome with CDS coordinates")
    proteins = {p.protein_id: p for p in dataset.proteins}
    ctd_hits, sprf_hits = tuple(ctd_hits), tuple(sprf_hits)
    hits = ctd_hits + sprf_hits
    if len({h.hit_id for h in hits}) != len(hits):
        raise InputError("Duplicate domain hit IDs")
    for group, models in ((ctd_hits, CTD_MODELS), (sprf_hits, SPRF_MODELS)):
        for h in group:
            p = proteins.get(h.protein_id)
            if p is None or h.protein_length != p.length or p.gene_id is None:
                raise InputError(f"Hit protein/genome mismatch: {h.hit_id}")
            if h.model.model_id not in models or h.method != "hmm" or h.hmm is None:
                raise InputError(f"Unexpected HMM model/statistics: {h.hit_id}")
            if h.hmm.ga is None:
                raise InputError(f"Missing GA evidence: {h.hit_id}")
            if not math.isclose(h.query_coverage, (h.query_end - h.query_start) / h.query_length):
                raise InputError(f"HMM coverage disagrees with coordinates: {h.hit_id}")
            if not math.isclose(h.target_coverage, (h.end - h.start) / h.protein_length):
                raise InputError(f"Protein coverage disagrees with coordinates: {h.hit_id}")
    by_ctd, by_sprf, gene_proteins = defaultdict(list), defaultdict(list), defaultdict(list)
    for h in ctd_hits:
        by_ctd[h.protein_id].append(h)
    for h in sprf_hits:
        by_sprf[h.protein_id].append(h)
    for p in dataset.proteins:
        gene_proteins[p.gene_id].append(p)
    genes = {g.gene_id: g for g in dataset.genes}
    index = NeighborhoodIndex(dataset)
    cargos = []
    for protein_id in sorted(by_ctd):
        p = proteins[protein_id]
        neighbors = []
        for relation in index.query(p.gene_id, window).neighbors:
            for neighbor in gene_proteins[relation.gene.gene_id]:
                if neighbor.protein_id in by_sprf:
                    neighbors.append(
                        LocalAssociation(
                            protein_id=neighbor.protein_id,
                            source_protein_id=neighbor.source_id,
                            protein_length=neighbor.length,
                            relationship=relation,
                            domain_hit_ids=tuple(h.hit_id for h in by_sprf[neighbor.protein_id]),
                        )
                    )
        cargos.append(
            CargoAssociation(
                gene=genes[p.gene_id],
                protein_id=protein_id,
                source_protein_id=p.source_id,
                protein_length=p.length,
                domain_hit_ids=tuple(h.hit_id for h in by_ctd[protein_id]),
                hit_models=tuple(sorted({h.model.model_id for h in by_ctd[protein_id]})),
                ga_models=tuple(
                    sorted(
                        {
                            h.model.model_id
                            for h in by_ctd[protein_id]
                            if h.hmm.ga.status == "passed_ga"
                        }
                    )
                ),
                status="local_hit" if neighbors else "no_local_hit",
                neighbors=tuple(neighbors),
            )
        )
    return AssociationResult(
        genome_id=dataset.genome.genome_id,
        window=window,
        cargos=tuple(cargos),
        domains=domain_evidence(hits),
    )
