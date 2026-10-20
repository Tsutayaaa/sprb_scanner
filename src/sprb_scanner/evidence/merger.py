"""Union retained evidence and query the existing neighborhood index once per candidate."""

import hashlib
from collections import defaultdict

from sprb_scanner.association.analyzer import SPRF_MODELS, domain_evidence
from sprb_scanner.association.models import LocalAssociation
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.evidence.models import (
    SOURCES,
    EvidenceProfile,
    EvidenceResult,
    NeighborhoodEvidence,
    ProteinEvidence,
)
from sprb_scanner.neighborhood import NeighborhoodIndex


def merge(dataset, runs, hmm_hits=(), alignments=(), window=5, neighborhood=True):
    """Any retained scanner hit includes a protein; never combine scores or coverage."""
    if type(window) is not int or window < 0:
        raise InputError("window must be a nonnegative integer")
    if set(runs) != set(SOURCES) or any(k != r.source_id for k, r in runs.items()):
        raise InputError("All six source run states must be provided")
    proteins = {p.protein_id: p for p in dataset.proteins}
    genes = {g.gene_id: g for g in dataset.genes}
    grouped, sprf = defaultdict(lambda: defaultdict(list)), defaultdict(list)
    seen = set()
    hmm_hits, alignments = tuple(hmm_hits), tuple(alignments)
    if any(h.method != "hmm" or h.hmm is None or h.hmm.ga is None for h in hmm_hits):
        raise InputError("HMM evidence requires HMM statistics and GA state")
    if any(h.method != "mmseqs2" for h in alignments):
        raise InputError("RN28 evidence requires sequence alignments")
    for h in (*hmm_hits, *alignments):
        source = h.model.model_id if h.method == "hmm" else "RN28"
        p = proteins.get(h.protein_id)
        if source not in runs or runs[source].state != "succeeded":
            raise InputError(f"Evidence attached to unsuccessful/unknown source: {source}")
        if p is None or h.protein_length != p.length:
            raise InputError(f"Evidence protein mismatch: {h.hit_id}")
        if (source, h.hit_id) in seen:
            raise InputError(f"Duplicate evidence ID: {source}/{h.hit_id}")
        seen.add((source, h.hit_id))
        grouped[p.protein_id][source].append(h)
        if source in SPRF_MODELS:
            sprf[p.protein_id].append(h)
    for source, run in runs.items():
        count = sum(len(e.get(source, ())) for e in grouped.values())
        if run.state == "succeeded" and run.retained_hits != count:
            raise InputError(f"Retained hit count disagrees with source: {source}")
    good_sprf = tuple(m for m in SPRF_MODELS if runs[m].state == "succeeded")
    index = None
    neighborhood_error = None
    if neighborhood and good_sprf and dataset.genes and dataset.contigs:
        try:
            index = NeighborhoodIndex(dataset)
        except InputError as exc:
            neighborhood_error = str(exc)
    if not neighborhood:
        context_state, context_message = "not_run", "Neighborhood analysis disabled"
    elif not dataset.genes or not dataset.contigs:
        context_state, context_message = "unavailable", "No genomic CDS coordinates"
    elif not good_sprf:
        context_state = (
            "failed" if any(runs[m].state == "failed" for m in SPRF_MODELS) else "not_run"
        )
        context_message = "No successful PorP/SprF source available"
    elif neighborhood_error:
        context_state, context_message = "failed", neighborhood_error
    else:
        context_state = "succeeded" if len(good_sprf) == len(SPRF_MODELS) else "partial"
        context_message = (
            None if context_state == "succeeded" else "Only one PorP/SprF model succeeded"
        )
    gene_proteins = defaultdict(list)
    for p in dataset.proteins:
        gene_proteins[p.gene_id].append(p)
    profiles = []
    for pid in sorted(grouped):
        p = proteins[pid]
        sources = {}
        for source in SOURCES:
            run, hits = runs[source], grouped[pid].get(source, ())
            state = (
                run.state
                if run.state != "succeeded"
                else ("succeeded_with_hits" if hits else "succeeded_no_hits")
            )
            domains = domain_evidence(hits) if source != "RN28" else ()
            sources[source] = ProteinEvidence(
                source_id=source,
                state=state,
                domains=tuple(
                    e.model_copy(update={"scan_directory": run.directory}) for e in domains
                ),
                alignments=tuple(hits) if source == "RN28" else (),
            )
        common = {"window": window, "sprf_sources": good_sprf}
        if context_state in {"not_run", "failed", "unavailable"}:
            context = NeighborhoodEvidence(
                state=context_state, outcome="unknown", **common, message=context_message
            )
        elif p.gene_id is None:
            context = NeighborhoodEvidence(
                state="unavailable",
                outcome="unknown",
                **common,
                message="No genomic CDS coordinates",
            )
        else:
            neighbors = []
            for relation in index.query(p.gene_id, window).neighbors:
                for neighbor in gene_proteins[relation.gene.gene_id]:
                    if neighbor.protein_id in sprf:
                        neighbors.append(
                            LocalAssociation(
                                protein_id=neighbor.protein_id,
                                source_protein_id=neighbor.source_id,
                                protein_length=neighbor.length,
                                relationship=relation,
                                domain_hit_ids=tuple(h.hit_id for h in sprf[neighbor.protein_id]),
                            )
                        )
            complete = len(good_sprf) == len(SPRF_MODELS)
            context = NeighborhoodEvidence(
                state="succeeded" if complete else "partial",
                outcome="local_hit" if neighbors else "no_local_hit" if complete else "unknown",
                neighbors=tuple(neighbors),
                **common,
                message=None if complete else "Only one PorP/SprF model succeeded",
            )
        profiles.append(
            EvidenceProfile(
                protein_id=pid,
                source_protein_id=p.source_id,
                protein_length=p.length,
                protein_sequence_sha256=hashlib.sha256(p.sequence.encode()).hexdigest(),
                description=p.description,
                gene=genes.get(p.gene_id),
                sources=sources,
                neighborhood=context,
            )
        )
    return EvidenceResult(
        genome_id=dataset.genome.genome_id,
        runs=runs,
        profiles=tuple(profiles),
        neighborhood_state=context_state,
        neighborhood_message=context_message,
    )
