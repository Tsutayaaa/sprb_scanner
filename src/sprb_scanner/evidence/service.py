"""Offline orchestration and atomic JSON/TSV export of unified evidence profiles."""

import shutil
import tempfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from sprb_scanner import __version__
from sprb_scanner.association.analyzer import CTD_MODELS, SPRF_MODELS
from sprb_scanner.association.inputs import check_unchanged, read_manifest
from sprb_scanner.association.service import write_rows
from sprb_scanner.core.exceptions import InputError
from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.evidence.inputs import load_hmm, load_pair, load_rn28
from sprb_scanner.evidence.merger import merge
from sprb_scanner.evidence.models import SOURCES, SourceRun
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared, write_json
from sprb_scanner.neighborhood.analyzer import GEOMETRY_POLICY

RECORD_FIELDS = [
    "model_id",
    "model_version",
    "model_sha256",
    "protein_length",
    "start",
    "end",
    "query_id",
    "query_length",
    "query_start",
    "query_end",
    "query_coverage",
    "target_coverage",
    "sequence_identity",
    "alignment_length",
    "evalue",
    "evalue_kind",
    "bit_score",
    "ga_status",
    "sequence_evalue",
    "sequence_bit_score",
    "independent_domain_evalue",
    "conditional_domain_evalue",
    "hmm_coverage",
]


def record_fields(hit):
    fields = {
        k: getattr(hit, k, None)
        for k in RECORD_FIELDS
        if k
        not in {
            "model_id",
            "model_version",
            "model_sha256",
            "ga_status",
            "sequence_evalue",
            "sequence_bit_score",
            "independent_domain_evalue",
            "conditional_domain_evalue",
            "hmm_coverage",
        }
    }
    fields.update(
        model_id=hit.model.model_id,
        model_version=hit.model.version,
        model_sha256=hit.model.sha256,
        ga_status=hit.hmm.ga.status if hit.hmm else None,
        hmm_coverage=hit.query_coverage if hit.hmm else None,
    )
    for name in (
        "sequence_evalue",
        "sequence_bit_score",
        "independent_domain_evalue",
        "conditional_domain_evalue",
    ):
        fields[name] = getattr(hit.hmm, name) if hit.hmm else None
    return fields


def tables(stage, result):
    rows, evidence_rows, relations = [], [], []
    for profile in result.profiles:
        row = {
            "protein_id": profile.protein_id,
            "source_protein_id": profile.source_protein_id,
            "protein_length": profile.protein_length,
            "gene_id": profile.gene.gene_id if profile.gene else None,
            "locus_tag": profile.gene.locus_tag if profile.gene else None,
            "description": profile.description,
        }
        for source, evidence in profile.sources.items():
            row[source + "_state"] = evidence.state
            row[source + "_hit_count"] = (
                len(evidence.domains) + len(evidence.alignments)
                if evidence.state.startswith("succeeded_")
                else None
            )
            row[source + "_ga_domain_count"] = (
                sum(e.hit.hmm.ga.status == "passed_ga" for e in evidence.domains)
                if source != "RN28" and evidence.state.startswith("succeeded_")
                else None
            )
            for domain in evidence.domains:
                evidence_rows.append(
                    {
                        "protein_id": profile.protein_id,
                        "source_id": source,
                        "hit_id": domain.hit.hit_id,
                        "evidence_kind": "hmm_domain",
                        "scan_directory": result.runs[source].directory,
                        "evidence": domain.model_dump(mode="json"),
                        **record_fields(domain.hit),
                    }
                )
            for alignment in evidence.alignments:
                evidence_rows.append(
                    {
                        "protein_id": profile.protein_id,
                        "source_id": source,
                        "hit_id": alignment.hit_id,
                        "evidence_kind": "alignment",
                        "scan_directory": result.runs[source].directory,
                        "evidence": alignment.model_dump(mode="json"),
                        **record_fields(alignment),
                    }
                )
        context = profile.neighborhood
        row.update(
            neighborhood_state=context.state,
            neighborhood_outcome=context.outcome,
            local_sprf_proteins=len(context.neighbors)
            if context.state in {"succeeded", "partial"}
            else None,
            directly_adjacent_sprf_proteins=sum(
                n.relationship.directly_adjacent for n in context.neighbors
            )
            if context.state in {"succeeded", "partial"}
            else None,
        )
        rows.append(row)
        for neighbor in context.neighbors:
            relations.append(
                {
                    "protein_id": profile.protein_id,
                    "neighbor_protein_id": neighbor.protein_id,
                    "neighbor_source_protein_id": neighbor.source_protein_id,
                    "domain_hit_ids": neighbor.domain_hit_ids,
                    **neighbor.relationship.model_dump(mode="json"),
                }
            )
    fields = [
        "protein_id",
        "source_protein_id",
        "protein_length",
        "gene_id",
        "locus_tag",
        "description",
        *[
            source + suffix
            for source in SOURCES
            for suffix in ("_state", "_hit_count", "_ga_domain_count")
        ],
        "neighborhood_state",
        "neighborhood_outcome",
        "local_sprf_proteins",
        "directly_adjacent_sprf_proteins",
    ]
    write_rows(stage / "evidence_profiles.tsv", rows, fields)
    write_rows(
        stage / "evidence_records.tsv",
        evidence_rows,
        [
            "protein_id",
            "source_id",
            "hit_id",
            "evidence_kind",
            "scan_directory",
            *RECORD_FIELDS,
            "evidence",
        ],
    )
    from sprb_scanner.neighborhood.models import NeighborRecord

    write_rows(
        stage / "neighborhood_associations.tsv",
        relations,
        [
            "protein_id",
            "neighbor_protein_id",
            "neighbor_source_protein_id",
            "domain_hit_ids",
            *NeighborRecord.model_fields,
        ],
    )


def run_merge(
    prepared,
    output,
    *,
    bctd_scan=None,
    sprf_scan=None,
    repeat_scan=None,
    rn28_scan=None,
    registry=Path("references/registry.json"),
    window=5,
    neighborhood=True,
    source_states=None,
    rn28_version=None,
    rn28_artifact=None,
):
    """Merge existing runs only. Failed sources are data states, corrupt inputs are errors."""
    prepared, output = Path(prepared), Path(output)
    if output.exists():
        raise InputError(f"Output already exists: {output}")
    fingerprints = {}
    runs = {source: SourceRun(source_id=source, state="not_run") for source in SOURCES}
    for source, run in (source_states or {}).items():
        if (
            source not in runs
            or not isinstance(run, SourceRun)
            or run.source_id != source
            or run.state not in {"not_run", "failed"}
        ):
            raise InputError("Explicit source states must be not_run/failed canonical sources")
        runs[source] = run
    hits, alignments = [], ()
    try:
        prepared_meta = read_manifest(prepared, fingerprints)
        if validate_prepared(prepared)["errors"]:
            raise InputError("Prepared dataset has quality errors")
        dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
        for directory, models in ((bctd_scan, CTD_MODELS), (sprf_scan, SPRF_MODELS)):
            if directory is not None:
                statuses, group = load_pair(directory, models, dataset, fingerprints)
                if set(statuses) & set(source_states or {}):
                    raise InputError("Source state conflicts with supplied scanner output")
                runs.update(statuses)
                hits.extend(group)
        if repeat_scan is not None:
            if "PF13573" in (source_states or {}):
                raise InputError("Source state conflicts with repeat output")
            run, group = load_hmm(repeat_scan, "PF13573", dataset, fingerprints)
            runs["PF13573"] = run
            hits.extend(group)
        if rn28_scan is not None:
            if "RN28" in (source_states or {}):
                raise InputError("Source state conflicts with RN28 output")
            run, alignments = load_rn28(
                rn28_scan, dataset, registry, fingerprints, rn28_version, rn28_artifact
            )
            runs["RN28"] = run
        result = merge(dataset, runs, hits, alignments, window, neighborhood)
        check_unchanged(fingerprints)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise InputError(f"Cannot merge evidence inputs: {exc}") from exc
    summary = {
        "genome_id": result.genome_id,
        "candidate_rule": result.candidate_rule,
        "profiles": len(result.profiles),
        "source_runs": {
            source: {"state": run.state, "retained_hits": run.retained_hits, "message": run.message}
            for source, run in runs.items()
        },
        "protein_states": {
            source: dict(Counter(p.sources[source].state for p in result.profiles))
            for source in SOURCES
        },
        "neighborhood_run_state": result.neighborhood_state,
        "neighborhood_run_message": result.neighborhood_message,
        "neighborhood_states": dict(Counter(p.neighborhood.state for p in result.profiles)),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".evidence-", dir=output.parent))
    try:
        write_json(stage / "evidence_profiles.json", result.model_dump(mode="json"))
        tables(stage, result)
        write_json(stage / "summary.json", summary)
        complete = all(r.state == "succeeded" for r in runs.values()) and (
            not neighborhood or result.neighborhood_state == "succeeded"
        )
        write_json(stage / "status.json", {"state": "succeeded", "evidence_complete": complete})
        check_unchanged(fingerprints)
        write_json(
            stage / "run_metadata.json",
            {
                "software_version": __version__,
                "created_at": datetime.now(UTC).isoformat(),
                "genome_id": result.genome_id,
                "parameters": {
                    "window": window,
                    "neighborhood": neighborhood,
                    "candidate_rule": result.candidate_rule,
                    "coverage_cutoff": None,
                },
                "coordinates": result.coordinates,
                "geometry_policy": GEOMETRY_POLICY,
                "input_files_sha256": fingerprints,
                "prepared_metadata": prepared_meta,
                "source_runs": {k: r.model_dump(mode="json") for k, r in runs.items()},
                "semantics": "no score fusion, filtering, ranking or functional classification; "
                "counts are retained domain/alignment records, not physical repeat counts; "
                "no_local_hit requires both PorP/SprF sources to have succeeded; "
                "incomplete-source local counts are observed lower bounds",
                "output_sha256": {p.name: sha256(p) for p in stage.iterdir()},
            },
        )
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return summary
