"""Generic single-profile protein HMM engine with lossless domain evidence."""

import csv
import hashlib
import io
import json
import shutil
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import pyhmmer
from pydantic import Field

from sprb_scanner.core.config import ProfileConfig
from sprb_scanner.core.exceptions import HMMError
from sprb_scanner.core.models import (
    DomainHit,
    GatheringEvidence,
    HMMStatistics,
    Interval,
    Positive,
    Record,
    ScanStatus,
    internal_id,
)
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import write_json, write_tsv


class HMMSearchConfig(Record):
    sequence_report_evalue: float = Field(default=10, gt=0, allow_inf_nan=False)
    domain_report_evalue: float = Field(default=10, gt=0, allow_inf_nan=False)
    threads: Positive = 1
    seed: int = Field(default=42, ge=1, le=2147483647, strict=True)
    max_sensitivity: bool = False


def gathering_evidence(sequence_score, domain_score, cutoffs):
    """Use the loaded model's GA; absent cutoffs are not inferred."""
    if cutoffs is None:
        return GatheringEvidence(status="unavailable")
    sequence_cutoff, domain_cutoff = cutoffs
    sequence_passed = sequence_score >= sequence_cutoff
    domain_passed = domain_score >= domain_cutoff
    return GatheringEvidence(
        status="passed_ga" if sequence_passed and domain_passed else "below_ga",
        sequence_cutoff=sequence_cutoff,
        domain_cutoff=domain_cutoff,
        sequence_margin=sequence_score - sequence_cutoff,
        domain_margin=domain_score - domain_cutoff,
        sequence_passed=sequence_passed,
        domain_passed=domain_passed,
    )


def evidence_summary(hits, ga_available):
    """Disjoint protein groups: any passed domain takes precedence over below-GA."""
    states = ("passed_ga", "below_ga", "unavailable")
    domains = {state: 0 for state in states}
    proteins = {}
    for hit in hits:
        state = hit.hmm.ga.status
        domains[state] += 1
        proteins.setdefault(hit.protein_id, set()).add(state)
    counts = {state: 0 for state in states}
    for values in proteins.values():
        state = next(state for state in states if state in values)
        counts[state] += 1
    return {
        "ga_available": ga_available,
        "domain_counts": domains,
        "protein_counts": counts,
        "total_domains": len(hits),
        "total_hit_proteins": len(proteins),
        "protein_policy": "any passed_ga domain wins; otherwise below_ga or unavailable",
        "scope": "all backend-returned domains, including unreported; not a final classification",
    }


def write_evidence_table(path, hits):
    fields = [
        "hit_id",
        "protein_id",
        "model_id",
        "model_version",
        "model_sha256",
        "ga_status",
        "sequence_ga",
        "domain_ga",
        "sequence_margin",
        "domain_margin",
        "sequence_passed_ga",
        "domain_passed_ga",
        "sequence_bit_score",
        "domain_bit_score",
        "sequence_evalue",
        "independent_domain_evalue",
        "conditional_domain_evalue",
        "start",
        "end",
        "query_start",
        "query_end",
        "ctd_distance",
        "query_coverage",
        "target_coverage",
        "sequence_reported",
        "domain_reported",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            stats, ga = hit.hmm, hit.hmm.ga
            writer.writerow(
                dict(
                    hit_id=hit.hit_id,
                    protein_id=hit.protein_id,
                    model_id=hit.model.model_id,
                    model_version=hit.model.version,
                    model_sha256=hit.model.sha256,
                    ga_status=ga.status,
                    sequence_ga=ga.sequence_cutoff,
                    domain_ga=ga.domain_cutoff,
                    sequence_margin=ga.sequence_margin,
                    domain_margin=ga.domain_margin,
                    sequence_passed_ga=ga.sequence_passed,
                    domain_passed_ga=ga.domain_passed,
                    sequence_bit_score=stats.sequence_bit_score,
                    domain_bit_score=hit.bit_score,
                    sequence_evalue=stats.sequence_evalue,
                    independent_domain_evalue=stats.independent_domain_evalue,
                    conditional_domain_evalue=stats.conditional_domain_evalue,
                    start=hit.start,
                    end=hit.end,
                    query_start=hit.query_start,
                    query_end=hit.query_end,
                    ctd_distance=hit.ctd_distance,
                    query_coverage=hit.query_coverage,
                    target_coverage=hit.target_coverage,
                    sequence_reported=stats.sequence_reported,
                    domain_reported=stats.domain_reported,
                )
            )


def text(value):
    return value.decode() if isinstance(value, bytes) else value


def load_hmm(profile: ProfileConfig):
    """Validate bytes, checksum, alphabet and identity; reject empty/multi-model files."""
    try:
        if profile.method != "hmm":
            raise ValueError("Profile method must be hmm")
        data = profile.path.read_bytes()
        if hashlib.sha256(data).hexdigest() != profile.provenance.sha256:
            raise ValueError("HMM checksum mismatch")
        with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as handle:
            models = list(handle)
        if len(models) != 1:
            raise ValueError("Expected exactly one HMM per profile file")
        hmm = models[0]
        hmm.validate()
        if not hmm.alphabet.is_amino():
            raise ValueError("Only amino-acid HMMs are supported")
        identity = text(hmm.accession or hmm.name)
        if identity.split(".")[0] != profile.provenance.model_id.split(".")[0]:
            raise ValueError(f"HMM identity mismatch: {identity}")
        if hmm.accession and "." in identity and profile.provenance.version != identity:
            raise ValueError(f"Version must match HMM accession: {identity}")
        return hmm
    except (OSError, ValueError, EOFError) as exc:
        raise HMMError(f"Cannot load HMM: {exc}") from exc


def scan(proteins, profile: ProfileConfig, config: HMMSearchConfig):
    """Return every backend domain, raw evidence and native tables; never length-filter."""
    hmm = load_hmm(profile)
    proteins = tuple(proteins)
    if len({p.protein_id for p in proteins}) != len(proteins):
        raise HMMError("Duplicate internal protein IDs")
    # Backend tokens avoid truncation/encoding ambiguity in arbitrary public identifiers.
    lookup = {f"p{i}": p for i, p in enumerate(proteins)}
    targets = pyhmmer.easel.DigitalSequenceBlock(
        hmm.alphabet,
        [
            pyhmmer.easel.TextSequence(name=key, sequence=p.sequence).digitize(hmm.alphabet)
            for key, p in lookup.items()
        ],
    )
    options = dict(
        E=config.sequence_report_evalue, domE=config.domain_report_evalue, seed=config.seed
    )
    if config.max_sensitivity:
        options.update(F1=1.0, F2=1.0, F3=1.0, bias_filter=False)
    hits, raw = [], []
    native = {"targets.tbl": b"", "domains.domtbl": b""}
    metadata = {
        "model_name": text(hmm.name),
        "accession": text(hmm.accession),
        "model_length": hmm.M,
        "target_count": len(proteins),
        "Z": 0,
        "domZ": 0,
        "backend": f"pyhmmer {pyhmmer.__version__}",
        "protein_ids": {k: p.protein_id for k, p in lookup.items()},
        "input_proteins_sha256": hashlib.sha256(
            json.dumps(
                [(p.protein_id, p.sequence) for p in proteins],
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode()
        ).hexdigest(),
        "trusted_cutoffs": hmm.cutoffs.trusted,
        "gathering_cutoffs": hmm.cutoffs.gathering,
    }
    if not proteins:
        return hits, raw, native, metadata
    try:
        top = next(pyhmmer.hmmsearch(hmm, targets, cpus=config.threads, **options))
        metadata.update(Z=top.Z, domZ=top.domZ)
        for name, fmt in (("targets.tbl", "targets"), ("domains.domtbl", "domains")):
            buf = io.BytesIO()
            top.write(buf, format=fmt)
            native[name] = buf.getvalue()
        for target in top:
            p = lookup[text(target.name)]
            for number, domain in enumerate(target.domains, 1):
                a = domain.alignment
                hit_id = internal_id(profile.provenance.model_id, p.protein_id, str(number))
                record = dict(
                    hit_id=hit_id,
                    protein_id=p.protein_id,
                    backend_target=text(target.name),
                    domain_number=number,
                    target_from=a.target_from,
                    target_to=a.target_to,
                    hmm_from=a.hmm_from,
                    hmm_to=a.hmm_to,
                    env_from=domain.env_from,
                    env_to=domain.env_to,
                    sequence_evalue=target.evalue,
                    sequence_bit_score=target.score,
                    independent_domain_evalue=domain.i_evalue,
                    conditional_domain_evalue=domain.c_evalue,
                    domain_bit_score=domain.score,
                    sequence_reported=target.reported,
                    sequence_included=target.included,
                    domain_reported=domain.reported,
                    domain_included=domain.included,
                    hmm_sequence=text(a.hmm_sequence),
                    target_sequence=text(a.target_sequence),
                )
                raw.append(record)
                stats = HMMStatistics(
                    **{
                        k: record[k]
                        for k in HMMStatistics.model_fields
                        if k not in {"envelope", "ga"}
                    },
                    envelope=Interval(start=domain.env_from - 1, end=domain.env_to),
                    ga=gathering_evidence(target.score, domain.score, hmm.cutoffs.gathering),
                )
                hits.append(
                    DomainHit(
                        hit_id=hit_id,
                        protein_id=p.protein_id,
                        method="hmm",
                        model=profile.provenance,
                        query_id=text(hmm.accession or hmm.name),
                        query_length=hmm.M,
                        query_start=a.hmm_from - 1,
                        query_end=a.hmm_to,
                        start=a.target_from - 1,
                        end=a.target_to,
                        protein_length=p.length,
                        evalue=domain.i_evalue,
                        evalue_kind="independent_domain",
                        bit_score=domain.score,
                        query_coverage=(a.hmm_to - a.hmm_from + 1) / hmm.M,
                        target_coverage=(a.target_to - a.target_from + 1) / p.length,
                        backend_version=metadata["backend"],
                        raw_artifact="raw_domains.json",
                        hmm=stats,
                    )
                )
        return hits, raw, native, metadata
    except (ValueError, RuntimeError, OverflowError) as exc:
        raise HMMError(f"HMM search failed: {exc}") from exc


def run_hmm(proteins, profile: ProfileConfig, output: Path, config=None):
    """Write a new run directory, including explicit failed status on engine errors."""
    config = config or HMMSearchConfig()
    if output.exists():
        raise HMMError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".hmm-", dir=output.parent))
    failure = None
    try:
        try:
            hits, raw, native, metadata = scan(proteins, profile, config)
            status = ScanStatus(model_id=profile.provenance.model_id, state="succeeded")
            write_json(stage / "hits.json", [h.model_dump(mode="json") for h in hits])
            write_tsv(stage / "hits.tsv", hits, list(DomainHit.model_fields))
            write_json(stage / "raw_domains.json", raw)
            evidence = evidence_summary(hits, metadata["gathering_cutoffs"] is not None)
            write_json(stage / "evidence_summary.json", evidence)
            write_evidence_table(stage / "hit_evidence.tsv", hits)
            for name, data in native.items():
                (stage / name).write_bytes(data)
            metadata["domain_count"] = len(hits)
        except HMMError as exc:
            failure = exc
            status = ScanStatus(
                model_id=profile.provenance.model_id, state="failed", message=str(exc)
            )
            metadata = {}
        write_json(stage / "status.json", status.model_dump())
        write_json(
            stage / "run_metadata.json",
            metadata
            | {
                "created_at": datetime.now(UTC).isoformat(),
                "configuration": config.model_dump(),
                "profile": profile.model_dump(mode="json"),
                "coordinates": "normalized: 0-based half-open; raw: HMMER 1-based inclusive",
                "retention": (
                    "all backend-returned domains, including unreported; "
                    "sequence reporting thresholds and heuristic filters may omit candidates"
                ),
                "output_sha256": {p.name: sha256(p) for p in stage.iterdir()},
            },
        )
        stage.rename(output)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    if failure:
        raise failure
    return {
        "state": status.state,
        "domain_count": len(hits),
        "output": str(output),
        "ga_available": evidence["ga_available"],
        "ga_domain_counts": evidence["domain_counts"],
        "ga_protein_counts": evidence["protein_counts"],
    }
