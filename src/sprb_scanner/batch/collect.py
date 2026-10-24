"""Stream existing scientific tables; validate artifacts and retain source availability."""

import csv
import json
from contextlib import ExitStack
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.evidence.models import SOURCES
from sprb_scanner.evidence.service import RECORD_FIELDS
from sprb_scanner.io.common import sha256
from sprb_scanner.neighborhood.models import NeighborRecord
from sprb_scanner.pipeline.storage import atomic_json, now

from .plan import load_manifest, read_json
from .status import batch_status, pipeline_report, task_record
from .task import disk_usage


def _writer(stack, path, fields):
    handle = stack.enter_context(path.open("w", newline="", encoding="utf-8"))
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
    writer.writeheader()
    return writer


def native_schemas():
    return {
        "evidence_profiles.tsv": [
            "protein_id",
            "source_protein_id",
            "protein_length",
            "gene_id",
            "locus_tag",
            "description",
            *[
                s + suffix
                for s in SOURCES
                for suffix in ("_state", "_hit_count", "_ga_domain_count")
            ],
            "neighborhood_state",
            "neighborhood_outcome",
            "local_sprf_proteins",
            "directly_adjacent_sprf_proteins",
        ],
        "evidence_records.tsv": [
            "protein_id",
            "source_id",
            "hit_id",
            "evidence_kind",
            "scan_directory",
            *RECORD_FIELDS,
            "evidence",
        ],
        "neighborhood_associations.tsv": [
            "protein_id",
            "neighbor_protein_id",
            "neighbor_source_protein_id",
            "domain_hit_ids",
            *NeighborRecord.model_fields,
        ],
    }


TABLES = {
    "evidence_profiles.tsv": "protein_summary.tsv",
    "evidence_records.tsv": "evidence_records.tsv",
    "neighborhood_associations.tsv": "association_summary.tsv",
}


def extra_fields(source_file):
    return (
        ["cargo_gene"]
        if source_file == "evidence_profiles.tsv"
        else ["source_protein_id", "cargo_gene"]
    )


def collect_batch(path, output):
    path, output = Path(path).resolve(), Path(output).resolve()
    manifest = load_manifest(path)
    status = batch_status(path, verify=True)
    if output.exists():
        raise InputError("Collection output already exists; choose a new snapshot directory")
    output.mkdir(parents=True)
    base = ["assembly_accession", "genome_id", "array_index"]
    genome_fields = [
        *base,
        "state",
        "error",
        "preparation_state",
        "pipeline_state",
        "proteins",
        "evidence_profiles",
        "evidence_complete",
        "typeb_counts_complete",
        "association_counts_complete",
        "typeb_retained_proteins",
        "typeb_ga_supported_proteins",
        "direct_adjacent_typeb_cargos",
        "local_hit_typeb_cargos",
        "no_local_hit_typeb_cargos",
        "rn28_proteins",
        *[
            s + suffix
            for s in SOURCES
            for suffix in (
                "_state",
                "_retained_records",
                "_retained_proteins",
                "_ga_supported_proteins",
            )
        ],
        "apparent_bytes",
        "allocated_bytes",
        "files",
    ]
    skipped = []
    native_writers = {}
    schemas = native_schemas()
    counters = {"collected_genomes": 0, "complete_evidence_genomes": 0, "profiles": 0}
    atomic_json(output / "status.json", {"state": "running"})
    try:
        with ExitStack() as stack:
            for source_file, output_file in TABLES.items():
                fields = [*base, *schemas[source_file], *extra_fields(source_file)]
                native_writers[output_file] = (fields, _writer(stack, output / output_file, fields))
            subsets = {}
            for filename, source_file in (
                ("cargo_summary.tsv", "evidence_profiles.tsv"),
                ("domain_summary.tsv", "evidence_records.tsv"),
                ("rn28_alignments.tsv", "evidence_records.tsv"),
                ("architecture_summary.tsv", "evidence_records.tsv"),
            ):
                fields = [*base, *schemas[source_file], *extra_fields(source_file)]
                subsets[filename] = _writer(stack, output / filename, fields)
            genome_writer = _writer(stack, output / "genome_summary.tsv", genome_fields)
            sources_writer = _writer(
                stack,
                output / "source_status.tsv",
                [*base, "source_id", "state", "retained_hits", "message"],
            )
            provenance = stack.enter_context((output / "provenance.jsonl").open("w"))
            for entry, row in zip(manifest.entries, status["entries"], strict=True):
                prefix = {
                    "assembly_accession": entry.assembly_accession,
                    "genome_id": entry.genome_id,
                    "array_index": entry.index,
                }
                summary = prefix | {
                    k: row.get(k) for k in ("state", "error", "preparation_state", "pipeline_state")
                }
                root = Path(manifest.output_root) / "genomes" / entry.genome_id
                summary.update(disk_usage(root))
                if row["state"] not in {"succeeded", "partial"} or not row["verified"]:
                    skipped.append(prefix | {"state": row["state"], "error": row.get("error")})
                    for source in SOURCES:
                        summary[source + "_state"] = "unavailable"
                        sources_writer.writerow(
                            prefix
                            | {
                                "source_id": source,
                                "state": "unavailable",
                                "message": "No verified merged source evidence; "
                                "inspect genome/task state",
                            }
                        )
                    genome_writer.writerow(summary)
                    continue
                record = task_record(manifest, entry)
                run = pipeline_report(manifest, entry, record)
                evidence = run["stages"]["evidence"]
                if not evidence.get("result_usable"):
                    summary.update(
                        state="missing_evidence", error="No usable final evidence output"
                    )
                    skipped.append(prefix | {"state": summary["state"], "error": summary["error"]})
                    genome_writer.writerow(summary)
                    continue
                native = Path(evidence["result_directory"])
                data = read_json(native / "evidence_profiles.json")
                if data["genome_id"] != entry.genome_id:
                    raise InputError("Collected evidence belongs to another genome")
                profiles = data["profiles"]
                by_protein = {p["protein_id"]: p for p in profiles}
                complete = read_json(native / "status.json")["evidence_complete"]
                summary.update(
                    proteins=record["preparation"]["proteins"],
                    evidence_profiles=len(profiles),
                    evidence_complete=complete,
                )
                for source, source_run in data["runs"].items():
                    sources_writer.writerow(
                        prefix
                        | {"source_id": source}
                        | {k: source_run.get(k) for k in ("state", "retained_hits", "message")}
                    )
                    summary[source + "_state"] = source_run["state"]
                    summary[source + "_retained_records"] = source_run["retained_hits"]
                    if source_run["state"] == "succeeded":
                        summary[source + "_retained_proteins"] = sum(
                            p["sources"][source]["state"] == "succeeded_with_hits" for p in profiles
                        )
                        if source != "RN28":
                            summary[source + "_ga_supported_proteins"] = sum(
                                any(
                                    d["hit"]["hmm"]["ga"]["status"] == "passed_ga"
                                    for d in p["sources"][source]["domains"]
                                )
                                for p in profiles
                            )
                cargos = [
                    p
                    for p in profiles
                    if any(p["sources"][s]["domains"] for s in ("TIGR04131", "PF13585"))
                ]
                # These are observed retained-evidence counts, never functional classifications.
                summary.update(
                    typeb_retained_proteins=len(cargos),
                    typeb_ga_supported_proteins=sum(
                        any(
                            d["hit"]["hmm"]["ga"]["status"] == "passed_ga"
                            for s in ("TIGR04131", "PF13585")
                            for d in p["sources"][s]["domains"]
                        )
                        for p in cargos
                    ),
                    local_hit_typeb_cargos=sum(
                        p["neighborhood"]["outcome"] == "local_hit" for p in cargos
                    ),
                    no_local_hit_typeb_cargos=sum(
                        p["neighborhood"]["outcome"] == "no_local_hit" for p in cargos
                    ),
                    direct_adjacent_typeb_cargos=sum(
                        any(
                            n["relationship"]["directly_adjacent"]
                            for n in p["neighborhood"]["neighbors"]
                        )
                        for p in cargos
                    ),
                    rn28_proteins=summary.get("RN28_retained_proteins"),
                )
                ctd_ready = all(
                    data["runs"][s]["state"] == "succeeded" for s in ("TIGR04131", "PF13585")
                )
                neighborhood_ready = ctd_ready and data["neighborhood_state"] == "succeeded"
                summary.update(
                    typeb_counts_complete=ctd_ready, association_counts_complete=neighborhood_ready
                )
                if not any(
                    data["runs"][s]["state"] == "succeeded" for s in ("TIGR04131", "PF13585")
                ):
                    summary["typeb_retained_proteins"] = None
                    summary["typeb_ga_supported_proteins"] = None
                if data["neighborhood_state"] not in {"succeeded", "partial"}:
                    for key in (
                        "direct_adjacent_typeb_cargos",
                        "local_hit_typeb_cargos",
                        "no_local_hit_typeb_cargos",
                    ):
                        summary[key] = None
                # Copy native rows verbatim, including scores, coverage, positions and raw JSON.
                for source_file, output_file in TABLES.items():
                    with (native / source_file).open(newline="") as handle:
                        reader = csv.DictReader(handle, delimiter="\t")
                        extra = extra_fields(source_file)
                        fields = [*base, *reader.fieldnames, *extra]
                        if output_file not in native_writers:
                            native_writers[output_file] = (
                                fields,
                                _writer(stack, output / output_file, fields),
                            )
                        expected, writer = native_writers[output_file]
                        if fields != expected:
                            raise InputError("Incompatible native scientific table schemas")
                        for native_row in reader:
                            profile = by_protein[native_row["protein_id"]]
                            extras = {"cargo_gene": json.dumps(profile["gene"], ensure_ascii=False)}
                            if "source_protein_id" in extra:
                                extras["source_protein_id"] = profile["source_protein_id"]
                            exported = prefix | native_row | extras
                            writer.writerow(exported)
                            if source_file == "evidence_profiles.tsv":
                                if any(
                                    profile["sources"][s]["domains"]
                                    for s in ("TIGR04131", "PF13585")
                                ):
                                    subsets["cargo_summary.tsv"].writerow(exported)
                            if source_file == "evidence_records.tsv":
                                destination = (
                                    "rn28_alignments.tsv"
                                    if native_row["evidence_kind"] == "alignment"
                                    else "domain_summary.tsv"
                                )
                                subsets[destination].writerow(exported)
                                if native_row["source_id"] in {"TIGR03519", "PF11751"}:
                                    subsets["architecture_summary.tsv"].writerow(exported)
                provenance.write(
                    json.dumps(
                        prefix
                        | {
                            "manifest_sha256": sha256(path),
                            "annotation_identity": entry.annotation_identity,
                            "input_files": {k: v.model_dump() for k, v in entry.files.items()},
                            "pipeline_run": run,
                            "evidence_metadata": read_json(native / "run_metadata.json"),
                        },
                        ensure_ascii=False,
                    )
                    + "\n"
                )
                # Recheck after reading native rows: mutation during collection is fatal,
                # leaving an explicitly failed collection snapshot rather than silent mixing.
                pipeline_report(manifest, entry, record)
                genome_writer.writerow(summary)
                counters["collected_genomes"] += 1
                counters["complete_evidence_genomes"] += complete
                counters["profiles"] += len(profiles)
        atomic_json(
            output / "collection_metadata.json",
            {
                "created_at": now(),
                "manifest": str(path),
                "manifest_sha256": sha256(path),
                "pinned_files": manifest.pinned_files,
                "software": manifest.software,
                "resources": manifest.resources.model_dump(),
                "counts": counters,
                "excluded": skipped,
                "semantics": "Evidence union; counts are observed retained records/proteins, "
                "not family "
                "confirmation. GA remains model-specific. no_local_hit is not orphan; "
                "PorP/SprF HMM hits do not establish shuttle function. Missing counts are unknown.",
                "output_sha256": {
                    p.name: sha256(p) for p in output.iterdir() if p.name != "status.json"
                },
            },
        )
        atomic_json(
            output / "status.json",
            {"state": "succeeded", **counters, "excluded_genomes": len(skipped)},
        )
    except BaseException as exc:
        atomic_json(output / "status.json", {"state": "failed", "error": str(exc)})
        raise
    return read_json(output / "status.json")
