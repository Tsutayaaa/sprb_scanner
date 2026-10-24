"""Batch execution uses real tiny scanner services, preserving every native evidence row."""

import csv
import json
from pathlib import Path

from Bio.Data.CodonTable import unambiguous_dna_by_id
from test_pipeline import small_pipeline as small_pipeline

from sprb_scanner.batch import batch_status, collect_batch, plan_batch, run_task
from sprb_scanner.batch.plan import load_manifest
from sprb_scanner.batch.status import pipeline_report


def test_real_batch_evidence_and_corrupt_output_isolation(small_pipeline, tmp_path):
    fixture = small_pipeline
    data = tmp_path / "ncbi"
    folder = data / "GCF_000000001.1"
    folder.mkdir(parents=True)
    codons = {v: k for k, v in unambiguous_dna_by_id[11].forward_table.items()}
    dna, gff, faa = "", ["##gff-version 3"], ""
    for i, protein in enumerate(fixture["dataset"].proteins):
        start = len(dna) + 1
        sequence = "M" + protein.sequence[1:]
        dna += "".join(codons[aa] for aa in sequence) + "TAA"
        gff.append(f"chr\ttest\tCDS\t{start}\t{len(dna)}\t.\t+\t0\tID=cds-{i};protein_id=p{i}")
        dna += "AAA" * 3
        faa += f">p{i}\n{sequence}\n"
    (folder / "genomic.fna").write_text(f">chr\n{dna}\n")
    (folder / "genomic.gff").write_text("\n".join(gff) + "\n")
    (folder / "protein.faa").write_text(faa)
    resources = tmp_path / "resources.json"
    resources.write_text(
        json.dumps(
            {
                "partition": "test",
                "account": "test",
                "qos": "test",
                "cpus": 1,
                "memory": "4G",
                "time": "00:10:00",
                "concurrency": 1,
                "max_array_tasks": 1,
                "max_array_size": 1001,
                "python": "/test/python",
            }
        )
    )
    out = tmp_path / "batch"
    plan_batch(data, fixture["config"], resources, out)
    path = out / "manifest.json"
    first = run_task(path, 0, local=True)
    assert first["state"] == "succeeded", first
    second = run_task(path, 0, local=True)
    plan = load_manifest(path)
    report = pipeline_report(plan, plan.entries[0], second)
    assert all(s["action"] == "reused" for s in report["stages"].values())
    collection = tmp_path / "collected"
    result = collect_batch(path, collection)
    assert result["complete_evidence_genomes"] == 1 and result["profiles"] >= 1
    with (collection / "evidence_records.tsv").open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["source_id"] for r in records} == {
        "TIGR04131",
        "PF13585",
        "TIGR03519",
        "PF11751",
        "PF13573",
        "RN28",
    }
    assert all(r["source_protein_id"] == "p0" for r in records)
    assert all(
        "protein_id" in json.loads(r["cargo_gene"]) or "gene_id" in json.loads(r["cargo_gene"])
        for r in records
    )
    native = Path(report["stages"]["evidence"]["result_directory"])
    with (native / "evidence_records.tsv").open() as handle:
        original = list(csv.DictReader(handle, delimiter="\t"))
    assert len(original) == len(records)
    assert all(
        all(record[k] == v for k, v in old.items())
        for record, old in zip(records, original, strict=True)
    )
    (native / "evidence_records.tsv").write_text("corrupt")
    assert batch_status(path, scheduler=False, verify=True)["counts"] == {"invalid_result": 1}
    excluded = collect_batch(path, tmp_path / "corrupt-collection")
    assert excluded["collected_genomes"] == 0 and excluded["excluded_genomes"] == 1
