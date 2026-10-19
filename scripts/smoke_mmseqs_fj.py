"""Opt-in real FJ proteome smoke test; no sequences/results are bundled or downloaded."""

import argparse
import csv
import json
from pathlib import Path

from sprb_scanner.core.models import GenomeDataset
from sprb_scanner.io.common import sha256
from sprb_scanner.io.exporters import validate_prepared, write_json
from sprb_scanner.scanners.mmseqs import MMseqsConfig, create_database, run_search

REFERENCE = "https://journals.asm.org/doi/10.1128/JB.00218-19"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New directory outside Git")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output must be a new directory")
    if validate_prepared(args.prepared)["errors"]:
        raise ValueError("FJ prepared data has quality errors")
    ds = GenomeDataset.model_validate_json((args.prepared / "dataset.json").read_text())
    sprb = next(p for p in ds.proteins if p.source_id == "WP_012023065.1")
    assert sprb.length == 6497
    args.output.mkdir(parents=True)
    # Native SprB segments corresponding to published C-terminal lengths; no GFP/linker tags.
    lengths = (6497, 218, 448, 1182)
    regions = [
        dict(
            query_id="SprB_full" if n == 6497 else f"SprB_C{n}",
            start=sprb.length - n,
            end=sprb.length,
            length=n,
            target_id=sprb.protein_id,
            target_source_id=sprb.source_id,
            reference=REFERENCE,
        )
        for n in lengths
    ]
    query = args.output / "sprb_queries.faa"
    query.write_text(
        "".join(
            f">{r['query_id']} source={sprb.source_id} [{r['start']},{r['end']})\n"
            f"{sprb.sequence[r['start'] : r['end']]}\n"
            for r in regions
        )
    )
    write_json(
        args.output / "query_regions.json",
        {
            "coordinates": "0-based half-open",
            "query_fasta_sha256": sha256(query),
            "regions": regions,
            "note": "native reference subsequences; not the experimental GFP fusion constructs",
        },
    )
    config = MMseqsConfig(threads=args.threads, sensitivity=7.5, evalue=0.001)
    database = args.output / "database"
    create_database(args.prepared, database, config)
    run_search(
        query, args.output / "search", prepared=args.prepared, database=database, config=config
    )
    hits = json.loads((args.output / "search/alignments.json").read_text())
    summary = []
    for region in regions:
        matching = [
            h
            for h in hits
            if h["query_id"] == region["query_id"]
            and h["protein_id"] == sprb.protein_id
            and (h["start"], h["end"]) == (region["start"], region["end"])
            and (h["query_start"], h["query_end"]) == (0, region["length"])
            and h["query_coverage"] == h["sequence_identity"] == 1
        ]
        assert matching, f"Expected exact SprB alignment missing: {region}"
        h = matching[0]
        summary.append(
            {
                k: h[k]
                for k in (
                    "query_id",
                    "protein_id",
                    "query_length",
                    "protein_length",
                    "start",
                    "end",
                    "query_start",
                    "query_end",
                    "query_coverage",
                    "target_coverage",
                    "sequence_identity",
                    "evalue",
                    "bit_score",
                )
            }
        )
    with (args.output / "expected_matches.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    write_json(
        args.output / "smoke_summary.json",
        {
            "state": "succeeded",
            "proteins": len(ds.proteins),
            "total_alignments": len(hits),
            "expected_exact_matches": summary,
            "source": REFERENCE,
            "artifact_sha256": {
                str(p.relative_to(args.output)): sha256(p)
                for p in (
                    query,
                    args.output / "query_regions.json",
                    args.output / "expected_matches.tsv",
                    args.output / "database/run_metadata.json",
                    args.output / "search/run_metadata.json",
                )
            },
        },
    )
    print(
        json.dumps(
            {
                "state": "succeeded",
                "proteins": len(ds.proteins),
                "total_alignments": len(hits),
                "expected_exact_matches": len(summary),
                "output": str(args.output.resolve()),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
