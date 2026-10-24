"""Combine validated collection snapshots without reopening legacy execution plans."""

import csv
import json
from contextlib import ExitStack
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256
from sprb_scanner.pipeline.storage import atomic_json, now

from .plan import load_manifest, read_json

SOURCE_FIELD = "collection_manifest_sha256"


def _rows(path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _source(root):
    root = Path(root).resolve()
    metadata = read_json(root / "collection_metadata.json")
    if read_json(root / "status.json")["state"] != "succeeded":
        raise InputError(f"Incomplete collection snapshot: {root}")
    for name, checksum in metadata["output_sha256"].items():
        target = root / name
        if Path(name).name != name or target.is_symlink() or sha256(target) != checksum:
            raise InputError(f"Collection artifact checksum mismatch: {target}")
    manifest_path = Path(metadata["manifest"])
    if sha256(manifest_path) != metadata["manifest_sha256"]:
        raise InputError("Collection input manifest changed")
    # Do not compare legacy software pins against the currently installed implementation.
    # The immutable manifest identity is preserved alongside its original software provenance.
    manifest = load_manifest(manifest_path, pins=False)
    if (
        metadata["software"] != manifest.software
        or metadata["pinned_files"] != manifest.pinned_files
    ):
        raise InputError("Collection provenance differs from immutable manifest")
    rows = list(_rows(root / "genome_summary.tsv"))
    by_genome = {r["genome_id"]: r for r in rows}
    if len(rows) != len(by_genome) or set(by_genome) != {e.genome_id for e in manifest.entries}:
        raise InputError("Collection genome inventory does not match its manifest")
    if any(r["assembly_accession"] != r["genome_id"] for r in rows):
        raise InputError("Collection assembly identity mismatch")
    return root, metadata, manifest, by_genome


def combine_collections(collections, output):
    """Oldest first; replace only unavailable genomes, never existing scientific success.

    This is a collection migration, not stage-cache reuse or a scientific rerun. Every
    row carries its original manifest identity, which resolves its array/software identity.
    """
    if len(collections) < 2:
        raise InputError("Provide at least two versioned collection snapshots")
    sources = [_source(p) for p in collections]
    identities = [m["manifest_sha256"] for _, m, _, _ in sources]
    if len(set(identities)) != len(identities):
        raise InputError("Duplicate collection manifest")
    outputs = set(sources[0][1]["output_sha256"])
    if not {"genome_summary.tsv", "source_status.tsv", "provenance.jsonl"} <= outputs:
        raise InputError("Incomplete scientific collection schema")
    if any(set(m["output_sha256"]) != outputs for _, m, _, _ in sources):
        raise InputError("Collection schemas differ")
    selected, input_files, superseded = {}, {}, []
    for i, (_, _metadata, manifest, _genome_rows) in enumerate(sources):
        for entry in manifest.entries:
            gid = entry.genome_id
            identity = {k: (f.path, f.size_bytes, f.sha256) for k, f in entry.files.items()}
            if gid in input_files and input_files[gid] != identity:
                raise InputError(f"Genome source inputs differ across batches: {gid}")
            input_files[gid] = identity
            if gid in selected:
                previous = sources[selected[gid]][3][gid]
                if previous["state"] in {"succeeded", "partial"}:
                    raise InputError(
                        f"Refusing to replace successful/partial scientific output: {gid}"
                    )
                superseded.append(
                    {
                        "genome_id": gid,
                        "previous": previous,
                        "previous_manifest_sha256": identities[selected[gid]],
                        "replacement_manifest_sha256": identities[i],
                    }
                )
            selected[gid] = i
    output = Path(output).resolve()
    if output.exists():
        raise InputError("Combined output exists; choose a new immutable snapshot")
    output.mkdir(parents=True)
    atomic_json(output / "status.json", {"state": "running"})
    try:
        for name in sorted(outputs):
            if name == "provenance.jsonl":
                with (output / name).open("w") as target:
                    for i, (root, _, _, _) in enumerate(sources):
                        with (root / name).open() as handle:
                            for line in handle:
                                row = json.loads(line)
                                gid = row["genome_id"]
                                if gid not in selected:
                                    raise InputError("Unknown genome in collection provenance")
                                if selected[gid] != i:
                                    raise InputError(
                                        "Superseded unavailable genome has scientific provenance"
                                    )
                                if SOURCE_FIELD in row:
                                    raise InputError("Recursive combined provenance is unsupported")
                                target.write(
                                    json.dumps({SOURCE_FIELD: identities[i], **row}) + "\n"
                                )
                continue
            if not name.endswith(".tsv"):
                raise InputError(f"Unsupported collection artifact: {name}")
            with ExitStack() as stack:
                readers = []
                for root, _, _, _ in sources:
                    readers.append(
                        csv.DictReader(
                            stack.enter_context((root / name).open(newline="")), delimiter="\t"
                        )
                    )
                fields = readers[0].fieldnames
                if (
                    not fields
                    or SOURCE_FIELD in fields
                    or any(reader.fieldnames != fields for reader in readers)
                ):
                    raise InputError(f"Incompatible table schema: {name}")
                writer = csv.DictWriter(
                    stack.enter_context((output / name).open("w", newline="")),
                    fieldnames=[SOURCE_FIELD, *fields],
                    delimiter="\t",
                )
                writer.writeheader()
                if name == "genome_summary.tsv":
                    for gid in sorted(selected):
                        i = selected[gid]
                        writer.writerow({SOURCE_FIELD: identities[i], **sources[i][3][gid]})
                    continue
                for i, reader in enumerate(readers):
                    for row in reader:
                        gid = row["genome_id"]
                        if gid not in selected:
                            raise InputError("Unknown genome in scientific table")
                        if selected[gid] == i:
                            writer.writerow({SOURCE_FIELD: identities[i], **row})
                        elif name != "source_status.tsv":
                            raise InputError("Superseded unavailable genome has scientific rows")
        # Revalidate snapshots after copying, detecting concurrent changes.
        for root, metadata, manifest, genome_rows in sources:
            current = _source(root)
            if current[1:] != (metadata, manifest, genome_rows):
                raise InputError("Collection snapshot changed concurrently during migration")
        counts = {
            state: sum(sources[i][3][gid]["state"] == state for gid, i in selected.items())
            for state in sorted({sources[i][3][gid]["state"] for gid, i in selected.items()})
        }
        atomic_json(
            output / "collection_metadata.json",
            {
                "schema_version": "combined-1.0.0",
                "created_at": now(),
                "migration_implementation_sha256": sha256(Path(__file__)),
                "sources": [
                    {
                        "collection": str(root),
                        "collection_metadata_sha256": sha256(root / "collection_metadata.json"),
                        "manifest_sha256": identities[i],
                        "manifest": str(manifest.output_root) + "/manifest.json",
                        "software": metadata["software"],
                        "pinned_files": metadata["pinned_files"],
                    }
                    for i, (root, metadata, manifest, _) in enumerate(sources)
                ],
                "superseded_unavailable": superseded,
                "genomes": len(selected),
                "states": counts,
                "semantics": (
                    "Native rows retained; collection_manifest_sha256 resolves original "
                    "array/software identity. No scientific rescans or classification."
                ),
                "output_sha256": {name: sha256(output / name) for name in sorted(outputs)},
            },
        )
        atomic_json(
            output / "status.json",
            {"state": "succeeded", "genomes": len(selected), "genome_states": counts},
        )
    except BaseException as exc:
        atomic_json(output / "status.json", {"state": "failed", "error": str(exc)})
        raise
    return read_json(output / "status.json")
