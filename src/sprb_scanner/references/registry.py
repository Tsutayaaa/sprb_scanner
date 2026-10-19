"""Offline artifact resolution, type inspection and adapters to existing scanners."""

import csv
import hashlib
import io
import json
from pathlib import Path

from Bio import AlignIO, SeqIO

from sprb_scanner.core.config import ProfileConfig
from sprb_scanner.core.exceptions import ReferenceRegistryError, ReferenceUnavailableError
from sprb_scanner.core.models import ModelProvenance, ProteinRecord, Record, internal_id
from sprb_scanner.references.models import Reference, RegistryDocument

AMINO = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def json_document(data):
    def invalid_constant(value):
        raise ValueError(f"Non-finite JSON constant: {value}")

    return json.loads(data, object_pairs_hook=unique_object, parse_constant=invalid_constant)


class ResolvedReference(Record):
    reference: Reference
    artifact_path: Path
    summary: dict

    def as_profile(self):
        """A validated HMM becomes the same ProfileConfig consumed by the HMM engine."""
        entry = self.reference
        if entry.resource_type != "profile_hmm":
            raise ReferenceRegistryError("Only profile_hmm resources can become ProfileConfig")
        model_id = (
            entry.biological_target.model_id
            or (self.summary["accession"] or self.summary["model_name"]).split(".")[0]
        )
        return ProfileConfig(
            method="hmm",
            path=self.artifact_path,
            provenance=ModelProvenance(
                model_id=model_id,
                # Resource release labels and HMM accession versions are distinct identities.
                version=(
                    self.summary["accession"]
                    if self.summary["accession"] and "." in self.summary["accession"]
                    else entry.version
                ),
                source=entry.source.locator,
                sha256=entry.sha256,
            ),
        )


def inspect_artifact(entry, data):
    """Inspect one immutable byte snapshot; never rewrite residues, alignments or GA."""
    fmt, kind, target = entry.artifact.format, entry.resource_type, entry.biological_target
    if kind == "profile_hmm":
        import pyhmmer

        with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as handle:
            models = list(handle)
        if len(models) != 1:
            raise ValueError("Expected exactly one profile HMM")
        hmm = models[0]
        hmm.validate()
        if not hmm.alphabet.is_amino():
            raise ValueError("Expected an amino-acid HMM")

        def text(value):
            return value.decode() if isinstance(value, bytes) else value

        accession, name = text(hmm.accession), text(hmm.name)
        if target.model_id and (accession or name).split(".")[0] != target.model_id:
            raise ValueError("HMM model identity mismatch")
        if target.model_accession and accession != target.model_accession:
            raise ValueError("HMM accession/version mismatch")
        return {
            "model_name": name,
            "accession": accession,
            "model_length": hmm.M,
            "gathering_cutoffs": hmm.cutoffs.gathering,
        }, None
    text = data.decode("utf-8")
    if kind == "annotation_catalog":
        if fmt == "json":
            obj = json_document(text)
            if not isinstance(obj, (dict, list)):
                raise ValueError("Annotation JSON must be an object or array")
            return {"json_type": type(obj).__name__, "entries": len(obj)}, None
        rows = list(csv.reader(io.StringIO(text), delimiter="\t", strict=True))
        if (
            not rows
            or len(rows[0]) < 2
            or any(not s for s in rows[0])
            or len(set(rows[0])) != len(rows[0])
            or any(len(row) != len(rows[0]) for row in rows[1:])
        ):
            raise ValueError("Malformed TSV catalog: expected unique header and rectangular rows")
        return {"columns": rows[0], "entries": len(rows) - 1}, None
    if not text.startswith(">") and fmt != "stockholm":
        raise ValueError("Protein FASTA/MSA must start with a header")
    if fmt == "stockholm":
        alignment = AlignIO.read(io.StringIO(text), "stockholm")
        rows = list(alignment)
    else:
        rows = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if not rows or any(not r.id or not r.seq for r in rows):
        raise ValueError("Empty sequence/alignment or record")
    ids = tuple(r.id for r in rows)
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate FASTA/MSA sequence IDs")
    if target.sequence_ids and ids != target.sequence_ids:
        raise ValueError("Sequence IDs/order disagree with biological_target")
    if kind in {"sequence", "sequence_set"}:
        if kind == "sequence" and len(rows) != 1:
            raise ValueError("A sequence resource requires exactly one protein")
        proteins = tuple(
            ProteinRecord(
                protein_id=rec.id,
                source_id=rec.id,
                description=rec.description,
                genome_id=internal_id("reference", entry.reference_id, entry.version),
                sequence=str(rec.seq),
            )
            for rec in rows
        )
        if kind == "sequence":
            sequence = proteins[0].sequence
            if (
                target.sequence_sha256
                and hashlib.sha256(sequence.encode()).hexdigest() != target.sequence_sha256
            ):
                raise ValueError("Single-sequence checksum mismatch")
            if (
                target.protein_interval
                and len(sequence) != target.protein_interval.end - target.protein_interval.start
            ):
                raise ValueError("Sequence length disagrees with parent interval")
        return {
            "sequence_count": len(rows),
            "sequence_ids": ids,
            "lengths": [p.length for p in proteins],
        }, proteins
    # No alignment or HMM construction: inspect gap/insert conventions in existing MSAs.
    sequences = [str(r.seq) for r in rows]
    if any(set(seq.upper()) - AMINO - {"-", "."} for seq in sequences):
        raise ValueError("Invalid protein residues in MSA")
    if any(not set(seq.upper()) & AMINO for seq in sequences):
        raise ValueError("MSA contains an all-gap sequence")
    match_sequences = (
        ["".join(c for c in seq if not c.islower() and c != ".") for seq in sequences]
        if fmt == "a3m"
        else sequences
    )
    lengths = {len(seq) for seq in match_sequences}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ValueError("MSA match-column lengths disagree or are empty")
    return {
        "sequence_count": len(rows),
        "sequence_ids": ids,
        "match_columns": next(iter(lengths)),
        "format": fmt,
    }, None


class ReferenceRegistry:
    def __init__(self, path, document, registry_sha256):
        self.path = Path(path).expanduser().resolve()
        self.root = self.path.parent
        self.document = document
        self.registry_sha256 = registry_sha256
        self._entries = {entry.key: entry for entry in document.references}

    @classmethod
    def load(cls, path):
        path = Path(path).expanduser()
        try:
            path = path.resolve()
            data = path.read_bytes()
            return cls(
                path,
                RegistryDocument.model_validate(json_document(data)),
                hashlib.sha256(data).hexdigest(),
            )
        except (OSError, ValueError, RuntimeError) as exc:
            raise ReferenceRegistryError(f"Cannot load reference registry {path}: {exc}") from exc

    def get(self, reference_id, version=None):
        matches = [
            entry
            for key, entry in self._entries.items()
            if key[0] == reference_id and (version is None or key[1] == version)
        ]
        if not matches:
            raise ReferenceRegistryError(f"Unknown reference: {reference_id}@{version or '*'}")
        if len(matches) != 1:
            raise ReferenceRegistryError(f"Ambiguous reference {reference_id}; specify --version")
        return matches[0]

    def lineage(self, reference_id, version=None):
        """Return the transitive pinned ancestors once, without requiring parent artifact bytes."""
        entry, ancestors, seen = self.get(reference_id, version), [], set()
        pending = list(reversed(entry.derived_from))
        while pending:
            parent = pending.pop()
            if parent.key not in seen:
                seen.add(parent.key)
                ancestors.append(self._entries[parent.key])
                pending.extend(reversed(self._entries[parent.key].derived_from))
        return tuple(ancestors)

    def _resolve(self, reference_id, version=None, artifact=None):
        entry = self.get(reference_id, version)
        label = f"{entry.reference_id}@{entry.version}"
        if artifact is not None and entry.artifact.storage != "external":
            raise ReferenceRegistryError(f"{label}: artifact override requires external storage")
        try:
            path = (
                Path(artifact).expanduser()
                if artifact is not None
                else self.root / entry.artifact.path
            ).resolve()
        except (OSError, RuntimeError) as exc:
            raise ReferenceRegistryError(f"{label}: cannot resolve artifact path: {exc}") from exc
        if entry.artifact.storage == "bundled" and not path.is_relative_to(self.root):
            raise ReferenceRegistryError(f"{label}: bundled artifact escapes registry directory")
        if not path.exists():
            raise ReferenceUnavailableError(
                f"{label}: missing {entry.artifact.storage} artifact: {path}; no automatic download"
            )
        try:
            if not path.is_file():
                raise ValueError("Artifact path is not a file")
            data = path.read_bytes()
            digest = hashlib.sha256(data).hexdigest()
            if digest != entry.sha256:
                raise ValueError(f"SHA-256 mismatch: expected {entry.sha256}, got {digest}")
            summary, proteins = inspect_artifact(entry, data)
            return ResolvedReference(reference=entry, artifact_path=path, summary=summary), proteins
        except (OSError, ValueError, EOFError, csv.Error, RuntimeError) as exc:
            raise ReferenceRegistryError(f"{label}: invalid artifact {path}: {exc}") from exc

    def resolve(self, reference_id, version=None, artifact=None):
        return self._resolve(reference_id, version, artifact)[0]

    def sequences(self, reference_id, version=None, artifact=None):
        entry = self.get(reference_id, version)
        if entry.resource_type not in {"sequence", "sequence_set"}:
            raise ReferenceRegistryError("ProteinRecord access requires sequence or sequence_set")
        return self._resolve(reference_id, version, artifact)[1]

    def to_profile(self, reference_id, version=None, artifact=None):
        return self.resolve(reference_id, version, artifact).as_profile()

    def validate(self, allow_missing_external=True, artifact_overrides=None):
        overrides = artifact_overrides or {}
        if set(overrides) - set(self._entries):
            raise ReferenceRegistryError("Artifact overrides contain unknown reference identities")
        results = []
        for entry in self.document.references:
            result = {
                "reference_id": entry.reference_id,
                "version": entry.version,
                "resource_type": entry.resource_type,
                "sha256": entry.sha256,
            }
            try:
                resolved = self.resolve(*entry.key, artifact=overrides.get(entry.key))
                result.update(
                    state="validated",
                    artifact_path=str(resolved.artifact_path),
                    summary=resolved.summary,
                )
            except ReferenceUnavailableError as exc:
                permitted = allow_missing_external and entry.artifact.storage == "external"
                result.update(state="missing_external" if permitted else "failed", message=str(exc))
            except ReferenceRegistryError as exc:
                result.update(state="failed", message=str(exc))
            results.append(result)
        states = {row["state"] for row in results}
        return {
            "state": "failed"
            if "failed" in states
            else "partial"
            if "missing_external" in states
            else "validated",
            "registry": str(self.path),
            "registry_sha256": self.registry_sha256,
            "resource_count": len(results),
            "resources": results,
            "missing_external_policy": "permitted" if allow_missing_external else "required",
        }
