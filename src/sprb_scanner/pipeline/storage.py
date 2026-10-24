"""Deterministic fingerprints, immutable attempt directories and atomic checkpoints."""

import hashlib
import json
import os
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

from sprb_scanner.core.exceptions import InputError
from sprb_scanner.io.common import sha256

PROTOCOL = "1.0.0"


def now():
    return datetime.now(UTC).isoformat()


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        ).encode()
    ).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


@contextmanager
def workspace_lock(root):
    import fcntl

    if (root / ".pipeline.lock").is_symlink():
        raise InputError("Pipeline lock must not be a symlink")
    with (root / ".pipeline.lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise InputError(f"Another pipeline runner holds this workspace: {root}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def inventory(root, exclude=()):
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise InputError(f"Symlinks are not allowed inside an attempt: {path}")
        if path.is_file() and str(path.relative_to(root)) not in exclude:
            files[str(path.relative_to(root))] = sha256(path)
    return files


def verify_files(root, checks):
    if not isinstance(checks, dict) or not checks:
        raise InputError("Missing attempt output checksums")
    if inventory(root, ("execution_manifest.json",)) != checks:
        raise InputError("Attempt artifacts changed or are incomplete")


def next_directory(parent, prefix):
    parent.mkdir(parents=True, exist_ok=True)
    numbers = [
        int(p.name[len(prefix) :])
        for p in parent.iterdir()
        if p.name.startswith(prefix) and p.name[len(prefix) :].isdigit()
    ]
    directory = parent / f"{prefix}{max(numbers, default=0) + 1:04d}"
    directory.mkdir()
    return directory
