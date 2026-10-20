"""Unified offline evidence aggregation without score fusion or new searches."""

from sprb_scanner.evidence.merger import merge
from sprb_scanner.evidence.models import EvidenceProfile, EvidenceResult, SourceRun
from sprb_scanner.evidence.service import run_merge

__all__ = ["EvidenceProfile", "EvidenceResult", "SourceRun", "merge", "run_merge"]
