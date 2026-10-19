"""Independent PorP/SprF family evidence via existing HMM engines."""

from sprb_scanner.scanners.pair import collect_evidence as collect_evidence
from sprb_scanner.scanners.pair import run_pair

MODEL_IDS = ("TIGR03519", "PF11751")


def run_sprf(proteins, profiles, output, config=None):
    return run_pair(
        proteins,
        profiles,
        output,
        MODEL_IDS,
        "family evidence only; no substrate pairing or final SprB classification",
        config,
    )
