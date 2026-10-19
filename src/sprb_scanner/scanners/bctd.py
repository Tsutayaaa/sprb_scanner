"""Independent Type B CTD HMM evidence; no score fusion or final classification."""

from sprb_scanner.scanners.pair import run_pair

MODEL_IDS = ("TIGR04131", "PF13585")


def run_bctd(proteins, profiles, output, config=None):
    return run_pair(
        proteins,
        profiles,
        output,
        MODEL_IDS,
        "Type B CTD HMM evidence; no automatic final SprB classification",
        config,
    )
