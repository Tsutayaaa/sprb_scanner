"""Submit bounded waves from a small SLURM allocation; never retry errors automatically."""

import argparse
import json
import time

from sprb_scanner.batch import batch_status, submit_batch

parser = argparse.ArgumentParser()
parser.add_argument("--manifest", required=True)
parser.add_argument("--prepare-only", action="store_true")
args = parser.parse_args()
while True:
    result = submit_batch(args.manifest, prepare_only=args.prepare_only, wait=True)
    print(json.dumps(result), flush=True)
    if result["state"] == "active_wave":
        time.sleep(30)
        continue
    if result["state"] == "nothing_to_submit":
        print(json.dumps({"final_counts": batch_status(args.manifest)["counts"]}), flush=True)
        break
