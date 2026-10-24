"""Create a version-aware union of immutable batch collection snapshots."""

import argparse
import json

from sprb_scanner.batch.combine import combine_collections
from sprb_scanner.core.exceptions import InputError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection", action="append", required=True, help="Oldest first")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(combine_collections(args.collection, args.output), indent=2))
    except (InputError, OSError, ValueError) as exc:
        parser.exit(2, f"Collection migration failed: {exc}\n")


if __name__ == "__main__":
    main()
