import argparse
import json
import pandas as pd
from pathlib import Path

def main(
    preseq_file: str | Path,
    reads_col: str,
    fractions_json: str | Path,
    key: str,
):
    # Load fraction
    with open(fractions_json) as f:
        fractions = json.load(f)

    fraction = fractions[key]

    # Read preseq output
    df = pd.read_csv(
        preseq_file,
        sep="\t",
    )

    # Correct read depth
    df["corrected_reads"] = df[reads_col] / fraction

    # Save
    df.to_csv(preseq_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Correct preseq read depth by UMI-bearing read fraction"
    )
    parser.add_argument("--preseq_file", required=True)
    parser.add_argument("--reads_col", required=True)
    parser.add_argument("--fractions_json", required=True)
    parser.add_argument(
        "--key",
        required=True,
        choices=["included", "included_tissue"],
        help="Which fraction to use",
    )

    args = parser.parse_args()
    main(
        args.preseq_file,
        args.reads_col,
        args.fractions_json,
        args.key,
    )
