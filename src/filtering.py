import argparse
import json

import pandas as pd

from parser import parse_vcf

DEFAULT_THRESHOLDS = {
    "qual_min": 100,
    "dp_min": 10,
    "af_min": 0.05,
}


def load_thresholds(config_path=None, overrides=None):
    thresholds = dict(DEFAULT_THRESHOLDS)
    if config_path:
        with open(config_path) as f:
            thresholds.update(json.load(f))
    if overrides:
        thresholds.update({k: v for k, v in overrides.items() if v is not None})
    return thresholds


def filter_variants(df, qual_min=30, dp_min=10, af_min=0.05):
    thresholds = {"qual": qual_min, "info_DP": dp_min, "info_AF": af_min}
    mask = pd.Series(True, index=df.index)
    skipped = []
    for column, threshold in thresholds.items():
        if column in df.columns:
            mask &= df[column] > threshold
        else:
            skipped.append(column)
    if skipped:
        print(f"Skipping filters for columns not present in this VCF: {', '.join(skipped)}")
    return df[mask]


def main():
    arg_parser = argparse.ArgumentParser(description="Quality-filter parsed VCF records.")
    arg_parser.add_argument("vcf_path", help="Path to the input VCF file")
    arg_parser.add_argument("--config", help="Path to a JSON file with threshold overrides, e.g. {\"qual_min\": 30}")
    arg_parser.add_argument("--qual-min", type=float, help=f"Minimum QUAL (default {DEFAULT_THRESHOLDS['qual_min']})")
    arg_parser.add_argument("--dp-min", type=float, help=f"Minimum INFO DP (default {DEFAULT_THRESHOLDS['dp_min']})")
    arg_parser.add_argument("--af-min", type=float, help=f"Minimum INFO AF (default {DEFAULT_THRESHOLDS['af_min']})")
    arg_parser.add_argument("--out", help="Write filtered records to this CSV path instead of printing them")
    args = arg_parser.parse_args()

    thresholds = load_thresholds(
        args.config,
        {"qual_min": args.qual_min, "dp_min": args.dp_min, "af_min": args.af_min},
    )

    df = parse_vcf(args.vcf_path)
    filtered = filter_variants(df, **thresholds)

    print(f"Variants before filtering: {len(df)}")
    print(
        f"Variants after filtering:  {len(filtered)} "
        f"(QUAL > {thresholds['qual_min']}, DP > {thresholds['dp_min']}, AF > {thresholds['af_min']})"
    )

    if args.out:
        filtered.to_csv(args.out, index=False)
        print(f"Wrote {len(filtered)} records to {args.out}")
    else:
        print(filtered)


if __name__ == "__main__":
    main()
