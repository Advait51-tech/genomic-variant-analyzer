import argparse

from parser import parse_vcf
from parser_cyvcf2 import parse_vcf_cyvcf2

# Known, documented differences between the two backends -- not bugs, just API
# limitations of cyvcf2 that make byte-for-byte parity impossible on these columns.
KNOWN_DIVERGENT_COLUMNS = {
    "filter": "cyvcf2 collapses both \"PASS\" and \".\" to None; the hand-rolled parser keeps the raw text",
}
FLOAT_PRECISION_NOTE = (
    "htslib stores VCF Float-typed fields as float32 internally, so cyvcf2 loses precision "
    "relative to the hand-rolled parser's float64 text parsing; differences only appear when "
    "that rounding lands on the opposite side of a .round(3) boundary (rare, not systemic)"
)


def compare(path):
    df_manual = parse_vcf(path)
    df_cyvcf2 = parse_vcf_cyvcf2(path)

    report = {"path": path, "row_count_match": len(df_manual) == len(df_cyvcf2)}
    print(f"Rows: hand-rolled={len(df_manual)}, cyvcf2={len(df_cyvcf2)}  match={report['row_count_match']}")

    manual_cols = set(df_manual.columns)
    cyvcf2_cols = set(df_cyvcf2.columns)
    report["columns_only_in_manual"] = sorted(manual_cols - cyvcf2_cols)
    report["columns_only_in_cyvcf2"] = sorted(cyvcf2_cols - manual_cols)
    if report["columns_only_in_manual"]:
        print(f"Columns only in hand-rolled parser: {report['columns_only_in_manual']}")
    if report["columns_only_in_cyvcf2"]:
        print(f"Columns only in cyvcf2 parser:       {report['columns_only_in_cyvcf2']}")

    shared_columns = sorted(manual_cols & cyvcf2_cols)
    mismatched_columns = []
    if report["row_count_match"]:
        for col in shared_columns:
            left, right = df_manual[col], df_cyvcf2[col]
            is_float = False
            if col == "alt":
                equal = (left.apply(tuple) == right.apply(tuple)).all()
                diff_count = (~(left.apply(tuple) == right.apply(tuple))).sum()
            elif left.dtype.kind in "fc" or right.dtype.kind in "fc":
                is_float = True
                left_r, right_r = left.astype(float).round(3), right.astype(float).round(3)
                equal = left_r.equals(right_r)
                diff_count = (~(left_r.isna() & right_r.isna()) & (left_r != right_r)).sum()
            else:
                left_s, right_s = left.astype(str), right.astype(str)
                equal = left_s.equals(right_s)
                diff_count = (left_s != right_s).sum()
            if not equal:
                mismatched_columns.append((col, diff_count, is_float))

    report["mismatched_columns"] = [col for col, _, _ in mismatched_columns]
    if mismatched_columns:
        print(f"\nColumns with differing values ({len(mismatched_columns)}):")
        for col, diff_count, is_float in mismatched_columns:
            if col in KNOWN_DIVERGENT_COLUMNS:
                note = KNOWN_DIVERGENT_COLUMNS[col]
            elif is_float:
                note = FLOAT_PRECISION_NOTE
            else:
                note = "unexplained -- worth investigating"
            print(f"  {col}: {diff_count}/{len(df_manual)} rows differ -- {note}")
    else:
        print("\nAll shared columns match (within float rounding tolerance).")

    return report


def main():
    arg_parser = argparse.ArgumentParser(
        description="Compare the hand-rolled parser (parser.py) against the cyvcf2-based parser for output parity."
    )
    arg_parser.add_argument("vcf_path", help="Path to the input VCF file")
    args = arg_parser.parse_args()
    compare(args.vcf_path)


if __name__ == "__main__":
    main()
