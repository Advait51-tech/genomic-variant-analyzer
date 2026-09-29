import argparse

import pandas as pd
import gzip


def _coerce_numeric(value):
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


# Original exploratory version — kept for reference, superseded by parse_vcf/parse_record below.
# def read_vcf(path):
#     with open(path) as f:
#         for line in f:
#             line = line.rstrip("\n")
#             if line.startswith("##"):
#                 continue  # metadata — skip for now
#             elif line.startswith("#CHROM"):
#                 columns = line.lstrip("#").split("\t")
#                 print("Columns:", columns)
#             else:
#                 fields = line.split("\t")
#                 chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
#                 print(chrom, pos, ref, alt, qual)
#                 break


def parse_info(info):
    parsed = {}
    for entry in info.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            parsed[f"info_{key}"] = _coerce_numeric(value)
        else:
            parsed[f"info_{entry}"] = True
    return parsed


def parse_record(line, columns):
    fields = line.split("\t")
    chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]

    record = {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt.split(","),
        "qual": float(qual) if qual != "." else None,
        "filter": filt,
    }
    record.update(parse_info(info))

    if len(fields) > 9:
        format_keys = fields[8].split(":")
        for sample_name, sample_str in zip(columns[9:], fields[9:]):
            values = sample_str.split(":")
            for key, value in zip(format_keys, values):
                record[f"{sample_name}_{key}"] = _coerce_numeric(value)

    return record


def parse_vcf(path):
    records = []
    columns = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                columns = line.lstrip("#").split("\t")
            else:
                records.append(parse_record(line, columns))
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Parse a VCF file into a tabular DataFrame.")
    parser.add_argument("vcf_path", help="Path to the input VCF file")
    parser.add_argument("--out", help="Write the parsed records to this CSV path instead of printing them")
    args = parser.parse_args()

    df = parse_vcf(args.vcf_path)

    if args.out:
        df.to_csv(args.out, index=False)
        print(f"Wrote {len(df)} records to {args.out}")
    else:
        print(df)

    # Sanity-check diagnostics — kept for reference, re-enable when investigating the data:
    # print(df.dtypes)                          # are numeric columns actually numeric, or stuck as text?
    # print(df.isnull().sum())                  # how many missing values per column, and where
    # print(df["filter"].unique())              # any unexpected FILTER values beyond PASS?
    # print(df.assign(alt=df["alt"].apply(tuple)).duplicated(subset=["chrom", "pos", "ref", "alt"]).sum())  # any duplicate variant records?


if __name__ == "__main__":
    main()