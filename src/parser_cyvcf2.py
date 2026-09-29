import argparse

import pandas as pd
from cyvcf2 import VCF


def _format_gt(genotype):
    a, b, phased = genotype
    allele = lambda x: "." if x == -1 else str(x)
    return f"{allele(a)}{'|' if phased else '/'}{allele(b)}"


def _coerce_numeric(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def parse_vcf_cyvcf2(path):
    vcf = VCF(path)
    samples = vcf.samples
    records = []

    for variant in vcf:
        record = {
            "chrom": variant.CHROM,
            "pos": variant.POS,
            "ref": variant.REF,
            "alt": list(variant.ALT),
            "qual": float(variant.QUAL) if variant.QUAL is not None else None,
            # cyvcf2 collapses both "PASS" and "." to None -- the raw FILTER text isn't
            # recoverable from its API, so we can't tell these two cases apart here.
            "filter": variant.FILTER,
        }

        for key, value in variant.INFO:
            record[f"info_{key}"] = value if not isinstance(value, str) else _coerce_numeric(value)

        if samples:
            genotypes = variant.genotypes
            for sample_idx, sample_name in enumerate(samples):
                record[f"{sample_name}_GT"] = _format_gt(genotypes[sample_idx])
                for fmt_key in variant.FORMAT:
                    if fmt_key == "GT":
                        continue
                    values = variant.format(fmt_key)
                    if values is None:
                        continue
                    row = values[sample_idx]
                    if hasattr(row, "__len__") and len(row) > 1:
                        record[f"{sample_name}_{fmt_key}"] = ",".join(str(v) for v in row)
                    else:
                        record[f"{sample_name}_{fmt_key}"] = row[0].item()

        records.append(record)

    return pd.DataFrame(records)


def main():
    arg_parser = argparse.ArgumentParser(description="Parse a VCF file into a DataFrame using cyvcf2.")
    arg_parser.add_argument("vcf_path", help="Path to the input VCF file")
    arg_parser.add_argument("--out", help="Write the parsed records to this CSV path instead of printing them")
    args = arg_parser.parse_args()

    df = parse_vcf_cyvcf2(args.vcf_path)

    if args.out:
        df.to_csv(args.out, index=False)
        print(f"Wrote {len(df)} records to {args.out}")
    else:
        print(df)


if __name__ == "__main__":
    main()
