import argparse

from parser import parse_vcf

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def classify_variant(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    elif len(ref) < len(alt):
        return "insertion"
    elif len(ref) > len(alt):
        return "deletion"
    else:
        return "MNP"


def classify_genotype(gt):
    alleles = gt.replace("|", "/").split("/")
    if "." in alleles:
        return "missing"
    return "homozygous" if len(set(alleles)) == 1 else "heterozygous"


def compute_statistics(df):
    variant_counts = {"SNP": 0, "insertion": 0, "deletion": 0, "MNP": 0}
    ts_count = 0
    tv_count = 0
    genotype_counts = {"homozygous": 0, "heterozygous": 0, "missing": 0}

    sample_gt_columns = [col for col in df.columns if col.endswith("_GT")]

    for _, row in df.iterrows():
        ref = row["ref"]
        for alt in row["alt"]:
            kind = classify_variant(ref, alt)
            variant_counts[kind] += 1
            if kind == "SNP":
                if (ref, alt) in TRANSITIONS:
                    ts_count += 1
                else:
                    tv_count += 1

        for col in sample_gt_columns:
            zygosity = classify_genotype(row[col])
            genotype_counts[zygosity] += 1

    ts_tv_ratio = ts_count / tv_count if tv_count else float("nan")

    return {
        "total_variants": len(df),
        "variant_counts": variant_counts,
        "ts_count": ts_count,
        "tv_count": tv_count,
        "ts_tv_ratio": ts_tv_ratio,
        "genotype_counts": genotype_counts,
    }


def print_summary(stats):
    width = 40
    print("=" * width)
    print("VARIANT SUMMARY")
    print("=" * width)
    print(f"{'Total variants':<28}{stats['total_variants']:>12}")
    print("-" * width)
    print(f"{'SNPs':<28}{stats['variant_counts']['SNP']:>12}")
    print(f"{'Insertions':<28}{stats['variant_counts']['insertion']:>12}")
    print(f"{'Deletions':<28}{stats['variant_counts']['deletion']:>12}")
    if stats["variant_counts"]["MNP"]:
        print(f"{'MNPs (other)':<28}{stats['variant_counts']['MNP']:>12}")
    print("-" * width)
    print(f"{'Transitions (Ts)':<28}{stats['ts_count']:>12}")
    print(f"{'Transversions (Tv)':<28}{stats['tv_count']:>12}")
    ratio = stats["ts_tv_ratio"]
    ratio_str = f"{ratio:.2f}" if ratio == ratio else "N/A"
    print(f"{'Ts/Tv ratio':<28}{ratio_str:>12}")
    print("-" * width)
    print(f"{'Homozygous genotypes':<28}{stats['genotype_counts']['homozygous']:>12}")
    print(f"{'Heterozygous genotypes':<28}{stats['genotype_counts']['heterozygous']:>12}")
    if stats["genotype_counts"]["missing"]:
        print(f"{'Missing genotypes':<28}{stats['genotype_counts']['missing']:>12}")
    print("=" * width)


def main():
    arg_parser = argparse.ArgumentParser(
        description="Compute SNP/indel, Ts/Tv, and zygosity statistics for a VCF file."
    )
    arg_parser.add_argument("vcf_path", help="Path to the input VCF file")
    args = arg_parser.parse_args()

    df = parse_vcf(args.vcf_path)
    stats = compute_statistics(df)
    print_summary(stats)


if __name__ == "__main__":
    main()
