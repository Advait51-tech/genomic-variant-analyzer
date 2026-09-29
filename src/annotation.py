import argparse
from collections import Counter, defaultdict

import pandas as pd

from parser import parse_vcf


def parse_geneinfo(geneinfo):
    """GENEINFO looks like 'SYMBOL:ENTREZID' or 'GENE1:ID1|GENE2:ID2' for overlapping genes."""
    if not isinstance(geneinfo, str) or not geneinfo:
        return []
    return [entry.partition(":")[0] for entry in geneinfo.split("|")]


def parse_consequences(mc):
    """MC looks like 'SO:0001583|missense_variant' or multiple such entries joined by ','."""
    if not isinstance(mc, str) or not mc:
        return []
    return [entry.partition("|")[2] or entry for entry in mc.split(",")]


def annotate_variants(df):
    rows = []
    for _, row in df.iterrows():
        rows.append({
            "chrom": row["chrom"],
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
            "genes": parse_geneinfo(row.get("info_GENEINFO")),
            "consequences": parse_consequences(row.get("info_MC")),
            "clinical_significance": row.get("info_CLNSIG"),
        })
    return pd.DataFrame(rows)


def gene_burden_summary(annotated_df):
    variant_counts = Counter()
    consequence_counts = defaultdict(Counter)
    clnsig_counts = defaultdict(Counter)

    for _, row in annotated_df.iterrows():
        for gene in row["genes"]:
            variant_counts[gene] += 1
            for consequence in row["consequences"]:
                consequence_counts[gene][consequence] += 1
            if pd.notna(row["clinical_significance"]):
                clnsig_counts[gene][row["clinical_significance"]] += 1

    summary_rows = []
    for gene, count in variant_counts.most_common():
        top = consequence_counts[gene].most_common(1)
        pathogenic = sum(n for sig, n in clnsig_counts[gene].items() if "Pathogenic" in sig)
        summary_rows.append({
            "gene": gene,
            "variant_count": count,
            "top_consequence": top[0][0] if top else None,
            "pathogenic_or_likely_pathogenic": pathogenic,
        })

    return pd.DataFrame(summary_rows)


def print_gene_burden(summary_df, top_n=20):
    print(f"{'Gene':<16} {'Variants':>9} {'Top consequence':<32} {'Pathogenic/LP':>13}")
    print("-" * 74)
    for _, row in summary_df.head(top_n).iterrows():
        print(f"{row['gene']:<16} {row['variant_count']:>9} {str(row['top_consequence']):<32} {row['pathogenic_or_likely_pathogenic']:>13}")


def main():
    arg_parser = argparse.ArgumentParser(
        description="Map variants to genes/consequences and summarize variant burden per gene "
                     "(requires ClinVar-style GENEINFO/MC/CLNSIG INFO fields)."
    )
    arg_parser.add_argument("vcf_path", help="Path to the input VCF file")
    arg_parser.add_argument("--top", type=int, default=20, help="Show the top N genes by variant count (default: 20)")
    arg_parser.add_argument("--out", help="Write the full gene burden summary to this CSV path")
    args = arg_parser.parse_args()

    df = parse_vcf(args.vcf_path)
    if "info_GENEINFO" not in df.columns:
        arg_parser.error("This VCF has no INFO/GENEINFO field -- gene annotation currently only supports ClinVar-style VCFs")

    annotated = annotate_variants(df)
    summary = gene_burden_summary(annotated)

    print(f"Annotated {len(annotated)} variants across {len(summary)} genes\n")
    print_gene_burden(summary, top_n=args.top)

    if args.out:
        summary.to_csv(args.out, index=False)
        print(f"\nWrote full gene burden summary to {args.out}")


if __name__ == "__main__":
    main()
