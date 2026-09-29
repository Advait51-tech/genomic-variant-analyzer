import argparse

from cyvcf2 import VCF

CLINVAR_GRCH38_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"


def fetch_region(url, region, out_path):
    vcf = VCF(url)
    with open(out_path, "w") as f:
        f.write(vcf.raw_header)
        count = 0
        for variant in vcf(region):
            f.write(str(variant))
            count += 1
    return count


def main():
    arg_parser = argparse.ArgumentParser(
        description="Download a genomic region from the public ClinVar VCF (GRCh38) into a local VCF file."
    )
    arg_parser.add_argument("--region", default="1:1000000-1200000", help="Region to fetch, e.g. '1:1000000-1200000' (default: %(default)s)")
    arg_parser.add_argument("--url", default=CLINVAR_GRCH38_URL, help="Remote tabix-indexed VCF URL (default: ClinVar GRCh38)")
    arg_parser.add_argument("--out", default="data/raw/clinvar_chr1_subset.vcf", help="Output VCF path")
    args = arg_parser.parse_args()

    count = fetch_region(args.url, args.region, args.out)
    print(f"Wrote {count} variants from {args.region} to {args.out}")


if __name__ == "__main__":
    main()
