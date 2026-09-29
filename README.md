# Genomic Variant Analyzer

A Python pipeline for parsing, filtering, and analyzing genomic variants from VCF files —
built from a hand-rolled parser up through real public clinical data, a validated
`cyvcf2`-based parser, and gene-level annotation.

## What it does

- **Parses** VCF files into pandas DataFrames, including INFO and per-sample FORMAT fields
- **Filters** variants by configurable QUAL / depth / allele-frequency thresholds
- **Classifies** variants (SNP / insertion / deletion), computes Ts/Tv ratio and zygosity
- **Visualizes** quality, chromosome, SNP/indel, Ts/Tv, allele-frequency, and depth distributions
- **Analyzes sequences**: GC content, composition, reverse complement, transcription,
  translation, ORF finding, and coding-sequence variant effect classification
  (missense / synonymous / nonsense)
- **Validates against real data**: pulls a real chromosome-1 slice of the public ClinVar
  dataset and re-runs the pipeline against it
- **Cross-checks parsing correctness** against `cyvcf2` (a standard htslib-backed VCF library)
- **Annotates variants by gene** using ClinVar's built-in gene/consequence/significance fields
  and summarizes variant burden per gene

## Project structure

```text
src/
  parser.py            Hand-rolled VCF parser -> pandas DataFrame
  parser_cyvcf2.py      Same output schema, built on cyvcf2
  compare_parsers.py    Parity check between the two parsers
  filtering.py           Quality filtering (QUAL / DP / AF thresholds)
  statistics.py           SNP/indel counts, Ts/Tv ratio, zygosity summary
  visualization.py       Saves summary plots as PNGs
  sequence_analysis.py    GC content, reverse complement, transcription/translation,
                          ORF finding, variant effect classification
  annotation.py           Gene-level annotation and variant burden summary
  fetch_real_data.py      Downloads a real region from the public ClinVar VCF
  prioritization.py       Not yet implemented
data/raw/
  ex.vcf                          Small hand-built toy VCF (5 variants, 3 samples)
  clinvar_chr1_subset.vcf         Real ClinVar data, chr1:1,000,000-1,200,000 (3,013 variants)
reports/figures/                 Plots generated from the toy VCF
reports/figures_real/            Plots generated from the real ClinVar subset
notebooks/, tests/, reports/analysis_report.html   Scaffolding, not yet filled in
```

## Setup

Requires Python 3.10+.

```bash
git clone <repo-url>
cd genomic_variant_analyzer
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage

All commands are run from the project root.

### Parse a VCF

```bash
python3 src/parser.py data/raw/ex.vcf                  # print parsed DataFrame
python3 src/parser.py data/raw/ex.vcf --out out.csv     # write to CSV
```

### Quality filtering

```bash
python3 src/filtering.py data/raw/ex.vcf
python3 src/filtering.py data/raw/ex.vcf --qual-min 50 --dp-min 15 --af-min 0.1
python3 src/filtering.py data/raw/ex.vcf --config thresholds.json --out filtered.csv
```

`thresholds.json` looks like `{"qual_min": 30, "dp_min": 10, "af_min": 0.05}`. CLI flags
override the config file, which overrides the defaults. On VCFs that lack a given field
(e.g. ClinVar has no per-sample DP), that threshold is skipped with a printed note instead
of failing.

### Classification & statistics

```bash
python3 src/statistics.py data/raw/ex.vcf
```

Prints SNP/insertion/deletion counts, transition/transversion ratio, and
homozygous/heterozygous genotype counts.

### Visualizations

```bash
python3 src/visualization.py data/raw/ex.vcf                    # saves 6 PNGs to reports/figures/
python3 src/visualization.py data/raw/ex.vcf --out-dir some/dir
```

### Sequence analysis

```bash
python3 src/sequence_analysis.py --seq ACCGGGTTTN                              # length, GC%, composition, reverse complement
python3 src/sequence_analysis.py --fasta path/to/file.fasta

python3 src/sequence_analysis.py --seq ATGGCCATTGTAATGGGCCGCTGA --transcribe
python3 src/sequence_analysis.py --seq ATGGCCATTGTAATGGGCCGCTGA --translate
python3 src/sequence_analysis.py --seq ATGGCCATTGTAATGGGCCGCTGA --find-orfs --min-orf-length 5

python3 src/sequence_analysis.py --seq AAAGCCTTT --variant-pos 1 --variant-alt G   # missense/synonymous/nonsense report
```

### Fetching real data

```bash
python3 src/fetch_real_data.py --region 1:1000000-1200000 --out data/raw/clinvar_chr1_subset.vcf
```

Streams a genomic region directly from the public ClinVar VCF (GRCh38) over a remote tabix
query, so it doesn't require downloading the full ~185 MB file.

### Robust parsing (cyvcf2) and parity check

```bash
python3 src/parser_cyvcf2.py data/raw/clinvar_chr1_subset.vcf --out out.csv
python3 src/compare_parsers.py data/raw/clinvar_chr1_subset.vcf
```

The parity check confirms the hand-rolled and `cyvcf2` parsers agree on every shared column,
with two documented, non-systemic exceptions: `cyvcf2` collapses both `"PASS"` and `"."`
FILTER values to `None` (an API limitation, not a bug), and it stores VCF `Float`-typed
fields as float32 internally, which occasionally lands on the other side of a rounding
boundary versus the hand-rolled parser's float64 text parsing.

### Gene-level annotation

```bash
python3 src/annotation.py data/raw/clinvar_chr1_subset.vcf --top 20 --out gene_burden.csv
```

Maps variants to genes and molecular consequence via ClinVar's `GENEINFO`/`MC`/`CLNSIG`
INFO fields and prints a per-gene variant burden table. Requires a ClinVar-style VCF (errors
clearly if `GENEINFO` isn't present).

## Known limitations

- `src/prioritization.py`, `tests/test_parser.py`, `notebooks/exploratory_analysis.ipynb`,
  and `reports/analysis_report.html` are placeholders and not yet implemented.
- Gene annotation currently only works on ClinVar-style VCFs (relies on ClinVar's INFO
  fields rather than a general-purpose GTF/consequence annotator).
- `sequence_analysis.py`'s variant-position argument is 0-based relative to whatever coding
  sequence you pass in, not a genomic VCF `POS` — mapping genomic coordinates to coding-sequence
  offsets would require exon/CDS boundaries, which this project doesn't currently model.
