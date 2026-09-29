import argparse
from collections import Counter

# IUPAC nucleotide complement table (covers ambiguity codes, not just ACGT).
COMPLEMENT = str.maketrans(
    "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
    "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
)

# Standard genetic code, keyed by RNA codon. "*" marks a stop codon.
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {"UAA", "UAG", "UGA"}
START_CODON = "AUG"

AMINO_ACID_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate",
    "C": "Cysteine", "E": "Glutamate", "Q": "Glutamine", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
    "*": "Stop", "X": "Unknown",
}


def sequence_length(seq):
    return len(seq)


def gc_content(seq):
    if not seq:
        return 0.0
    gc = sum(1 for base in seq.upper() if base in "GC")
    return gc / len(seq) * 100


def nucleotide_composition(seq):
    return dict(Counter(seq.upper()))


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def transcribe(dna_seq):
    return dna_seq.replace("T", "U").replace("t", "u")


def translate(seq, to_stop=True):
    rna = seq.upper().replace("T", "U")
    protein = []
    for i in range(0, len(rna) - 2, 3):
        amino_acid = CODON_TABLE.get(rna[i:i + 3], "X")
        if amino_acid == "*":
            if to_stop:
                break
            protein.append("*")
        else:
            protein.append(amino_acid)
    return "".join(protein)


def find_orfs(seq, min_protein_length=1):
    """Scan all 3 forward and 3 reverse-complement reading frames for ATG...stop ORFs.

    Each ATG start is paired with the first in-frame stop codon that follows it
    (the standard longest-ORF-per-start convention); starts with no downstream
    in-frame stop are not reported.
    """
    orfs = []
    for strand_label, strand_seq in (("+", seq.upper()), ("-", reverse_complement(seq.upper()))):
        rna = strand_seq.replace("T", "U")
        for frame in range(3):
            i = frame
            while i + 3 <= len(rna):
                if rna[i:i + 3] == START_CODON:
                    protein_codons = []
                    j = i
                    while j + 3 <= len(rna):
                        codon = rna[j:j + 3]
                        if codon in STOP_CODONS:
                            protein = "".join(protein_codons)
                            if len(protein) >= min_protein_length:
                                orfs.append({
                                    "strand": strand_label,
                                    "frame": frame + 1,
                                    "start": i,
                                    "end": j + 3,
                                    "length_nt": j + 3 - i,
                                    "protein": protein,
                                })
                            break
                        protein_codons.append(CODON_TABLE.get(codon, "X"))
                        j += 3
                i += 3
    orfs.sort(key=lambda orf: orf["length_nt"], reverse=True)
    return orfs


def get_codon_at(seq, pos):
    codon_start = (pos // 3) * 3
    codon = seq[codon_start:codon_start + 3]
    if len(codon) != 3:
        raise ValueError(f"Position {pos} falls in an incomplete codon (sequence too short to translate it)")
    return codon, pos // 3, pos - codon_start


def classify_mutation(ref_aa, alt_aa):
    if ref_aa == alt_aa:
        return "synonymous"
    if alt_aa == "*":
        return "nonsense"
    if ref_aa == "*":
        return "stop-loss"
    return "missense"


def analyze_variant_effect(ref_seq, pos, alt_base):
    """Classify the effect of substituting alt_base at 0-based position pos in ref_seq."""
    ref_codon, codon_index, offset_in_codon = get_codon_at(ref_seq, pos)
    if ref_codon[offset_in_codon].upper() == alt_base.upper():
        raise ValueError(f"Alt base {alt_base!r} matches the reference base at position {pos}")

    alt_codon = ref_codon[:offset_in_codon] + alt_base.upper() + ref_codon[offset_in_codon + 1:]

    ref_aa = translate(ref_codon, to_stop=False)
    alt_aa = translate(alt_codon, to_stop=False)

    return {
        "position": pos,
        "codon_index": codon_index,
        "offset_in_codon": offset_in_codon,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "mutation_type": classify_mutation(ref_aa, alt_aa),
    }


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def summarize_sequence(seq):
    return {
        "length": sequence_length(seq),
        "gc_content": gc_content(seq),
        "composition": nucleotide_composition(seq),
        "reverse_complement": reverse_complement(seq),
    }


def print_summary(name, summary):
    print(f"=== {name} ===")
    print(f"Length:              {summary['length']}")
    print(f"GC content:          {summary['gc_content']:.2f}%")
    print(f"Nucleotide counts:   {summary['composition']}")
    print(f"Reverse complement:  {summary['reverse_complement']}")
    print()


def print_orfs(name, orfs):
    print(f"=== {name}: {len(orfs)} ORF(s) found ===")
    for orf in orfs:
        print(
            f"  strand={orf['strand']} frame={orf['frame']} "
            f"start={orf['start']} end={orf['end']} length_nt={orf['length_nt']} "
            f"protein={orf['protein']}"
        )
    print()


def print_variant_effect(name, effect):
    ref_name = AMINO_ACID_NAMES.get(effect["ref_aa"], effect["ref_aa"])
    alt_name = AMINO_ACID_NAMES.get(effect["alt_aa"], effect["alt_aa"])
    label = effect["mutation_type"].capitalize()
    print(f"=== {name}: position {effect['position']} ===")
    print(f"Reference codon: {effect['ref_codon']} ({ref_name})")
    print(f"Variant codon:   {effect['alt_codon']} ({alt_name})")
    print(f"{label} mutation: {ref_name} → {alt_name}")
    print()


def main():
    arg_parser = argparse.ArgumentParser(
        description="Sequence statistics and central dogma operations (transcription, translation, ORF finding)."
    )
    group = arg_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--seq", help="A raw nucleotide sequence to analyze")
    group.add_argument("--fasta", help="Path to a FASTA file; every record is analyzed")
    arg_parser.add_argument("--transcribe", action="store_true", help="Print the DNA -> RNA transcript")
    arg_parser.add_argument("--translate", action="store_true", help="Print the RNA -> protein translation")
    arg_parser.add_argument("--find-orfs", action="store_true", help="Print ORFs found across all 6 reading frames")
    arg_parser.add_argument("--min-orf-length", type=int, default=1, help="Minimum ORF protein length in amino acids (default: 1)")
    arg_parser.add_argument("--variant-pos", type=int, help="0-based position of a substitution to analyze against the reference sequence")
    arg_parser.add_argument("--variant-alt", help="Alternate base substituted at --variant-pos")
    args = arg_parser.parse_args()

    if bool(args.variant_pos is not None) != bool(args.variant_alt):
        arg_parser.error("--variant-pos and --variant-alt must be given together")

    records = {"sequence": args.seq} if args.seq else read_fasta(args.fasta)
    any_operation = args.transcribe or args.translate or args.find_orfs or args.variant_alt

    for name, seq in records.items():
        if not any_operation:
            print_summary(name, summarize_sequence(seq))
            continue
        if args.transcribe:
            print(f"=== {name}: transcript ===\n{transcribe(seq)}\n")
        if args.translate:
            print(f"=== {name}: protein ===\n{translate(seq)}\n")
        if args.find_orfs:
            print_orfs(name, find_orfs(seq, min_protein_length=args.min_orf_length))
        if args.variant_alt:
            print_variant_effect(name, analyze_variant_effect(seq, args.variant_pos, args.variant_alt))


if __name__ == "__main__":
    main()
