import pandas as pd

def sum_table():
    # Step 1: Load the data
    normalized_counts = pd.read_csv("codon_profile.csv", index_col=0)  
    raw_counts = pd.read_csv("codon_profile_rawcounts.csv", index_col=0) 

    # Step 2: Map codons to amino acids
    codon_to_aa = {
        "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
        "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
        "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
        "TGT": "Cys", "TGC": "Cys",
        "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
        "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
        "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
        "ATG": "Met",
        "TGG": "Trp",
        "TTT": "Phe", "TTC": "Phe",
        "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
        "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
        "TAT": "Tyr", "TAC": "Tyr",
        "AAT": "Asn", "AAC": "Asn",
        "CAA": "Gln", "CAG": "Gln",
        "AAA": "Lys", "AAG": "Lys",
        "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
        "CAT": "His", "CAC": "His",
        "GAT": "Asp", "GAC": "Asp",
        "GAA": "Glu", "GAG": "Glu",
        "TAA": "*", "TAG": "*", "TGA": "*",
}

    # Step 3: Aggregate raw counts and normalized counts
    raw_totals = raw_counts.sum(axis=0)  # Total counts for each codon
    normalized_totals = normalized_counts.mean(axis=0).round(4)  # Mean normalized frequency

    # Step 4: Create a summary table
    summary_data = []
    for codon, aa in codon_to_aa.items():
        aa_codons = [c for c, a in codon_to_aa.items() if a == aa]  # Codons for this amino acid
        total_for_aa = raw_totals[aa_codons].sum()  # Total raw counts for this amino acid
        fraction = (raw_totals[codon] / total_for_aa).round(4) if total_for_aa > 0 else 0
        summary_data.append({
            "Cod": codon,
            "AA": aa,
            "Fract": fraction,
            "Frequ": normalized_totals[codon],
            "Number": raw_totals[codon],
        })

    # Convert to DataFrame
    summary_df = pd.DataFrame(summary_data)

    # Step 5: Export to CSV or TXT
    summary_df.to_csv("codon_summary.csv", index=False)
    summary_df.to_csv("codon_summary.txt", index=False, sep="\t")

if __name__ == "__main__":
    sum_table()