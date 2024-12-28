import numpy as np
import pandas as pd
from Bio import SeqIO

codons = [
    # Nonpolar
    "GGT", "GGC", "GGA", "GGG",  # Glycine
    "GCT", "GCC", "GCA", "GCG",  # Alanine
    "GTT", "GTC", "GTA", "GTG",  # Valine
    "TGT", "TGC",                # Cysteine
    "CCT", "CCC", "CCA", "CCG",  # Proline
    "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",  # Leucin
    "ATT", "ATC", "ATA",         # Isoleucine
    "ATG",                       # Methionine
    "TGG",                       # Tryptophan
    "TTT", "TTC",                # Phenylalanine
    # Polar
    "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",  # Serine
    "ACT", "ACC", "ACA", "ACG",  # Threonine
    "TAT", "TAC",                # Tyrosine
    "AAT", "AAC",                # Asparagine
    "CAA", "CAG",                # Glutamine
    # Charged (positive)
    "AAA", "AAG",                # Lysine
    "CGT", "CGC", "CGA", "CGG", "AGA", "AGG",  # Arginine
    "CAT", "CAC",                # Histidine
    # Charged (negative)
    "GAT", "GAC",                # Aspartic Acid
    "GAA", "GAG",                # Glutaminc Acid
    # Stop codons
    "TAA", "TAG", "TGA"
]

def create_dataframe(fasta_file):
    # Initialize a list for all transcripts
    transcript_ids = []
    # Go through all the records and collect unique transcript IDs
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcript_id = record.id.split("|")[0]
            if transcript_id not in transcript_ids:
                transcript_ids.append(transcript_id)

    # Create a mapping from transcript ID to row index
    transcript_id_map = {transcript_id: idx for idx, transcript_id in enumerate(transcript_ids)}
    
    # Initialize an array to store codon counts
    codon_count_array = np.zeros((len(transcript_ids), len(codons)), dtype=int)
    
    return transcript_id_map, codon_count_array, transcript_ids

def codon_profile(fasta_file):
    transcript_id_map, codon_count_array, transcript_ids = create_dataframe(fasta_file)

    i = 0
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            i += 1
            transcript_id = record.id.split("|")[0]
            sequence = str(record.seq)

            # Extract CDS region
            cds_header = record.id.split("CDS:")[1].split("|")[0]
            start, end = map(int, cds_header.split("-"))
            cds_sequence = sequence[start-1:end]

            # Check for frameshift mutations
            if len(cds_sequence) % 3 == 0:
                # Find row index for this transcript
                transcript_idx = transcript_id_map[transcript_id]
                
                # Count codons
                for j in range(0, len(cds_sequence), 3):
                    codon = cds_sequence[j:j+3]
                    if codon in codons:
                        codon_idx = codons.index(codon)
                        codon_count_array[transcript_idx, codon_idx] += 1

            if i % 1000 == 0:
                print(f"Processed {i} sequences")

    # Convert the array to a DataFrame
    codon_count_df = pd.DataFrame(codon_count_array, columns=codons)
    codon_count_df.insert(0, "Transcript ID", transcript_ids)

    return codon_count_df

if __name__ == "__main__":
    input = "ref_transcript_IDs.fa"
    codon_profile_df = codon_profile(input)
    output = "codon_profile.csv"
    codon_profile_df.to_csv(output, index=False)
    print(f"DataFrame exported to {output}")
