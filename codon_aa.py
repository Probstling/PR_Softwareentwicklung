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

amino_acids = [
    # Nonpolar
    "G", "A", "V", "C", "P", "L", "I", "M", "W", "F",
    #Polar
    "S", "T", "Y", "N", "Q", 
    # Charged (positive)
    "K", "R", "H",
    # Charged (negative)
    "D", "E"
]


def dataframe(fasta_file, columns):
    # Initialize a list for all transcripts
    transcript_ids = []
    # Go through all the records and collect unique transcript IDs
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcript_id = record.id.split("|")[0]
            transcript_ids.append(transcript_id)

    # Create a mapping from transcript ID to row index
    transcript_id_map = {transcript_id: idx for idx, transcript_id in enumerate(transcript_ids)}
    
    # Initialize an array to store codon counts
    count_array = np.zeros((len(transcript_ids), len(columns)), dtype=float)
    
    return transcript_id_map, count_array, transcript_ids

def codon_profile(fasta_file):

    transcript_id_map, codon_count_array, transcript_ids = dataframe(fasta_file, codons)

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
                
                # Normalize codon counts by sequence length
                total_codons = len(cds_sequence)/3
                codon_count_array[transcript_idx] /= total_codons
                codon_count_array[transcript_idx] = codon_count_array[transcript_idx].round(4)

            if i % 1000 == 0:
                print(f"Processed {i} sequences")

    # Convert the array to a DataFrame
    codon_count_df = pd.DataFrame(codon_count_array, columns=codons)
    codon_count_df.insert(0, "Transcript ID", transcript_ids)

    return codon_count_df

def aa_profile(fasta_file):

    transcript_id_map, aa_count_array, transcript_ids = dataframe(fasta_file, amino_acids)

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcript_id = record.id.split("|")[0]

            # Extract CDS region
            sequence = record.seq
            cds_header = record.id.split("CDS:")[1].split("|")[0]
            start, end = map(int, cds_header.split("-"))
            cds_sequence = sequence[start-1:end]

            # Exclude frameshift mutations, count aa and normalize
            if len(cds_sequence) % 3 == 0:
                aa_sequence = cds_sequence.translate(to_stop=True)
                # Find row index for this transcript
                transcript_idx = transcript_id_map[transcript_id]

                # Count aa
                for aa in aa_sequence:
                    if aa in amino_acids:
                        aa_idx = amino_acids.index(aa)
                        aa_count_array[transcript_idx, aa_idx] += 1

                # Normalize codon counts by sequence length
                total_aa = len(aa_sequence)
                aa_count_array[transcript_idx] /= total_aa
                aa_count_array[transcript_idx] = aa_count_array[transcript_idx].round(4)

    # Convert the array to a DataFrame
    aa_count_df = pd.DataFrame(aa_count_array, columns=amino_acids)
    aa_count_df.insert(0, "Transcript ID", transcript_ids)                    

    return aa_count_df


if __name__ == "__main__":

    input = "ref_transcript_IDs.fa"

    # Codon profile
    #codon_profile_df = codon_profile(input)
    #codon_output = "codon_profile.csv"
    #codon_profile_df.to_csv(codon_output, index=False)
    #print(f"Codon profile exported to {codon_output}")

    # Amino acid profile
    aa_profile_df = aa_profile(input)
    aa_output = "aa_profile.csv"
    aa_profile_df.to_csv(aa_output, index=False)
    print(f"Amino acid profile exported to {aa_output}")
    
