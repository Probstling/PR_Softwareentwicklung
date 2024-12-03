from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Length of entire transcript
def length():
    array_length = []
    sequences = SeqIO.parse('ref_transcript_IDs.fa', 'fasta')
    for record in sequences:
        array_length.append(record.id.split('|')[6])
    array_length_np = np.array(array_length)
    #print(len(array_length_np))
    return array_length_np

# Length of the 5' UTR region (not every transcript has one defined)
def fiveprime():
    array_fiveprime = []
    sequences = SeqIO.parse('ref_transcript_IDs.fa', 'fasta')
    for record in sequences:
        id_part = record.id.split('|')[7]
        if id_part.startswith("UTR5:"):     
            utr5_length = int(id_part.split("-")[1]) 
            array_fiveprime.append(utr5_length)
    array_fiveprime_np = np.array(array_fiveprime)
    #print(len(array_fiveprime_np))
    return array_fiveprime_np

#Length of the 3' UTR region (not every transcript has one defined)
def threeprime():
    array_threeprime = []
    sequences = SeqIO.parse('ref_transcript_IDs.fa', 'fasta')
    for record in sequences:
        id_parts = record.id.split('|')
        if len(id_parts) > 9: 
            utr3 = record.id.split('|')[9]
            if utr3.startswith("UTR3:"):     
                utr3_range = utr3.split(":")[1]
                utr3_start, utr3_end = map(int, utr3_range.split("-"))
                utr3_length = utr3_end - utr3_start + 1 
                array_threeprime.append(utr3_length)
    array_threeprime_np = np.array(array_threeprime)
    #print(len(array_threeprime_np))
    return array_threeprime_np

#Plotting the length arrays
def plot_length(array, title, xlabel):
    plt.figure(figsize=(8, 5))
    sns.histplot(array, kde=True, color="skyblue", bins=30)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
   total_length = length()
   utr5_length = fiveprime()
   utr3_length = threeprime()
   plot_length(total_length, "Distribution of Total Transcript Lengths", "Transcript Length (bp)")
   #plot_length(utr5_length, "Distribution of 5' UTR Lengths", "5' UTR Length (bp)")
   #plot_length(utr3_length, "Distribution of 3' UTR Lengths", "3' UTR Length (bp)")