from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Length of entire transcript
def length():
    array_length = []
    sequences = SeqIO.parse('ref_transcript_IDs.fa', 'fasta')
    for record in sequences:
        array_length.append(int(record.id.split('|')[6]))
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

#Length of the coding sequence
def cds_length():
    array_cds_length = []
    sequences = SeqIO.parse('ref_transcript_IDs.fa', 'fasta')
    for record in sequences:
        id_part = record.id.split('|')[7]
        if id_part.startswith("CDS:"):
            cds_range = id_part.split(':')[1]
            cds_start, cds_end = map(int, cds_range.split('-'))
            cds_length = cds_end - cds_start + 1
            array_cds_length.append(cds_length)
        else:
            id_part = record.id.split('|')[8]
            if id_part.startswith("CDS:"):
                cds_range = id_part.split(':')[1]
                cds_start, cds_end = map(int, cds_range.split('-'))
                cds_length = cds_end - cds_start + 1
                array_cds_length.append(cds_length)
    array_cds_length_np = np.array(array_cds_length)
    #print(len(array_cds_length_np))
    return array_cds_length_np

#Plotting the length arrays in a violin plot
def violin_all_lengths(total, utr5, cds, utr3):
    # Create category labels
    total_labels = np.array(["Total"] * len(total))
    utr5_labels = np.array(["5' UTR"] * len(utr5))
    cds_labels = np.array(["CDS"] * len(cds))
    utr3_labels = np.array(["3' UTR"] * len(utr3))

    # Concatenate data and labels
    all_lengths = np.concatenate([total, utr5, cds, utr3])
    all_labels = np.concatenate([total_labels, utr5_labels, cds_labels, utr3_labels])

    # Create a DataFrame
    data = pd.DataFrame({"Region": all_labels, "Length (bp)": all_lengths})

    # Function to remove outliers based on IQR
    def remove_outliers(df, column):
        Q1 = df[column].quantile(0.25)
        Q3 = df[column].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

    # Get the number of data points before removing outliers
    total_data_points = len(data)

    # Remove outliers from each region
    filtered_data = remove_outliers(data, "Length (bp)")

    # Get the number of data points after removing outliers
    filtered_data_points = len(filtered_data)

    # Calculate the percentage of data excluded
    excluded_percentage = ((total_data_points - filtered_data_points) / total_data_points) * 100

    # Print the percentage of data excluded
    print(f"Percentage of data excluded: {excluded_percentage:.2f}%")

    # Plot
    sns.violinplot(data = filtered_data, x="Region", y="Length (bp)", palette="muted")
    plt.title("Violin Plot of Transcript Lengths by Region")
    plt.show()

if __name__ == "__main__":
   total = length()
   utr5 = fiveprime()
   utr3 = threeprime()
   cds = cds_length()
   violin_all_lengths(total, utr5, cds, utr3)



