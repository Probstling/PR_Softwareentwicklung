"""
This script is generating violin plots from transcript data.

It includes:
- A function to create a violin plot for transcript lengths.
- A function to remove outliers based on the interquartile range.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def violinplot(dataframe, columns_to_plot, title, y_label, remove_outliers_flag=False):
    """
    Create a violin plot for specified properties for each region.

    Args:
        dataframe (pandas.DataFrame): contains all properties
        columns_to_plot (list of str): List of column names to plot.
        title (str): Title for the violin plot.
        y_label (str): Y-axis label
        remove_outliers_flag (bool): Wether to remove outliers or not. 

    Returns:
        None: The function displays the plot but does not return any value.
    """
    dataframe = dataframe[columns_to_plot]
    dataframe.replace("N/A", pd.NA, inplace=True)

    # Convert columns to numeric, ignoring errors (e.g., N/A)
    for column in dataframe.columns[1:]:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")

    # Melt the DataFrame into long format for plotting
    melted_df = dataframe.melt(
        id_vars=["Transcript ID"], var_name="Region", value_name=y_label
    )

    # Optionally remove outliers
    if remove_outliers_flag:
        melted_df = melted_df.groupby(
            "Region", group_keys=False
        ).apply(remove_outliers)

    # Plot with seaborn
    plt.figure(figsize=(10, 6))
    sns.violinplot(
        data=melted_df,
        x="Region",
        y=y_label,
        palette="muted",
        order=columns_to_plot[1:],
    )
    plt.title(title)
    plt.xlabel("Region")
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.show()


def remove_outliers(group):
    """
    Remove outliers from a DataFrame group based on the interquartile range.

    Args:
        group (pandas.DataFrame): Data grouped by a column, such as 'Region'.

    Returns:
        pandas.DataFrame: The filtered group with outliers removed.
    """
    Q1 = group.iloc[:, -1].quantile(0.25)
    Q3 = group.iloc[:, -1].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return group[
        (group.iloc[:, -1] >= lower_bound)
        & (group.iloc[:, -1] <= upper_bound)
    ]  


def barplot_codon(dataframe, title, y_label):

    codon_groups = {
        'Glycine': ["GGT", "GGC", "GGA", "GGG"],
        'Alanine': ["GCT", "GCC", "GCA", "GCG"],
        'Valine': ["GTT", "GTC", "GTA", "GTG"],
        'Cysteine': ["TGT", "TGC"],
        'Proline': ["CCT", "CCC", "CCA", "CCG"],
        'Leucine': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        'Isoleucine': ["ATT", "ATC", "ATA"],
        'Methionine': ["ATG"],
        'Tryptophan': ["TGG"],
        'Phenylalanine': ["TTT", "TTC"],
        'Serine': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        'Threonine': ["ACT", "ACC", "ACA", "ACG"],
        'Tyrosine': ["TAT", "TAC"],
        'Asparagine': ["AAT", "AAC"],
        'Glutamine': ["CAA", "CAG"],
        'Lysine': ["AAA", "AAG"],
        'Arginine': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        'Histidine': ["CAT", "CAC"],
        'Aspartic Acid': ["GAT", "GAC"],
        'Glutamic Acid': ["GAA", "GAG"],
        'Stop': ["TAA", "TAG", "TGA"]
    }

    # Collect codons from groups
    codons = [codon for group in codon_groups.values() for codon in group]

    # Select only the relevant columns (Transcript ID and codons)
    dataframe = dataframe[["Transcript ID"] + codons]

    # Calculate the mean and error SEM
    codon_means = dataframe[codons].mean()
    codon_errors = dataframe[codons].sem()

    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        "Codon": codon_means.index,
        y_label: codon_means.values,
        "Error": codon_errors.values
    })

    # Assign amino acid group to each codon
    plot_df['Amino Acid Group'] = plot_df['Codon'].apply(
        lambda codon: next(group for group, codons in codon_groups.items() if codon in codons)
    )

    # Ensure Codon column is categorical with the specified order
    plot_df["Codon"] = pd.Categorical(plot_df["Codon"], categories=codons, ordered=True)

    # Plot using Seaborn with grouped colors
    plt.figure(figsize=(16, 10))
    ax = sns.barplot(
        data=plot_df,
        x="Codon",
        y=y_label,
        hue="Amino Acid Group",
        ci=None,
        palette="Set2",
        dodge=False  # Single group columns
    )

    # Add error bars
    for i, row in plot_df.iterrows():
        ax.errorbar(
            x=i, y=row[y_label], yerr=row["Error"], fmt='none', c='black', capsize=5
        )

    # Adjust legend and layout
    plt.legend(
        title="Amino Acid Group",
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.
    )

    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel(y_label)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.subplots_adjust(right=0.8, bottom=0.1)  # Add space for the legend
    plt.show()


def barplot_aa(dataframe, title, y_label):

    # Convert one-letter amino acid codes to three-leter codes
    aa_mapping = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln',
    'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'
    }

    aa_groups = {
        'Nonpolar': ["Gly", "Ala", "Val", "Cys", "Pro", "Leu", "Ile", "Met", "Trp", "Phe"],
        'Polar': ["Ser", "Thr", "Tyr", "Asn", "Gln"],
        'Charged Positive': ["Lys", "Arg", "His"],
        'Charged Negative': ["Asp", "Glu"]
    }

    # Map columns to three-letter amino acid codes
    dataframe.rename(columns=aa_mapping, inplace=True)

    # Collect amino acids from groups
    amino_acids = [aa for group in aa_groups.values() for aa in group]

    # Select only the relevant columns (Transcript ID and amino acids)
    dataframe = dataframe[["Transcript ID"] + amino_acids]

    # Calculate the mean and error SEM
    aa_means = dataframe[amino_acids].mean()
    aa_errors = dataframe[amino_acids].sem()

    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        "Amino Acid": aa_means.index,
        y_label: aa_means.values,
        "Error": aa_errors.values
    })

    # Assign amino acid group to each amino acid
    plot_df['Amino Acid Group'] = plot_df['Amino Acid'].apply(
        lambda aa: next(group for group, aas in aa_groups.items() if aa in aas)
    )

    # Ensure Amino Acid column is categorical with the specified order
    plot_df["Amino Acid"] = pd.Categorical(plot_df["Amino Acid"], categories=amino_acids, ordered=True)

    # Plot using Seaborn with grouped colors
    plt.figure(figsize=(12, 8))
    ax = sns.barplot(
        data=plot_df,
        x="Amino Acid",
        y=y_label,
        hue="Amino Acid Group",
        ci=None,
        palette="Set2",
        dodge=False  # Single group columns
    )

    # Add error bars
    for i, row in plot_df.iterrows():
        ax.errorbar(
            x=i, y=row[y_label], yerr=row["Error"], fmt='none', c='black', capsize=5
        )

    # Adjust legend and layout
    plt.legend(
        title="Amino Acid Group",
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.
    )

    plt.title(title)
    plt.xlabel("Amino Acid")
    plt.ylabel(y_label)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.8, bottom=0.1) 
    plt.show()

if __name__ == "__main__":
    """
    Main script execution.
    Reads CSV file as input and passes on
    certain properties to create a plot.
    """
    #file_path = "output.csv"
    #dataframe = pd.read_csv(file_path)

    #codon_profile = "codon_profile.csv"
    #codon_profile_df = pd.read_csv(codon_profile)

    aa_profile = "aa_profile.csv"
    aa_profile_df = pd.read_csv(aa_profile)

    #lengths = [
    #    "Transcript ID",
    #    "Total Length",
    #    "5' UTR Length",
    #    "CDS Length",
    #    "3' UTR Length",
    #]

    #violinplot(
    #    dataframe, 
    #    lengths, 
    #    title = "Transcription Length Distribution by Region", 
    #    y_label = "Length (bp)",
    #    remove_outliers_flag=True,
    #)

    #cg_contents = [
    #    "Transcript ID",
    #    "Total CG Content",
    #    "5' UTR CG Content",
    #    "CDS CG Content",
    #    "3' UTR CG Content",
    #]

    #violinplot(
    #    dataframe,
    #    cg_contents,
    #    title="CG Content Distribution by Region",
    #    y_label="CG Content", 
    #)

    #barplot_codon(
    #    codon_profile_df,
    #    title="Codon Frequency Distribution by Amino Acid",
    #    y_label="Frequency (normalized)",
    #)

    barplot_aa(
        aa_profile_df,
        title = "Amino Acid Frequency Profile",
        y_label = "Frequency (normalized)",
    )