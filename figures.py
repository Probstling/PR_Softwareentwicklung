"""
This script is generating violin plots from transcript data.

It includes:
- A function to create a violin plot for transcript lengths.
- A function to remove outliers based on the interquartile range.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def violinplot_length(dataframe, lengths):
    """
    Create a violin plot for specified transcript length regions.

    Args:
        dataframe (pandas.DataFrame): contains properties
        columns_to_plot (list of str): List of column names to plot.

    Returns:
        None: The function displays the plot but does not return any value.
    """
    dataframe = dataframe[lengths]
    dataframe.replace("N/A", pd.NA, inplace=True)

    # Convert columns to numeric, ignoring errors (e.g., N/A)
    for column in dataframe.columns[1:]:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")

    # Melt the DataFrame into long format for plotting
    melted_df = dataframe.melt(
        id_vars=["Transcript ID"], var_name="Region", value_name="Length (bp)"
    )

    melted_df = melted_df.groupby(
        "Region", group_keys=False
    ).apply(remove_outliers)
    custom_order = [
        "Total Length",
        "5' UTR Length",
        "CDS Length",
        "3' UTR Length"
    ]

    plt.figure(figsize=(10, 6))
    sns.violinplot(
        data=melted_df,
        x="Region",
        y="Length (bp)",
        palette="muted",
        order=custom_order
    )
    plt.title("Transcript Length Distribution by Region")
    plt.xlabel("Region")
    plt.ylabel("Length (bp)")
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
    Q1 = group["Length (bp)"].quantile(0.25)
    Q3 = group["Length (bp)"].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return group[
        (group["Length (bp)"] >= lower_bound)
        & (group["Length (bp)"] <= upper_bound)
    ]


if __name__ == "__main__":
    """
    Main script execution.
    Reads CSV file as input and passes on
    certain properties to create a plot.
    """
    file_path = "output.csv"
    dataframe = pd.read_csv(file_path)
    lengths = [
        "Transcript ID",
        "Total Length",
        "5' UTR Length",
        "CDS Length",
        "3' UTR Length",
    ]
    violinplot_length(dataframe, lengths)
