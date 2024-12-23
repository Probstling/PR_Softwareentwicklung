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

    violinplot(
        dataframe, 
        lengths, 
        title = "Transcription Length Distribution by Region", 
        y_label = "Length (bp)",
        remove_outliers_flag=True,
    )

    cg_contents = [
        "Transcript ID",
        "Total CG Content",
        "5' UTR CG Content",
        "CDS CG Content",
        "3' UTR CG Content",
    ]

    violinplot(
        dataframe,
        cg_contents,
        title="CG Content Distribution by Region",
        y_label="CG Content", 
    )
