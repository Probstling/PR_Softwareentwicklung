import pandas as pd

def calculate_statistics(dataframe, columns):
    """
    Calculate statistics for specified columns in a DataFrame.

    Args:
        dataframe (pandas.DataFrame): Input data.
        columns (list of str): List of column names to analyze.
    Returns:
        pandas.DataFrame: Summary statistics for each column.
    """    
    stats = []
    for column in columns:
        if column in dataframe.columns:
            col_data = pd.to_numeric(dataframe[column], errors='coerce').dropna()
            stats.append({
                "Property": column,
                "Mean": col_data.mean().round(2),
                "Median": col_data.median().round(2),
                "Standard Deviation": col_data.std().round(2),
                "Minimum": col_data.min().round(2),
                "Maximum": col_data.max().round(2),
                "25th Percentile": col_data.quantile(0.25).round(2),
                "75th Percentile": col_data.quantile(0.75).round(2),
                "IQR": (col_data.quantile(0.75) - col_data.quantile(0.25)).round(2),
                "Valid Count": col_data.count(),
            })
    return pd.DataFrame(stats)


def save_statistics_table(stats_df, output_file):
    """
    Save the statistics summary as a CSV or text file.

    Args:
        stats_df (pandas.DataFrame): DataFrame containing statistics.
        output_file (str): Path to save the file.
    """
    if output_file.endswith(".csv"):
        stats_df.to_csv(output_file, sep = '\t', index=False)
    else:
        with open(output_file, "w") as file:
            file.write(stats_df.to_string(index=False))


if __name__ == "__main__":
    input_file = "output.csv"
    output_file = "statistics_summary.txt"

    dataframe = pd.read_csv(input_file)

    columns_to_analyze = [
        "Total Length", "5' UTR Length", "CDS Length", "3' UTR Length",
        "Total CG Content", "5' UTR CG Content", "CDS CG Content", "3' UTR CG Content"
    ]

    statistics_df = calculate_statistics(dataframe, columns_to_analyze)

    save_statistics_table(statistics_df, output_file)
    print(f"Statistics summary saved to {output_file}.")