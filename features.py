"""
Extract sequence-specific properties.

Create a cvs file as output.
"""

from Bio import SeqIO
import pandas as pd
import numpy as np


def create_dataframe(fasta_file):
    """
    Create a pandas Dataframe from the reference transcript IDs.

    Args:
        fasta_file (.fa): holds all reference Transcript IDs
    Returns:
        dataframe: Transcript ID as identifier
    """
    transcript_ids = []

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcript_id = record.id.split("|")[0]
            transcript_ids.append(transcript_id)

    dataframe = pd.DataFrame(transcript_ids, columns=["Transcript ID"])
    return dataframe


# Extract sequence specific properties:
def length(fasta_file):
    """
    Extract the length of four regions.

    Args:
        fasta_file (.fa): holds all reference Transcript IDs
    Returns:
        dictionary: Holds the lengths of all regions coupled to Transcript ID
    """
    properties = {}

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcript_id = record.id.split("|")[0]

            # Length of entire transcript
            total = int(record.id.split("|")[6])

            # Length of 5'UTR region
            utr5 = "N/A"
            if "UTR5:" in record.id:
                utr5_header = record.id.split("UTR5:")[1].split("|")[0]
                start, end = map(int, utr5_header.split("-"))
                utr5 = end - start + 1

            # Length of CDS
            if "CDS:" in record.id:
                cds_header = record.id.split("CDS:")[1].split("|")[0]
                start, end = map(int, cds_header.split("-"))
                cds = end - start + 1

            # Length of 3' UTR region
            utr3 = "N/A"
            if "UTR3:" in record.id:
                utr3_header = record.id.split("UTR3:")[1].split("|")[0]
                start, end = map(int, utr3_header.split("-"))
                utr3 = end - start + 1

            properties[transcript_id] = {
                "Total Length": total,
                "5' UTR Length": utr5,
                "CDS Length": cds,
                "3' UTR Length": utr3,
            }

    return properties


def map_features(dataframe, properties):
    """
    Map extracted properties to the DataFrame according to Transcript ID.

    Args:
        dataframe (pandas Dataframe): holds Transcript IDs as identifieres
        properties (dictionary): holds properties coupled to transcript ID
    Returns:
        dataframe: properties coupled to correct Transcript IDs
    """
    dataframe["Total Length"] = dataframe["Transcript ID"].map(
        lambda x: properties[x]["Total Length"] if x in properties else np.nan
    )
    dataframe["5' UTR Length"] = dataframe["Transcript ID"].map(
        lambda x: properties[x]["5' UTR Length"] if x in properties else np.nan
    )
    dataframe["CDS Length"] = dataframe["Transcript ID"].map(
        lambda x: properties[x]["CDS Length"] if x in properties else np.nan
    )
    dataframe["3' UTR Length"] = dataframe["Transcript ID"].map(
        lambda x: properties[x]["3' UTR Length"] if x in properties else np.nan
    )

    return dataframe


def export(dataframe, output):
    """
    Export Dataframe to CSV file.

    Args:
        dataframe (pandas Dataframe): holds properties coupled to Transcript ID
        output (csv file): lists properties to Transcript ID in a csv file
    """
    dataframe.to_csv(output, index=False)
    print(f"DataFrame exportet to {output}")


# Call funcitons
if __name__ == "__main__":
    """
    Main script execution.
    Extracts sequence specific features,
    stores it in a pandas Dataframe
    and creates a cvs file as output for further plotting.
    """
    input = "ref_transcript_IDs.fa"
    dataframe = create_dataframe(input)
    properties = length(input)
    dataframe = map_features(dataframe, properties)
    export(dataframe, "output.csv")
