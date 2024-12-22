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


def length(fasta_file):
    """
    Extract the length of four regions.

    Args:
        fasta_file (.fa): holds all reference Transcript IDs
    Returns:
        dictionary: Holds the lengths of all regions coupled to Transcript ID
    """
    length_properties = {}

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

            length_properties[transcript_id] = {
                "Total Length": total,
                "5' UTR Length": utr5,
                "CDS Length": cds,
                "3' UTR Length": utr3,
            }

    return length_properties


def cg_content(fasta_file):

    cg_properties = {}

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            transcript_id = record.id.split("|")[0]
            sequence = str(record.seq)

            # Function to calculate CG content for each region
            def calculate_cg_content(seq_slice):
                if seq_slice:
                    cg_count = seq_slice.count("C") + seq_slice.count("G")
                    return round(cg_count / len(seq_slice), 4)
                return "N/A"
            
            # Entire transcript
            total_content = calculate_cg_content(sequence)            

            # 5' UTR
            if "UTR5:" in record.id:
                utr5_header = record.id.split("UTR5:")[1].split("|")[0]
                start, end = map(int, utr5_header.split("-"))
                utr5_sequence = sequence[start -1:end]
                utr5_content = calculate_cg_content(utr5_sequence)
            else:
                utr5_content = "N/A"

            # CDS
            cds_header = record.id.split("CDS:")[1]. split("|")[0]
            start, end = map(int, cds_header.split("-"))
            cds_sequence = sequence[start -1:end]
            cds_content = calculate_cg_content(cds_sequence)

            # 3' UTR 
            if "UTR3:" in record.id:
                utr3_header = record.id.split("UTR3:")[1].split("|")[0]
                start, end = map(int, utr3_header.split("-"))
                utr3_sequence = sequence[start -1:end]
                utr3_content = calculate_cg_content(utr3_sequence)
            else:
                utr3_content = "N/A"

            # Add results to dictionary
            cg_properties[transcript_id] = {
                "Total CG Content": total_content,
                "5' UTR CG Content": utr5_content,
                "CDS CG Content": cds_content,
                "3' UTR CG Content": utr3_content,
            }

    return cg_properties


def map_features(dataframe, length_properties, cg_properties):
    """
    Map extracted properties to the DataFrame according to Transcript ID.

    Args:
        dataframe (pandas Dataframe): holds Transcript IDs as identifieres
        length_properties (dictionary): holds length coupled to transcript ID
        cg_properties (dictionary): holds CG content coupled to transcript ID
    Returns:
        dataframe: properties coupled to correct Transcript IDs
    """
    # Mapping of Length
    dataframe["Total Length"] = dataframe["Transcript ID"].map(
        lambda x: length_properties[x]["Total Length"] if x in length_properties else np.nan
    )
    dataframe["5' UTR Length"] = dataframe["Transcript ID"].map(
        lambda x: length_properties[x]["5' UTR Length"] if x in length_properties else np.nan
    )
    dataframe["CDS Length"] = dataframe["Transcript ID"].map(
        lambda x: length_properties[x]["CDS Length"] if x in length_properties else np.nan
    )
    dataframe["3' UTR Length"] = dataframe["Transcript ID"].map(
        lambda x: length_properties[x]["3' UTR Length"] if x in length_properties else np.nan
    )

    # Mapping of CG content
    dataframe["Total CG Content"] = dataframe["Transcript ID"].map(
        lambda x: cg_properties[x]["Total CG Content"] if x in cg_properties else np.nan
    )
    dataframe["5' UTR CG Content"] = dataframe["Transcript ID"].map(
        lambda x: cg_properties[x]["5' UTR CG Content"] if x in cg_properties else np.nan
    )
    dataframe["CDS CG Content"] = dataframe["Transcript ID"].map(
        lambda x: cg_properties[x]["CDS CG Content"] if x in cg_properties else np.nan
    )
    dataframe["3' UTR CG Content"] = dataframe["Transcript ID"].map(
        lambda x: cg_properties[x]["3' UTR CG Content"] if x in cg_properties else np.nan
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


if __name__ == "__main__":
    """
    Main script execution.
    Extracts sequence specific features,
    stores it in a pandas Dataframe
    and creates a cvs file as output for further plotting.
    """
    input = "ref_transcript_IDs.fa"
    dataframe = create_dataframe(input)

    length_properties = length(input)
    cg_properties = cg_content(input)

    dataframe = map_features(dataframe, length_properties, cg_properties)

    export(dataframe, "output.csv")
