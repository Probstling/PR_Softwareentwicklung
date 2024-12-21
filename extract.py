"""
Filter Gencode Transcript IDs from MANE select.

Creates a new fasta file as output.
"""

from Bio import SeqIO
import re


def counting(fasta_file):
    """
    Check total number of transcripts in fasta files.

    Args:
        fasta_file (.fa): all transcripts before and after filtering
    """
    sequences = SeqIO.parse(fasta_file, "fasta")
    count = 0
    for record in sequences:
        count += 1
    print(str(count) + " records are in file " + fasta_file)


def MANE_select_transcript_IDs():
    """
    Extract transcript entries from gencode matching MANE select.

    All reference transcript IDs from MANE select are listed in an array.
    Entries from Gencode fasta file are listed in a new fasta file
    if they are listed in the MANE selection.

    """
    # extracting all transcript ID's from the MANE selection
    mane_select = "MANE.GRCh38.v1.4.summary.txt"
    with open(mane_select, "r") as mane:
        lines = mane.readlines()

    # creating regular expression pattern for extracting transcript ID
    pattern = r"ENST\d+\.\d+"
    mane_ids = []
    for line in lines:
        transcript_ID = re.findall(pattern, line)
        mane_ids.extend(transcript_ID)

    print(
        f"{len(mane_ids) } IDs listed in MANE select."
      )

    # Parsing through gencode fasta file
    gencode = "gencode.v47.pc_transcripts.fa"
    sequences = SeqIO.parse(gencode, "fasta")
    output_file = "ref_transcript_IDs.fa"

    # writing selection into a new fasta file
    with open(output_file, "w") as f:
        for record in sequences:
            transcript_ID = record.id.split("|")[0]
            if transcript_ID in mane_ids:
                SeqIO.write(record, f, "fasta")


if __name__ == "__main__":
    """
    Main script execution.
    Counts before and after filtering.
    """
    counting("gencode.v47.pc_transcripts.fa")
    MANE_select_transcript_IDs()
    counting("ref_transcript_IDs.fa")
