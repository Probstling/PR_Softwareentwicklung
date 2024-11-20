from Bio import SeqIO
import re

#check total number of transcripts in the GENCODE fasta file
gencode = 'gencode.v47.pc_transcripts.fa'
sequences = SeqIO.parse(gencode, 'fasta')
count = 0
for record in sequences:
   count += 1
print("There were " + str(count) + " records in file " + gencode)

#extracting all transcript ID's from the MANE selection
mane_select = "MANE.GRCh38.v1.4.summary.txt"
with open(mane_select, 'r') as mane:
   lines = mane.readlines()

pattern = r"ENST\d+\.\d+"
all_transcripts_ID = []
for line in lines:
   transcript_ID = re.findall(pattern, line)
   all_transcripts_ID.extend(transcript_ID)

output_file = "MANE_transcript_IDs.txt"
with open(output_file, 'w') as f:
   for match in all_transcripts_ID:
      f.write(str(match) + '\n')
