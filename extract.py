from Bio import SeqIO
import re

#check total number of transcripts in the GENCODE fasta file
def counting(fasta_file):
   sequences = SeqIO.parse(fasta_file, 'fasta')
   count = 0
   for record in sequences:
      count += 1
   print(str(count) + " records are in file " + fasta_file)

#extract all transcript entries from gencode which are listet in the MANE selection
def MANE_select_transcript_IDs():

   #extracting all transcript ID's from the MANE selection
   mane_select = "MANE.GRCh38.v1.4.summary.txt"
   with open(mane_select, 'r') as mane:
      lines = mane.readlines()

   #creating regular expression pattern to extract only transcript IDs into the array
   pattern = r"ENST\d+\.\d+"
   mane_ids = []
   for line in lines:
      transcript_ID = re.findall(pattern, line)
      mane_ids.extend(transcript_ID)

   print(f"{len(mane_ids) } reference transcripts are listed in the MANE seleciton.")

   #Parsing through gencode fasta file
   gencode = 'gencode.v47.pc_transcripts.fa'
   sequences = SeqIO.parse(gencode, 'fasta')
   output_file = 'ref_transcript_IDs.fa'

   #writing selection into a new fasta file
   with open (output_file, 'w') as f:
      for record in sequences:
         #strip after | and keep first entry of the id (ENST...)
         transcript_ID = record.id.split('|')[0] 
         if transcript_ID in mane_ids: 
            #write matching record to new output fasta file
            SeqIO.write(record, f, 'fasta')

if __name__ == "__main__":
   counting('gencode.v47.pc_transcripts.fa')
   MANE_select_transcript_IDs()
   counting('ref_transcript_IDs.fa')