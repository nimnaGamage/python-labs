'''
Transcribe the sequence in the cds_seq FASTA file and save the transcribed mRNA sequence in another FASTA file (mRNA_seq.fasta)
Input: cds_seq.fasta file
Output: mRNA_seq.fasta file with the transcribed mRNA sequence
21/02/2023
Nimna Gamage
s14682
Lab 13-Question2_Sub-question2
'''

# Import Biopython sub-modules
from Bio.Seq import Seq

# Read the input DNA sequence
with open("cds_seq.fasta", 'r') as file:
    # for each line in file
    for line in file:
        # if the line is not empty
        if line != '\n':
            # Removing unwanted characters
            line = line.strip()
            # header
            if '>' in line:
                # concatenate
                header = line + " transcribed"
            # sequence
            else:
                sequence = Seq(line)
                # transcribe the sequence
                messenger_rna = sequence.transcribe()


# open the fasta file
with open("mRNA_seq.fasta", 'w') as mRNA_file:
    # write the header and the sequence of the fasta file
    mRNA_file.write("{}\n{}".format(header, messenger_rna))



