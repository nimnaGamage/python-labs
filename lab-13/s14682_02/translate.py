'''
Translate the sequence in the mRNA_seq.fasta file and save the translated amino acid sequence in another FASTA file (aa_seq.FASTA)
Input: mRNA_seq.fasta file
Output: aa_seq.fasta file with the translated amino acid sequence
21/02/2023
Nimna Gamage
s14682
Lab 13-Question2_Sub-question3
'''

# Import Biopython sub-modules
from Bio.Seq import Seq

# Read the input mRNA sequence
with open("mRNA_seq.fasta", 'r') as file_mRNA:
    # for each line in file
    for line in file_mRNA:
        # if the line is not empty
        if line != '\n':
            # Removing unwanted characters
            line = line.strip()
            # header
            if '>' in line:
                # concatenate
                header = line.replace(" transcribed", " translated")
            # sequence
            else:
                sequence = Seq(line)
                # transcribe the sequence
                aa_seq = sequence.translate(to_stop=True)


# open the fasta file
with open("aa_seq.fasta", 'w') as aa_file:
    # write the header and the sequence of the fasta file
    aa_file.write("{}\n{}".format(header, aa_seq))



