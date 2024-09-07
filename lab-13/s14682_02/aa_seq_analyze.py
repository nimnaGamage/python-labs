'''
Analyze the aa_seq.fasta file and calculate the length, molecular weight, alanine percentage, and glycine percentage of the sequence
Input: aa_seq.fasta file
Output: aa_stats.txt including the length, molecular weight, alanine percentage, and glycine percentage of the amino acid sequence
21/02/2023
Nimna Gamage
s14682
Lab 13-Question2_Sub-question4
'''

# Import Biopython sub-modules
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Read the input mRNA sequence
with open("aa_seq.fasta", 'r') as file_aa:
    # for each line in file
    for line in file_aa:
        # if the line is not empty
        if line != '\n':
            # Removing unwanted characters
            line = line.strip()
            # removing the header
            if '>' not in line:
                prot_analysis = ProteinAnalysis(line)


# Assign statistics to a variable
pro_length = len(line)
molarW = prot_analysis.molecular_weight()
alanine_p = prot_analysis.get_amino_acids_percent()['A']*100
glycine_p = prot_analysis.get_amino_acids_percent()['G']*100

# open the text file
with open("aa_stats.txt", 'w') as stat_file:
    # write the statistics of the protein sequence
    stat_file.write("Statistical analysis of the protein sequence ; \n")
    # length of the amino acid sequence
    stat_file.write("  The length of the sequence : {} aa \n".format(pro_length))
    # molecular weigth of the amino acid sequence
    stat_file.write("  The molecular weight of the sequence : %.2f\n" % molarW)
    # The alanine percentage of the amino acid sequence
    stat_file.write("  The alanine percentage of the sequence : %.2f" % alanine_p + "%\n")
    # The glycine percentage of the amino acid sequence
    stat_file.write("  The glycine percentage of the sequence : %.2f" % glycine_p + "%\n")





