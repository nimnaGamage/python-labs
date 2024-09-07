'''
Retrieve the GenBank records for the accession number input by the user
Input: user input [the accession number (with version)]
Output: The GenBank record for an accession number (with version) of a protein coding DNA sequence or a reverse-transcribed mRNA complement
21/02/2023
Nimna Gamage
s14682
Lab 13-Question2_Sub-question1
'''

# Import Biopython sub-modules
from Bio import Entrez
from Bio import SeqIO

# get the input from the user
accession = input("Enter the accession number (with version) : ")

# provide E-mail
Entrez.email = "nimnagamage65@gmail.com"
# retrieve the fasta sequence from Entrez
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
# close the handle
handle.close()
# open the fasta file
file = open("cds_seq.fasta", 'w')
# write the header and the sequence of the fasta file
file.write(">{}\n{}".format(record.description, record.seq))
# close the file
file.close()

