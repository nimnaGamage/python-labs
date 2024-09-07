'''
Retrieve the GenBank records for the "AAK43967.1", "AED90870.1", "NP_567720.1", "AAK59861.1" accessions
Input: "AAK43967.1", "AED90870.1", "NP_567720.1", "AAK59861.1" accessions
Output: The FASTA files named by their respective accession number with the fasta sequence written inside
21/02/2023
Nimna Gamage
Lab 13-Question1_Sub-question1
'''

# Import Biopython sub-modules
from Bio import Entrez
from Bio import SeqIO

# Create a list of accession numbers
acc_list = ["AAK43967.1", "AED90870.1", "NP_567720.1", "AAK59861.1"]

# loop through each accession number in the list
for acc_no in acc_list:
    # provide E-mail
    Entrez.email = "nimnagamage65@gmail.com"
    # retrieve full records from entrez
    handle = Entrez.efetch(db="protein", id=acc_no, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    # close the handle
    handle.close()
    # open the fasta file named with the respective accession number
    file = open("{}.fasta".format(acc_no), 'w')
    # write the respective fasta sequence in each file
    file.write(">{}\n{}".format(record.description, record.seq))
    # close the file
    file.close()


