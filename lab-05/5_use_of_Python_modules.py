'''
lab 5
Using Biopython and re modules
Input: Load the downloaded FASTA file as a sequence record object
Output: print the attributes(sequence ID, description, sequence, and sequence length) of the record
Nimna Gamage
'''

# Step 1: import Biopython
from Bio import SeqIO

'''without using a for loop'''

# Step 2: Read the fasta file and store it in a variable

seq_record = SeqIO.read("ATdreb2a.fasta", "fasta")

# Step 3: Extract the relevant information from the fasta file

print("Sequence ID: ", seq_record.id)
print("Description: ", seq_record.description)
print("Sequence: ", seq_record.seq)
print("Sequence in a shorter form: ", repr(seq_record.seq))
print("Sequence Length: ", len(seq_record.seq), "bp")


'''using a for loop'''

#parsing FASTA file to extract the information from the sequence
for seq_record in SeqIO.parse("ATdreb2a.fasta", "fasta"):
    print("Sequence ID: ", seq_record.id)
    print("Description: ", seq_record.description)
    print("Sequence: ", seq_record.seq)
    print("Sequence in a shorter form: ", repr(seq_record.seq))
    print("Sequence Length: ", len(seq_record.seq), "bp")



'''
Using Biopython and re modules
Input: “ATdreb2a.fasta” FASTA file
Output: The blast output of the web-based nucleotide blast program on the ATDREB2A sequence
Nimna Gamage
'''

# Step 1: import Biopython

from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Step 2: Read the fasta file and store it in a variable

record = SeqIO.read("ATdreb2a.fasta", format="fasta")

# Step 3: Run BLAST over the Internet and store the output in a variable

result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

# Step 4: Open the .xml file and write the blast output

with open("dreb2a_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

# Step 5: Close the .xml file

result_handle.close()

### Currently qblast only works with blastn, blastp, blastx, tblast and tblastx ###


'''
Using Biopython and re modules
Input: “dreb2a_blast.xml” XML file
Output: Print the attributes of each blast hit, which are below the threshold
        (blast hit title, alignment length, E-value, score, hit/subject sequence, and hit sequence length)
Nimna Gamage
'''

# Step 1: import Biopython

from Bio.Blast import NCBIXML

# Step 2: Open and read the xml file and store it in a variable

result_handle = open("dreb2a_blast.xml")
blast_records = NCBIXML.read(result_handle)

# Step 3: Define an E-value threshold of 0.05

E_VALUE_THRESH = 0.05

# Step 4: Select the hits that are below the threshold

for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print("Blast Hit Title: ", alignment.title)
            print("Alignment Length: ", alignment.length)
            print("E-value: ", hsp.expect)
            print("Score: ", hsp.score)
            print("Hit/ Subject Sequence: ", hsp.sbjct)
            print("Hit Sequence Length: ", len(hsp.sbjct), "bp")
            print("*****************\n")






'''
Using Biopython and re modules
Input: “dreb2a_blast.xml” XML file
Output: check each blast hit sequence for the presence of the ABRE element
Nimna Gamage
'''

# Step 1: import Biopython and re modules

from Bio.Blast import NCBIXML
import re

# Step 2: Open and read the xml file and store it in a variable

result_handle = open("dreb2a_blast.xml")
blast_records = NCBIXML.read(result_handle)

# Step 3: Define an E-value threshold of 0.05

E_VALUE_THRESH = 0.05

# Step 4: Assign a variable to zero to count the number of blast hits with ABRE element present in the sequence

abre_hit_count = 0

# Step 4: Select the blast hits with ABRE element present in the sequence

print("The blast hits with ABRE element present in the sequence;\n")

# parsing through the blast hits
for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        # the hits which have an e-value less than 0.05 is selected
        if hsp.expect < E_VALUE_THRESH:
            # blast hit title
            print("\tBlast Hit Title(sequence): ", alignment.title)
            # blast hit length
            print("\tAlignment Length: ", alignment.length)
            # blast hit E-value
            print("\tE-value: ", hsp.expect)
            # blast hit score
            print("\tScore: ", hsp.score)
            # subject sequence
            print("\tHit/ Subject Sequence: ", hsp.sbjct)
            # length of the subject sequence
            print("\tHit/ Subject Sequence length: ", len(hsp.sbjct), "\n")
            sbjt_seq = hsp.sbjct
            # re compile the element that needs to be found
            pattern = re.compile("[TC]ACGT[GT]C")
            # find the element in each blast hit
            mos = re.finditer(pattern, sbjt_seq)
            for hit in mos:
                print("\t\tThe ABRE Element is present in this sequence!")
                # the element
                print("\t\tThe ABRE Element: ", hit.group())
                # location of the element
                print("\t\tThe Location of the Element in the Sequence: ", hit.span(), "\n")
                # if the element is found, count is incremented by 1
                abre_hit_count += 1

# Step 5: Return the value of the counter

print("\nThe number of blast hits with ABRE element present in the sequence: ", abre_hit_count, "hits")

