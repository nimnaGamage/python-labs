'''
Translating mRNA sequences using a Python script
Input: Fasta file with the mRNA sequence (with Uracil bases instead of Thymine)
       The text file of the amino acid codons table
Output: A new fasta file with the translated mRNA sequences/ amino acid sequence
Nimna Gamage
'''

# Step 1: Create an empty dictionary (dict_codon) to store codons to amino acid mappings.
# Step 2: Open the codon_table text file and identify the three letter codons as the key of the dictionary.
# Step 3: identify the amino acid letter as the value of the corresponding key codon of the dictionary.

#Create the empty dictionary
dict_codon = {}

#Read the codon table text file
with open("codon_table.txt", "r") as file_codon:
    for line in file_codon.readlines():
        if "#" not in line and line != "\n":
            newline = line.strip().split("\t")

            codon = newline[0]
            aaletter = newline[2]
            dict_codon[codon] = aaletter

#Printing the codon dictionary
print("The dictionary created using codon table : \n", dict_codon)

# Step 4: Create a new empty string to hold the sequence of the mRNA fasta file
# Step 5: Open the mRNA sequence and extract the sequence without the header to the newly created string('seq')

seq = ""

# read the mRNA sequence
with open("OSDREB1A_mRNA.fasta", 'r') as filename_mRNA:
    for line in filename_mRNA:
        if line != '\n':
            # remove unwanted characters
            line = line.strip()
            # remove the header
            if '>' not in line:
                seq += line

# Step 6: Create a new empty string to hold the newly translated amino acid sequence
# Step 7: Translate the three letter codons into amino acids using the dict_codon, while stop codon is met
# Step 8: Assign the new amino acid sequence to the newly created string

#empty string for amino acids
aminoacid_seq = ""

#while loop
i = 0
while i < len(seq):
    triple = seq[i:i+3]
    i += 3
    # break the while loop when hit a stop codon
    if triple in ["UAA", "UAG", "UGA"]:
        break
    aminoacid_seq += (dict_codon[triple])

# Step 9: Create and Write the translated sequence in a new output file in FASTA format

with open("translated_seq.fasta", "w") as file_protein:
    file_protein.write(aminoacid_seq)

#printing the amino acid sequence and the length of the amino acid sequence
print("The translated Amino Acid sequence is : ", aminoacid_seq)
print("The length of the translated Amino Acid sequence is : ", len(aminoacid_seq), "Amino acids.")

