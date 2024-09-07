'''
Using control structures to differentiate between DNA and amino acid sequences
Input: Two FASTA records (the amino acid sequence and the mRNA sequence) in OSDREB1A.fasta file
Output: A new FASTA file (“OSDREB1A_mRNA.fasta”) with the transcribed mRNA sequence of DREB1A protein
Nimna Gamage
'''

# Step 1: Create an empty dictionary(dict_fasta) for the sequence information

dict_fasta = {}
#header = ""
#sequence = ""

# Step 2: Open the fasta file with the amino acid sequence and mRNA sequence
# Step 3: Identify the headers of the fasta sequences using the '>' sign and make them as keys of the dict_fasta
# Step 4: Identify the sequences and make the sequences as the corresponding value of the key identifier

# read the input fasta sequence
with open("OSDREB1A.fasta", 'r') as file:
    for line in file:
        if line != '\n':
            # Removing unwanted characters
            line = line.strip()
            #header
            if '>' in line:
                header = line
                sequence = ""
            #sequence
            else:
                # concatenate different lines
                sequence += line
                dict_fasta[header] = sequence

#Printing the dictionary
print("Dictionary : ", dict_fasta)

# Step 5:Create a list of amino acids without 'A', 'C', 'G', 'T' amino acids

list_amino_acids = ['K', 'N', 'R', 'I', 'Q', 'M', 'H', 'P', 'L', 'E', 'D', 'V', 'Y', 'S', 'W', 'F']

# Step 6:Differentiate the two fasta sequences in the dictionary(values of the dictionary) using the list of amino acids
# If the elements of the above list_amino_acids are not be seen in the dictionary values, that dictionary value/ sequence is identified as a mRNA sequence.
# Unless it is an amino acid sequence.
# Step 7:After differentiating the mRNA record from the amino acid sequence, replace its Thymine (“T”) bases with Uracil (“U”) bases and store it in a variable

for head in dict_fasta:
    for amino_acid in list_amino_acids:
        if amino_acid not in dict_fasta[head]:
            # replace 'T' with 'U'
            output_seq = dict_fasta[head].replace("T", "U")

#print(output_seq)

# Step 8: Create and open the new OSDREB1A_mRNA.fasta file
# Step 9: Write the header with the word “transcribed” and the output sequence (the sequence which 'T' bases replaced with 'U' bases) in the new OSDREB1A_mRNA.fasta file

with open("OSDREB1A_mRNA.fasta", "w") as file_w:
    file_w.write(head + " transcribed" + "\n")
    file_w.write(output_seq)

#printing the content of the file
print("The OSDREB1A_mRNA.fasta file contains: \n", head, " transcribed \n")
print(output_seq)

