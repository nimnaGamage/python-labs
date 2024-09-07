'''
Lab 1-2
Calculating the nucleotide base counts of a given sequence
Input: DNA sequence in FASTA format
Output: The nucleotide base counts of a given DNA sequence
Nimna Gamage
'''

# Step 1: Assign four variables to zero to get the count of each base

countBaseA = 0
countBaseT = 0
countBaseG = 0
countBaseC = 0

# Step 2: Remove the identifier line (heading) and other escaping characters in the fasta file
# Step 3: Increment the value of the relevant counter variable by one for each particular base in the given DNA sequence

# read the fasta sequence
with open("BRCA1_sequence.fasta", 'r') as file:
    for line in file:
        if line != '\n':
            # removing unwanted characters
            line = line.strip()
            # removing the header
            if '>' not in line:
                for base in line:
                    if base == 'A':
                        countBaseA += 1
                    if base == 'T':
                        countBaseT += 1
                    if base == 'G':
                        countBaseG += 1
                    if base == 'C':
                        countBaseC += 1

# Step 4: Store the total base count in a variable and return the values of the variables

total = countBaseA + countBaseT + countBaseG + countBaseC

print(f"The nucleotide base counts of this DNA sequence;\nNo. of Adenine bases: {countBaseA} bp\n"
      f"No. of Thymine bases: {countBaseT} bp\nNo. of Guanine bases: {countBaseG} bp\nNo. of Cytosine bases: {countBaseC} bp\n"
      f"The total nucleotide base count of this DNA sequence: {total} bp")

## Console Output
# The nucleotide base counts of this DNA sequence;
# No. of Adenine bases: 22779 bp
# No. of Thymine bases: 23556 bp
# No. of Guanine bases: 17899 bp
# No. of Cytosine bases: 16955 bp
# The total nucleotide base count of this DNA sequence: 81189 bp