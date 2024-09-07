'''
Lab 1-1-1
Calculating the length of a given DNA sequence
Input: DNA sequence in FASTA format
Output: The length of a given DNA sequence
Nimna Gamage
'''

# Step 1: Assign a variable to zero to count the number of bases

baseCount = 0

# Step 2: Remove the identifier line (heading) and other escaping characters in the fasta file
# Step 3: Increment the counter by one for each base in the given DNA sequence

# read the input fasta sequence
with open("BRCA1_sequence.fasta", 'r') as file:
    # check each line in the file
    for line in file:
        # checking if the line is not empty
        if line != '\n':
            # removing unwanted characters
            line = line.strip()
            # removing the header
            if '>' not in line:
                for base in line:
                    # increment the baseCount
                    baseCount += 1

# Step 4: Return the value of the counter

print(f"This DNA sequence is {baseCount} base pairs long.")


## Console Output
# This DNA sequence is 81189 base pairs long.

'''
Lab 1-1-2
Calculating the length of a given DNA sequence using len() 
Input: DNA sequence in FASTA format
Output: The length of a given DNA sequence
Nimna Gamage
'''

# Step 1: Assign a variable to zero to count number of bases

base_count = 0

# Step 2: Remove the identifier line (heading) and other escaping characters in the fasta file
# Step 3: Increment the counter by adding length of each line

with open("BRCA1_sequence.fasta", 'r') as file:
    for line in file:
        if line != '\n':
            # removing unwanted characters
            line = line.strip()
            # removing the header
            if '>' not in line:
                base_count += len(line)

# Step 4: Return the value of the counter

print(f"The length of this DNA sequence is {base_count} base pairs.")


## Console Output
# The length of this DNA sequence is 81189 base pairs.