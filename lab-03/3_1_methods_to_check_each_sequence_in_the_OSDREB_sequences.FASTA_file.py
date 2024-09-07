'''
lab 3-1
Custom methods to check each sequence in the 'OSDREB_sequences.FASTA' file and print the AT content
Input: 'OSDREB_sequences.FASTA' file
Output: check each sequence in the 'OSDREB_sequences.FASTA' file and print the AT content
30/11/2022
Nimna Alupotha Gamage
Index:s14682
'''

#custom Python method to calculate the AT content of a given DNA sequence
def total_A_T_Content(inputseq):

    #Assign a variable to zero to count the number of 'A' and 'T' bases
    baseCountAT = 0

    #Assign a variable to len(inputseq) to get the total base count
    baseTotal = len(inputseq)

    #for loop to go through each base
    for base in inputseq:
        #Increment the 'baseCountAT' by one only if base A or base T encountered
        if (base == 'A' or base == "T"):
            baseCountAT += 1

    #get the AT Content as a percentage value from the total base count
    baseATContent = (baseCountAT / baseTotal) * 100

    #Return the AT Content as a percentage
    return baseATContent

#As an alternative way default function 'count' also can be used to get the AT Content

    # #count the Adenine bases
    # baseCountA = inputseq.count("A")

    # #count the Thymine bases
    # baseCountT = inputseq.count("T")
    #
    # #total A and T count
    # baseCountAT = baseCountA + baseCountT

    # #total base count
    # total = len(inputseq)
    #
    # base_AT_Content = (baseCountAT/total)*100
    # return base_AT_Content


#custom Python method to split multiple FASTA sequences in a single text file and return a dictionary containing sequence headers as keys
#and the sequences as values.
def fastaDictionary(fasta_file):

    #Assign an empty dictionary
    dictFasta = {}

    #Assign two empty strings
    header = ""
    sequence = ""

    #Open the fasta file
    with open(fasta_file, 'r') as file:
        #check each line in the file
        for line in file:
            #checking if the line is not empty
            if line != '\n':
                #strip the line to remove spaces
                line = line.strip()
                #header
                #if '>' is in the line, the line is assign to a header in the dictionary
                if '>' in line:
                    header = line
                    sequence = ""
                #sequence
                #else the lines are added to the sequence as values in the dictionary
                else:
                    sequence += line
                    dictFasta[header] = sequence


    #return the dictionary with fasta headers and sequences
    return dictFasta


#custom Python method to check whether the given sequence type is DNA, mRNA or amino acid sequence
def checkingSeqType(seq):

    #list of amino acid letters except the letter ('A', 'T', 'G', 'C')
    listAminoAcids = ['K', 'N', 'R', 'S', 'I', 'Q', 'M', 'H', 'P', 'L', 'E', 'D', 'V', 'O', 'Y', 'S', 'W', 'F']

    #For loop to go through each amino acid of the amino acid list
    for amino_acid in listAminoAcids:
        #sequence is recognized as a protein, if an item in amino acid list is encountered in the sequence
        if amino_acid in seq:
            seqType = 'Amino Acid Sequence'
            break
        #else the sequence is DNA or mRNA
        else:
            #the sequence is recognized as a mRNA, if 'U' is encountered in the sequence
            if 'U' in seq:
                seqType = 'mRNA Sequence'
                break
            #else the sequence is a DNA sequence
            elif 'U' not in seq:
                seqType = 'DNA Sequence'

    #Return the type of the sequence
    return seqType


#Use the above written methods to check each sequence in the OSDREB_sequences.FASTA file and print the AT content.

#Read the fasta file and store it in a variable
fastaFile = "OSDREB_sequences.FASTA"

#Give the above fasta file contained variable as an argument to the fastaDictionary method/ function
fastFileDict = fastaDictionary(fastaFile)

#For loop to go through each header in the dictionary
for header in fastFileDict:
    #Print each header in the dictionary
    print(header)
    #values/sequences of the dictionary are assigned to 'seq'
    seq = fastFileDict[header]
    # print seq
    print("Sequence: ", seq)
    #sequence type is checked by the checkingSeqType method
    seqType = checkingSeqType(seq)
    #print the sequence type
    print("Sequence Type: ", seqType)
    #if the sequence is a DNA sequence, the AT content is printed
    if seqType == "DNA Sequence":
        print("The AT Content of the DNA sequence(as a percentage): ", total_A_T_Content(seq), "%")

