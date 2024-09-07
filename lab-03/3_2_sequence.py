'''
lab 3-2
Python Sequence class to store any biological sequence (DNA, mRNA, and amino acid sequences)
Input : Fasta file with mRNA/nucleotide sequences and protein sequences, codon table text file
Output : Specific Details for the particular sequences and the total number of sequence objects
Nimna Alupotha Gamage
Index: s14682
'''

'''Main Class - Sequence Class'''

class Sequence:

    #Attributes/Variables

    #Class variables
    #Sequence count (to count the number of sequences created by the sequence class)
    sequence_count = 0

    # Constructor method with input parameters
    def __init__(self, gene_name, gene_id, species_name, subspecies_name, sequence):
        # Instance variables
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.sequence = sequence
        self.sequence_type = self.get_Seq_Type()
        self.sequence_length = len(sequence)
        self.species_name = species_name
        self.subspecies_name = subspecies_name

        #When an object is created, the constructor is initiated and sequence_count is incremented by 1 for each object
        Sequence.sequence_count += 1


    #static method
    #a method to split multiple FASTA sequences in a single text file and return as a dictionary
    @staticmethod
    def fasta_Split(fasta_file):
    # Give the input fasta file as an argument to the method/ function

        # Assign an empty dictionary
        dict_fasta = {}

        # Assign an empty strings
        # key = ""
        # value = ""

        # Open the fasta file
        with open(fasta_file, 'r') as file:
            # check each line in the file
            for line in file:
                # checking if the line is not empty
                if line != '\n':
                    # strip the line to remove spaces
                    line = line.strip()
                    # fasta header
                    # if '>' is in the line
                    if '>' in line:
                        #line is stripped to remove '>' and splitted using '-'
                        line = (line.strip(">")).split("-")
                        # first item of the 'line' list is the key of the dictionary
                        key = line[0]
                        # 'line' list is assigned as the values of the dictionary
                        value = line[0:]

                    # fasta sequence
                    else:
                        # lines are appended to the 'value' list
                        value.append(line)

                        #dictionary
                        dict_fasta[key] = value
                        #seperated sequence letters are merged according to the sequence type
                        for header in dict_fasta:
                            geneDetails = dict_fasta[header]
                            #for coding/mRNA sequences
                            if "_CDS" in header:
                                geneDetails[4:] = [''.join(geneDetails[4:])]
                            #for protein sequences
                            elif "_P" in header:
                                geneDetails[6:] = [''.join(geneDetails[6:])]


        # return the dictionary
        return dict_fasta


    #instance method
    #A method to check the sequence type
    def get_Seq_Type(self):

        # list of amino acid letters except the letter ('A', 'T', 'G', 'C')
        listAminoAcids = ['K', 'N', 'R', 'I', 'Q', 'M', 'H', 'P', 'L', 'E', 'D', 'V', 'Y', 'S', 'W', 'F']

        # For loop to go through each amino acid of the amino acid list
        for amino_acid in listAminoAcids:
            # sequence is recognized as a protein, if an item in amino acid list is encountered in the sequence
            if amino_acid in self.sequence:
                seqType = 'Amino Acid Sequence'
                break
            # else the sequence is DNA or mRNA
            else:
                # the sequence is recognized as a mRNA, if 'U' is encountered in the sequence
                if 'U' in self.sequence:
                    seqType = 'mRNA Sequence'
                    break
                # else the sequence is a DNA sequence
                elif 'U' not in self.sequence:
                    seqType = 'DNA Sequence'

        # Return the type of the sequence
        return seqType

    #instance method
    #a method to return a dictionary of character counts with each character as the key and count as the value
    def get_Character_Count(self):

        #assign an empty dictionary
        dict_character = {}

        #for loop to go through each character in a sequence
        for character in self.sequence:
            #each character is the key of the dictionary
            characterBase = character
            #count of each character is the value of dictionary
            baseCount = self.sequence.count(character)
            dict_character[characterBase] = baseCount

        #return the dictionary
        return dict_character


#main method
if __name__ == "__main__":

    '''Creating Object sequences from a fasta file-lab 4'''

    # Read the file and store it in a variable
    file = "OSDREB_sequences.FASTA"

    # Give the above fastafile contained variable as an argument to the fasta_Split method/ function
    fastaSeqDict = Sequence.fasta_Split(file)

    # #printing the fasta sequence dictionary
    # print(fastaSeqDict)

    # for loop to go through each value of the dictionary
    for seqDetailList in fastaSeqDict.values():
        # #Printing the sequence detailed list
        # print(seqDetailList)

        # Create objects for each sequence detailed list
        seqObjcts = Sequence(*[seqDetailList[0], seqDetailList[1], seqDetailList[2], seqDetailList[3], seqDetailList[-1]])

        # Extract and print the details of the gene where the gene name equals to "DREB1A_CDS"
        if seqDetailList[0] == "DREB1A_CDS":
            print("The DREB1A DNA sequence; ")
            # print gene id
            print("\tGene ID: ", seqObjcts.gene_id)
            # print sequence length
            print("\tSequence Length: ", seqObjcts.sequence_length, "bp")
            # print sequence type
            print("\tSequence Type: ", seqObjcts.get_Seq_Type())

            # print character count using the 'get_Character_Count' method
            print("\nThe base count of the four bases of the DREB1A coding sequence: \n\t",
                  seqObjcts.get_Character_Count())


