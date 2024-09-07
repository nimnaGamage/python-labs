'''
lab 4
Python Sequence class and subclasses to store any biological sequence (DNA, mRNA, and amino acid sequences)
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


#Subclasses of the Sequence class

# A subclass for DNA sequences
#Subclass-"DNAseq"
class DNAseq(Sequence):

    #Constructer method to create DNA sequence objects
    def __init__(self, gene_name, gene_id, species_name, subspecies_name, sequence):

        # Using the super method to pass object variables to the parent class
        super().__init__(gene_name, gene_id, species_name, subspecies_name, sequence)

        #unique/additional attributes
        self.at_content = self.get_ATcontent()
        self.transcribed_sequence = self.transcribe_Sequence()


    #unique/additional methods

    #transcribe the given DNA sequence into its mRNA sequence
    def transcribe_Sequence(self):



        # replace 'T' with 'U'
        # Store the new sequence in the transcribed_Sequence instance variable
        transcribed_sequence = self.sequence.replace("T", "U")

        # return the transcribed_Sequence instance variable
        return transcribed_sequence


    #return the AT content of the given sequence
    def get_ATcontent(self):

        # Assign a variable to len(input_sequence) to get the total base count
        baseTotal = len(self.sequence)

        # # Assign a variable to zero to count the number of 'A' and 'T' bases
        # baseCountAT = 0
        #
        # # for loop to go through each base
        # for base in self.sequence:
        #     # Increment the 'baseCountAT' by one only if base A or base T encountered
        #     if (base == 'A' or base == "T"):
        #         baseCountAT += 1


        # Count the number of 'A' and 'T' bases and assign it to a variable
        A_count = self.sequence.count('A')
        T_count = self.sequence.count('T')
        baseCountAT = A_count + T_count

        # get the AT Content as a percentage value from the total base count
        # update the at_content instance variable
        at_content = (baseCountAT / baseTotal)

        # Return the at_content instance variable
        return at_content


#Subclass-"MRNAseq"
class MRNAseq(Sequence):

    # Attributes/Variables

    # Class variable
    #hidden/private variable/dictionary to store codon-amino acid pairs
    # hide the dict for protect the variable from getting overriden in subclasses
    __amino_acid_codons = {}

    # Constructer method to create mRNA sequence objects
    def __init__(self, gene_name, gene_id, species_name, subspecies_name, sequence):

        super().__init__(gene_name, gene_id, species_name, subspecies_name, sequence)

        # unique/additional attributes
        self.at_content = self.get_ATcontent()
        self.translated_sequence = self.translate_Sequence()


    #unique/additional methods

    # return the AT content of the given sequence
    def get_ATcontent(self):

        # #replace 'U' with 'T'
        # mRNASeq = sequence.replace("U", "T")
        # count the Adenine bases
        baseCountA = self.sequence.count("A")
        # count the Thymine bases
        baseCountT = self.sequence.count("U")

        # total A and T count
        baseCountAT = baseCountA + baseCountT
        # total base count
        total = len(self.sequence)

        # update the at_content object variable
        at_content = (baseCountAT / total)

        # return the at_content object variable
        return at_content

    #class method
    @classmethod
    # store codon-amino acid pairs from a text file
    # into the Amino_acid_codons variable
    def upload_Codons(cls, codonTableFile):

        #codonTableFile = "codon_table.txt"

        # Read the 'codonTableFile'
        with open(codonTableFile, "r") as file_codon:
            for line in file_codon.readlines():
                if "#" not in line and line != "\n":
                    #all lines are stripped and splitted and form a list
                    newline = line.strip().split("\t")

                    #dictionary key
                    codon = newline[0]
                    #dictionary value
                    aaletter = newline[2]
                    #'__amino_acid_codons' dictionary of codon-amino acid (key-value) pairs
                    cls.__amino_acid_codons[codon] = aaletter

        #return the '__amino_acid_codons' private variable
        return cls.__amino_acid_codons


    #translate a given mRNA sequence into its amino acid sequence
    def translate_Sequence(self):

        #aadict = MRNAseq.upload_Codons(codonTableFile)

        #assign an empty string to the 'translated_sequence' variable
        self.translated_sequence = ""

        # read the mRNA sequence in 3 letter codons
        for i in range(0, len(self.sequence), 3):
            triple_seq = self.sequence[i:i + 3]
            # break when hit a stop codon
            if triple_seq in ["UAA", "UAG", "UGA"]:
                break
            # update the 'translated_sequence' variable
            self.translated_sequence += (aadict[triple_seq])


        # return the 'translated_sequence' variable
        return self.translated_sequence



#Subclass-"Proteinseq"
class Proteinseq(Sequence):

    # Constructor method to create Protein sequence objects
    def __init__(self, gene_name, gene_id, species_name, subspecies_name, uniprot_id, reviewed_status, sequence):
        # Attributes
        super().__init__(gene_name, gene_id, species_name, subspecies_name, sequence)

        # unique/additional attributes
        self.uniprot_id = uniprot_id
        self.reviewed_status = reviewed_status
        self.hydrophobicity = self.get_Hydrophobicity()


    # Unique/additional methods

    # return the percentage of the total hydrophobic amino acid residues
    def get_Hydrophobicity(self):

        #sequence = self.sequence

        # # get the total hydrophobic amino acid residues (A, I, L, M, F, W, Y, V) in the sequence
        # hydrophobicity_Content = sequence.count("A") + sequence.count("I") + sequence.count("L") + sequence.count("M") + \
        #                          sequence.count("F") + sequence.count("W") + sequence.count("Y") + sequence.count("V")

        hydrophobic_aa_list = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']

        protein_length = len(self.sequence)
        total = 0
        for hydro_aa in hydrophobic_aa_list:
            count = self.sequence.count(hydro_aa)
            total += count


        # update the 'hydrophobicity' variable with the percentage hydrophobocity
        hydrophobicity = (total / protein_length) * 100

        # return the 'hydrophobicity' variable
        return hydrophobicity


#main method
if __name__ == "__main__":

    '''Creating sequence objects from a fasta file'''


    #Read the file and store it in a variable
    file = "OSDREB_sequences.FASTA"

    #Give the above fastafile contained variable as an argument to the fasta_Split method/ function
    #Store the resulting dictionary in the 'fastaSeqDict' variable
    fastaSeqDict = Sequence.fasta_Split(file)
    #printing the fasta sequence dictionary
    #print(fastaSeqDict)

    # for loop to go through each value of the dictionary
    for seqData in fastaSeqDict.values():
        #print(seqData[-1])

        # define a set of amino acid letters excluding A, T, G, C amino acids; to improve this better remove 'N'
        listAminoAcids = ['K', 'N', 'R', 'I', 'Q', 'M', 'H', 'P', 'L', 'E', 'D', 'V', 'Y', 'S', 'W', 'F']

        if (1 in [letter in seqData[-1] for letter in listAminoAcids]):
            # Creating objects using 'Proteinseq' subclass
            protein_Objcts = Proteinseq(seqData[0], seqData[1], seqData[2], seqData[3], seqData[4], seqData[5],
                                        seqData[-1])
            # Question_02-subquestion_04
            # printing the relevent detailes of the 'DREB2A_P' protein
            if seqData[0] == "DREB2A_P":
                print("The DREB2A Protein Sequence; ")
                print("\tThe Uniprot ID: ", protein_Objcts.uniprot_id)
                print("\tThe Reviewed Status: ", protein_Objcts.reviewed_status)
                print("\tThe Sequence Type: ", protein_Objcts.sequence_type)
                print("\tThe Amino Acid Composition: ", protein_Objcts.get_Character_Count())
                print("\tThe Hydrophobicity: ", protein_Objcts.get_Hydrophobicity(), "%")
        # elif 'U' in seqData[-1]:
        #     #return 'mRNA Sequence'
        else:
            # Creating objects using 'DNAseq' subclass
            dna_Objcts = DNAseq(seqData[0], seqData[1], seqData[2], seqData[3], seqData[-1])
            # Question_02-subquestion_01
            # printing the relevent detailes of the 'OSDREB1A' DNA sequence
            if seqData[0] == "DREB1A_CDS":
                print("The DREB1A DNA sequence; ")
                print("\tGene ID: ", dna_Objcts.gene_id)
                print("\tSequence Length: ", dna_Objcts.sequence_length, "bp")
                print("\tSequence Type: ", dna_Objcts.sequence_type)
                print("\tAT Content: ", dna_Objcts.get_ATcontent())

            # Question_02-subquestion_02
            if seqData[0] == "DREB2B.isoform 1_CDS":
                # transcribe the sequence
                transcribed_mRNA = dna_Objcts.transcribe_Sequence()

                #defined the codon table text file
                codonTableFile = "codon_table.txt"
                #defined the 'aadict' dictionary
                aadict = MRNAseq.upload_Codons(codonTableFile)

                # Creating a new object for the resulting transcribed_mRNA sequence
                mRNA_Objct = MRNAseq(seqData[0], seqData[1], seqData[2], seqData[3], transcribed_mRNA)

                # printing the relevent detailes of the transcribed_mRNA sequence
                print("The DREB2B.isoform 1 mRNA sequence; ")
                print("\tSequence Length: ", mRNA_Objct.sequence_length, "bp")
                print("\tSequence Type: ", mRNA_Objct.sequence_type)
                print("\tAT Content: ", mRNA_Objct.at_content)
                print("\tTranscribed mRNA Sequence: ", transcribed_mRNA)

                # Question_02-subquestion_03
                print("\tAmino Acid Sequence: ", mRNA_Objct.translate_Sequence())
                print("\tLength of the Amino Acid Sequence: ", len(mRNA_Objct.translate_Sequence()), "aa")

    # Question_02-subquestion_05
    # number of sequences created
    print("No. of Sequence Objects created: ", Sequence.sequence_count)


