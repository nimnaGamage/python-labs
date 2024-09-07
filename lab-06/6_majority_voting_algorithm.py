'''
Lab 6
Implementing the majority voting network-based candidate protein prediction algorithm
Input : A list of known proteins annotated to the particular function as a text file (“AT_stress_proteins.txt”)
        Protein-protein interaction data in tabular format as a TSV file (“DREB2A_short_tabular.tsv”)
Output : The list of unknown proteins with the predicted majority voting score
Nimna Alupotha Gamage
Index: s14682
'''

# Step 01: Import the networkx python package and OrderedDict submodule from the Collections Python package

import networkx as nx
from collections import OrderedDict

# Step 02: Create an empty graph with no nodes and no edges

dreb2a_network = nx.Graph()

# Step 03: Read the lines in TSV file and grow the graph by adding edges
# Create a network of Proteins using the Protein - protein interaction data obtained from the string database

# read the tsv file
with open("DREB2A_short_tabular.tsv", 'r') as file_network:
    # for each line in file
    for line in file_network:
        # remove empty lines and the lines with '#' sign
        if line != '\n' and "#" not in line:
            # remove unwanted characters
            line = line.strip().split('\t')
            node1 = line[0].upper()
            node2 = line[1].upper()
            dreb2a_network.add_edge(node1, node2, weight=line[12])

# Step 04: Create a list of nodes and remove duplicate nodes

network_proteins = list(set(dreb2a_network.nodes))

# Step 05: Read the text file and Create a list of function-known proteins using the text file given

known_proteins_list = []

# read the text file
with open("AT_stress_proteins.txt", 'r') as protein_file_known:
    # for each line in file
    for line in protein_file_known:
        # remove unwanted characters
        line = line.strip().split('\t')
        known_p = line[1].upper()
        known_proteins_list.append(known_p)

# Step 06: Remove duplicates from the 'known_proteins_list'

known_proteins = list(set(known_proteins_list))

# Step 07: Create a list of function-unknown proteins of the network by comparing both 'network_proteins' list and 'known_proteins' list.
# Elements of the 'known_proteins' list is substracted from the 'network_proteins' list and get the difference between two sets to a list ('unknown_proteins' list)

unknown_proteins = list(set(network_proteins) - set(known_proteins))

# Step 08: For each unknown protein,
# Count the number of function-known proteins interacting with each function-unknown protein and consider that count as majority voting score of each function-unknown protein
# Step 09: Create a dictionary by considering each function-unknown protein as a key and the majority voting score/ count as their respective value

unknown_P_dict = {}

for unknown_node in unknown_proteins:
    score = 0
    # if neighbours of unknown proteins are in the list of known proteins, score is incremented by one
    for neighbor in dreb2a_network.neighbors(unknown_node):
        if neighbor in known_proteins:
            score += 1
    unknown_P_dict[unknown_node] = score

# Step 10: The majority voting scores should be sorted in descending order based on the scores,
# with proteins with high scores at the top

sorted_up_dict = OrderedDict(sorted(unknown_P_dict.items(), key=lambda t: t[1], reverse=True))

# Step 11: write the ordered list to an output file

mvFile = open("AT_stress_proteins_majority_voting.txt", "w")
mvFile.write("Protein" + "\t" + "The majority voting score" + "\n")
for protein, score in sorted_up_dict.items():
    mvFile.write(protein + "\t" + str(score) + "\n")
mvFile.close()

# Assign the 'DREB2A' protein to a variable

protein = "DREB2A"

# Console Output

print("The degree of the 'ATDREB2A' protein is: ", dreb2a_network.degree(protein))
print("The number of unknown proteins in the network for stress tolerance is: ", len(unknown_proteins), "proteins")
print("The unknown protein with highest majority vote score is: ", list(sorted_up_dict.items())[0])



