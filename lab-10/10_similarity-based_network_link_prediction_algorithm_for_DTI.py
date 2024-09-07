'''
Lab 10
Implementing a similarity-based network link prediction algorithm for Drug-Target Interactions (DTI)
Input: The "DTIsubset.tsv" file
Output: CN scores for novel drug-target interactions in descending order
Name: Nimna Alupotha Gamage
Date: 15/02/2023
Index No.: s14682
'''

# import python packages

import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict

# Constructing the empty DTI network using NetworkX

dti_network = nx.Graph()

# Create lists

drug_list = []
target_list = []

# read the tsv file and construct the DTI network

# read the tsv file
with open("DTIsubset.tsv", "r") as file:
    # for each line in file
    for line in file:
        # if line is noe empty
        if line != '\n':
            # remove unwanted characters and split
            line = line.strip().split("\t")
            if "Drug_name" not in line:
                # add drugs in to a separate list
                drug_list.append(line[0])
                # remove duplicates
                drug_list = list(set(drug_list))
                # add proteins/targets in to a separate list
                target_list.append(line[1])
                # remove duplicates
                target_list = list(set(target_list))
                # add nodes and edges to the network
                dti_network.add_nodes_from(drug_list, bipartite=0)
                dti_network.add_nodes_from(target_list, bipartite=1)
                dti_network.add_edge(line[0], line[1], weight=int(line[2]))


# plot the graph/network

plt.figure(figsize=(15, 10))
nx.draw_networkx(dti_network, with_labels=True, node_size=200, width=2)
plt.savefig("DTI_Network.png")
plt.show()

plt.figure(figsize=(15, 10))
nx.draw_networkx(dti_network, pos=nx.drawing.layout.bipartite_layout(dti_network, target_list), width=2)
plt.savefig("DTI_Network_bipartite.png")
plt.show()

# print(dti_network.degree)

# Create the dictionary of protein neighbors of drug nodes

drug_neighbours = {}

for drug in drug_list:
    drug_n = list(dti_network.neighbors(drug))
    # get the neighbours of drugs into a dictionary
    drug_neighbours[drug] = drug_n

#print(drug_neighbours)

# create the dictionary of protein neighbors of protein nodes

target_neighbours = {}

for target in target_list:
    nbr_nghbors = []
    target_n = list(dti_network.neighbors(target))
    for drug1 in target_n:
        drug1_n = list(dti_network.neighbors(drug1))
        nbr_nghbors.extend(drug1_n)
        # get the neighbours of target neighbours into a dictionary
        target_neighbours[target] = nbr_nghbors

#print(target_neighbours)

# create the dictionary of common neighbors interaction scores

dictCN = {}

for dg in drug_list:
    for tg in target_list:
        # remove duplicates of the drug_neighbour list
        drug_nb = list(drug_neighbours[dg])
        if tg not in drug_nb:
            # get the targets which are in the target_neighbours list and not in the drug_neighbours list
            # remove duplicates
            target_nb = list(target_neighbours[tg])
            # get the common neighbours
            commonNeighbors = set(drug_nb).intersection(target_nb)
            # count the common neighbours
            score = len(commonNeighbors)
            #print(score)
            dictCN_key = dg + " and " + tg
            dictCN[dictCN_key] = score

# get the scores in descending order

sorted_CN_dict = OrderedDict(sorted(dictCN.items(), key=lambda t: t[1], reverse=True))
#print(sorted_CN_dict)

print("CN scores for novel drug-target interactions in descending order ; ")
for dg_tg_pair in sorted_CN_dict:
    print(dg_tg_pair, " : ", sorted_CN_dict[dg_tg_pair])

# write the results into a file
with open ("CN.txt", 'w') as file:
    # write the header
    file.write("drug_target_combination\t\tCN Score\n")
    # writing the vote count for all unknown proteins
    for dg_tg_pair in sorted_CN_dict:
        file.write("%s\t%s\n" %(dg_tg_pair, sorted_CN_dict[dg_tg_pair]))



