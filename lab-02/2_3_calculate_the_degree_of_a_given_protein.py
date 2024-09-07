'''
Lab 2-3
Using the graph data structure in Python to calculate the degree of a given protein
Input: TSV file
Output: The degree of the rice DREB1A protein
Nimna Gamage
'''

# Step 1: Import the networkx package

import networkx as nx

# Step 2: Create an empty graph with no nodes and no edges.

Graph_network = nx.Graph()

# Step 3: Create a variable to get the count of the neighbors with more than 0.7 combined score and assign it to zero

neighbour_Combined_Score = 0

# Step 3: Read lines in tsv file
# Step 4: Grow the graph by adding edges

# read the input tsv file
with open("DREB1A_protein.tsv",'r') as file:
    for line in file:
        # remove the commented line
        if "#" not in line:
            # remove the unwanted characters and split
            line = line.strip().split('\t')
            # add edges to the newly created graph
            Graph_network.add_edge(line[0], line[1])
            # add edges to the newly created graph with their respective weights
            Graph_network.add_edge(line[0], line[1], weight=line[12])


# Step 5: Write the graph to a file

nx.write_gml(Graph_network, "DREB1A_network.gml")

# Step 6: Output the degree of the 'ERF24' protein

print("The degree of the rice DREB1A protein is : ", Graph_network.degree("ERF24"), "\n")

# Step 7: Create a list and assign all the neighbours to the list

n_list = Graph_network.neighbors("ERF24")

# Step 8: Increment the counter by one, for the neighbors with more than 0.7 combined score for the DREB1A protein

for neighbour in n_list:
    if Graph_network.edges[("ERF24", neighbour)]['weight'] > str(0.7):
        neighbour_Combined_Score += 1
        print("The combined scores which are greater than 0.7 ; ", Graph_network.edges[("ERF24", neighbour)]['weight'])

print("The number of neighbors with more than 0.7 combined score for the DREB1A protein is : ", neighbour_Combined_Score)

