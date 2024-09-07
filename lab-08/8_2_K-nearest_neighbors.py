'''
lab 8 - Question 02
K-nearest neighbors (KNN) algorithm
Input : The dataset (“iris.csv” file)
Output : The results/ predictions of the K-nearest neighbors (KNN) algorithm
Date : 01/02/2023
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# SubQuestion 01
# Import packages/sub modules

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn import datasets

# Load the bundled iris data set

iris = datasets.load_iris()
#print(iris)

# Remove null/missing values from the data set

#iris = irisData.dropna(inplace=False)

# Load the iris data

irData = iris.data
#print(irData)

# Extract the labels

irLabels = iris.target
#print(irLabels)

# standardize the data

scaler = StandardScaler().fit(irData)
stdIrdata = scaler.transform(irData)
#print(stdIrdata)

# Read the Plant 1 data record as a test data record

testdata = [[4.6, 3.0, 1.5, 0.2]]

# standardize the test data

stdTestdata= scaler.transform(testdata)
#print(stdTestdata)

# Create the 2-d array for the two data records (Plant 1 and Plant 2)

twoDArray = np.array([[4.6, 3.0, 1.5, 0.2], [6.2, 3.0,4.1,1.2]])

# Standardize the array

stdArray= scaler.transform(twoDArray)


'''Train the KNN model using only the sepal measurements'''


# SubQuestion 02-01
# Train the KNN model using only the sepal measurements

# Extract the standardized sepal data

sepalData = stdIrdata[:, [0, 1]]
#print(sepalData)

# Train the KNN model to find the two nearest neighbors of the sepal data of plant 1

nn2_sepal = NearestNeighbors(n_neighbors=2).fit(sepalData, iris["target"])
#print(nn2_sepal)

# print the indices, data values, labels, and species of the two-nearest neighbors

testSepal = stdTestdata[:, [0, 1]]
nn2_distances_s, nn2_indices_s = nn2_sepal.kneighbors(testSepal)
# #print(iris)
# #print(nn2_distances_s)
# #print(nn2_indices_s)
# #lbl = iris['target'][nn2_indices_s]
# #print(iris['target_names'][lbl])

print("The data of the two-nearest neighbors of sepal data of plant 1 ;")
# Indices
print("Indices ; \tnn1 :", nn2_indices_s[0][0], "\tnn2 :", nn2_indices_s[0][1])

#Data values
print("Data values ; \tnn1 :", iris['data'][nn2_indices_s][0][0], "\tnn2 :", iris['data'][nn2_indices_s][0][1])

#labels
nn2_lbl_s = iris['target'][nn2_indices_s]
print("Labels ; \tnn1 :", nn2_lbl_s[0][0], "\tnn2 :", nn2_lbl_s[0][1])

#Species
species = iris['target_names']
nn2_sp_s = species[nn2_lbl_s]
print("Species ; \tnn1 :", nn2_sp_s[0][0], "\tnn2 :", nn2_sp_s[0][1])

# SubQuestion 02-02
# Train the KNN model to find the 5-nearest neighbors of the sepal dat of plant 1

nn5_sepal = KNeighborsClassifier(n_neighbors=5).fit(sepalData, iris["target"])
#print(nn5_sepal)

# predict the probability of each plant categorizing into each species

testSepalData = stdArray[:, [0, 1]]
nn5_distances_s, nn5_indices_s = nn5_sepal.kneighbors(testSepalData)

testSepalData_prob = nn5_sepal.predict_proba(testSepalData)
#print(testSepalData_prob)

prob_sepal_df = pd.DataFrame(testSepalData_prob, index=['Plant 01', 'Plant 02'], columns=species)
print("\nThe probability of each plant categorizing into each species are : \n", prob_sepal_df)

# predict the class labels for the two plants

prd_nn5_lbl_s = nn5_sepal.predict(testSepalData)
print("Predicted Class Labels ; \tPlant 01 :", prd_nn5_lbl_s[0], "\tPlant 02 :", prd_nn5_lbl_s[1])

# The species predicted from the KNN algorithm

prd_species_s = species[prd_nn5_lbl_s]
print("Predicted Species ; \tPlant 01 :", prd_species_s[0], "\tPlant 02 :", prd_species_s[1])


'''Train the KNN model using only the petal measurements'''


# SubQuestion 03-01
# Train the KNN model using only the petal measurements

# Extract the standardized petal data

petalData = stdIrdata[:, [2, 3]]
#print(petalData)

# Train the KNN model to find the two nearest neighbors of the sepal data of plant 1

nn2_petal = NearestNeighbors(n_neighbors=2).fit(petalData, iris["target"])
#print(nn2_petal)

# print the indices, data values, labels, and species of the two-nearest neighbors

testPetal = stdTestdata[:, [2, 3]]
nn2_distances_p, nn2_indices_p = nn2_petal.kneighbors(testPetal)

print("\nThe data of the two-nearest neighbors of petal data of plant 1 ;")
# Indices
print("Indices ; \tnn1 :", nn2_indices_p[0][0], "\tnn2 :", nn2_indices_p[0][1])

#Data values
print("Data values ; \tnn1 :", iris['data'][nn2_indices_p][0][0], "\tnn2 :", iris['data'][nn2_indices_p][0][1])

#labels
nn2_lbl_p = iris['target'][nn2_indices_p]
print("Labels ; \tnn1 :", nn2_lbl_p[0][0], "\tnn2 :", nn2_lbl_p[0][1])

#Species
species = iris['target_names']
nn2_sp_p = species[nn2_lbl_p]
print("Species ; \tnn1 :", nn2_sp_p[0][0], "\tnn2 :", nn2_sp_p[0][1])

# SubQuestion 03-02
# Train the KNN model to find the 5-nearest neighbors of the petal dat of plant 1

nn5_petal = KNeighborsClassifier(n_neighbors=5).fit(petalData, iris["target"])
#print(nn5_petal)

# predict the probability of each plant categorizing into each species

testPetalData = stdArray[:, [2, 3]]
nn5_distances_p, nn5_indices_p = nn5_petal.kneighbors(testPetalData)

testPetalData_prob = nn5_petal.predict_proba(testPetalData)
#print(testSepalData_prob)

prob_petal_df = pd.DataFrame(testPetalData_prob, index=['Plant 01', 'Plant 02'], columns=species)
print("\nThe probability of each plant categorizing into each species are : \n", prob_petal_df)

# predict the class labels for the two plants

prd_nn5_lbl_p = nn5_petal.predict(testPetalData)
print("Predicted Class Labels ; \tPlant 01 :", prd_nn5_lbl_p[0], "\tPlant 02 :", prd_nn5_lbl_p[1])

# The species predicted from the KNN algorithm

prd_species_p = species[prd_nn5_lbl_p]
print("Predicted Species ; \tPlant 01 :", prd_species_p[0], "\tPlant 02 :", prd_species_p[1])


'''Train the KNN model using both the sepal and petal measurements'''


# SubQuestion 04-01
# Train the KNN model for all the standardized measurement data to find the two nearest neighbors of the test data record/ Plant 1

nn2_test = NearestNeighbors(n_neighbors=2).fit(stdIrdata, iris["target"])
#print(nn2_test)

# print the indices, data values, labels, and species of the two-nearest neighbors

nn2_distances_test, nn2_indices_test = nn2_test.kneighbors(stdTestdata)

print("\nThe data of the two-nearest neighbors of sepal and petal data of plant 1 ;")
# Indices
print("Indices ; \tnn1 :", nn2_indices_test[0][0], "\tnn2 :", nn2_indices_test[0][1])

#Data values
print("Data values ; \tnn1 :", iris['data'][nn2_indices_test][0][0], "\tnn2 :", iris['data'][nn2_indices_test][0][1])

#labels
nn2_lbl_test = iris['target'][nn2_indices_test]
print("Labels ; \tnn1 :", nn2_lbl_test[0][0], "\tnn2 :", nn2_lbl_test[0][1])

#Species
species = iris['target_names']
nn2_sp_test = species[nn2_lbl_test]
print("Species ; \tnn1 :", nn2_sp_test[0][0], "\tnn2 :", nn2_sp_test[0][1])

# SubQuestion 04-02
# Train the KNN model to find the 5-nearest neighbors for all measurement data of plant 1

nn5_test = KNeighborsClassifier(n_neighbors=5).fit(stdIrdata, iris["target"])
#print(nn5_test)

# predict the probability of each plant categorizing into each species

nn5_distances_test, nn5_indices_test = nn5_test.kneighbors(stdArray)

testData_prob = nn5_test.predict_proba(stdArray)
#print(testData_prob)

prob_test_df = pd.DataFrame(testData_prob, index=['Plant 01', 'Plant 02'], columns=species)
print("\nThe probability of each plant categorizing into each species are : \n", prob_test_df)

# predict the class labels for the two plants

prd_nn5_lbl_test = nn5_test.predict(stdArray)
print("Predicted Class Labels ; \tPlant 01 :", prd_nn5_lbl_test[0], "\tPlant 02 :", prd_nn5_lbl_test[1])

# The species predicted from the KNN algorithm

prd_species_test = species[prd_nn5_lbl_test]
print("Predicted Species ; \tPlant 01 :", prd_species_test[0], "\tPlant 02 :", prd_species_test[1])


