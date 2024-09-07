'''
lab 8 - Question 01
Implement the K-means clustering algorithm on the famous iris data set
Input : The iris data (“iris.csv” file)
Output : The results/ predictions of the K-means clustering algorithm
Date : 01/02/2023
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# Sub-question 01
# Import packages/sub modules

import numpy as np
from sklearn.cluster import KMeans
import matplotlib. pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler

# read the iris data file into a Pandas DataFrame

irisData = pd.read_csv("iris.csv")
#print(irisData)

# Remove null/missing values from the data set

iris = irisData.dropna(inplace=False)

# sub-question 02
# Implement the K-means algorithm on sepal data

# Extract sepal data

sepal_iris = iris.iloc[:,0:2]
#print(sepal_iris)

#Standardize sepal data

scaler = StandardScaler().fit(sepal_iris)
trans_sepal = scaler.transform(sepal_iris)
#print(trans_sepal)

# Pick the K value for the algorithm and train the model

kmeans = KMeans(n_clusters=3).fit(trans_sepal)
#print(kmeans)

# Centroids

centroids = kmeans.cluster_centers_
#print(centroids)
centroids_df = pd.DataFrame(centroids, columns=['sepal length', 'sepal width'])
print("The sepal length and sepal width values for the centroids: \n", centroids_df)

# sub-question 03
# Predicted Labels

testdata_sepal = np.array([[4.6, 3.0], [6.2, 3.0]])
trans_testdata = scaler.transform(testdata_sepal)
predicted = kmeans.predict(trans_testdata)
predicted_df = pd.DataFrame(predicted, index=['Plant 1', 'Plant 2'], columns=['Predicted label'])
print("\nThe predicted labels : \n", predicted_df)

# sub-question 04
# scatter plot for the sepal measurement data

labels = kmeans.labels_
variety = np.array(iris["v_short"])

plt.scatter(trans_sepal[:,0], trans_sepal[:,1], c=labels.astype(float), s=50, alpha=0.5)
plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=100)
plt.scatter(trans_testdata[:, 0], trans_testdata[:, 1], c='green', s=50)
#label data points
for i, txt in enumerate(labels):
    plt.annotate(txt, (trans_sepal[:, 0][i], trans_sepal[:, 1][i]))

for i, txt in enumerate(variety):
    plt.annotate(txt, (trans_sepal[:, 0][i], trans_sepal[:, 1][i]), xytext=(-18, 18), arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=0"),
                 textcoords='offset points')

for i, txt in enumerate(predicted):
    plt.annotate(txt, (trans_testdata[:, 0][i], trans_testdata[:, 1][i]))
plt.title("The scatter plot for the sepal measurement data")
plt.savefig("sepal_scatter.jpg")
plt.show()



# sub-question 05
# Implement the K-means algorithm on petal data

# Extract petal data

petal_iris = iris.iloc[:,2:4]
#print(petal_iris)

#Standardize petal data

scaler = StandardScaler().fit(petal_iris)
trans_petal = scaler.transform(petal_iris)
#print(trans_petal)

# Pick the K value for the algorithm and train the model

kmeans_p = KMeans(n_clusters=3).fit(trans_petal)
#print(kmeans_p)

# Centroids

centroids_p = kmeans_p.cluster_centers_
#print(centroids_p)
centroids_p_df = pd.DataFrame(centroids_p, columns=['petal length', 'petal width'])
print("The petal length and petal width values for the centroids: \n", centroids_p_df)

# Predicted Labels

testdata_petal = np.array([[1.5, 0.2], [4.1,1.2]])
trans_testdata_p = scaler.transform(testdata_petal)
predicted_p = kmeans.predict(trans_testdata_p)
predicted_p_df = pd.DataFrame(predicted_p, index=['Plant 1', 'Plant 2'], columns=['Predicted label'])
print("\nThe predicted labels : \n", predicted_p_df)

# sub-question 06
# scatter plot for the petal measurement data

labels_p = kmeans_p.labels_
variety = np.array(iris["v_short"])

plt.scatter(trans_petal[:,0], trans_petal[:,1], c=labels_p.astype(float), s=50, alpha=0.5)
plt.scatter(centroids_p[:, 0], centroids_p[:, 1], c='red', s=100)
plt.scatter(trans_testdata_p[:, 0], trans_testdata_p[:, 1], c='green', s=50)
#label data points
for i, txt in enumerate(labels_p):
    plt.annotate(txt, (trans_petal[:, 0][i], trans_petal[:, 1][i]))

for i, txt in enumerate(variety):
    plt.annotate(txt, (trans_petal[:, 0][i], trans_petal[:, 1][i]), xytext=(-18, 18), arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=0"),
                 textcoords='offset points')

for i, txt in enumerate(predicted_p):
    plt.annotate(txt, (trans_testdata_p[:, 0][i], trans_testdata_p[:, 1][i]))
plt.title("The scatter plot for the petal measurement data")
plt.savefig("petal_scatter.jpg")
plt.show()
