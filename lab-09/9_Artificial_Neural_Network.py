'''
Lab 9
Implementing the Artificial Neural Network (ANN)
Input : The iris data (“iris.csv” file)
Output : The results/ predictions of the Artificial Neural Network (ANN)
Date : 07/02/2023
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# subQuestion 1.1
# Import packages/sub modules

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report, confusion_matrix

# Pre-processing the Iris dataset

# Load the Iris dataset into a Pandas DataFrame from the given CSV file

iris_data = pd.read_csv("iris.csv")
# print(iris_data)

# Remove null/ missing values from the dataset

iris_Data = iris_data.dropna(inplace=False)

# split the measurement data and class labels (species) to X and Y variables

# X
X = iris_Data.iloc[:, 0:4]

# y
Y = iris_Data.iloc[:, 5]

# subQuestion 1.2
# Output the species category values of the Y variable

print("The species category values of the Y variable : \n", Y.unique())

# replace the category values with integer numbers

le_Y = LabelEncoder().fit(Y)
y = le_Y.transform(Y)

# output the Y variable values

print("\nValues of the 'Y' variable after replacing with integer numbers : \n", pd.unique(y))

# subQuestion 1.3
# split the Iris dataset (X and Y variables) to train and test data

# Split to train and test data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, train_size=0.80)
# print(X_train)
# print(X_test)

# subQuestion 1.4
# Standardize the training and testing data of the input variable

scaler = StandardScaler()
scaler.fit(X_train)

X_train_std = scaler.transform(X_train)
X_test_std = scaler.transform(X_test)
# print(X_train_std)
# print(X_test_std)

# subQuestion 02
# Train an ANN that contains 3 hidden layers

# model
mlp = MLPClassifier(hidden_layer_sizes=(10, 10, 10), max_iter=1000)
# fit the model
mlp. fit(X_train_std, y_train.ravel())

# use the trained ANN for predicting the species labels of the test data

predictions = mlp.predict(X_test_std)
#print(predictions)

# subQuestion 03
# Evaluate the ANN using the classification report and the confusion matrix

# The classification report
print("\nThe classification report : \n", classification_report(y_test, predictions))
# The confusion matrix
print("\nThe confusion matrix : \n", confusion_matrix(y_test, predictions))

# subQuestion 04
# predict the species of the 3 plants

test_data = np.array([[5.9, 3.0, 7.0, 5.0], [4.6, 3.0, 1.5, 0.2], [6.2, 3.0, 4.1, 1.2]])
testData_std = scaler.transform(test_data)
prdTestData = mlp.predict(testData_std)
print("\nPredicted species of the 3 plants; \n Plant 01 - %s \n Plant 02 - %s \n Plant 03 - %s  " % (prdTestData[0], prdTestData[1], prdTestData[2]))





'''
Implementing the Artificial Neural Network (ANN)
Input : The iris data (“iris.csv” file)
Output : The results/ predictions of the Artificial Neural Network (ANN)
Date : 07/02/2023
Author : Nimna Alupotha Gamage
Index No.: s14682
lab 11 - Question 01 - SubQuestion 05
'''

# Import packages/sub modules

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report, confusion_matrix

# Pre-processing the Iris dataset

# Load the Iris dataset into a Pandas DataFrame from the given CSV file

iris_data = pd.read_csv("iris.csv")

# Remove null/ missing values from the dataset

iris_Data = iris_data.dropna(inplace=False)

# split the measurement data and class labels (species) to X and Y variables

# X
X = iris_Data.iloc[:, 0:4]

# y
Y = iris_Data.iloc[:, 5]

# replace the category values with integer numbers

le_Y = LabelEncoder().fit(Y)
y = le_Y.transform(Y)

# split the Iris dataset (X and Y variables) to train and test data

# Split to train and test data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, train_size=0.80)

# Standardize the training and testing data of the input variable

scaler = StandardScaler()
scaler.fit(X_train)

X_train_std = scaler.transform(X_train)
X_test_std = scaler.transform(X_test)

# subQuestion 02
# Train an ANN that contains 3 hidden layers

# model
mlp2 = MLPClassifier(hidden_layer_sizes=(2, 2, 2), max_iter=1000)
# fit the model
mlp2. fit(X_train_std, y_train.ravel())

# use the trained ANN for predicting the species labels of the test data

predictions2 = mlp2.predict(X_test_std)

# subQuestion 03
# Evaluate the ANN using the classification report and the confusion matrix

# The classification report
print("\nThe classification report : \n", classification_report(y_test, predictions2))
# The confusion matrix
print("\nThe confusion matrix : \n", confusion_matrix(y_test, predictions2))

