'''
lab 7-1
One sample t-test
Input : The dataset (“Temperature.csv” file)
Output : The results of the Normality test and the One sample t-test
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# Import packages/sub modules

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from statsmodels.graphics.gofplots import qqplot
from scipy import stats

# read the csv file and import the dataset into a Pandas DataFrame

temp = pd.read_csv("Temperature.csv")

# print(temp)
print("The statistics for the human temperature variable: \n", temp.describe())

# Testing the normality assumption for the temperature variable

# Histogram for the temperature variable

fig, axs = plt.subplots(figsize=(8, 4))
sns.histplot(data=temp, ax=axs, kde=True)
axs.set_title("Histogram for the temperature variable")
#plt.savefig("q1_histogram.jpg")
plt.show()

# QQplot for the temperature variable

qqplot(temp['temperature'], line='s')
plt.title("QQplot for the temperature variable")
#plt.savefig("q1_qqplot.jpg")
plt.show()

# The Shapiro-Wilk test for normality on the variable

stat, p = stats.shapiro(temp)
print("\nThe Shapiro-Wilk test: \n Statistic=%.3f  p-value=%.3f" % (stat, p))

# One sample t-test

checkValue = 98.6
t, p1 = stats.ttest_1samp(temp, checkValue)
print("\nThe One sample t-test: \n Statistic=%.3f  p-value=%.3f" % (t, p1))



