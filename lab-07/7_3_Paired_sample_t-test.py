'''
lab 7-3
Paired sample t-test
Input : The dataset (“BlackbirdTestosterone.csv” file)
Output : The results of the Normality test and the Paired sample t-test
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# Import packages/sub modules

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from statsmodels.graphics.gofplots import qqplot
from scipy import stats

# Read the csv file and import the data set into a Pandas DataFrame

dataset_Q3 = pd.read_csv("BlackbirdTestosterone.csv")

# Remove the missing values in the data sets

df = dataset_Q3.dropna(inplace=False)

# Statistics for each variable

# log before
log_bfr = df['log before']
print("statistics for log before: \n ", log_bfr.describe())

# log after
log_after = df['log after']
print("statistics for log after: \n", log_after.describe())

# log difference
log_difference = df['dif in logs']
print("statistics for difference in log: \n", log_difference.describe())

# Testing the normality assumption for the log difference variable

# histogram

fig1, axs1 = plt.subplots(figsize=(8, 4))
sns.histplot(data=log_difference, ax=axs1, kde=True, bins=30)
plt.title("Histogram of log_difference")
plt.xlabel("log difference")
# plt.savefig("q3_histogram.jpg")
plt.show()

# qqplot

qqplot(log_difference, line='s')
plt.title("QQplot of log_difference")
plt.savefig("q3_qqplot.jpg")
plt.show()

# Shapiro-Wilk test for normality

stat, p = stats.shapiro(log_difference)
print("\nThe normality test for log_difference: \n Statistic=%.3f  p-value=%.3f" % (stat, p))

# Comparison of means

# Box plot

fig2, axs2 = plt.subplots(figsize=(5, 5))
sns.boxplot(data=df[['log before', 'log after']], ax=axs2, showmeans=True)
axs2.set_title("Boxplot of the data")
plt.ylabel("Rate of antibody production")
plt.savefig("q3_boxplot.jpg")
plt.show()

# Violin plot
fig3, axs3 = plt.subplots(figsize=(5, 5))
sns.violinplot(data=df[['log before', 'log after']], ax=axs3, linewidth=1)
axs3.set_title("Violin plot of data")
plt.ylabel("Rate of antibody production")
plt.savefig("q3_violinplot.jpg")
plt.show()

# Paired sample t-test

stat1, p1 = stats.ttest_rel(log_bfr, log_after, alternative='greater')
print("\nThe paired sample t-test: \n statistic=%.3f  p-value=%.6f" % (stat1, p1))


