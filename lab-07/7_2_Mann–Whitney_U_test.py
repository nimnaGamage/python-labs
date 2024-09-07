'''
lab 7-2
Mann–Whitney U test
Input : The dataset (“HornedLizards.csv” file)
Output : The results of the Normality test and the Mann–Whitney U test
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# Import packages/sub modules

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.graphics.gofplots import qqplot
from scipy import stats

# Read the csv file and import the data set into a Pandas DataFrame

dataset = pd.read_csv("HornedLizards.csv")

# to identify number of null values

print(dataset.info())

# Remove the missing values in the data sets

df = dataset.dropna(inplace=False)

# Statistics for 'survived' variable

data_survived = df.loc[df['Survive'] == 'survived', 'Squamosal horn length']
print("statistics for horn lengths of survived lizards: \n", data_survived.describe())

# Statistics for 'dead' variable

data_dead = df.loc[df['Survive'] == 'dead', 'Squamosal horn length']
print("statistics for horn lengths of dead lizards: \n", data_dead.describe())

# Testing the normality assumption for the two independent samples

# Histograms for both variables in the same plot

fig1, axes = plt.subplots()
for plot in [data_survived, data_dead]:
    sns.histplot(plot, bins=30, kde=True, ax=axes)
axes.set_title("Histogram of horn lengths")
plt.legend(['Survived', 'Dead'])
#plt.savefig("q2_histogram_samePlot.jpg")
plt.show()

# Histograms for both variables in the same figure

fig, axs = plt.subplots(1, 2, figsize=(12, 4))

sns.histplot(data=data_survived, bins=30, kde=True, ax=axs[0])
axs[0].set_title("Histogram of horn lengths of survived lizards")

sns.histplot(data=data_dead, bins=30, kde=True, ax=axs[1])
axs[1].set_title("Histogram of horn lengths of dead lizards")

#plt.savefig("q2_histogram_sameFigure.jpg")
plt.show()

# QQplots for variables

# survived
qqplot(data=data_survived, line='s')
plt.title("QQplot of horn lengths of survived lizards")
#plt.savefig("q2_qqplot_survived.jpg")
plt.show()

# dead
qqplot(data=data_dead, line='s')
plt.title("QQplot for horn lengths of dead lizards")
#plt.savefig("q2_qqplot_dead.jpg")
plt.show()

# Shapiro-Wilk test for normality

# survived
stat_s, p_s = stats.shapiro(data_survived)
print("\nThe normality test for horn lengths of survived lizards: \n statistic=%.3f p-value=%.3f" % (stat_s, p_s))

# dead
stat_d, p_d = stats.shapiro(data_dead)
print("\nThe normality test for horn lengths of dead lizards: \n statistic=%.3f p-value=%.3f" % (stat_d, p_d))

# Comparison of means

# Box plot

bp = sns.boxplot(data=df, x='Survive', y='Squamosal horn length', showmeans=True)
plt.title("Box plot comparison of squamosal horn length between dead and survived lizards")
plt.ylabel("squamosal horn length")
#plt.savefig("q2_boxPlot.jpg")
plt.show()

# Violin plot

vp = sns.violinplot(data=df, x='Survive', y='Squamosal horn length', linewidth=1)
plt.title("Violin plot comparison of squamosal horn length between dead and survived lizards")
plt.ylabel("squamosal horn length")
#plt.savefig("q2_violinPlot.jpg")
plt.show()

# Mann–Whitney U test

U_m, p_m = stats.mannwhitneyu(data_dead, data_survived, alternative='less')
print("\nThe Mann-Whitney U test: \n statistic=%.3f  p-value=%.6f" % (U_m, p_m))




