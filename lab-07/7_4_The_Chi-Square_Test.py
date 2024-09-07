'''
lab 7-4
The Chi-Square Test
Input : The dataset (Observed frequencies of fish eaten or not eaten by birds according to trematode infection level is given)
Output : The results of the Chi-Square Test
Author : Nimna Alupotha Gamage
Index No.: s14682
'''

# Import packages/sub modules

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from statsmodels.graphics.mosaicplot import mosaic

# Create the contingency table in a Pandas DataFrame

data_ct = pd.DataFrame([[1, 10, 37], [49, 35, 9]], index=["Eaten by birds", "Not eaten by birds"], columns=["Uninfected", "Lightly infected", "Highly infected"])
print("Contingency Table: \n", data_ct)

# Mosaic Plot

fig, axs = plt.subplots(figsize=(9, 6))
mos_plot = mosaic(data_ct.stack(), title="Mosaic Plot of data", axes_label=True, ax=axs)
#plt.savefig("q4_mosaicPlot.jpg")
plt.show()

# Chi-square contingency test

chiSq_results = stats.chi2_contingency(data_ct, correction=True, lambda_=None)
print("\nThe Chi-square test results: \n", chiSq_results)
print("The Chi-square statistic: \n", chiSq_results.statistic)
print("The P-value: \n", chiSq_results.pvalue)
print("The degree of freedom value: \n", chiSq_results.dof)

# Output the expected value table

ex_value_df = pd.DataFrame(chiSq_results.expected_freq, index=["Eaten by birds", "Not eaten by birds"], columns=["Uninfected", "Lightly infected", "Highly infected"])
print("\nThe expected value table: \n", ex_value_df)




