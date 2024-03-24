# read in from Trubetskoy2022broadcoding_significant_rows_sorted.csv

# we want only positive fl,
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv("Ora/Enriched/Trubetskoy2022broadcoding_significant_rows_sorted.csv")

# Filter the data for rows where FL is True
data = data[data['FL'] == True]

# Combine 'alg' and 'cl' into one column
data['alg_cl'] = data['alg'] + '_' + data['cl'].astype(str)

# Split 'overlapGenes' into separate rows
data = data.assign(overlapGenes=data['overlapGenes'].str.split(', ')).explode('overlapGenes')

# Create a pivot table to prepare data for heatmap
pivot_data = data.pivot_table(index='alg_cl', columns='overlapGenes', values='pval', aggfunc='mean')

# Disregard rows and columns with only one value
# pivot_data = pivot_data.loc[:, pivot_data.nunique() > 1]
# #pivot_data = pivot_data.loc[pivot_data.nunique(axis=1) > 1]

pivot_data = pivot_data.dropna(axis=0, how='all')
pivot_data = pivot_data.dropna(axis=1, how='all')

print(pivot_data)

# Create the heatmap
plt.figure(figsize=(18, 10))
heatmap = sns.heatmap(pivot_data, cmap="viridis", annot=False, fmt=".2f", xticklabels=True, yticklabels=True, cbar=True, linecolor='grey', linewidth=0.5)

# Set labels and title
plt.title('Heatmap of pval value for each enriched algorithm cluster\'s SCH genes', fontsize=16)
plt.ylabel('Algorithm and Cluster Number', fontsize=14)
plt.xlabel('SCH Enriched Genes', fontsize=14)

# Rotate the x-axis labels for better readability
plt.xticks(rotation=90, fontsize=12)
plt.yticks(fontsize=12)

# Show the plot
plt.savefig("BroadCodingEnrichmentHeatmap.png")
