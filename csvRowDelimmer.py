import csv

# We read in csv, we get those with list as value and add in blanks for those that don't

import pandas as pd

# Read the CSV file
df = pd.read_csv("tempCSV/Enriched Clusters - Full.csv")

split_values = ['","', "','"]


# Function to split rows based on values
def split_rows(row):
    values = row["Values"].split(",")
    for value in values:
        if value in split_values:
            new_row = row.copy()
            new_row["Values"] = value
            return new_row
    return row


# Apply the split_rows function to each row
df = df.apply(split_rows, axis=1)

# Reset index to clean up the resulting DataFrame
df = df.reset_index(drop=True)

# Save the result to a new CSV file
df.to_csv("tempCSV/EnrichedClustersFullFixed", index=False)
