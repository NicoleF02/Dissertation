import pandas as pd

# Read the CSV file into a DataFrame
df = pd.read_csv('Ora/Overlapped/overlappedEnrichedGenesTrubetskoy.csv', header=None)

# Set the first row as the column names
df.columns = df.iloc[0]

# Drop the first row (column names) and reset the index
df = df.drop(0).reset_index(drop=True)

# Initialize an empty dictionary to store frequency counts
frequency_table = {}

# Iterate through the DataFrame
for column in df.columns:
    values = df[column].dropna().values
    for val in values:
        if ',' in str(val):
            for num in map(int, val.split(',')):
                frequency_table[num] = frequency_table.get(num, 0) + 1

# Print the frequency table
print("Number\tFrequency")
for num, freq in frequency_table.items():
    print(f"{num}\t{freq}")