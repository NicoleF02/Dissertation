import pandas as pd

df = pd.read_csv('Ora/Overlapped/overlappedEnrichedGenesTrubetskoy.csv', header=None)

df.columns = df.iloc[0]

df = df.drop(0).reset_index(drop=True)

frequency_table = {}

for column in df.columns:
    values = df[column].dropna().values
    for val in values:
        if ',' in str(val):
            for num in map(int, val.split(',')):
                frequency_table[num] = frequency_table.get(num, 0) + 1

print("Number\tFrequency")
for num, freq in frequency_table.items():
    print(f"{num}\t{freq}")