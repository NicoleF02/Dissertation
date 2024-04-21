import pandas as pd
import os


def process_csv_file(file_path):
    # Read CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)

    max_cl = df['cl'].max()

    # Create a DataFrame with all integers from 1 to the maximum 'cl' value
    all_cl_values = pd.DataFrame({'cl': range(1, max_cl + 1)})

    # Merge the original DataFrame with the new DataFrame to fill in missing values
    df = pd.merge(all_cl_values, df, on='cl', how='left')

    # Replace NaN values in the 'OR' column with 'NA'
    df['OR'].fillna('NA', inplace=True)

    or_values = ", ".join(map(str, df['OR']))
    print(or_values)


algs = ['wt', 'fc', 'infomap', 'louvain',
        'sgG1', 'sgG2', 'sgG5', 'spectral']

for alg in algs:
    csv_file = f"FullDBClustered/{alg}.csv"
    print(f"Algorithm {alg}")
    process_csv_file(csv_file)
    print("\n")
