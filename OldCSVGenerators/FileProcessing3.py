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

    # Print the values for the 'OR' column in the desired format
    or_values = ", ".join(map(str, df['OR']))
    print(or_values)

# Specify the directory where your CSV files are located
csv_directory = "/path/to/your/csv/files"

algs = ['wt', 'fc', 'infomap', 'louvain',
        'sgG1', 'sgG2', 'sgG5', 'spectral']

# Process each CSV file
for alg in algs:
    csv_file = f"FullDBClustered/{alg}.csv"
    print(f"Algorithm {alg}")
    process_csv_file(csv_file)
    print("\n")

