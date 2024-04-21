# read in from Trubetskoy2022broadcoding_significant_rows_sorted.csv

# we want only positive fl,
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Entrez


def get_gene_info(entrez_id):
    Entrez.email = 's1958199@ed.ac.uk'
    handle = Entrez.esummary(db='gene', id=entrez_id)
    record = Entrez.read(handle)
    official_symbol = record['DocumentSummarySet']['DocumentSummary'][0]['Name']
    summary = record['DocumentSummarySet']['DocumentSummary'][0]['Summary']
    return official_symbol, summary


def gene_frequency_table(pivot_data):
    frequency_table = {}

    for column in pivot_data.columns:
        values = pivot_data[column].dropna().values
        for pivot_data[column] in values:
            frequency_table[column] = frequency_table.get(column, 0) + 1


def get_algorithms(pivot_data, gene):
    algorithms = []
    # Get the index of non NaN for pivot_data[gene]
    for i in pivot_data[gene].index[pivot_data[gene].notna()]:
        algorithms.append(i)

    return algorithms


def gen_latex_table(pivot_data):
    # We want to generate a latex table, we want a column for gene id, gene symbol, algorithms and summary
    # We want to get the gene symbol and summary from the gene id using get_gene_info

    latexPandasList = []

    for column in pivot_data.columns:
        gene_id = column
        gene_symbol, summary = get_gene_info(gene_id)

        dictGene = {'Entrez Gene ID': gene_id, 'Gene Symbol': gene_symbol,
                    'Algorithm & Cluster Associated': get_algorithms(pivot_data, gene_id), 'Summary': summary}
        latexPandasList.append(dictGene)

    sortedLatex = sorted(latexPandasList, key=lambda x: len(x['Algorithm & Cluster Associated']), reverse=True)

    latexPandas = pd.DataFrame(sortedLatex)
    print(latexPandas.to_latex(index=False))


def visualise(pivot_data, filename, reduced=True):
    # Disregard rows and columns with only one value
    if reduced:
        pivot_data = pivot_data.loc[:, pivot_data.nunique() > 1]
    pivot_data = pivot_data.dropna(axis=0, how='all')
    pivot_data = pivot_data.dropna(axis=1, how='all')

    print(pivot_data)

    plt.figure(figsize=(18, 10))
    heatmap = sns.heatmap(pivot_data, cmap="viridis", annot=False, fmt=".2f", xticklabels=True, yticklabels=True,
                          cbar=True, linecolor='grey', linewidth=0.5)

    plt.title('Heatmap of pval value for each enriched algorithm cluster\'s SCH genes', fontsize=16)
    plt.ylabel('Algorithm and Cluster Number', fontsize=14)
    plt.xlabel('SCH Enriched Genes', fontsize=14)

    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(fontsize=12)

    plt.savefig(filename)


csvs = ["Trubetskoy2022broadcoding_significant_rows_sorted.csv",
        "Trubetskoy2022priortisedcoding_significant_rows_sorted.csv"]
types = ["broad", "prioritised"]

networks = ["Consensus", "SynGO"]
filepaths = ["Consensus/", "SynGO/Ora/Enriched/"]
reduced = [True, False]

for i, network in enumerate(networks):
    for j, csv in enumerate(csvs):
        for reduce in reduced:
            print("Network: ", network, "Type: ", types[j], "Reduced: ", reduce)
            data = pd.read_csv(f"{filepaths[i]}{csv}")

            # Only positive enrichment
            data = data[data['FL'] == True]

            data['alg_cl'] = data['alg'] + '_' + data['cl'].astype(str)

            # Makes each gene into own row
            data = data.assign(overlapGenes=data['overlapGenes'].str.split(', ')).explode('overlapGenes')

            # Create a pivot table to prepare data for heatmap
            pivot_data = data.pivot_table(index='alg_cl', columns='overlapGenes', values='pval', aggfunc='mean')

            gen_latex_table(pivot_data)
            visualise(pivot_data, f"{network}{types[j]}CodingEnrichmentHeatmap{reduce}.png", reduced=reduce)
            print("-------------")