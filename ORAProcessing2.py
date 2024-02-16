import os
import pandas as pd
import sqlite3
import pronto
import csv

ontologycsvDict = {
    "Trubetskoy2022broadcoding_significant_rows_sorted.csv": "Trubetskoy2022broadcoding",
    "Trubetskoy2022priortisedcoding_significant_rows_sorted.csv": "Trubetskoy2022priortisedcoding",
    "GOBPID_significant_rows_sorted.csv": "GOBPID",
    "SynapseLocations_significant_rows_sorted.csv": "SynapseLocations",
    "syngo_significant_rows_sorted.csv": "SynGO",
    "TopOntoOVGHDOID_significant_rows_sorted.csv": "TopOntoOVGHOID",
    "GOMFID_significant_rows_sorted.csv": "GOMFID"
}

goTerms = pronto.Ontology('go-basic.obo')


def generate_dataframe(path):
    folder_path = path
    all_csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

    dfs = []
    for csv_file in all_csv_files:
        df = pd.read_csv(os.path.join(folder_path, csv_file))
        df['ontology_csv'] = csv_file
        dfs.append(df)

    # Concatenate all dataframes into one
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df


def generate_table(dataframe, ontology1, ontology2):
    # We want to go through every ontology with another, we don't want to compare those in the same csv

    df1 = dataframe.loc[dataframe['ontology_csv'] == ontology1]
    df2 = dataframe.loc[dataframe['ontology_csv'] == ontology2]

    # We want to compare where alg and cl are the same, compare geneoverlap

    cols = ["alg",
            "cluster",
            "Ontology1",
            "FL1",
            "Ontology2",
            "FL2",
            "OverlappedGenes",
            "NoOverlappedGenes"]

    complied_table = []

    for index, row in df1.iterrows():
        result = df2[(df2['alg'] == row['alg']) & (df2['cl'] == row['cl'])]
        # We then want to extract number of gene overlap
        if result.empty:
            continue

        genes1 = set(row['overlapGenes'].split(','))

        for index2, row2 in result.iterrows():
            genes2 = set(row2['overlapGenes'].split(','))

            matchingGenes = genes1.intersection(genes2)

            if len(matchingGenes) == 0:
                continue

            # Now we want to add in a thing, create a dict for what we want

            dict = {
                "alg": row['alg'],
                "cluster": row['cl'],
                "Ontology1": ontologycsvDict[row['ontology_csv']],
                "FL1": row['FL'],
                "Ontology2": ontologycsvDict[row2['ontology_csv']],
                "FL2": row2['FL'],
                "OverlappedGenes": str(matchingGenes),
                "NoOverlappedGenes": len(matchingGenes)
            }

            complied_table.append(dict)
    return complied_table


def generate_alg_table(dataframe, alg):
    # We want to go through dataframe and find stuff, firstly we get all algorithm and then sort by cluster.
    # We want cluster size, no. trub etc.. I will run twice, once for enriched and once w/o

    complied_table = []
    # We want to connect to the database to get Disease HDOID -> Description
    # GOID -> Description

    colnames = ['alg', 'clNo', 'clustSize', 'no. Trubetskoy', 'GO Terms']

    diseaseDict = dict(get_ids("Disease", "HDOID"))

    broken = []

    goDict = dict(get_ids("GO", "GOID"))

    df = dataframe.loc[dataframe['alg'] == alg]

    clusterNumbers = df['cl'].unique()

    print(clusterNumbers)

    # iterating through cluster numbers
    for cl in clusterNumbers:
        reduceddf = df.loc[df['cl'] == cl]

        gobpids = []
        gobpidadj = []

        gomfids = []
        gomfidspadj = []

        synapselocations = []

        trubGenes = 0
        trubPadj = 0
        trubPval = 0

        dieases = []
        diseasespadj = []

        falseTrub = False
        trubGenesPrior = 0
        trubPriorpadj = 0
        trubPriorPval = 0

        for index2, row2 in reduceddf.iterrows():

            if row2['ontology_csv'] == "Trubetskoy2022broadcoding_significant_rows_sorted.csv":
                # Else it is a false enrichment
                if row2[('FL')] != False:
                    if trubGenes != 0:
                        print("error poss")

                    trubGenes = trubGenes + row2['Mu']
                    trubPval = row2['pval']
                    trubPadj = row2['padj']

                else:
                    trubGenes = row2['Mu']
                    falseTrub = True
                    trubPval = row2['pval']
                    trubPadj = row2['padj']



            elif row2['ontology_csv'] == "GOBPID_significant_rows_sorted.csv":
                try:
                    if row2['padj'] <= 0.05:
                        description = goTerms[row2["FL"]]
                        gobpids.append(description.name)
                        gobpidadj.append(row2['padj'])

                except:
                    broken.append(row2["FL"])

            elif row2["ontology_csv"] == "GOMFID_significant_rows_sorted.csv":
                try:
                    if row2['padj'] <= 0.05:
                        description = goTerms[row2["FL"]]
                        gomfids.append(description.name)
                        gomfidspadj.append(row2['padj'])
                except:
                    broken.append(row2["FL"])

            elif row2["ontology_csv"] == "TopOntoOVGHDOID_significant_rows_sorted.csv":
                try:
                    if row2['padj'] <= 0.05:
                        description = diseaseDict[row2["FL"]]
                        dieases.append(description)
                        diseasespadj.append(row2['padj'])
                except:
                    broken.append(row2["FL"])

            elif row2['ontology_csv'] == "Trubetskoy2022priortisedcoding_significant_rows_sorted.csv":
                trubGenesPrior = row2['Mu']
                trubPriorPval = row2['pval']
                trubPriorpadj = row2['padj']

        # don't count if no trub enrichment
        if (trubGenes == 0 and trubGenesPrior == 0) and not (falseTrub):
            continue

        if falseTrub:
            enrichment = "Negative"
        else:
            enrichment = "Positive"

        dictonary = {
            'alg': row2['alg'],
            'clNo': row2['cl'],
            'clustsize': row2['Cn'],
            'no.Trubetskoy Broad': trubGenes,
            'Enrichment': enrichment,
            'Broad pval': trubPval,
            'Broad padj': trubPadj,
            'no. Trubetskoy Priortised': trubGenesPrior,
            'Prior pval': trubPriorPval,
            'Prior padj': trubPriorpadj,
            'Molecular Function': gomfids,
            'Molecular padj': gomfidspadj,
            'Biological Function': gobpids,
            'Biological padj': gomfidspadj,
            'Diseases': dieases,
            'Diseases padj': diseasespadj,

        }

        complied_table.append(dictonary)

    return complied_table


def expand_lists(row):
    expanded_rows = [row[:]]  # Create a copy of the original row
    for i, item in enumerate(row):
        if isinstance(item, list):
            new_rows = []
            for existing_row in expanded_rows:
                for subitem in item:
                    new_row = existing_row.copy()
                    new_row[i] = subitem
                    new_rows.append(new_row)
            expanded_rows = new_rows
    return expanded_rows

def fix_csv(input_file, output_file):
    with open(input_file, 'r') as csv_file, open(output_file, 'w') as output:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)
        output.write(','.join(header) + '\n')

        for row in csv_reader:
            expanded_rows = expand_lists(row)
            for expanded_row in expanded_rows:
                output.write(','.join(map(str, expanded_row)) + '\n')


def get_ids(table, primary_key):
    connection = sqlite3.connect("Database/DS_10283_3877/synaptic.proteome_SR_20210408.db.sqlite")
    cursor = connection.cursor()

    # This gets the list of IDs and descriptions from db
    sql_query = f"SELECT {primary_key}, description FROM {table}"

    cursor.execute(sql_query)

    result = cursor.fetchall()

    connection.close()

    return result


if __name__ == "__main__":
    path = "Ora/Enriched"

    combined_df = generate_dataframe(path)

    print("Start")

    ont1 = "Trubetskoy2022broadcoding_significant_rows_sorted.csv"

    otherOnts = ['GOBPID_significant_rows_sorted.csv', 'SynapseLocations_significant_rows_sorted.csv',
                 'syngo_significant_rows_sorted.csv', 'TopOntoOVGHOID_significant_rows_sorted.csv',
                 'GOMFID_significant_rows_sorted.csv', "Trubetskoy2022priortisedcoding_significant_rows_sorted.csv"]

    algs = ["spectral", "sgG1", "sgG2", "sgG5", "infomap"]

    print()

    for alg in algs:
        print(alg)
        fix_csv(f"{path}/algorithmSummary/enriched{alg}.csv", f"{path}/algorithmSummary/fixedEnriched{alg}.csv")

        # print(alg)
        # algTable = generate_alg_table(combined_df, alg)
        #
        # algdf = pd.DataFrame(algTable)
        #
        # columns_to_explode = ['alg', 'clNo', 'clustsize', 'no.Trubetskoy Broad', 'Enrichment', 'Broad pval', 'Broad padj', 'no. Trubetskoy Priortised', 'Prior pval', 'Prior padj', 'Molecular Function', 'Molecular padj', 'Biological Function', 'Biological padj', 'Diseases', 'Diseases padj']
        # df_expanded = algdf.copy()

        # for col in columns_to_explode:
        #     print(col)
        #     df_expanded[col] = df_expanded[col].apply(lambda x: [x] if isinstance(x, str) else x)
        #     df_expanded = df_expanded.explode(col)
        #
        #
        #     # Replace repeated values with blank for the exploded column
        #     df_expanded[col] = df_expanded[col].mask(df_expanded[col].duplicated(), '')

        # df_expanded.to_csv(f"{path}/algorithmSummary/enriched{alg}.csv", index=False)

    # now we want make the table with a few things: alg

    # do it for all, then have enriched turn on or off after, shouldn't be too difficult

# Check enriched cluster simularity indivudually
# Check biological terms correlations between clusters
# see whats happen when you add schizophrenia + genetic data together

# 1. fix spreadsheets
# checking enrichment cluster simularity
# check biology term correlation
# consider adding schizophrenia to genetic data -- extra time only
# come to my own conclusion
# then compare

# Look at trubetskoy paper again which conclusions they got from syngo etc. check their conlcusions, see if they
# conclude something that I can't find,
