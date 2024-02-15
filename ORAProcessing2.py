import os
import pandas as pd
import sqlite3
import pronto

ontologycsvDict = {
    "Trubetskoy2022broadcoding_significant_rows_sorted.csv": "Trubetskoy2022broadcoding",
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

    goDict = dict(get_ids("GO","GOID"))


    df = dataframe.loc[dataframe['alg'] == alg]

    clusterNumbers = df['cl'].unique()

    print(clusterNumbers)

    # iterating through cluster numbers
    for cl in clusterNumbers:
        reduceddf = df.loc[df['cl'] == cl]


        gobpids = []
        gomfids = []
        synapselocations = []
        trubGenes = 0
        dieases = []
        falseTrub = False

        for index2, row2 in reduceddf.iterrows():

            if row2['ontology_csv'] == "Trubetskoy2022broadcoding_significant_rows_sorted.csv":
                # Else it is a false enrichment
                if row2[('FL')] != False:
                    trubGenes = trubGenes + row2['Mu']
                else:
                    trubGenes = row2['Mu']
                    falseTrub = True

                if trubGenes != 0:
                    print("error poss")

            elif row2['ontology_csv'] == "GOBPID_significant_rows_sorted.csv":
                try:
                    description = goDict[row2["FL"]]
                    gobpids.append(description)
                except:
                    broken.append(row2["FL"])

            elif row2["ontology_csv"] == "GOMFID_significant_rows_sorted.csv":
                try:
                    description = goTerms[row2["FL"]]
                    gomfids.append(description.name)
                except:
                    broken.append(row2["FL"])

            elif row2["ontology_csv"] == "TopOntoOVGHDOID_significant_rows_sorted.csv":
                try:
                    description = diseaseDict[row2["FL"]]
                    dieases.append(description)
                except:
                    broken.append(row2["FL"])


        # don't count if no trub enrichment
        if trubGenes <= 0 and not falseTrub:
            continue



        dictonary = {
            'alg':row2['alg'],
            'clNo':row2['cl'],
            'clustsize':row2['Cn'],
            'no.Trubetskoy': trubGenes,
            'Molecular Function': str(gomfids),
            'Biological Function': str(gobpids),
            'Diseases': str(dieases),
            'False Trub Enrichment': falseTrub
        }
        complied_table.append(dictonary)

    return complied_table




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
    combined_df = generate_dataframe("Ora/Enriched")

    print("Start")


    ont1 = "Trubetskoy2022broadcoding_significant_rows_sorted.csv"

    otherOnts = ['GOBPID_significant_rows_sorted.csv', 'SynapseLocations_significant_rows_sorted.csv',
                 'syngo_significant_rows_sorted.csv', 'TopOntoOVGHOID_significant_rows_sorted.csv',
                 'GOMFID_significant_rows_sorted.csv']

    algs = ["spectral","sgG1","sgG2","sgG5","infomap"]

    print()

    for alg in algs:
        algTable = generate_alg_table(combined_df, alg)

        algdf = pd.DataFrame(algTable)

        algdf.to_csv(f"Ora/Enriched/overlapped2electricboogaloo/reduced{alg}.csv", index=False)








    # now we want make the table with a few things: alg

    # do it for all, then have enriched turn on or off after, shouldn't be too difficult
