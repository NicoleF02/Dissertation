import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles

from tabulate import tabulate


# We want to do a seperate csv for each ontology to allow for easier comparision, we can combine by algorithm though.
def read_csv_files(folder_path):
    combined_df = pd.DataFrame()
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path)
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df


def saveORAEnrichment(dataframe, ontology):
    significant_rows = dataframe[dataframe['padj'] <= 0.05 or dataframe["pval"] <= 0.05]

    # Specify columns to print
    columns_to_print = ['alg', 'cl', 'Cn', 'Mu', 'FL', 'padj', 'overlapGenes']

    # Sort the DataFrame by 'alg' and 'padj'
    df = significant_rows_sorted = significant_rows.sort_values(by=['alg', 'padj'])

    significant_rows_sorted[columns_to_print].to_csv(f'Ora/Enriched/{ontology}_significant_rows_sorted.csv',
                                                     index=False)

    return df


def returnGeneSet(df, algorithm = None):
    gene_set = set()
    if not (algorithm == None):
        df = df[df["alg"] == algorithm].copy()

    for _, row in df.iterrows():
        genes_list = [gene.strip() for gene in row['overlapGenes'].split(',')]
        gene_set.update(genes_list)

    return gene_set



def generateVenn(df1, df2 = None, algs = [], ontologyName = "Trubetskoy2022broadcoding"):
    print(df1)
    if df2 == None:
        noAlgs = len(algs)
        print(algs)
        geneAlg1 = returnGeneSet(df1, algs[0])
        geneAlg2 = returnGeneSet(df1, algs[1])
        overlap12 = len(geneAlg1.intersection(geneAlg2))
        plt.title(ontologyName)
        if noAlgs == 3:
            geneAlg3 = returnGeneSet(df1, algs[2])
            overlap13 = len(geneAlg1.intersection(geneAlg3))
            overlap23 = len(geneAlg2.intersection(geneAlg3))
            overlap123 = len((geneAlg2.intersection(geneAlg3)).intersection(geneAlg1))

            # (Abc, aBc, ABc, abC, AbC, aBC, ABC)
            venn3(subsets=[len(geneAlg1),len(geneAlg2),overlap12,len(geneAlg3),overlap13, overlap23, overlap123], set_labels=algs)
            plt.savefig(f"Ora/Enriched/{ontologyName}{algs[0]}{algs[1]}{algs[2]}.png")
        else:
            # (Ab, aB, AB) for two
            venn2(subsets=[len(geneAlg1),len(geneAlg2),overlap12], set_labels=algs)
            plt.savefig(f"Ora/Enriched/{ontologyName}{algs[0]}{algs[1]}.png")
        plt.show()
    else:
        noAlgs = len(algs)
        print(algs)
        geneAlg1 = returnGeneSet(df1, algs[0])
        geneAlg2 = returnGeneSet(df2, algs[1])
        overlap12 = len(geneAlg1.intersection(geneAlg2))
        plt.title(ontologyName)
        if noAlgs == 3:
            geneAlg3 = returnGeneSet(df1, algs[2])
            overlap13 = len(geneAlg1.intersection(geneAlg3))
            overlap23 = len(geneAlg2.intersection(geneAlg3))
            overlap123 = len((geneAlg2.intersection(geneAlg3)).intersection(geneAlg1))

            # (Abc, aBc, ABc, abC, AbC, aBC, ABC)
            venn3(subsets=[len(geneAlg1),len(geneAlg2),overlap12,len(geneAlg3),overlap13, overlap23, overlap123], set_labels=algs)
            plt.savefig(f"Ora/Enriched/{ontologyName}{algs[0]}{algs[1]}{algs[2]}.png")
        else:
            # (Ab, aB, AB) for two
            venn2(subsets=[len(geneAlg1),len(geneAlg2),overlap12], set_labels=algs)
            plt.savefig(f"Ora/Enriched/{ontologyName}{algs[0]}{algs[1]}.png")
        plt.show()


if __name__ == "__main__":
    list_ontologies = ["Trubetskoy2022broadcoding", "SynapseLocations", "GOMFID", "syngo", "TopOntoOVGHDOID", "GOBPID"]

    trubetskoydf = read_csv_files("Ora/" + list_ontologies[0])
    synapseLocationdf = read_csv_files("Ora/" + list_ontologies[1])
    gomfiddf = read_csv_files("Ora/" + list_ontologies[2])
    syngodf = read_csv_files("Ora/" + list_ontologies[3])
    topontoovghdoiddf = read_csv_files("Ora/" + list_ontologies[4])
    gobpid = read_csv_files("Ora/" + list_ontologies[5])

    trubetskoyEnriched = saveORAEnrichment(trubetskoydf, "Trubetskoy2022broadcoding")
    synapseEnriched = saveORAEnrichment(synapseLocationdf,list_ontologies[1])
    gomfidEnriched = saveORAEnrichment(gomfiddf,list_ontologies[2])
    syngoEnriched = saveORAEnrichment(syngodf,list_ontologies[3])
    topontoovgdoidEnriched = saveORAEnrichment(topontoovghdoiddf,list_ontologies[4])
    gobpidEnriched = saveORAEnrichment(gobpid,list_ontologies[5])

    algs = ["infomap","wt","louvian"]
    generateVenn(df1 = trubetskoyEnriched, df2= None, algs=algs)


# Look at what are associated and not, consider the false enrichment as well, generate, generate graph for trubetskoy in other ontologies,
# Bridgeness do best two, justify why is best algorithm, combine with sigmoid functions and compare sigmoid and with manual selection
# Flag up issues I've identified.

# Compare infomap and spectral communties, infomap as bigger one and the smaller ones for spectral, upset diagram