from sklearn.metrics import fowlkes_mallows_score
from itertools import combinations
import pandas as pd
import networkx as nx
import os
def get_enriched(filepath):
    # Read in all of the files in the folder and gets alg and clusternum
    enrichedClusters = []
    fulldf = pd.DataFrame()

    for filename in os.listdir(filepath):
        if filename.endswith(".csv") and filename.startswith("reduced"):
            file_path = os.path.join(filepath, filename)
            df = pd.read_csv(file_path)
            fulldf = pd.concat([fulldf, df])

            for index, row in df.iterrows():
                tupleAlgClNo = (row['alg'], row['clNo'])
                enrichedClusters.append(tupleAlgClNo)


    return enrichedClusters, fulldf

def fowlkes_mallows_scr(algcl1, algcl2):

    # Calculate Rand Index
    rand_index_value = fowlkes_mallows_score(algcl1, algcl2)

    return rand_index_value

def compare_cluster_similarity(cluster1, cluster2):
    # Jaccard Similarity
    intersection = len(set(cluster1) & set(cluster2))
    union = len(set(cluster1) | set(cluster2))
    similarity = intersection / union if union != 0 else 0
    return similarity

def list_alg_cl(alg, clNo, df):
    # Get the cluster from the graph
    df[alg] = pd.to_numeric(df[alg], errors='coerce')
    clusterdf = df[df[alg] == int(clNo)]
    cluster = clusterdf["name"].tolist()
    return cluster


def print_latex(fulltabledf, jaccardDf):
    start_table = r"""\begin{table}[h]
    \begin{tabularx}{\linewidth}{|XXX|XXX|c|} \hline
    \textbf{Alg, Cluster No.} &
      \textbf{Cluster Size} &
      \textbf{No. SCZ} &
      \textbf{Alg, Cluster No.} &
      \textbf{Cluster Size} &
      \textbf{No. SCZ} &
      \textbf{Jaccard Index} \\ \hline"""

    end_table =r"""\end{tabularx} 
    \caption{Table to compare enriched clusters between different algorithms for the SynGO Synaptic Graph. This compares their similarity using the Jaccard Index, explained in \ref{Jaccard}. Only those with a Jaccard Index greater than 0 were included.}
    \label{table:enriched_cluster_comparison_syngo}
    \end{table}"""

    latexRows = ""

    jaccardDf = jaccardDf.sort_values(by='Jaccard', ascending=False)
    for index, row in jaccardDf.iterrows():
        #get the algorithm and cluster number from the Cluster1 and Cluster2 cols
        alg1 = row['Cluster1'][0]
        cl1 = row['Cluster1'][1]

        alg2 = row['Cluster2'][0]
        cl2 = row['Cluster2'][1]
        jaccard = row['Jaccard']


        alg1SczGenes = fulltabledf[(fulltabledf['alg'] == str(alg1)) & (fulltabledf['clNo'] == int(cl1))]['no.Trubetskoy Broad']
        alg2SczGenes = fulltabledf[(fulltabledf['alg'] == str(alg2)) & (fulltabledf['clNo'] == int(cl2))]['no.Trubetskoy Broad']
        cl1Size = fulltabledf[(fulltabledf['alg'] == str(alg1)) & (fulltabledf['clNo'] == int(cl1))]['clustsize'].iloc[0]
        cl2Size = fulltabledf[(fulltabledf['alg'] == str(alg2)) & (fulltabledf['clNo'] == int(cl2))]['clustsize'].iloc[0]

        if alg1SczGenes.empty:
            print(f"Empty: {alg1} {cl1}")
        else:
            alg1SczGenes = alg1SczGenes.iloc[0]

        if alg2SczGenes.empty:
            print(f"Empty: {alg2} {cl2}")
        else:
            alg2SczGenes = alg2SczGenes.iloc[0]

        print(fulltabledf.head(10))




        #latex table inner
        latexrow = f"{alg1} {cl1} & {cl1Size} & {alg1SczGenes} & {alg2} {cl2} & {alg2SczGenes} & {cl2Size} & {jaccard} \\\\ \n"

        latexRows += latexrow

    print(start_table)
    print(latexRows)
    print(end_table)









file_path = "PostsynapticNetwork/SynGoNetwork.gml"

graph = nx.read_gml(file_path, label='id')
df = pd.DataFrame.from_dict(graph.nodes, orient='index')


algorithms = ['infomap', 'sgG1', 'sgG2', 'sgG5', 'spectral']

clusters, fulldf = get_enriched("SynGO/Ora/Enriched/algorithmSummary/reduced/")

clusterDf = pd.DataFrame()

clustersList = {}

for cluster in clusters:
    # We have tuple here of (alg, clNo), this be be used as a map to generate the lists
    listNodes = list_alg_cl(cluster[0], cluster[1], df)
    clustersList[cluster] = listNodes



clustersDf = []
for cluster1, cluster2 in combinations(clustersList.keys(), 2):
    jaccard = compare_cluster_similarity(clustersList[cluster1], clustersList[cluster2])
    print(f"Comparing {cluster1} and {cluster2}")
    print(f"Jaccard Similarity: {jaccard}")

    recordCluster = {"Cluster1": cluster1, "Cluster2": cluster2, "Jaccard": jaccard}
    clustersDf.append(recordCluster)

clustersDf = pd.DataFrame(clustersDf)
clustersDf.to_csv("SynGO/Ora/Enriched/algorithmSummary/reduced/clusterSimilarity.csv")



print("rows of fulldf", fulldf.shape)

print_latex(fulldf, clustersDf)





# compare by cluster numbers using rand index

# Read in the enriched cluster numbers from Ora




# WE know all of the enriched clusters from Ora/,,, , we can get the geneoverlap from FullDBClustered and their cluster.
# We want to compare clustering accross different algorithms enrichment, maybe consider it especially clusters

# Basically just compare trub clusters to each other, rand index again




# Read in a csv from the Ora folder, we want to compare the enriched clusters from the different algorithms
