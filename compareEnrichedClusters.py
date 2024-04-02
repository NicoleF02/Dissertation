from sklearn.metrics import fowlkes_mallows_score
from itertools import combinations
import pandas as pd
import networkx as nx
import os
def get_enriched(filepath):
    # Read in all of the files in the folder and gets alg and clusternum
    enrichedClusters = []

    for filename in os.listdir(filepath):
        if filename.endswith(".csv"):
            file_path = os.path.join(filepath, filename)
            df = pd.read_csv(file_path)

            for index, row in df.iterrows():
                tupleAlgClNo = (row['alg'], row['clNo'])
                enrichedClusters.append(tupleAlgClNo)


    return enrichedClusters

def fowlkes_mallows_scr(algcl1, algcl2):

    # Calculate Rand Index
    rand_index_value = fowlkes_mallows_score(algcl1, algcl2)

    return rand_index_value

def compare_cluster_similarity(cluster1, cluster2):
    """
    Compare the Jaccard similarity between two clusters.
    """
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


file_path = "PostsynapticNetwork/FullPSDDBNetwork.gml"

graph = nx.read_gml(file_path, label='id')
df = pd.DataFrame.from_dict(graph.nodes, orient='index')
print(df.head())

algorithms = ['infomap', 'sgG1', 'sgG2', 'sgG5', 'spectral']

clusters = get_enriched("Ora/Enriched/algorithmSummary/reduced")

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

clustersDf.to_csv("Ora/Enriched/algorithmSummary/reduced/clusterSimilarity.csv")






# compare by cluster numbers using rand index

# Read in the enriched cluster numbers from Ora




# WE know all of the enriched clusters from Ora/,,, , we can get the geneoverlap from FullDBClustered and their cluster.
# We want to compare clustering accross different algorithms enrichment, maybe consider it especially clusters

# Basically just compare trub clusters to each other, rand index again




# Read in a csv from the Ora folder, we want to compare the enriched clusters from the different algorithms
