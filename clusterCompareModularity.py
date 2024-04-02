import networkx as nx
import pandas as pd
from itertools import combinations

file_path = "PostsynapticNetwork/FullPSDDBNetwork.gml"

graph = nx.read_gml(file_path, label='name')

df = pd.DataFrame.from_dict(graph.nodes, orient='index')



# Now we want to rank modularity using nx

# Add edges based on the clustering results

algorithms = ['infomap', 'sgG1', 'sgG2', 'sgG5', 'spectral','louvain','wt']


for algo in algorithms:
    clustersDict = {}
    # Makes sure clear from prev run

    # Gets clusters in dictionary for modularity score

    for alg_value in df[algo].unique():
        indices = df.index[df[algo] == alg_value].tolist()
        clustersDict[str(alg_value)] = indices


    modScore = nx.community.modularity(graph, list(clustersDict.values()))

    print("Modularity score for ", algo, " is ", modScore)
    clustersDict.clear()













