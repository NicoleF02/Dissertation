# We will be using this for rand index (RI) etc...
import networkx as nx
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from itertools import combinations
from multiprocessing import Pool, cpu_count


def jaccard(cluster1, cluster2):
    intersection = len(set(cluster1) & set(cluster2))
    union = len(set(cluster1) | set(cluster2))
    similarity = intersection / union if union != 0 else 0
    return similarity

def community_rename_sort(df, algo1, algo2):
    # This does jaccard index between clusters, and basically renames clusters with most simularity

    dfSimilarities = pd.DataFrame()

    # 1. get the clusters

    algo1Dict = get_clusters(df, algo1)
    algo2Dict = get_clusters(df, algo2)

    # 2. Compare every cluster using jaccard

    for cluster, genes in algo1Dict.items():
        for cluster2, genes2 in algo2Dict.items():
            jaccard_simularity = jaccard(set(genes), set(genes2))
            dfSimilarities["Alg1Cl"] = cluster
            dfSimilarities["Alg2Cl"] = cluster2
            dfSimilarities["Jaccard"] = jaccard_simularity

    # 3. Rank by simularity (most to least)

    dfSimilarities = dfSimilarities.sort_values(by='Jaccard', ascending=False)

    # We also want a list of all possible clusters, remove them from

    alg1ClusterList = df[algo1].unique()
    alg2ClusterList = df[algo2].unique()

    # 4. remove from both as we re-name them to arbitary ones in the dataframe, go by id then new name

    clusterNum = 0
    while not dfSimilarities.empty:
        most_simular = dfSimilarities.iloc[0]

        cluster1 = most_simular["Alg1Cl"]
        cluster2 = most_simular["Alg2Cl"]

        # Removes the rows with these clusters to avoid counting same one twice
        dfSimilarities = dfSimilarities[(dfSimilarities['Alg1Cl'] == cluster1) | (dfSimilarities['Alg2Cl'] == cluster2)]

        alg1ClusterList.remove(cluster1)
        alg2ClusterList.remove(cluster2)

        # This renames cluster values
        for gene in algo1Dict[cluster1]:
            df.loc[df.index == gene, algo1] = clusterNum

        for gene in algo1Dict[cluster2]:
            df.loc[df.index == gene, algo2] = clusterNum

        clusterNum = clusterNum + 1

    # Any left have no simularity at all, this sorts them
    if alg1ClusterList.any():
        for cluster in alg1ClusterList:
            for gene in algo1Dict[cluster]:
                df.loc[df.index == gene, algo1] = clusterNum
            clusterNum = clusterNum + 1
    elif alg2ClusterList.any():
        for cluster in alg2ClusterList:
            for gene in algo2Dict[cluster]:
                df.loc[df.index == gene, algo2] = clusterNum
            clusterNum = clusterNum + 1

    return df


def get_clusters(df, algo):
    clusters = {}

    for index, row in df.iterrows():
        if row[algo] not in clusters:

            clusters[row[algo]] = [index]
        else:
            clusters[row[algo]].append(index)

    return clusters


def calculate_rand_index(pair):
    # rand is symmetric so true and predicted don't matter.
    algo1, algo2, df = pair

    df = community_rename_sort(df, algo1, algo2)

    true_labels = df[algo1].tolist()
    predicted_labels = df[algo2].tolist()
    rand_index_value = adjusted_rand_score(true_labels, predicted_labels)
    return algo1, algo2, rand_index_value


def main():
    file_path = "PostsynapticNetwork/FullPSDDBNetwork.gml"
    graph = nx.read_gml(file_path, label='name')
    df = pd.DataFrame.from_dict(graph.nodes, orient='index')

    algorithms = ['infomap', 'sgG1', 'sgG2', 'sgG5', 'spectral', 'wt', 'louvain']

    matrix_size = len(algorithms)
    rand_index_matrix = pd.DataFrame(index=algorithms, columns=algorithms)
    pairs = [(algo1, algo2, df) for algo1, algo2 in combinations(algorithms, 2)]

    # # Multiprocess, pool will only use what you give it.
    # with Pool(cpu_count()) as pool:
    #     results = pool.map(calculate_rand_index, pairs)
    results = []

    for pair in pairs:
        results.append(calculate_rand_index(pair))


    for algo1, algo2, rand_index_value in results:
        rand_index_matrix.at[algo1, algo2] = rand_index_value
        rand_index_matrix.at[algo2, algo1] = rand_index_value

        print(f"Adjusted Rand Index between {algo1} and {algo2}: {rand_index_value}")

    print("Adjusted Rand Index Matrix:")
    print(rand_index_matrix)

    rand_index_matrix.to_csv('rand_index_matrix_full.csv')


if __name__ == "__main__":
    main()

# Using manuel inspection
