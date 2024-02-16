# We will be using this for rand index (RI) etc...
import networkx as nx
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from itertools import combinations



if __name__ == "__main__":
    file_path = "PostsynapticNetwork/SynGoNetwork.gml"
    graph = nx.read_gml(file_path, label = 'id')
    df = pd.DataFrame.from_dict(graph.nodes, orient='index')

    algorithms = ['infomap', 'sgG1', 'sgG2', 'sgG5', 'spectral']

    matrix_size = len(algorithms)
    rand_index_matrix = pd.DataFrame(index=algorithms, columns=algorithms)



    for algo1, algo2 in combinations(algorithms, 2):
        true_labels = df[algo1].tolist()
        predicted_labels = df[algo2].tolist()

        # Calculate Rand Index
        rand_index_value = adjusted_rand_score(true_labels, predicted_labels)
        # Fill the matrix symmetrically
        rand_index_matrix.at[algo1, algo2] = rand_index_value
        rand_index_matrix.at[algo2, algo1] = rand_index_value

        print(f"Adjusted Rand Index between {algo1} and {algo2}: {rand_index_value}")

    print("Adjusted Rand Index Matrix:")
    print(rand_index_matrix)

    rand_index_matrix.to_csv('SynGO/rand_index_matrix.csv')
    print("Matrix saved as rand_index_matrix.csv")



# Using manuel inspection