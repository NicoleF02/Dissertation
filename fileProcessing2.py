import networkx as nx
import pandas as pd

def read_gml_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    G = nx.parse_gml(data, label="name")  # Specify the attribute used as the 'name' category
    return G


def read_csv_file(csv_path):
    # Assuming that your CSV files have columns named 'alg', 'cl', 'pval', 'padj', etc.
    df = pd.read_csv(csv_path)
    return df

def find_max_value_for_category(G, category):
    return max((int(node_data[category]) for node_data in G.nodes.values() if category in node_data), default=None)


def find_no_trubetskoy(G, algorithm, community):
    totalTrubetskoy = 0
    #print(community)

    for node_id, node_data in G.nodes.items():
        if (node_data.get(algorithm) == str(community) and node_data.get("Trubetskoy2022broadcoding") == "TRUE"):
            totalTrubetskoy += 1

    return totalTrubetskoy

def main():
    file_path = "PostsynapticNetwork/FullPSDDBNetwork.gml"
    G = read_gml_file(file_path)
    algorithmsConsensus = ['lec', 'wt', 'fc', 'infomap', 'louvain','sgG1', 'sgG2', 'sgG5', 'spectral']
    algorithmsFull = ['wt', 'fc', 'infomap', 'louvain','sgG1', 'sgG2', 'sgG5', 'spectral']

    fullAlgoDict = {}
    algoTrubetskoy = {}

    for algorithm in algorithmsFull:
        numberCommunities = find_max_value_for_category(G, algorithm)

        fullAlgoDict[algorithm] = numberCommunities

        print(algorithm)
        print("Community, NoTrubetskoy, pval, padj")
        algoTrubetskoy[algorithm] = 0
        totalTrubetskoy = 0

        csv_path = f"FullDBClustered/{algorithm}.csv"
        df = read_csv_file(csv_path)



        for i in range(1, numberCommunities + 1):
            noTrubetskoy = find_no_trubetskoy(G,algorithm,i)

            matching_rows = df[(df['alg'] == algorithm) & (df['cl'] == i)]

            if not matching_rows.empty:
                row = matching_rows.iloc[0]
                pval = row['pval']
                padj = row['padj']
            else:
                pval = "N/A"
                padj = "N/A"


            print(f"{i}, {noTrubetskoy}, {pval}, {padj}")
            totalTrubetskoy += noTrubetskoy

        algoTrubetskoy[algorithm] += totalTrubetskoy



    print("Algorithms and communities")
    print("Algorithm, No. Communities, No. Trubetskoy")
    for algorithm in algorithmsFull:
        print(algorithm,",",fullAlgoDict[algorithm], ",",algoTrubetskoy[algorithm])



if __name__ == "__main__":
    main()
