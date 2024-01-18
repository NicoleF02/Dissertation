import networkx as nx


def read_gml_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    G = nx.parse_gml(data, label="name")  # Specify the attribute used as the 'name' category
    return G


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
        print("Community, NoTrubetskoy")
        algoTrubetskoy[algorithm] = 0
        totalTrubetskoy = 0
        for i in range(1, numberCommunities):
            noTrubetskoy = find_no_trubetskoy(G,algorithm,i)
            if noTrubetskoy != 0:
                print(f"{i}, {noTrubetskoy}")
            totalTrubetskoy += noTrubetskoy

        algoTrubetskoy[algorithm] += totalTrubetskoy



    print("Algorithms and communities")
    print("Algorithm, No. Communities, No. Trubetskoy")
    for algorithm in algorithmsFull:
        print(algorithm,",",fullAlgoDict[algorithm], ",",algoTrubetskoy[algorithm])



if __name__ == "__main__":
    main()
