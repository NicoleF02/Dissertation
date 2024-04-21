import pandas as pd

graphFile = []
trubetskoyBroad = []
trubetskoyBroadCross = []
result_list = []


def readTrubetskoy(file_path):
    trubetskoyList = []
    with open(file_path, 'r') as file:
        next(file)

        for line in file:
            try:
                value = int(line.split('\t')[0])
                trubetskoyList.append(value)
            except (ValueError, IndexError):
                pass

    return trubetskoyList


def readGraph(file_path):
    graphNodes = []
    with open(file_path, 'r') as file:
        for line in file:
            if 'name' in line:
                try:
                    value = int(line.split('"')[1])
                    graphNodes.append(value)
                except (ValueError, IndexError):
                    pass

    return graphNodes


graphFile2 = []
file_path = "../PostsynapticNetwork/Useless/PSD_FULL_clustered.gml"
fullGraphNodesWorking = readGraph(file_path)

file_path = "../Files/Trubetskoy_2022_broad_coding.txt"
trubetskoyBroadList = readTrubetskoy(file_path)

file_path = "../Files/Trubetskoy_broad_db_cross.txt"
trubetskoyBroadCrossList = readTrubetskoy(file_path)

file_path = "../PostsynapticNetwork/FullPSDDBNetwork.gml"
fullGraphNodes = readGraph(file_path)


file_path = "../Files/geneTable.txt"
geneTable = pd.read_csv(file_path, sep=' ')
geneTable = geneTable.dropna(subset=['HumanEntrez'])

entrezGeneTableList = [int(value) for value in geneTable["HumanEntrez"].tolist()]


fullDB = pd.read_csv("../Files/Full_DB_Rat_Aut22.txt", sep="\t")


print("Set of nodes for broken network ", len(set(fullGraphNodes)))
print("Set of nodes for unbroken network ", len(set(fullGraphNodesWorking)))

print("Intersection of sets ", len(set(fullGraphNodes).intersection(set(fullGraphNodesWorking))))

uniqueFull = sorted(list(set(fullGraphNodes) - set(fullGraphNodesWorking)))
uniqueFullWorking = sorted(list(set(fullGraphNodesWorking) - set(fullGraphNodes)))




uniqueFullDFCross = fullDB[fullDB["HUMAN.ENTREZ.ID"].isin(uniqueFullWorking)]
print(len(uniqueFullDFCross))



