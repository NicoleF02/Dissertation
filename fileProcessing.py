import pandas as pd

# Initialize lists
graphFile = []
trubetskoyBroad = []
trubetskoyBroadCross = []
result_list = []

file_path = "Files/Trubetskoy_2022_broad_coding.txt"
with open(file_path, 'r') as file:
    next(file)

    for line in file:
        try:
            value = int(line.split('\t')[0])
            trubetskoyBroad.append(value)
        except (ValueError, IndexError):
            pass


file_path = "Files/Trubetskoy_broad_db_cross.txt"
with open(file_path, 'r') as file:
    next(file)

    for line in file:
        try:
            value = int(line.split('\t')[0])
            trubetskoyBroadCross.append(value)
        except (ValueError, IndexError):
            pass


file_path = "PostsynapticNetwork/FullPSDDBNetwork.gml"
with open(file_path, 'r') as file:
    for line in file:
        if 'name' in line:
            try:
                value = int(line.split('"')[1])
                graphFile.append(value)
            except (ValueError, IndexError):
                pass


graphFile2 = []
file_path = "PostsynapticNetwork/PSD_FULL_clustered.gml"
with open(file_path, 'r') as file:
    for line in file:
        if 'name' in line:
            try:
                value = int(line.split('"')[1])
                graphFile2.append(value)
            except (ValueError, IndexError):
                pass

file_path = "Files/geneTable.txt"
geneTable = pd.read_csv(file_path, sep=' ')
geneTable = geneTable.dropna(subset=['HumanEntrez'])

entrezGeneTableList = [int(value) for value in geneTable["HumanEntrez"].tolist()]

missingFromGraph = [item for item in entrezGeneTableList if item not in graphFile]

temp = []
for node in graphFile:
    if node in entrezGeneTableList:
        temp.append(node)

temp2 = []

for value in entrezGeneTableList:
    if value in graphFile:
        temp2.append(value)

# Find the differences

# Print the differences

myGraphSet = set(graphFile)
otherGraphSet = set(graphFile2)

trubetskoyBroadSet = set(trubetskoyBroad)


print("No unique nodes of my graph ", len(myGraphSet))
print("No unique nodes of your graph ", len(otherGraphSet))
print("No. same nodes", len(myGraphSet.intersection(otherGraphSet)))
print("My graph trubetskoy matching ", len(myGraphSet.intersection(trubetskoyBroadSet)))
print("your graph trubetskoy matching ", len(otherGraphSet.intersection(trubetskoyBroadSet)))




#%%
