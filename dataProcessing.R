
graph <- read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml", format = "gml")



# Get the node attributes
node_attrs <- data.frame(
  trubetskoy = as.logical(V(graph)$Trubetskoy2022broadcoding),
  louvain = as.numeric(V(graph)$louvain)
)

# Create a table counting the occurrences of trubetskoy true for each louvian value
result_table <- table(node_attrs$louvain[node_attrs$Trubetskoy2022broadcoding])

# Print the result_table
print(result_table)