library(BioNAR)

gFull <- igraph::read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml",format="gml") #graph from gml

algs<-c('infomap','sgG1', 'sgG2', 'sgG5', 'spectral','louvain','fc','wt')

list_ontologies <- c("Trubetskoy2022broadcoding","Trubetskoy2022priortisedcoding","SynapseLocations","GOMFID","syngo","TopOntoOVGHDOID","GOBPID")

runORA <- function(ontology){
  print("Starting for ontology:")
  print(ontology)
  dir.create(output_folder, showWarnings = FALSE)
  print(1)

  ora<-lapply(algs, function(alg){clusterORA(gFull, alg, name = ontology,
                                             vid = "name",alpha = 1, col = COLLAPSE)})
  names(ora)<-algs
  print("Ora generated")
  return(ora)

}

# Loop through the ontologies and perform an ORA analysis
for (ontology in list_ontologies) {
  output_folder <- ""
  output_folder <- paste0("Ora/",ontology)
  print(output_folder)

  ora <- runORA(ontology)
  # Save the results to a CSV file
  k <- 1
  for (j in seq_along(ora)) {

    column_name <- paste0("column", j)
    file_name <- paste0(output_folder, "/", ontology, "_", algs[k], ".csv")
    write.csv(ora[[j]], file = file_name, row.names = FALSE)

    k <- k + 1
  }
  print("Done with one")
}
