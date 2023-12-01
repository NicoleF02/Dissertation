library(knitr)
library(BioNAR)
library(synaptome.db)
library(ggplot2)
library(pander)
library(ggrepel)
library(randomcoloR)
library(dplyr)
require(BioNAR)
library(org.Hs.eg.db)
library(randomcoloR)
# We need to split this up
# No. 1. get all schizophrenia genes
# get them for postsynaptic
# get them for only postsynaptic only
# then do networks
# do refined network

initGeneTable <- function(count=2, localisation="Postsynaptic", diseasehdoid="DOID:5419"){
  narrowGeneList <- findGeneByCompartmentPaperCnt(cnt = count)
  tableSynaptic <- filter(narrowGeneList, Localisation == localisation)
  # then we want to get all genediease, then filter. Maybe check own db rather then their own?

  # if its been identified 2 or more times, consider maybe specifiying diff papers?
  if (!is.null(diseasehdoid)){
    columnGenes <- tableSynaptic$HumanEntrez
    tableDisease <- getGeneDiseaseByEntres(columnGenes)
    tableNarrowDisease <- filter(tableDisease, HDOID == diseasehdoid)
    return(tableNarrowDisease)
  }
  return(tableSynaptic)
}

networkSetup <- function(PPIList){
  network <- buildNetwork(PPIList)
  network <- annotateGeneNames(network)
  network <- annotateGOont(network)

  afile<-system.file("extdata", "flatfile_human_gene2HDO.csv", package = "BioNAR")
  dis <- read.table(afile,sep="\t",skip=1,header=FALSE,strip.white=TRUE,quote="")
  network <- annotateTopOntoOVG(network, dis)

  sg <- read.table("Files/SynGO.txt", sep = "\t", header = T, stringsAsFactors = F)
  network <- annotateVertex(network, name = "syngo", values = sg, idatt = "name")

  trubetskoyBroad<- read.table("Files/Trubetskoy_2022_broad_coding.txt", sep = "\t", header = T, stringsAsFactors = F)
  network <- annotateVertex(network, name="Trubetskoy_2022_broad_coding", values = trubetskoyBroad, idatt = "name")

  trubetskoyPriortised<- read.table("Files/Trubetskoy_2022_priortised_coding.txt", sep = "\t", header = T, stringsAsFactors = F)
  network <- annotateVertex(network, name="Trubetskoy_2022_priortised_coding", values = trubetskoyPriortised, idatt = "name")

  reducedRatDB <- read.table("Files/ReducedRatDB.txt", sep="\t", header= T, stringsAsFactors = F)
  network <- annotateVertex(network, name="Synapse Locations", values =reducedRatDB , idatt="name")

  priortisedSchizophreniaDB <- read.table("Files/SchizphreniaPriortisedDB.txt", sep="\t", header= T, stringsAsFactors = F)
  network <- annotateVertex(network, name="SchizophreniaGeneCount", values =priortisedSchizophreniaDB , idatt="name")

  return(network)
}

generateGraph <- function(count=2, localisation="Postsynaptic", diseasehdoid="DOID:5419", filename){
  geneTable <- initGeneTable(count, localisation, diseasehdoid)
  ppiGene <- getPPIbyEntrez(geneTable$HumanEntrez, type = "limited")
  networkGene <- networkSetup(ppiGene)
  write_graph(networkGene, file=filename, format="gml")
  return (networkGene)
}

generateGraphNew <- function(count=1, localisation="Postsynaptic", diseasehdoid=NULL, filename){
  geneTable <- initGeneTable(count, localisation, diseasehdoid)
  ppiGene <- getPPIbyEntrez(geneTable$HumanEntrez, type = "induced")
  networkGene <- networkSetup(ppiGene)
  write_graph(networkGene, file=filename, format="gml")
}

networkProperties <- function(network){
  # Nsim will make slower, 100 is default, 1000 will take a few mins
  pFit <- fitDegree(as.vector(igraph::degree(graph=g)), Nsim=100, plot=TRUE, WIDTH=2480, HEIGHT=2480)
  pwr <- slot(pFit, 'alpha')

}

compareNetwork <- function(network){


  lrgnp<-list()
  alphaGNP<-c()
  for(i in 1:5000){
    rgnp<-BioNAR:::getGNP(network)
    pFit <- fitDegree( as.vector(igraph::degree(graph=rgnp)),
                       Nsim=10, plot=FALSE,threads=5,
                       p <- slot(pFit,'alpha'),
                       lrgnp[[i]]<-rgnp,
                       alphaGNP[i]<-p)
  }
  qplot(alphaGNP)+geom_vline(xintercept = pwr)


}



# Postsynaptic Graphs

# # Graph for Schziphrenia Genes that appear >= 2 that are Postsynaptic
# generateGraph(file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml")
#
# # Graph for all Schziphrenia Genes that appear that are Postsynaptic
# generateGraph(count=1, file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml")
#
# # Graph for Genes that appear >= 2 that are Postsynaptic
consensusGraph <- generateGraph(count = 2,diseasehdoid = NULL, filename="PostsynapticNetwork/ConsensusPSDDBNetwork.gml")
#
# # Graph for all Genes that appear that are Postsynaptic
fullGraph <- generateGraph(count = 1, diseasehdoid = NULL, filename="PostsynapticNetwork/FullPSDDBNetwork.gml")


consensusGraph <- calcAllClustering(consensusGraph)
summary(consensusGraph)
write_graph(consensusGraph, file="PostsynapticNetwork/ConsensusPSDDBNetwork.gml", format="gml")

fullGraph <- calcAllClustering(fullGraph)
summary(fullGraph)
write_graph(fullGraph, file="PostsynapticNetwork/FullPSDDBNetwork.gml", format="gml")

consensusM <- clusteringSummary(consensusGraph)
fullM <- clusteringSummary(fullGraph)

View(consensusM)
View(fullM)


# geneList <- findGeneByCompartmentPaperCnt(cnt = 1)
# fullGraph <- read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml", format = "gml")
# node_info <- as_data_frame(id = V(fullGraph)$id, gene = V(fullGraph)$name)
#
#
# trubetskoyPriortisedDBCross <- read.table(Trubetskoy_broad_db_crossfile="Files/Trubetskoy_broad_db_cross.txt", sep="\t", header=TRUE)
#
# matchingList <- node_info[node_info$name %in% trubetskoyPriortisedDBCross]


# Presynaptic Graphs

# # Graph for Schziphrenia Genes that appear >= 2 that are Presynaptic
# generateGraph(file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml", localisation="Presynaptic")
#
# # Graph for all Schziphrenia Genes that appear that are Presynaptic
# generateGraph(count=1, file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml", localisation="Presynaptic")
#
# # Graph for Genes that appear >= 2 that are Presynaptic
# generateGraph(diseasehdoid = NULL, filename="PostsynapticNetwork/ConsensusPSDDBNetwork.gml", localisation="Presynaptic")
#
# # Graph for all Genes that appear that are Presynaptic
# generateGraph(count = 1, diseasehdoid = NULL, filename="PostsynapticNetwork/FullPSDDBNetwork.gml", localisation="Presynaptic")





