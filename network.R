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

# Table: each cluster, cluster size, how many genes are schizophrenic, test for fisher test to check for overpresentation
# annotation will be trubetskoy to check enrichment, it will show which ones are enriched, do it for all algorithms.
# also do the same for presynaptic.
# Keep p value and adjusted p value using ORA
#

initGeneTable <- function(count, localisation="Postsynaptic"){
  if (count != 1){
    gp <- findGeneByCompartmentPaperCnt(count)
    t <- filter(gp, Localisation == localisation)

  }else{
    cid<-match(localisation, getCompartments()$Name)
    t<-getAllGenes4Compartment(cid)
  }

  return(t)
}

networkSetup <- function(network){
  network <- annotateGeneNames(network)
  network <- annotateGOont(network)

  afile<-system.file("extdata", "flatfile_human_gene2HDO.csv", package = "BioNAR")
  dis <- read.table(afile,sep="\t",skip=1,header=FALSE,strip.white=TRUE,quote="")
  network <- annotateTopOntoOVG(network, dis)

  sg <- read.table("Files/SynGO.txt", sep = "\t", header = T, stringsAsFactors = F)
  network <- annotateVertex(network, name = "syngo", values = sg, idatt = "name")

  trubetskoyBroad<- read.table("Files/Trubetskoy_2022_broad_coding.txt", sep = "\t", header = T, stringsAsFactors = F)
  network <- annotateVertex(network, name="Trubetskoy_2022_broad_coding", trubetskoyBroad[,c(1,2)])

  trubetskoyPriortised<- read.table("Files/Trubetskoy_2022_priortised_coding.txt", sep = "\t", header = T, stringsAsFactors = F)
  network <- annotateVertex(network, "Trubetskoy_2022_priortised_coding", trubetskoyPriortised[,c(1,2)])

  reducedRatDB <- read.table("Files/ReducedRatDB.txt", sep="\t", header= T, stringsAsFactors = F)
  network <- annotateVertex(network, name="Synapse Locations", values =reducedRatDB , idatt="name")

  priortisedSchizophreniaDB <- read.table("Files/SchizphreniaPriortisedDB.txt", sep="\t", header= T, stringsAsFactors = F)
  network <- annotateVertex(network, "SchizophreniaGeneCount", priortisedSchizophreniaDB[,c(1,2)])

  return(network)
}

generateGraph <- function(count=1, localisation="Postsynaptic", filename, synGoOnly = FALSE){
  if (synGoOnly){
    print("SynGO Only")
    geneTable <- synGoOnly()
  }else{
    geneTable <- initGeneTable(count=count, localisation=localisation)

  }

  g <- graphFromSynaptomeByEntrez(geneTable$HumanEntrez)
  networkGene <- networkSetup(g)
  network <- findLCC(networkGene)
  write_graph(network, file=filename, format="gml")
  return (network)
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


synGoOnly <- function(){
  data <- read.table("Files/SynGO.txt")
  colnames(data) <- c("HumanEntrez", "synGo")
  return(data)
}


# SynGO graph
graph <- generateGraph(count = 1, file='PostsynapticNetwork/SynGoNetwork.gml', synGoOnly=TRUE)





# Postsynaptic Graphs

# # Graph for Schziphrenia Genes that appear >= 2 that are Postsynaptic
# generateGraph(file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml")
#
# # Graph for all Schziphrenia Genes that appear that are Postsynaptic
# generateGraph(count=1, file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml")
#

#
# Warning in igraph::leading.eigenvector.community(ugg, weights = weights) :
# At core/linalg/arpack.c:805 : ARPACK solver failed to converge (10001 iterations, 0/1 eigenvectors converged).
# Warning in getClustering(gg, alg, weights = weights) :
# Clustering calculations for algorithm "lec" failed. NULL is returned


# # Graph for all Genes that appear that are Postsynaptic
# fullGraph <- generateGraph(count = 1, filename="PostsynapticNetwork/FullPSDDBNetwork.gml")
# fullGraph <- findLCC(fullGraph)
# # Graph for Genes that appear >= 2 that are Postsynaptic
# consensusGraph <- generateGraph(count = 2, filename="PostsynapticNetwork/ConsensusPSDDBNetwork.gml")
# consensusGraph <- findLCC(consensusGraph)
# # consensusGraph <- read.graph("PostsynapticNetwork/ConsensusPSDDBNetwork.gml", type="gml")
# #
# consensusGraph <- calcAllClustering(consensusGraph)
# #Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "sgG1" failed. NULL is returned
# # Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "sgG2" failed. NULL is returned
# # Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "sgG5" failed. NULL is returned
#
#
# summary(consensusGraph)
# write_graph(consensusGraph, file="PostsynapticNetwork/ConsensusPSDDBNetwork.gml", format="gml")
#
# fullGraph <- calcAllClustering(fullGraph)
# summary(fullGraph)
# write_graph(fullGraph, file="PostsynapticNetwork/FullPSDDBNetwork.gml", format="gml")
#
# # Warning in igraph::leading.eigenvector.community(ugg, weights = weights) :
# # At core/linalg/arpack.c:805 : ARPACK solver failed to converge (10001 iterations, 0/1 eigenvectors converged).
# # Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "lec" failed. NULL is returned
# # Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "sgG1" failed. NULL is returned
# # Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "sgG2" failed. NULL is returned
# # Warning in getClustering(gg, alg, weights = weights) :
# # Clustering calculations for algorithm "sgG5" failed. NULL is returned


#
# consensusM <- clusteringSummary(consensusGraph)
# fullM <- clusteringSummary(fullGraph)

# View(consensusM)
# View(fullM)

# cid<-match("Postsynaptic", getCompartments()$Name)
# t<-getAllGenes4Compartment(cid)
# ggb<-buildNetwork(getPPIbyEntrez(t$HumanEntrez, type = "limited"))
#
# papercnt <- findGeneByPaperCnt(1)
# papercptcnt <- findGeneByCompartmentPaperCnt(1)
#
# papercntG <- buildNetwork(getPPIbyEntrez(papercnt$HumanEntrez, type = "limited"))
# papercptcntG <- buildNetwork(getPPIbyEntrez(papercptcnt$HumanEntrez, type= "limited"))
#
# m <- read.table("Files/Trub_broad.txt", sep = "\t", header = T, stringsAsFactors = F)
#
# papercntG <- applpMatrixToGraph(papercntG, m)
# papercptcntG <- applpMatrixToGraph(papercptcntG, m)
# gg <- applpMatrixToGraph(ggb, m)

# g <- igraph::read.graph("PostsynapticNetwork/PSD_FULL_clustered.gml",format="gml")
# m <- read.table("Files/Trub_broad.txt", sep = "\t", header = T, stringsAsFactors = F)
# network <- annotateVertex(network, "TrubBroad", trubetskoyBroad[,c(1,2)])



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





