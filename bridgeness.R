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


calculateBridgeness <- function (graph, algorithm){
  # Calculate the consensus Community structure

  conmat <- makeConsensusMatrix(graph, N=5, alg=algorithm, type = 2, mask = 10, reclust = FALSE, Cnmax = 2000)

  # Calculate robustness

  clrob <-getRobustness(graph, alg = "louvain", conmat)
  pander(clrob)

  # Calculate the bridgeness
  br <- getBridgeness(graph, alg = "louvain", conmat)
  head(br)

  graph <- calcBridgeness(graph, alg = "louvain", conmat)
  vertex_attr_names(graph)



  sfile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
  shan<- read.table(sfile,sep="\t",skip=1,header=FALSE,strip.white=TRUE,quote="")
  head(shan)

  table(shan$V2)
  shan[shan$V2 =="Protein_cluster",] -> prCl
  dim(prCl)

  plotBridgeness(graph,alg = algorithm,
                 VIPs=prCl$V3,
                 Xatt='SL',
                 Xlab = "Semilocal Centrality (SL)",
                 Ylab = "Bridgeness (B)",
                 bsize = 3,
                 spsize =7,
                 MainDivSize = 0.8,
                 xmin = 0,
                 xmax = 1,
                 ymin = 0,
                 ymax = 1,
                 baseColor="royalblue2",
                 SPColor="royalblue2")
}


# Due to earlier results, we will be using louvain, sgG2 and spectral

gFull <- igraph::read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml",format="gml") #graph from gml
gConsensus <- igraph::read.graph("PostsynapticNetwork/ConsensusPSDDBNetwork.gml",format="gml")

gFullAlgorithms <- c("louvain", "sgG2", "spectral")

