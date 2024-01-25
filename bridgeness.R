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
  print(algorithm)
  print("Starting for algorithm")
  conmat <- makeConsensusMatrix(graph, N=5, alg=algorithm, type = 2, mask = 10, reclust = FALSE, Cnmax = 2000)

  # Calculate robustness
  print("Consensus Matrix created")
  clrob <-getRobustness(graph, alg = algorithm, conmat)
  pander(clrob)

  # Calculate the bridgeness
  br <- getBridgeness(graph, alg = algorithm, conmat)
  head(br)

  graph <- calcBridgeness(graph, alg = algorithm, conmat)
  vertex_attr_names(graph)

  print("Bridgeness calculated")

  sfile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
  shan<- read.table(sfile,sep="\t",skip=1,header=FALSE,strip.white=TRUE,quote="")
  head(shan)

  table(shan$V2)
  shan[shan$V2 =="Protein_cluster",] -> prCl
  dim(prCl)

  print("Plotting Bridgeness")

  plotName <- paste0("Bridgeness/",paste0(algorithm,"Bridgeness.pdf"))
  g <- plotBridgeness(graph,alg = algorithm,
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

  ggsave(plotName, g)



}


# Due to earlier results, we will be using louvain, sgG2 and spectral

gFull <- igraph::read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml",format="gml") #graph from gml
gConsensus <- igraph::read.graph("PostsynapticNetwork/ConsensusPSDDBNetwork.gml",format="gml")

gFull <- calcCentrality(gFull)
gFullLCC <- findLCC(gFull)
algsFull<-c('wt', 'fc', 'infomap', 'louvain',
            'sgG1', 'sgG2', 'sgG5', 'spectral')

#
# for (algorithm in algsFull){
#   calculateBridgeness(gFull,algorithm)
#
# }

