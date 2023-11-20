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
# We need to split this up
# No. 1. get all schizophrenia genes
# get them for postsynaptic
# get them for only postsynaptic only
# then do networks
# do refined network



initGeneTable <- function(count=2, localisation="Postsynaptic", diseasehdoid="DOID:5419"){
  narrowGeneList <- findGeneByCompartmentPaperCnt(cnt = count)

  tablePostsynaptic <- filter(narrowGeneList, Localisation == localisation)

  columnGenes <- tablePostsynaptic$HumanEntrez

  # then we want to get all genediease, then filter. Maybe check own db rather then their own?

  tableDisease <- getGeneDiseaseByEntres(columnGenes)

  # if its been identified 2 or more, times, not nescessarily for schizophrenia !!!! double check!!

  tableNarrowDisease <- filter(tableDisease, HDOID == diseasehdoid)

  return(tableNarrowDisease)
}

tableSchizophrenicNarrow <- initGeneTable()
tableSchizophrenicBroad <- initGeneTable(count=1)


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

  return(network)
}

ppiNarrow <- getPPIbyEntrez(tableSchizophrenicNarrow$HumanEntrez, type = "limited")
networkNarrow <- networkSetup(ppiNarrow)


ppiBroad <- getPPIbyEntrez(tableSchizophrenicBroad$HumanEntrez, type = "limited")
networkBroad <- networkSetup(ppiBroad)


write_graph(networkNarrow, file="PostsynapticNetwork/NarrowPSDSchizphreniaNetwork.gml", format="gml")

write_graph(networkBroad, file="PostsynapticNetwork/BroadPSDSchizphreniaNetwork.gml", format="gml")


