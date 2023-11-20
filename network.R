library(knitr)
library(BioNAR)
library(synaptome.db)
library(ggplot2)
library(pander)
library(ggrepel)
library(randomcoloR)
library(dplyr)

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

#consider running both with other compartments and removing overlap.

system.file(package = "BioNAR")

# ggb <- buildFromSynaptomeGeneTable(tableSchizophrenicNarrow, type="limited")
# write_graph(ggb, file="PostsynapticNetworkBroad/PSDSchizphreniaNetwork.gml", format="gml")