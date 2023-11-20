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

#
# cid <- match("Postsynaptic",getCompartments()$Name)
#
# #\"DOID:5419\""
#
# tablePostsynaptic <- getAllGenes4Compartment(cid)


narrowGeneList <- findGeneByCompartmentPaperCnt(cnt = 2)

tablePostsynaptic <- filter(narrowGeneList, Localisation == "Presynaptic")

columnGenes <- tablePostsynaptic$HumanEntrez

# then we want to get all genediease, then filter. Maybe check own db rather then their own?

tableDisease <- getGeneDiseaseByEntres(columnGenes)

# if its been identified 2 or more, times, not nescessarily for schizophrenia !!!!
tableCountNarrow <- findGeneByPaperCnt(cnt = 2)


tableSchizophrenicNarrow <- filter(tableDisease, HDOID == "DOID:5419")


