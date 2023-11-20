library(knitr)
library(BioNAR)
library(synaptome.db)
library(ggplot2)
library(pander)
library(ggrepel)
library(randomcoloR)
library(RSQLite)
library(synaptome.db)


setwd("~/Documents/Dissertation")
query <- "SELECT HumanEntrez
FROM FullGeneFullDisease
WHERE HDOID = \"DOID:5419\""

db <- dbConnect(SQLite(), dbname = "Database/DS_10283_3877/synaptic.proteome_SR_20210408.db.sqlite")
dbData <- dbGetQuery(db, query)

gt <- buildNetwork(db)
summary(gt)

