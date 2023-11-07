# This is for reading

library(VennDiagram)
library(dplyr)
library(magrittr)
library(eulerr)


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("synaptome.db")
# BiocManager::install("BioNAR")

filePathFullDB <- "Files"

dbFile <- file.path("Files/Full_DB_Rat_Aut22.txt")


synapseDB <- read.table(file=dbFile, sep="\t", header=TRUE)


# For python, I'd make the 7 categories and append, but R must be easier

# Clears the grid for new page
grid.newpage()


psdpres = sum(subset(synapseDB, psd == 1 & pres == 1)$psd)
pressyn = sum(subset(synapseDB, syn == 1 & pres == 1)$syn)
psdsyn = sum(subset(synapseDB, psd == 1 & syn == 1)$psd)
psdpressyn = sum(subset(synapseDB, psd == 1 & pres == 1 & syn == 1)$psd)
psd = sum(synapseDB$psd)
pres = sum(synapseDB$pres)
syn = sum(synapseDB$syn)

# draw.triple.venn(area1 = psd, area2 = pres, area3= syn,
#                  n12 = psdpres, n23= pressyn, n13 = psdsyn, n123=psdpressyn,
#                  fill = c('yellow','brown','blue'),
#                  category = c("Post-Synaptic (psd)","Pre-Synaptic (pres)","Synaptosome (syn)"))
#
#TODO check synaptosome name


draw.triple.venn(area1 = psd, area2 = pres, area3= syn,
                 n12 = psdpres, n23= pressyn, n13 = psdsyn, n123=psdpressyn,
                 fill = c('yellow','brown','blue'),
                 category = c("PSD","PRES","SYN"))

