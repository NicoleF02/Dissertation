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
# Get summary table with the enrichment for each communities, how many communities are enriched, what are the p-values for this enrichment.

# How many clusters, how many have enrichment, what the p-value is, total number and enriched.

# Look at different ontology terms and find associations betweeen enirched communities and their functions
# See if things pop up only in certian clustering algorithms, and see if they depend on algorithms.

gFull <- igraph::read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml",format="gml") #graph from gml
gConsensus <- igraph::read.graph("PostsynapticNetwork/ConsensusPSDDBNetwork.gml",format="gml")


oraLouvian <- clusterORA(gFull, 'louvain', name = 'Trubetskoy2022broadcoding',
                         vid = "name", alpha = 1, col = COLLAPSE)

#Focus on what the data is, what my method is, how do you interpret it and how I actually
# process the data and what is actually new. what is actually novel

#
# algsConsensus<-c('lec', 'wt', 'fc', 'infomap', 'louvain',
#         'sgG1', 'sgG2', 'sgG5', 'spectral')
#
# algsFull<-c('wt', 'fc', 'infomap', 'louvain',
#         'sgG1', 'sgG2', 'sgG5', 'spectral')
#
#
# oraFull<-lapply(algsFull, function(alg){clusterORA(gFull, alg, name = 'Trubetskoy2022broadcoding',
#                                            vid = "name",alpha = 1, col = COLLAPSE)})
# names(oraFull)<-algsFull
# FeMaxFull<-log2(max(sapply(oraFull,function(d){max(d$Fe)})))
# FcMaxFull<-log2(max(sapply(oraFull,function(d){max(d$Fc)})))
#
#
# oraConsensus<-lapply(algsConsensus, function(alg){clusterORA(gConsensus, alg, name = 'Trubetskoy2022broadcoding',
#                                                vid = "name",alpha = 1, col = COLLAPSE)})
# names(oraConsensus)<-algsConsensus
# FeMaxConsensus<-log2(max(sapply(oraConsensus,function(d){max(d$Fe)})))
# FcMaxConsensus<-log2(max(sapply(oraConsensus,function(d){max(d$Fc)})))
#
#
# # Summary for Full one
# statsR1Full <- summaryStats(oraFull, 0.1, usePadj=FALSE, FeMAX=FeMaxFull, FcMAX=FcMaxFull)
# names(statsR1Full)
# View(head(statsR1Full$CAN))
#
#
# # Summary for Consensus one
# statsR1Consensus <- summaryStats(oraConsensus, 0.1, usePadj=FALSE, FeMAX=FeMaxConsensus, FcMAX=FcMaxConsensus)
# names(statsR1Consensus)
# View(head(statsR1Consensus$CAN))

#
# plots<-plotRatio(x=statsR1Full, desc="p.values",LEGtextSize=0.75, LEGlineSize=2)
# View(plots$ranktable)
# print(plots$p3)



#
# plots<-plotRatio(x=statsR1Consensus, desc="p.values",LEGtextSize=0.75, LEGlineSize=2)
# View(plots$ranktable)
# print(plots$p3)



