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

sfile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
shan<- read.table(sfile,sep="\t",skip=1,header=FALSE,strip.white=TRUE,quote="")
head(shan)

table(shan$V2)
shan[shan$V2 =="Protein_cluster",] -> prCl
dim(prCl)

plotBridgeness2 <- function(gg, alg, VIPs, Xatt = "SL", Xlab = "Semilocal Centrality (SL)", Ylab = "Bridgeness (B)", bsize = 3, spsize = 7, MainDivSize = 0.8, xmin = 0, xmax = 1, ymin = 0, ymax = 1, baseColor = "royalblue2", SPColor = "royalblue2") {
  indx <- match(V(gg)$name, VIPs)
  group <- ifelse(is.na(indx), 0, 1)
  X <- as.numeric(get.vertex.attribute(gg, Xatt, V(gg)))
  if (length(X) == 0) {
    stop("Graph vertices have no numerical attribute \"", Xatt, "\"\n")
  }
  X <- scale(X)
  Y <- as.numeric(get.vertex.attribute(gg, sprintf("BRIDGENESS.%s", alg), V(gg)))
  if (length(Y) == 0) {
    stop("Graph vertices have no numerical attribute \"", sprintf("BRIDGENESS.%s", alg), "\"\n", "Check that you\'ve calculated bridginess.\n")
  }


  trubetskoybroadcoding <- V(gg)$Trubetskoy2022broadcoding


  if ("GeneName" %in% vertex_attr_names(gg)) {
    lbls <- ifelse(!is.na(indx), V(gg)$GeneName, "")
    name <- V(gg)$GeneName
    dt <- data.frame(X = X, Y = Y, vips = group, entres = V(gg)$name, name = name, trubetskoybroadcoding = trubetskoybroadcoding)
  } else {
    lbls <- ifelse(!is.na(indx), V(gg)$name, "")
    name <- V(gg)$name
    dt <- data.frame(X = X, Y = Y, vips = group, entres = V(gg)$name, name = name, trubetskoybroadcoding = trubetskoybroadcoding)
  }
  dt_vips <- dt[dt$vips == 1, ]
  dt_res <- dt[dt$vips == 0, ]
  g <- ggplot(dt, aes(x = X, y = Y, label = name)) +
    geom_point(data = dt_vips, aes(x = X, y = Y), colour = baseColor, size = spsize, shape = 15, show.legend = FALSE) +

    geom_point(data = dt, aes(x = X, y = Y, color = trubetskoybroadcoding, alpha =0.1), size = bsize, show.legend = TRUE) +
    geom_label_repel(aes(label = as.vector(lbls)), fontface = "bold", color = "black", fill = "white", box.padding = 0.1,
                     point.padding = NA, label.padding = 0.15, segment.color = "black", force = 1, size = rel(3.8),
                     show.legend = FALSE, max.overlaps = 800) +
    labs(x = Xlab, y = Ylab, title = sprintf("%s", alg)) +
    scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
    scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red")) +
    theme(
      axis.title.x = element_text(face = "bold", size = rel(2.5)),
      axis.title.y = element_text(face = "bold", size = rel(2.5)),
      legend.title = element_text(face = "bold", size = rel(1.5)),
      legend.text = element_text(face = "bold", size = rel(1.5)),
      legend.key = element_blank(),
      panel.grid.major = element_line(colour = "grey40", linewidth = 0.2),
      panel.grid.minor = element_line(colour = "grey40", linewidth = 0.1),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(linetype = "solid", fill = NA)
    ) +
    geom_vline(xintercept = 0.5, colour = "grey40", linewidth = MainDivSize, linetype = 2, show.legend = FALSE) +
    geom_hline(yintercept = 0.5, colour = "grey40", linewidth = MainDivSize, linetype = 2, show.legend = FALSE)



  return(g)
}


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



  print("Plotting Bridgeness")

  plotName <- paste0("Bridgeness/",paste0(algorithm,"Bridgeness.png"))
  g <- plotBridgeness2(graph,alg = algorithm,
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
algsFull<-c('infomap','spectral')


for (algorithm in algsFull){
  calculateBridgeness(gFull,algorithm)

}

