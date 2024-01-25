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


fitSigmoid2 <- function(stat, SDv = c(0, 0.05, 0.1, 0.5)) {
  df <- stat$SUM3
  # x.range.value was changed from 6.0 to X4.0.1
  x.range.value <- "6.0"
  #Xmax <- match(x.range.value, colnames(df))
  Xmax <- "6.0"
  tt <- df[, seq_len(Xmax)]
  colnames(tt) <- colnames(df)[seq_len(Xmax)]
  N <- length(colnames(tt))
  x <- as.numeric(colnames(tt)[3 : N])
  SDlab <- as.character(SDv)
  rates <- c(-10, -5, -2, -1, -0.5)
  resout <- list()
  for (s in seq_along(SDv)) {
    models <- list()
    GPLOTS <- list()
    gof <- list()
    CNfit <- c("alg", "isConv", "finTol", "stopCode", "stopMessage")
    fitInfo <- data.frame()
    parInfo <- data.frame()
    for (i in seq_along(tt[, 1])) {
      y <- as.numeric(tt[i, 3 : N]) / as.numeric(tt[i, 2])
      y <- BioNAR:::addNoise(y, SD = SDv[s])
      pp <- list(a = 0, b = round(max(y)), c = -2, d = round(median(x)))
      m.s <- minpack.lm::nlsLM(`~`(y, a + ((b - a) / (1 + exp(-c * (x - d))))), start = pp, trace = FALSE, control = minpack.lm::nls.lm.control(maxiter = 150))
      models[[i]] <- m.s
      names(models)[i] <- as.character(tt[i, 1])
      rm(m.s)
      fitInfo <- rbind(fitInfo, storeFitInfo(names(models)[i], unlist(models[[i]]$convInfo)))
      parInfo <- rbind(parInfo, storeParInfo(names(models)[i], unlist(summary(models[[i]])$parameters[, c(1, 2)])))
    }
    colnames(fitInfo) <- c("alg", "isConv", "finIter", "finTol", "stopCode", "stopMessage")
    fitInfo$finIter <- as.numeric(fitInfo$finIter)
    fitInfo$finTol <- as.numeric(fitInfo$finTol)
    colnames(parInfo) <- c("alg", c(rbind(names(models[[1]]$m$getPars()), sprintf("sd_%s", names(models[[1]]$m$getPars())))))
    CN <- c("alg", sprintf("Rate_%.1f", rates))
    oo <- matrix("", ncol = length(CN), nrow = length(names(models)))
    colnames(oo) <- CN
    oo[, 1] <- names(models)
    for (i in seq_along(names(models))) {
      ks <- BioNAR:::gofs(x, rates, models[[i]])
      indx <- BioNAR:::highlightRate(rates = rates, val = -2)
      PV <- as.numeric(ks[[1]]$p.value)
      if (indx != -1) {
        PV <- as.numeric(ks[[indx[1]]]$p.value)
      }
      tmp <- BioNAR:::plotSigmoid(x = x, rates = rates, model = models[[i]], alg = names(models)[i], pv = PV)
      gof[[i]] <- ks
      GPLOTS[[i]] <- tmp$gplot
      names(gof)[i] <- names(models)[i]
      names(GPLOTS)[i] <- names(models)[i]
      for (j in seq_along(rates)) {
        oo[i, (j + 1)] <- ks[[j]]$p.value
      }
    }
    p <- cowplot::plot_grid(plotlist = GPLOTS, labels = "AUTO", label_size = 20, label_fontface = "bold")
    ooks <- as.data.frame(oo)
    ooks[-1] <- lapply(ooks[-1], as.numeric)
    pinf <- as.data.frame(parInfo)
    pinf[-1] <- lapply(pinf[-1], as.numeric)
    resout[[SDlab[s]]] <- list(gridplot = p, plots = GPLOTS, fitInfo = fitInfo, parInfo = pinf, ks = ooks)
    rm(models, GPLOTS, gof)
  }
  return(resout)
}



# Reduced here means those without any null values / infinity
gFull <- igraph::read.graph("PostsynapticNetwork/FullPSDDBNetwork.gml",format="gml") #graph from gml
gFull <- calcCentrality(gFull)
gFullRedcuded <- findLCC(gFull)

# Spectral and wt have been removed
algsFullReduced <-c('fc', 'infomap', 'louvain',
            'sgG1', 'sgG2', 'sgG5')

oraFull<-lapply(algsFullReduced, function(alg){clusterORA(gFullRedcuded, alg, name = 'Trubetskoy2022broadcoding',
                                                   vid = "name",alpha = 1, col = COLLAPSE)})
names(oraFull)<-algsFullReduced

FeMaxFull<-log2(max(sapply(oraFull,function(d){max(d$Fe)})))
FcMaxFull<-log2(max(sapply(oraFull,function(d){max(d$Fc)})))

statsR1Full <- summaryStats(oraFull, 0.1, usePadj=FALSE, FeMAX=FeMaxFull, FcMAX=FcMaxFull)
names(statsR1Full)


alg <- "louvain"
conmat <- makeConsensusMatrix(gFullRedcuded, N=5,
                              alg = alg, type = 2,
                              mask = 10,reclust = FALSE,
                              Cnmax = 10)


clrob<-getRobustness(gFullRedcuded, alg = alg, conmat)
pander(clrob)


alg <- "fc"
conmat <- makeConsensusMatrix(gFullRedcuded, N=5,
                              alg = alg, type = 2,
                              mask = 10,reclust = FALSE,
                              Cnmax = 10)


clrob<-getRobustness(gFullRedcuded, alg = alg, conmat)
pander(clrob)



alg <- "infomap"
conmat <- makeConsensusMatrix(gFullRedcuded, N=5,
                              alg = alg, type = 2,
                              mask = 10,reclust = FALSE,
                              Cnmax = 10)


clrob<-getRobustness(gFullRedcuded, alg = alg, conmat)
pander(clrob)



alg <- "sgG1"
conmat <- makeConsensusMatrix(gFullRedcuded, N=5,
                              alg = alg, type = 2,
                              mask = 10,reclust = FALSE,
                              Cnmax = 10)


clrob<-getRobustness(gFullRedcuded, alg = alg, conmat)
pander(clrob)
# fitres<-fitSigmoid2(statsR1Full)
# print(fitres[['0']]$gridplot)

