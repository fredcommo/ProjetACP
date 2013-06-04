# Filtering on Kim Lung data

require(splines)
source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/pcaSelect.R')
source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/pcaInfo.R')
source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/pcaTrace.R')
source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/heatmap.3.R')

op <- par(no.readonly = TRUE)

.computeBounds <- function(M, S){
  mCuts <- cut(M, breaks = seq(min(M, na.rm = TRUE), max(M, na.rm = TRUE), by = 0.2)) #
  labs <- levels(mCuts)
  mBounds <- cbind.data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                              upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  mBounds$lower[1] <- min(M, na.rm = TRUE)
  mBounds$upper[nrow(mBounds)] <- max(M, na.rm = TRUE)  
  mBounds <- cbind.data.frame(lower = mBounds$lower, med = (mBounds$upper + mBounds$lower)/2, upper = mBounds$upper)
  sBounds <- lapply(seq(1, nrow(mBounds)), function(x){
    index <- as.numeric(which(M >= mBounds$lower[x] & M < mBounds$upper[x]))
    if(length(index)<1) cbind(S = NA, M = mBounds$med[x], Mlower = mBounds$lower[x], Mupper = mBounds$upper[x])
    else cbind(S = S[index], M = mBounds$med[x], Mlower = mBounds$lower[x], Mupper = mBounds$upper[x])
  })
  sBounds <- do.call(rbind, sBounds)
  rownames(sBounds) <- seq(1, nrow(sBounds))
  return(as.data.frame(sBounds))
}

.pickM <- function(dM){
  m = as.numeric(sample(dM$x, 1, prob = dM$y))
  return(m)
}

.pickS <- function(MStable, m){
  j <- which(MStable$Mlower <= m & MStable$Mupper > m)
  #  if(length(j)<1) {s <- min(dS$x, na.rm = TRUE); cat('min S used for m:', m, '\n')}
  if(length(j)<2) {s <- min(MStable$S, na.rm = TRUE); cat('\nmin S used for m:', m, '\n')}
  else{
    tmpS <- MStable$S[j]
    if(length(tmpS)>2)
      s <- sample(tmpS, 1, prob = pchisq(tmpS, df = 1, ncp = m, lower.tail = FALSE))
    else s <- sample(tmpS, 1)
  }
  return(s)
}

.generateRandom <- function(Mvector, Svector, n, p){
  # provide a vector of means and a vector of corresponding Sdev
  # n samples, p probes
  # returns a (p, n) matrix of random probes
  dM <- density(Mvector)
  MStable <- .computeBounds(Mvector, Svector)
  output <- lapply(seq(1, p),
                   function(x){ if(x%%ceiling(p/10) == 0) cat(x, '\t')
                                m <- max(min(MStable$M),.pickM(dM))
                                if(is.na(m)) stop('m is NA')
                                s <- .pickS(MStable, m)
                                values <- rnorm(n, m, s)
                                if(any(is.na(values))) stop('s', s, 'is NA for m:', m)
                                return(values)
                   })
  output <- do.call(rbind, output)
  rownames(output) <- paste0('random', seq(1, nrow(output)))
  colnames(output) <- paste0('sample', seq(1, ncol(output)))
  return(output)
}

.generateGrps <- function(Mvector, Svector, n, p, grps = NULL, nGrp = 2, minP = 0.3, maxP = 0.8){
  # Provide vector of means and coresspondong Sdev
  # n, p : number of samples (provides if grps if not provided), number of probes, resp.
  # minP, maxP : min and max proportion of samples in a grp for which a specif probe is generated.
  # Returns a (p, n) matrix of probes specific to grps and the vector of grps as a list.
  
  # if grps not provided, generate random grps according to nGrps
  if(is.null(grps))
    grps <- factor(rbinom(n, nGrp - 1, 0.5), labels = LETTERS[1:nGrp])
  else if(!is.factor(grps)) grps <- as.factor(grps)
  if(length(grps) != n) n <- length(grps)
  nGrps <- nlevels(grps)
  dM <- density(Mvector)
  MStable <- .computeBounds(Mvector, Svector)
  output <- lapply(seq(1, p),
                   function(x){ if(x%%ceiling(p/10) == 0) cat(x, '\t')
                                m <- max(min(MStable$M),.pickM(dM))
                                if(is.na(m)) stop('m is NA')
                                grp <- sample(levels(grps), 1)
                                s <- .pickS(MStable, m)
                                values <- rnorm(n, m, s)
                                # with respect to grps but some samples only, according to minP/maxP
                                idx <- sample(which(grps == grp))
                                rbi <- rbinom(length(idx), 1, prob = sample(seq(minP, maxP, by = 0.05), 1))
                                idx <- idx[rbi == 1]
                                change <- sample(c('inc', 'dec'), 1)
                                if(change == 'inc')
                                  newmu <- min(max(MStable$M)*0.95, m*sample(seq(1.5, 2.5, by = 0.2), 1))
                                else
                                  newmu <- max(min(MStable$M)*1.05, m/sample(seq(1.5, 2.5, by = 0.2), 1))
                                values[idx] <- rnorm(length(idx), newmu, .pickS(MStable, newmu))
                                if(any(is.na(values))) stop('s', s, 'is NA for m:', m)                                      
                                return(values)
                   })
  output <- do.call(rbind, output)
  rownames(output) <- paste0('signif', seq(1, nrow(output)))
  colnames(output) <- paste0('sample', seq(1, ncol(output)))
  return(list(Data = output, Grps = grps))
}

################################################################
# dir.create('/Users/fredcommo/Documents/MyProjects/Simulations/ALLsimulation')

saveDir <- '/Users/fredcommo/Documents/MyProjects/Simulations/ALLsimulation/'
require(ALL)
data(ALL)
eset <- exprs(ALL)
samples <- pData(ALL)

M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)
trueDist <- .computeBounds(M, S)
factM <- factor(trueDist$M)

# Estimate row means, and row Sd
ns1 <- ns(M, df = 3)
lm1 <- lm(log10(S) ~ ns1)
png(file = paste0(saveDir, 'ALL_real_SdevMeans_distrib.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plot(M, log10(S), pch = 19, cex = 0.2, col = 'grey', xlab = 'Means', ylab = 'Log10(Sdev)')
lines(sort(M), lm1$fitted[order(M)], col = 'red3', lwd = 4)
dev.off()

png(file = paste0(saveDir, 'ALL_real_SdevMeans_boxplot.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
boxplot(log10(trueDist$S)~factor(trueDist$M), names = round(unique(trueDist$M), 2), 
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)')
dev.off()

# Compute white probes
whiteProbes <- .generateRandom(M, S, ncol(eset), ceiling(nrow(eset)*0.1))
  
whiteM <- apply(whiteProbes, 1, mean, na.rm = TRUE)
whiteS <- apply(whiteProbes, 1, sd, na.rm = TRUE)
whiteProbesTable <- .computeBounds(whiteM, whiteS)

png(file = paste0(saveDir, 'ALL_random_SdevMeans_boxplot.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
boxplot(log10(whiteProbesTable$S) ~ factor(whiteProbesTable$M), names = round(unique(whiteProbesTable$M), 2), 
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)')
dev.off()

# initialize the random samples in case of choice3
allType <- as.factor(ifelse(grepl('B', samples$BT), 'B', 'T'))
simulateGrps <- .generateGrps(M, S, n = ncol(eset), p = ceiling(nrow(eset)*0.01), grps = allType)
signifProbes <- simulateGrps$Data

signifM <- apply(signifProbes, 1, mean, na.rm = TRUE)
signifS <- apply(signifProbes, 1, sd, na.rm = TRUE)
signifProbesTable <- .computeBounds(signifM, signifS)

png(file = paste0(saveDir, 'ALL_simulated_SdevMeans_distrib.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
smoothScatter(M, log10(S))
points(whiteM, log10(whiteS), pch = 19, cex = 0.5, col = 'red3')
points(signifM, log10(signifS), pch = 19, cex = 0.75, col = 'pink')
dev.off()

boxplot(log10(signifProbesTable$S) ~ signifProbesTable$M, 
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)')

# Add the white probes into the eset matrix
colnames(whiteProbes) <- colnames(signifProbes) <- colnames(eset)
#eset <- .insertRow(eset, whiteProbes, sample(1:nrow(eset), nrow(whiteProbes))) 
#eset = kimData$eset
eset <- rbind(eset, whiteProbes, signifProbes)
randomIdx <- grep('random', rownames(eset))
signifIdx <- grep('signif', rownames(eset))


# Show the white probes
pcaProbes <- prcomp(eset)
png(file = paste0(saveDir, 'ALL_pca_simulProbesLocation.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
pairs(pcaProbes$x[,1:3], cex = 0.5,
       col = ifelse(grepl('random', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'darkblue', 'lightblue')))
dev.off()

# Use multtest to find the diff expr. genes between tumorTypes, before and after filtering
require(multtest)
fullTest <- mt.maxT(eset, classlabel = factor(allType), B = 5000)
bestFull <- fullTest$index[fullTest$adjp<0.01 & !grepl('signif', rownames(fullTest))]
length(bestFull)
#smoothScatter(pcaProbes$x[,2:3])

png(file = paste0(saveDir, 'ALL_pca_realProbesLocation.png'), width = 800, height = 800)
Cols <- ifelse(grepl('random', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'darkblue', 'grey'))
Cols[bestFull] <- 'cyan'
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
pairs(pcaProbes$x[,1:3], cex = 0.5, pch = 19,
     col = Cols)
#points(pcaProbes$x[bestFull,2:3], pch = 19, cex = 0.5, col = 'cyan')
dev.off()

png(file = paste0(saveDir, 'ALL_informationCurve.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
score <- pcaTrace(eset, pcaProbes, main = 'Information curve')
dev.off()

infos <- pcaInfo(score)
redDot = 10^(1/score$lModel$x.intercept)
select <- pcaSelect(score, redDot)
length(select)

filtTest <- mt.maxT(eset[select,], classlabel = factor(allType), B = 5000)
filtIndex <- filtTest$adjp<0.01 & !grepl('signif', rownames(filtTest))
any(grepl('random', rownames(filtTest)[filtIndex]))
bestFilt <- rownames(filtTest)[filtTest$adjp<0.01 & !grepl('signif', rownames(filtTest))]
bestFilt <- which(rownames(eset) %in% bestFilt)
length(bestFilt)

png(file = paste0(saveDir, 'ALL_pca_filteredProbesLocation.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
Cols <- ifelse(grepl('random', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'darkblue', 'grey'))
Cols[bestFull] <- 'cyan'
Cols[bestFilt] <- 'steelblue3'
plot(pcaProbes$x[,2:3], cex = 0.75, pch = 19, asp = 1, col = Cols)
#points(pcaProbes$x[bestFull,2:3], cex = 0.5, col = 'cyan')
#points(pcaProbes$x[bestFilt,2:3], cex = 0.5, col = 'steelblue3')
dev.off()

# Heatmap signature
Values <- t(eset[bestFilt,])
png(file = paste0(saveDir, 'ALL_heatmap_filteredProbes_Significant.png'), width = 800, height = 800)
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          RowSideColors = c('grey25', 'grey75')[factor(allType)],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

# Cut at a higher value and try to find clusters
infos
select <- pcaSelect(score, 0.50)
length(select)
#Values <- t(eset[which(rownames(eset) %in% bestFilt),])
Values <- eset[select,]
dim(Values)

png(file = paste0(saveDir, 'ALL_heatmap_filteredProbes_070.png'), width = 800, height = 800)
heatmap.3(t(Values),
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          #ColSideColors = ifelse(grepl('signif', colnames(Values)), 'steelblue', 'white'),
          RowSideColors = c('grey25', 'grey75')[factor(allType)],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

clust <- hclust(dist(t(Values)), 'ward')
plot(clust)
K = 7
cutClust <- cutree(clust, k = K)
png(file = paste0(saveDir, ,'ALL_Simulations/heatmap_filteredProbes_070_clusters.png'), width = 800, height = 800)
heatmap.3(t(Values),
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          RowSideColors = c('steelblue','seagreen','violet','indianred','purple','lightblue')[factor(cutClust)],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

clustTest <- mt.maxT(eset, classlabel = c(0:(K-1))[cutClust], test = 'f', B = 5000)
sum(clustTest$adjp < 0.001)/nrow(clustTest)
bestClust <- rownames(clustTest)[clustTest$adjp < 0.001]
bestClust <- which(rownames(eset) %in% bestClust)
length(bestClust)

par(op)
png(file = paste0(saveDir, , 'ALL_pca_subGrpsProbes.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
Cols <- ifelse(grepl('random', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'blue', 'grey'))
Cols[bestFull] <- 'cyan'
Cols[bestFilt] <- 'steelblue3'
Cols[setdiff(bestClust, c(bestFull, bestFilt))] <- 'seagreen4'
pairs(pcaProbes$x[,2:4], pch = 19, cex = 0.5, col = Cols)
points(pcaProbes$x[bestFull,2:3], cex = 0.25, col = 'cyan')
points(pcaProbes$x[which(rownames(eset) %in% bestFilt),2:3], cex = 0.25, col = 'steelblue3')
points(pcaProbes$x[bestClust,2:3], pch = 19, cex = 0.5, col = 'seagreen4')
dev.off()

png(file = paste0('/gluster/home/jcommo/Simulations/heatmap_filteredProbes_070_clusters.png'), width = 800, height = 800)
heatmap.3(t(eset[bestClust,]),
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          RowSideColors = c('steelblue','seagreen','violet','indianred','purple','lightblue')[factor(cutClust)],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

# Testing covariates
NAs <- which(is.na(samples$citog))
citog <- factor(samples$citog[-NAs])
citogTest <- mt.maxT(Values[-NAs,], classlabel = c(0:(nlevels(citog)-1))[citog], test = 'f', B = 5000)
sum(citogTest$adjp < 0.01)/nrow(citogTest)
bestCitog <- rownames(citogTest)[citogTest$adjp < 0.1]
bestCitog <- which(rownames(eset) %in% bestCitog)
length(bestCitog)

BT <- factor(samples$BT)
btTest <- mt.maxT(Values, classlabel = c(0:(nlevels(BT)-1))[BT], test = 'f', B = 5000)
sum(btTest$adjp < 0.01)/nrow(btTest)
bestBT <- rownames(btTest)[btTest$adjp < 0.01]
bestBT <- which(rownames(eset) %in% bestBT)
length(bestBT)

par(op)
png(file = paste0(saveDir, 'pca_subGrpsProbes.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plot(pcaProbes$x[,2:3], cex = 0.1,
     col = ifelse(grepl('random', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'blue', 'grey')))
points(pcaProbes$x[bestFull,2:3], cex = 0.25, col = 'cyan')
points(pcaProbes$x[bestFilt,2:3], cex = 0.25, col = 'steelblue3')
points(pcaProbes$x[bestBT,2:3], cex = 1, col = 'seagreen4')
dev.off()


# Try to find clusters and compare them
png(file = paste0(saveDir, 'ALL_heatmap_citog.png'), width = 800, height = 800)
heatmap.3(t(eset[bestCitog,]),
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          #ColSideColors = ifelse(grepl('signif', colnames(Values)), 'steelblue', 'white'),
          #RowSideColors = c(1:10)[factor(samples$citog)],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

# Test the cluster groups on the filtered probes
clustGrp <- mt.maxT(t(Values), classlabel = samples$TumourType, B = 5000)
clustTest <- mt.maxT(t(Values), classlabel = c(0:(K-1))[cutClust], test = 'f', B = 5000)
sum(clustTest$adjp < 0.01)/nrow(clustTest)
bestClust <- rownames(clustTest)[clustTest$adjp < 0.01]
bestClust <- which(rownames(eset) %in% bestClust)
length(bestClust)

par(op)
png(file = paste0(saveDir, 'pca_subGrpsProbes.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plot(pcaProbes$x[,2:3], cex = 0.1,
     col = ifelse(grepl('white', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'blue', 'grey')))
points(pcaProbes$x[bestFull,2:3], cex = 0.25, col = 'cyan')
points(pcaProbes$x[bestFilt,2:3], cex = 0.25, col = 'steelblue3')
points(pcaProbes$x[select,2:3], cex = 0.5, col = 'seagreen4')
dev.off()

# Test the cluster groups on the full matrix
clustTest <- mt.maxT(eset, classlabel = c(0:(K-1))[cutClust], test = 'f', B = 5000)
sum(clustTest$adjp < 0.01)/nrow(clustTest)
bestClust <- clustTest$index[clustTest$adjp < 0.01 & !grepl('signif', rownames(clustTest))]
length(bestClust)

par(op)
png(file = paste0(saveDir, 'pca_subGrpsProbes.png'), width = 1000, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5, mar = c(5.1, 5.1, 4.1, 2.1))
plot(pcaProbes$x[,2:3], cex = 0.1, col = 'grey80', main = 'Probes location')
abline(h = 0, v = 0, col = 'grey75', lty = 2, lwd = 3)
points(pcaProbes$x[whiteIdx,2:3], pch = 19, cex = 0.5, col = 'red3')
points(pcaProbes$x[bestFilt,2:3], pch = 19, cex = 0.5, col = 'steelblue3')
points(pcaProbes$x[setdiff(bestClust, bestFilt),2:3], pch = 19, cex = 0.5, col = 'seagreen4')
legend('topleft',
       legend = c('real probes', 'generated random probes', 'Ad Vs. Sq signif', 'sub-groups signif'),
       pch = 19, cex = 1.5, col = c('grey', 'red', 'steelblue3', 'seagreen4'),
       bty = 'n')
dev.off()

# Heatmap on subGrp probes
Values <- t(eset[setdiff(bestClust, bestFull),])
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
#        ColSideColors = ifelse(grepl('signif', colnames(Values)), 'steelblue', 'white'),
          RowSideColors = c('steelblue','seagreen','violet','indianred','purple','lightblue')[factor(cutClust)],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")

dev.off()


# Testting survival according to sub-grps
require(survival)
Time <- samples$RecFreeSurv_month
status <- samples$Recurrence
st <- Surv(Time, status)
K = 6
cutClust <- cutree(clust, k = K)
model <- coxph(st ~ cutClust)
summary(model)
km <- survfit(st ~ cutClust)
plot(km)

require(rbsurv)
time <- as.integer(samples$RecFreeSurv_month)
status <- samples$Recurrence
z <- cbind(samples$Age, samples$TumourType)
z <- cbind(factor(cutClust), samples$TumourType)
fit <- rbsurv(time = time, status = status, x = t(Values), z = z, method = "efron", max.n.genes = 10, n.iter = 2)
fit$model

#######################################################################

X = scale(PCA$x[,2:3])
Dist = rowSums(X^2)
whiteDist <- Dist[whiteProbesIndex]
qDist = quantile(whiteDist, 0.95)
theta <- seq(0, 90, by = 1)
x1 <- sqrt(qDist)*cos(theta)
y1 <- sqrt(qDist)*sin(theta)
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/whiteProbes_boundary.ps'), 
           width = 7, height = 7)
plot(X, cex = 0.2, col = ifelse(grepl('whiteProbe', rownames(eset)), 'red', 'grey80'))
points(y1~x1, pch = 19, col = "darkblue", cex = 0.5)
dev.off()

good <- which(Dist > qDist)
original <- prcomp(t(eset))
pcaGood <- prcomp(t(eset[good,]))
pcaBad <- prcomp(t(eset[-good,]))
#Cols = ifelse(samples$TumourType == 'Adk', rgb(0.15, 0, 0.75, 0.8), rgb(0.75, 0, 0.25, 0.8))
Cols = ifelse(samples$TumourType == 'Adk', 'orangered', 'darkblue')

postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/whiteProbes_based_selection.ps'), 
           width = 7, height = 7)
par(mfrow = c(2, 2))
plotPCA(original, pch = 19, col = Cols, main = 'Full dataset')
plot(X, cex = 0.2, col = ifelse(grepl('whiteProbe', rownames(eset)), 'red', 'grey80'))
points(y1~x1, pch = 19, col = "darkblue", cex = 0.5)
legend('topleft', legend = c('white probes', 'q95(dist)'), pch = 19,
       col = c('red', 'blue'), cex = 1.25, bty = 'n')
plotPCA(pcaGood, pch = 19, col = Cols, main = paste('Probes over white distance\nn =', length(good)))
plotPCA(pcaBad, xlim = range(pcaGood$x[,1]), ylim = range(pcaGood$x[,2]), pch = 19,
        col = Cols, main = paste('Probes inside white distance\nn =', nrow(eset) - length(good)))
par(op)
dev.off()

# where are the filtering probes
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/Sd_Means_whiteProbesLoc.ps'), 
           width = 7, height = 7)
plot(M, log10(S), pch = 19, cex = 0.2, col = ifelse(Dist[-whiteProbesIndex]>qDist, 'grey', 'red3'),
     xlab = 'Means', ylab = 'Log10(Sdev)')
dev.off()

# Compute the loadings
SVD <- svd(t(eset))
V <- SVD$v[,2:3]^2

# the filtering method
PCA <- prcomp(eset)
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/informationCurve.ps'), 
           width = 7, height = 7)
score <- pcaTrace(eset, PCA, lwd = 5)
dev.off()
Info <- pcaInfo(score); Info

postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/PCAonSelections.ps'), 
           width = 8, height = 8)
par(mfrow = c(3, 2))
redDot = 10^(1/score$lModel$x.intercept)
for(p in c(redDot, 0.1, 0.5)){
  select <- pcaSelect(score, p); length(select)
  final <- prcomp(t(eset[select,]))
  unselected <- prcomp(t(eset[-select,]))
  plotPCA(final, pch = 19, col = Cols, main = paste(length(select),'informative probes'))
  plotPCA(unselected, pch = 19, col = Cols, xlim = range(final$x[,1]), ylim = range(final$x[,2]),
          main = paste(nrow(eset)-length(select),'rejected probes'))  
}
par(op)
dev.off()


PC3 <- final$x[,3]
mp = min(PC3)
Mp = max(PC3)
pcex = 2*(PC3-mp)/(Mp - mp) + 0.5
plot(final$x, cex = pcex,
     pch = 19, col = rgb(0.5, 0, 0.5, 0.5), main = paste(length(select),'informative probes'))

PC3 <- unselected$x[,3]
mp = min(PC3)
Mp = max(PC3)
pcex = 2*(PC3-mp)/(Mp - mp) + 0.5
plot(unselected$x, cex = pcex, xlim = range(final$x[,1]), ylim = range(final$x[,2]),
     pch = 19, col = rgb(0.5, 0, 0.5, 0.5), main = paste(nrow(eset)-length(select),'rejected probes'))

plot(density(final$x[,1]))
pairs(final$x[,1:5], cex = pcex, col = Cols,
      pch = c(1, 19)[samples$Recurrence+1])

graph3D.8(final, class1 = samples$TumourType, class2 = samples$Recurrence)

require(survival)
Time <- samples$RecFreeSurv_month
status <- samples$Recurrence
st <- Surv(Time, status)
select <- pcaSelect(score, 0.05)
length(select)
subEset <- eset[select,]
dim(subEset)

require(rbsurv)

x <- as.matrix(eset[1:100,])
x <- matrix(rnorm(100*138, 10, 3), 100, 138); dim(x)
rownames(x) <- paste0('g', seq(1, nrow(x)))
colnames(x) <- paste0('Chip', seq(1, ncol(x)))
time <- as.integer(samples$RecFreeSurv_month)
status <- samples$Recurrence
z <- cbind(samples$Age, samples$TumourType)
fit <- rbsurv(time = time, status = status, x = x, method = "efron", max.n.genes = 10, n.iter = 5)
fit$model

x <- eset[select, ]
reord <- sample(1:ncol(x))
x <- x[,reord]
time <- samples$RecFreeSurv_month[reord]
status <- samples$Recurrence[reord]
st <- Surv(Time, status)

P <- c()
for(i in 1:nrow(x)){
  g <- as.numeric(x[i,])
  g <- ifelse(g > median(g), 'high', 'low')
  model <- coxph(st ~ factor(g))
  P <- rbind(P, cbind(index = i, geneName = rownames(x)[i], p = summary(model)$logtest[3]))
}
P <- as.data.frame(P)
P$index <- as.numeric(as.vector(P$index))
P$p <- as.numeric(as.vector(P$p))
P <- P[order(P$p, decreasing = FALSE),]
head(P)
length(which(P$p<1e-3))

Values <- t(x[P$index[P$p<1e-3],])
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "red3"),
          breaks = seq(-3, 3, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          #RowSideColors = c('darkblue', 'orangered')[factor(cutClust)],
          key = TRUE, keysize = 1, trace = "none", density.info = "none")

for(i in P$index[1:5]){
  g <- as.numeric(x[i,])
  g <- factor(ifelse(g > median(g), 'high', 'low'))
  model <- coxph(st ~ relevel(g, ref = 'low'))
  km <- survfit(st ~ factor(g))
  plot(km, main = rownames(x)[i], col = c('darkblue', 'orangered'), lwd = 4)
  legend('topright', legend = c('low expr.', 'high expr.'), lwd = 3, col = c('darkblue', 'orangered'))
  legend('bottomleft', legend = paste('OneByOne', 'p =', signif(summary(model)$logtest[3], 3)))
}


fit <- rbsurv(time=time, status=status, x=x, method = "efron", max.n.genes = 50, n.iter = 10, n.fold = 1, n.seq=1)
fit$mode
bestGenes <- as.numeric(fit$mode$Gene[grepl('\\*', fit$mode$Selected)])

Values <- t(x[bestGenes,])
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "red3"),
          breaks = seq(-3, 3, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          #RowSideColors = c('darkblue', 'orangered')[factor(cutClust)],
          key = TRUE, keysize = 1, trace = "none", density.info = "none")

for(i in bestGenes[1:5]){
  g <- as.numeric(x[i,])
  g <- factor(ifelse(g > median(g), 'high', 'low'))
  model <- coxph(st ~ relevel(g, ref = 'low'))
  km <- survfit(st ~ factor(g))
  plot(km, main = rownames(x)[i], col = c('darkblue', 'orangered'), lwd = 4)
  legend('topright', legend = c('low expr.', 'high expr.'), lwd = 3, col = c('darkblue', 'orangered'))
  legend('bottomleft', legend = paste('rbsurv', 'p =', signif(summary(model)$logtest[3], 3)))
}

clust <- hclust(dist(Values), 'ward')
cutClust <- cutree(clust, k = 2)
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "red3"),
          breaks = seq(-3, 3, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          RowSideColors = c('darkblue', 'orangered')[factor(cutClust)],
          #cexCol = 1.25, labCol = paste(annotId, annotType, sep = "-"),
          #ColSideColors = colSide[factor(colAnnot)],
          key = TRUE, keysize = 1, trace = "none", density.info = "none")

model <- coxph(st ~ factor(cutClust))
km <- survfit(st ~ factor(cutClust))
plot(km, main = paste(length(bestGenes), 'genes from rbsurv'), col = c('darkblue', 'orangered'), lwd = 4)
legend('topright', legend = c('grp1', 'grp2'), lwd = 3, col = c('darkblue', 'orangered'))
legend('bottomleft', legend = paste('rbsurv clusters', 'p =', signif(summary(model)$logtest[3], 3)))

x <- exprs(gliomaSet)
x <- log2(x)
time <- gliomaSet$Time
status <- gliomaSet$Status
z <- cbind(gliomaSet$Age, gliomaSet$Gender)
fit <- rbsurv(time=time, status=status, x=x, method="efron", max.n.genes=20)
fit$mode
bestGenes <- as.numeric(fit$mode$Gene[grepl('\\*', fit$mode$Selected)][-1])


newProbes <- rownames(eset)[select]
keepProbes = c()
iter = 1
while (abs(length(newProbes) - length(keepProbes)) > 10){
  keepProbes = newProbes
  newProbes <- .iterSurv(eset[which(rownames(eset) %in% keepProbes),])
  cat(iter, length(newProbes), '\n')
  iter = iter + 1
}


.insertRow <- function(DF, newrow, r) {
  while(length(r)>1){
    cat(r, '\t', dim(DF), '\n')
    DF <- rbind(DF[1:(r[1]-1),], newrow[1,], DF[r[1]:nrow(DF),])
    newrow  = newrow[-1,]
    r = r[-1]
  }
  DF <- rbind(DF[1:(r[1]-1),], as.numeric(newrow), DF[r[1]:nrow(DF),])
  return(DF)  
}

.iterSurv <- function(Data){
  k = 1
  checkedProbes = keepProbes = c()
  while(length(checkedProbes) < nrow(Data)){
    geneSamp <- sample(1:nrow(Data), ceiling(ncol(Data)/3))
    checkedProbes <- unique(c(checkedProbes, rownames(Data)[geneSamp]))
    tmp = as.data.frame(t(Data[geneSamp,]))
    model <- coxph(st ~ ., data = tmp)
    tmpCoefs <- summary(model)$coefficients
    keepRows <- which(tmpCoefs[,5] < 0.2)
    keepProbes <- c(keepProbes, gsub('`', '', names(keepRows)))
    keepProbes = unique(keepProbes)
    # cat(k, 'iteration(s), nProbes', length(checkedProbes), 'of', nrow(Data),'\n')
    k = k+1
  }
  return(keepProbes)  
}

keepProbes <- .iterSurv(subEset); length(keepProbes)
keepProbes <- .iterSurv(eset[which(rownames(eset) %in% keepProbes),]); length(keepProbes)
