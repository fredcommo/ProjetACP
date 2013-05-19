# Filtering on Kim Lung data
# Generates random probes by picking in the real probes, and randomly switch the samples

require(splines)
require(glmnet)
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaSelect.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaInfo.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace1.1.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/simulationFunctions.R')

op <- par(no.readonly = TRUE)


################################################################

kimData = readRDS('/Users/fredcommo/Documents/MyProjects/Kim_Lung/kimData.rds')
eset = kimData$eset
samples <- kimData$samples

M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)

whiteProbes <- .generateRandom1(eset, 2500)
#whiteProbes <- .generateRandom2(M, S, ncol(eset), 2500)
rownames(whiteProbes) <- paste0('random', seq(1, nrow(whiteProbes)))
colnames(whiteProbes) <- colnames(eset)
eset <- rbind(eset, whiteProbes)

eset = kimData$eset

pcaProbes <- prcomp(eset)
score <- pcaTrace1.1(eset, pcaProbes)
info <- pcaInfo(score);info
select <- pcaSelect(score, 0.05); length(select)
pcaS <- prcomp(t(eset[select,]))
pcaR <- prcomp(t(eset[-select]))

par(mfrow = c(1,2))
X <- abs(max(pcaS$x[,1]))
Y <- abs(max(pcaS$x[,2]))
plotPCA(pcaS)
plotPCA(pcaR, xlim = range(-X, X), ylim = range(-Y, Y))
par(op)

# Generate significant probes
simulateGrps <- .generateGrps(M, S, n = ncol(eset), p = ceiling(nrow(eset)*0.01), grps = samples$TumourType)
signifProbes <- simulateGrps$Data

# Add the white probes into the eset matrix
colnames(whiteProbes) <- colnames(signifProbes) <- colnames(eset)
eset <- rbind(eset, whiteProbes, signifProbes)

# Compute PCAs according to S cutoff
# Count the remaining white probes
par(mfrow = c(3, 2))
for(s in c(-0.5, -0.25, 0)){
  select <- as.numeric(which(log10(S) > s))
  nWhite <- length(intersect(whiteIdx, select))
  pca1 <- prcomp(t(eset[select, ]))
  pca2 <- prcomp(t(eset[-select,]))
  plotPCA(pca1, pch = 19,
          main = paste(length(select), 'selected probes, cutoff log10(Sdev) =', s,
                       '\n', nWhite, ' remaining whiteProbes'))
  plotPCA(pca2, pch = 19, xlim = range(pca1$x[,1]), ylim = range(pca1$x[,2]),
          main = paste(nrow(eset)-length(select), 'rejected probes, cutoff log10(Sdev) =', s))
}
par(op)

score <- pcaTrace1.1(eset, pcaProbes, main = 'Information curve')
infos <- pcaInfo(score)
redDot = 10^(1/score$lModel$x.intercept)

# Plot PCA filtering
par(mfrow = c(3,2), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
for(p in c(0.05, 0.1, 0.25)){
  select <- pcaSelect(score, p)
  n <- length(select)
  pcaOnSelect <- prcomp(t(eset[select,]))
  pcaOnReject <- prcomp(t(eset[-select,]))
  X <- max(abs(pcaOnSelect$x[,1]))
  Y <- max(abs(pcaOnSelect$x[,2]))
  plotPCA(pcaOnSelect, pch = 19, col = c('grey25', 'grey75')[samples$TumourType],
          main = paste('PCA on', n, 'selected probes\nfilter:', p))
  plotPCA(pcaOnReject, 19, col = c('grey25', 'grey75')[samples$TumourType],
          xlim = range(-X, X), ylim = range(-Y, Y),
          main = paste('PCA on', nrow(eset)-n, 'rejected probes\nfilter:', p))
}
par(op)




filtTest <- mt.maxT(eset[select,], classlabel = samples$TumourType, B = 5000)
filtIndex <- filtTest$adjp<0.01 & !grepl('signif', rownames(filtTest))
any(grepl('white', rownames(filtTest)[filtIndex]))
bestFilt <- rownames(filtTest)[filtTest$adjp<0.01 & !grepl('signif', rownames(filtTest))]
length(bestFilt)

png(file = paste0('/gluster/home/jcommo/Simulations/pca_filteredProbesLocation.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plot(pcaProbes$x[,2:3], cex = 0.1,
     col = ifelse(grepl('white', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'blue', 'grey')))
points(pcaProbes$x[bestFull,2:3], cex = 0.25, col = 'cyan')
points(pcaProbes$x[which(rownames(eset) %in% bestFilt),2:3], cex = 0.25, col = 'steelblue3')
dev.off()

# Heatmap signature
Values <- t(eset[rownames(eset) %in% bestFilt,])
png(file = paste0('/gluster/home/jcommo/Simulations/heatmap_filteredProbes_Significant.png'), width = 800, height = 800)
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          RowSideColors = c('grey25', 'grey75')[samples$TumourType],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

# Cut at a higher value and try to find clusters
source('/gluster/home/jcommo/Fred_Scripts/heatmap.3.R')
infos
select <- pcaSelect(score, 0.70)
length(select)
#Values <- t(eset[which(rownames(eset) %in% bestFilt),])
Values <- t(eset[select,])

png(file = paste0('/gluster/home/jcommo/Simulations/heatmap_filteredProbes_070.png'), width = 800, height = 800)
heatmap.3(Values,
          # Colv = NA,
          # dendrogram = "column",
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          #rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2), sepcolor = "grey85", sepwidth=c(0.001,0.001),
          #cexRow = 1, labRow = NA,
          ColSideColors = ifelse(grepl('signif', colnames(Values)), 'steelblue', 'white'),
          RowSideColors = c('grey25', 'grey75')[samples$TumourType],
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")
dev.off()

# Try to find clusters and compare them
clust <- hclust(dist(Values), 'ward')
plot(clust)
K = 6
cutClust <- cutree(clust, k = K)
png(file = paste0('/gluster/home/jcommo/Simulations/heatmap_filteredProbes_070_clusters.png'), width = 800, height = 800)
heatmap.3(Values,
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

clustGrp <- mt.maxT(t(Values), classlabel = samples$TumourType, B = 5000)
clustTest <- mt.maxT(t(Values), classlabel = c(0:(K-1))[cutClust], test = 'f', B = 5000)
sum(clustTest$adjp < 0.01)/nrow(clustTest)
bestClust <- rownames(clustTest)[clustTest$adjp < 0.01]
length(bestClust)

par(op)
png(file = paste0('/gluster/home/jcommo/Simulations/pca_subGrpsProbes.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plot(pcaProbes$x[,2:3], cex = 0.1,
     col = ifelse(grepl('white', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'blue', 'grey')))
points(pcaProbes$x[bestFull,2:3], cex = 0.25, col = 'cyan')
points(pcaProbes$x[which(rownames(eset) %in% bestFilt),2:3], cex = 0.25, col = 'steelblue3')
points(pcaProbes$x[select,2:3], cex = 0.5, col = 'seagreen4')
dev.off()

# Test the clsuter groups on the full matrix
clustTest <- mt.maxT(eset, classlabel = c(0:(K-1))[cutClust], test = 'f', B = 5000)
sum(clustTest$adjp < 0.01)/nrow(clustTest)
bestClust <- clustTest$index[clustTest$adjp < 0.01]
length(bestClust)

par(op)
png(file = paste0('/gluster/home/jcommo/Simulations/pca_subGrpsProbes.png'), width = 800, height = 800)
par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plot(pcaProbes$x[,2:3], cex = 0.1,
     col = ifelse(grepl('white', rownames(eset)), 'red', ifelse(grepl('signif', rownames(eset)), 'blue', 'grey')))
points(pcaProbes$x[bestFull,2:3], cex = 0.25, col = 'cyan')
points(pcaProbes$x[which(rownames(eset) %in% bestFilt),2:3], cex = 0.25, col = 'steelblue3')
points(pcaProbes$x[bestClust,2:3], pch = 19, cex = 0.5, col = 'seagreen4')
dev.off()

# Testing classification


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
