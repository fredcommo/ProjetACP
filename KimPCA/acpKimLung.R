# Filtering on Kim Lung data
op <- par(no.readonly = TRUE); dev.off()
kimData = readRDS('/gluster/home/jcommo/Kim_Lung/kimData.rds')
eset = kimData$eset
samples = kimData$samples

require(splines)
source('~/Fred_Scripts/pcaSelect.R')
source('~/Fred_Scripts/pcaInfo.R')
source('~/Fred_Scripts/pcaTrace.R')
source('~/Fred_Scripts/plotPCA.R')

# Estimate row means, and row Sd
M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)
ns1 <- ns(M, df = 3)
lm1 <- lm(log10(S) ~ ns1)
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/Sd_Means_distrib.ps'), 
           width = 7, height = 7)
plot(M, log10(S), pch = 19, cex = 0.2, col = 'grey', xlab = 'Means', ylab = 'Log10(Sdev)')
lines(sort(M), lm1$fitted[order(M)], col = 'red3', lwd = 4)
dev.off()

mCuts <- cut(M, breaks = seq(min(M), max(M), by = 0.2))
labs <- levels(mCuts)
mBounds <- cbind.data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
      upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
mBounds <- cbind.data.frame(lower = mBounds$lower, med = (mBounds$upper + mBounds$lower)/2, upper = mBounds$upper)
sBounds <- c()
for(i in 1:nrow(mBounds)){
  index <- which(M >= mBounds$lower[i] & M < mBounds$upper[i])
  sBounds <- rbind(sBounds, cbind(S = S[index], M = rep(mBounds$med[i], length(index))))
}
sBounds <- as.data.frame(sBounds)
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/Sd_Means_boxplot.ps'), 
           width = 7, height = 7)
boxplot(log10(sBounds$S) ~ factor(sBounds$M), names = unique(round(sBounds$M, 1)),
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)')
dev.off()

# Compute white probes
pickedValues <- c()
whiteProbes <- c()
for(i in 1:1000){
  m = as.numeric(sample(M, 1))
  j <- which(mBounds$l<=m & mBounds$upper>m)
  sIndex <- which(M >= mBounds$lower[j] & M < mBounds$upper[j])
  s = as.numeric(sample(S[sIndex], 1))
  tmp <- rnorm(ncol(eset), m, s)
  whiteProbes <- rbind(whiteProbes, tmp)
  rownames(whiteProbes)[nrow(whiteProbes)] = paste0('whiteProbe',i)
  pickedValues = rbind(pickedValues, c(m = m, s = s))
}
pickedValues <- as.data.frame(pickedValues)
ns2 <- ns(pickedValues$m, df = 3)
lm2 <- lm(log10(pickedValues$s) ~ ns2)
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/Sd_Means_PickedValues.ps'), 
           width = 7, height = 7)
plot(pickedValues$m, log10(pickedValues$s), pch = 19, cex = 0.75, col = 'grey')
lines(sort(pickedValues$m), lm2$fitted[order(pickedValues$m)], col = 'blue3', lwd = 4)
lines(sort(M), lm1$fitted[order(M)], col = 'red3', lwd = 4, lty = 2)
dev.off()

# Add the white probes into the eset matrix
colnames(whiteProbes) = colnames(eset)
#eset <- .insertRow(eset, whiteProbes, sample(1:nrow(eset), nrow(whiteProbes))) 
eset <- rbind(eset, whiteProbes)
whiteIdx <- grep('whiteProbe', rownames(eset))

# Recompute the Sdev vector, including white probes
S <- apply(eset, 1, sd, na.rm = TRUE)

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

# Show the white probes
PCA <- prcomp(eset)
whiteProbesIndex = which(grepl('whiteProbe', rownames(eset)))
postscript(file = paste0('/gluster/home/jcommo/Kim_Lung/whiteProbesLoc.ps'), 
           width = 7, height = 7)
pairs(PCA$x[,1:3], cex = 0.1, col = ifelse(grepl('whiteProbe', rownames(eset)), 'red', 'grey80'))
dev.off()

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

# Use multtest to find the diff expr. genes between tumorTypes, before and after filtering
require(multtest)
eset <- kimData$eset
samples <- kimData$samples
fullTest <- mt.maxT(eset, classlabel = samples$TumourType, B = 5000)

# PCA filtering
pcaProbes <- prcomp(eset)
score <- pcaTrace(eset, pcaProbes, main = 'Information curve')
infos <- pcaInfo(score)
redDot = 10^(1/score$lModel$x.intercept)
select <- pcaSelect(score, redDot)
filtTest <- mt.maxT(eset[select,], classlabel = samples$TumourType, B = 5000)

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
