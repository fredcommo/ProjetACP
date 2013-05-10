# Filtering ALL
source('~/Fred_Scripts/pcaSelect.R')
source('~/Fred_Scripts/pcaInfo.R')
source('~/Fred_Scripts/pcaTrace.R')

require(ALL)
data(ALL)
eset <- exprs(ALL)
samples <- pData(ALL)

# Estimate row means, and row Sd
M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)
postscript(file = paste0('/gluster/home/jcommo/ALL/Sd_Means_distrib.ps'), 
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
postscript(file = paste0('/gluster/home/jcommo/ALL/Sd_Means_boxplot.ps'), 
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
postscript(file = paste0('/gluster/home/jcommo/ALL/Sd_Means_PickedValues.ps'), 
           width = 7, height = 7)
plot(pickedValues$m, log10(pickedValues$s), pch = 19, cex = 0.75, col = 'grey')
lines(sort(pickedValues$m), lm2$fitted[order(pickedValues$m)], col = 'blue3', lwd = 4)
lines(sort(M), lm1$fitted[order(M)], col = 'red3', lwd = 4, lty = 2)
dev.off()

# Add the white probes into the eset matrix
colnames(whiteProbes) = colnames(eset)
# for(i in 1:nrow(whiteProbes)){
#   r <- sample(1:nrow(whiteProbes), 1)
#   eset <- .insertRow(eset, whiteProbes[r,], sample(1:nrow(eset), 1))
#   whiteProbes <- whiteProbes[-r,]
# }

eset <- rbind(eset, whiteProbes)

# Show the white probes
PCA <- prcomp(eset)
whiteProbesIndex = which(grepl('whiteProbe', rownames(eset)))
postscript(file = paste0('/gluster/home/jcommo/ALL/whiteProbesLoc.ps'), 
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
postscript(file = paste0('/gluster/home/jcommo/ALL/whiteProbes_boundary.ps'), 
           width = 7, height = 7)
plot(X, cex = 0.2, col = ifelse(grepl('whiteProbe', rownames(eset)), 'red', 'grey80'))
points(y1~x1, pch = 19, col = "darkblue", cex = 0.5)
dev.off()

good <- which(Dist > qDist)

original <- prcomp(t(eset))
pcaGood <- prcomp(t(eset[good,]))
pcaBad <- prcomp(t(eset[-good,]))
Cols = ifelse(grepl('T',samples$BT), rgb(0.15, 0, 0.75, 0.8), rgb(0.75, 0, 0.25, 0.8))
#Cols = c(2:(nlevels(samples$BT)+1))[samples$BT]

postscript(file = paste0('/gluster/home/jcommo/ALL/whiteProbes_based_selection.ps'), 
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
postscript(file = paste0('/gluster/home/jcommo/ALL/Sd_Means_whiteProbesLoc.ps'), 
           width = 7, height = 7)
plot(M, log10(S), pch = 19, cex = 0.2, col = ifelse(Dist[-whiteProbesIndex]>qDist, 'grey', 'red3'),
     xlab = 'Means', ylab = 'Log10(Sdev)')
dev.off()

# the filtering method
PCA <- prcomp(eset)
postscript(file = paste0('/gluster/home/jcommo/ALL/informationCurve.ps'), 
           width = 7, height = 7)
score <- pcaTrace(eset, PCA, lwd = 5)
dev.off()
Info <- pcaInfo(score); Info

BT <- factor(ifelse(grepl('B', samples$BT), 'B', 'T'))
#Cols = c(2:(nlevels(samples$BT)+1))[samples$BT]
Cols = c(2:(nlevels(BT)+1))[BT]

postscript(file = paste0('/gluster/home/jcommo/ALL/PCAonSelections.ps'), 
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
  legend('bottomleft', legend = levels(BT), col = 2:3, pch = 19)
}
par(op)
dev.off()
