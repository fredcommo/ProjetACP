source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaSelect.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaInfo.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace1.1.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/simulationFunctions.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/graph3D.8.R')
require(multtest)
require(glmnet)

op <- par(no.readonly = TRUE)

kimData = readRDS('/Users/fredcommo/Documents/MyProjects/Kim_Lung/kimData.rds')
eset = kimData$eset
samples <- kimData$samples

# Add nonspecific probes
Random <- .generateRandom1(eset, floor(nrow(eset)*.1))
colnames(Random) <- colnames(eset)
eset <- rbind(eset, Random)

# Add specific probes
M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)
signif <- .generateGrps(M, S, ncol(eset), nrow(eset)*.01,
                                      grps = samples$TumourType, minP = 0.5, maxP = 0.8)
signif <- signif$Data
colnames(signif) <- colnames(eset)  
eset <- rbind(eset, signif)

# PCA Filtering
pcaProbes <- prcomp(eset)
score <- pcaTrace1.1(eset, pcaProbes, main = 'Information curve')
Info <- pcaInfo(score); Info

# Visualize filtered PCA at different cutoffs
Cols = ifelse(samples$TumourType == 'Adk', 'orangered', 'darkblue')
par(mfrow = c(3, 2))
redDot = 10^(1/score$lModel$x.intercept)
for(p in c(redDot, 0.1, 0.5)){
     select <- pcaSelect(score, p)
     n <- length(select)
     pcaS <- prcomp(t(eset[select,]))
     pcaR <- prcomp(t(eset[-select,]))
     plotPCA(pcaS, pch = 19, col = Cols, main = paste(n,'informative probes'))
     plotPCA(pcaR, pch = 19, col = Cols, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
                         main = paste(nrow(eset)-n,'rejected probes'))  
}
par(op)

# Visualize the construction
  # S vs. M distribution
trueDist <- .computeBounds(M, S)
boxplot(log10(trueDist$S)~trueDist$M, names = round(unique(trueDist$M), 2), 
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)')

  # Visualize the added probes
randomM <- apply(Random, 1, mean, na.rm = TRUE)
randomS <- apply(Random, 1, sd, na.rm = TRUE)
signifM <- apply(signif, 1, mean, na.rm = TRUE)
signifS <- apply(signif, 1, sd, na.rm = TRUE)

par(cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
smoothScatter(M, log10(S))
points(randomM, log10(randomS), pch = 19, cex = 0.5, col = 'grey75')
points(signifM, log10(signifS), pch = 19, cex = 0.5, col = 'red3')
par(op)

probeClass <- ifelse(grepl('random',rownames(eset)), 'grey75',
                      ifelse(grepl('signif',rownames(eset)), 'red3', 'steelblue1'))
pairs(pcaProbes$x[,1:3], col = probeClass, cex = 0.2)

# Compute PCAs according to S cutoff: check out to filter on variances
# Count the remaining white probes
par(mfrow = c(3, 2))
for(s in c(-1, -0.5, 0)){
  select <- as.numeric(which(log10(S) > s))
  nWhite <- length(intersect(grep('random', rownames(eset)), select))
  pca1 <- prcomp(t(eset[select, ]))
  pca2 <- prcomp(t(eset[-select,]))
  plotPCA(pca1, pch = 19, col = Cols,
          main = paste(length(select), 'selected probes, cutoff log10(Sdev) =', s,
                       '\n', nWhite, ' remaining whiteProbes'))
  plotPCA(pca2, pch = 19, xlim = range(pca1$x[,1]), ylim = range(pca1$x[,2]),
          col = Cols,
          main = paste(nrow(eset)-length(select), 'rejected probes, cutoff log10(Sdev) =', s))
}
par(op)

# Perform multtest on S-filtering: test the remaining probes
select <- as.numeric(which(log10(S) > (-1)))
Stest <- mt.maxT(eset[-select,], classlabel = samples$TumourType, B = 5000)
bestS <- Stest$index[Stest$adjp<0.001]
length(bestS); nrow(eset) - length(select)
length(intersect(grep('random', rownames(eset)[-select]), bestS))
length(intersect(grep('signif', rownames(eset)[-select]), bestS))
                     
# Compare multtest Adk Vs. Sq before and after PCA filtering
fullTest <- mt.maxT(eset, classlabel = samples$TumourType, B = 5000)
bestFull <- fullTest$index[fullTest$adjp<0.001]
length(bestFull)
length(intersect(grep('random', rownames(eset)), bestFull))
length(intersect(grep('signif', rownames(eset)), bestFull))


select <- pcaSelect(score, 0.05)
filtTest <- mt.maxT(eset[select,], classlabel = samples$TumourType, B = 5000)
bestFilt <- filtTest$index[filtTest$adjp<0.001]
length(bestFilt)
length(intersect(grep('random', rownames(eset[select,])), bestFilt))
length(intersect(grep('signif', rownames(eset[select,])), bestFilt))

# Compare classifiers
select <- pcaSelect(score, 0.05)
idx <- sample(1:ncol(eset), ncol(eset)*.7)
train <- eset[,idx]
grpTrain <- samples$TumourType[idx]
test <- eset[,-idx]
grpTest <- samples$TumourType[-idx]

  # Model on full matrix
fullModel <- cv.glmnet(x = t(train), y = grpTrain, family = 'binomial', alpha = 0.01)
fullL <- fullModel$lambda.min
fullB <- fullModel$glmnet.fit$beta[,fullModel$lambda == fullL]
sum(fullB != 0)
fullFit <- predict(fullModel, newx = t(test), s = fullL, type = 'class')
table(fullFit, grpTest)

    # Model on filtered matrix
filtModel <- cv.glmnet(x = t(train[select,]), y = grpTrain, family = 'binomial', alpha = 0.01)
filtModel$cvm
filtL <- filtModel$lambda.min
filtB <- filtModel$glmnet.fit$beta[,filtModel$lambda == filtL]
sum(filtB != 0)
filtFit <- predict(filtModel, newx = t(test[select,]), s = filtL, type = 'class')
table(filtFit, grpTest)
