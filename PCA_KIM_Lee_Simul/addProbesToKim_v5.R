source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaSelect.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaInfo.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace1.1.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/simulationFunctions.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/graph3D.8.R')
require(multtest)
require(glmnet)
require(VennDiagram)
require(mclust)
require(graphics)

op <- par(no.readonly = TRUE)

ent <- synGet('syn1898672')
source(ent@filePath)

ent <- synGet('syn1898673')
kimData <- readRDS(ent@filePath)
eset = kimData$eset
samples <- kimData$samples

setwd('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_KIM_Lee_Simul/')

# Original M & S values
oriM <- apply(eset, 1, mean, na.rm = TRUE)
oriS <- apply(eset, 1, sd, na.rm = TRUE)

# Add nonspecific probes
Random <- .generateRandom1(eset, floor(nrow(eset)*.1))
colnames(Random) <- colnames(eset)
eset <- rbind(eset, Random)

# Add specific probes
signif <- .generateGrps(oriM, oriS, ncol(eset), nrow(eset)*.01,
                                      grps = samples$TumourType, minP = 0.5, maxP = 0.8)
signif <- signif$Data
colnames(signif) <- colnames(eset)  
eset <- rbind(eset, signif)

# final M & S values
finalM <- apply(eset, 1, mean, na.rm = TRUE)
finalS <- apply(eset, 1, sd, na.rm = TRUE)

# PCA Filtering
pcaProbes <- prcomp(eset)

png('Lee_informationCurve.png', width = 600, height = 400)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
score <- pcaTrace1.1(eset, pcaProbes, main = 'Information curve', lwd = 4)
par(op)
dev.off()

Info <- pcaInfo(score); Info

# Visualize filtered PCA at different cutoffs
Cols = ifelse(samples$TumourType == 'Adk', 'orangered', 'darkblue')
png('Lee_Filtering_effect_diffCuts.png', width = 600, height = 600)
par(mfrow = c(3, 2), mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
redDot = 10^(1/score$lModel$x.intercept)
for(p in c(0.05, 0.1, 0.25)){
     select <- pcaSelect(score, p)
     n <- length(select)
     pcaS <- prcomp(t(eset[select,]))
     pcaR <- prcomp(t(eset[-select,]))
     plotPCA(pcaS, pch = 19, col = Cols, main = paste(n,'informative probes'))
     plotPCA(pcaR, pch = 19, col = Cols, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
                         main = paste(nrow(eset)-n,'rejected probes'))  
}
par(op)
dev.off()

# Visualize the construction
  # S vs. M distribution
trueDist <- .computeBounds(oriM, oriS)
png('Lee_MSplot.png', width = 800, height = 600)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
boxplot(log10(trueDist$S)~trueDist$M, names = round(unique(trueDist$M), 2), 
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)', main = 'MS plot')
par(op)
dev.off()

# Visualize the added probes
randomM <- apply(Random, 1, mean, na.rm = TRUE)
randomS <- apply(Random, 1, sd, na.rm = TRUE)
signifM <- apply(signif, 1, mean, na.rm = TRUE)
signifS <- apply(signif, 1, sd, na.rm = TRUE)

png('Lee_ScatPLot_addedProbes.png', width = 800, height = 600)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
smoothScatter(oriM, log10(oriS), colramp = colorRampPalette(c("white", "blue4")),
              xlab = 'Means', ylab = 'Log10(sDev)', main = 'Added probes')
points(randomM, log10(randomS), pch = 19, cex = 0.5, col = 'grey75')
points(signifM, log10(signifS), pch = 19, cex = 0.5, col = 'red3')
legend('bottomright', legend = c('original', 'random added', 'specific added'),
       pch = 19, col = c('purple', 'grey70', 'red3'), cex = 1.5, bty = 'n')
par(op)
dev.off()

probeClass <- ifelse(grepl('random',rownames(eset)), 'grey75',
                      ifelse(grepl('signif',rownames(eset)), 'red3', rgb(0.5, 0.6, 0.9, 1)))
png('Lee_PairsPlot_addedProbes.png', width = 800, height = 600)
par(mar = c(5, 5, 4, 2.5), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.5)
pairs(pcaProbes$x[,1:3], col = probeClass, cex = 0.2)
par(op)
dev.off()

# Compute PCAs according to S cutoff: check how to filter on variances
# Count the remaining white probes
png('Lee_VarSelectionEffect.png', width = 800, height = 600)
par(mfrow = c(3, 3), mar = c(5, 5, 4, 2.5), cex.main = 1.75, cex.lab = 1.5, cex.axis = 1.25)
for(s in c(-1, -0.5, 0)){
  select <- which(log10(finalS) > s)
  cat('Rejected:', range(finalS[-select]), 'Selected:', range(finalS[select]), '\n')
  nWhite <- length(intersect(grep('random', rownames(eset)), select))
  pca1 <- prcomp(t(eset[select, ]))
  pca2 <- prcomp(t(eset[-select,]))
  smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "steelblue")))
  points(finalM[select], log10(finalS[select]), cex = 0.1, col ='steelblue3')
  points(finalM[-select], log10(finalS[-select]), cex = 0.1, col ='red3')
  legend('bottomright', legend = c('selected', 'rejected'),
         pch = 19, col =c('steelblue3', 'red3'), bty = 'n')
  plotPCA(pca1, pch = 19, col = Cols,
          main = paste(length(select), 'selected probes\ncutoff log10(Sdev) =', s))
  plotPCA(pca2, pch = 19, xlim = range(pca1$x[,1]), ylim = range(pca1$x[,2]),
          col = Cols,
          main = paste(nrow(eset)-length(select), 'rejected probes'))
}
par(op)
dev.off()


# distribution of distances, rejected probes: PCA (p = 0.05) Vs. random probes (same size)
# select <- pcaSelect(score, 0.05)
# n <- length(select)
# pcaR <- prcomp(t(eset[-select,]))
# pcaDist <- sum(pcaR$sdev^2)
# randomReject <- lapply(1:1000, function(i){
#   if(i%%10 == 0) cat(i, '\t')
#   Samp <- sample(1:nrow(eset), ceiling(nrow(eset)-n))
#   pcaR <- prcomp(t(eset[Samp,]))
#   pcaDist <- sum(pcaR$sdev^2)
# })
# randomReject <- do.call(c, randomReject)
# histogram(log(randomReject), nclass = 50)

  # Illustration
p = 0.05
select <- pcaSelect(score, p)
n <- length(select)
pcaS <- prcomp(t(eset[select,]))
pcaR <- prcomp(t(eset[-select,]))

png('Lee_SelectedByPCA.png', width = 800, height = 400)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              main = 'PCA-filtered probes')
legend('bottomright', legend = 'Filtered', pch = 19, col = 'red3', cex = 1.5, bty = 'n')
points(finalM[-select], log10(finalS[-select]), pch = 19, cex = 0.2, col = 'red')
plotPCA(pcaR, pch = 19, col = Cols, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
        main = paste(nrow(eset)-n,'rejected probes'))  
plotPCA(pcaS, pch = 19, col = Cols, main = paste(n,'informative probes'))
par(op)
dev.off()


# Random
png('Lee_SelectedRandom.png', width = 800, height = 400)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
Samp <- sample(1:nrow(eset), ceiling(nrow(eset)-n))
pcaRejectRand <- prcomp(t(eset[Samp,]))
pcaSelectRand <- prcomp(t(eset[-Samp,]))
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              main = 'Randomly-filtered probes')
points(finalM[Samp], log10(finalS[Samp]), pch = 19, cex = 0.2, col = 'red')
legend('bottomright', legend = 'Filtered', pch = 19, col = 'red3', cex = 1.5, bty = 'n')
plotPCA(pcaRejectRand, pch = 19, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
        col = Cols, main = paste(nrow(eset)-n, 'randomly rejected probes'))
plotPCA(pcaSelectRand, pch = 19, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
        col = Cols, main = paste(n, 'randomly selected probes'))
par(op)
dev.off()

  # on vars
png('Lee_SelectedByVar.png', width = 800, height = 400)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
q <- quantile(finalS, probs = (nrow(eset)-n)/nrow(eset))
lowVar <- which(finalS <= q)
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              main = 'Variance-filtered probes')
points(finalM[lowVar], log10(finalS[lowVar]), pch = 19, cex = 0.2, col = 'red')
pcaSelectVar <- prcomp(t(eset[-lowVar,]))
pcaRejectVar <- prcomp(t(eset[lowVar,]))
plotPCA(pcaRejectVar, pch = 19, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
        col = Cols, main = paste(length(lowVar), 'variance-rejected probes'))
plotPCA(pcaSelectVar, pch = 19, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
        col = Cols, main = paste(nrow(eset)-length(lowVar), 'variance-selected probes'))
par(op)
dev.off()

# distribution of distances, rejected probes: PCA (p = 0.05) Vs. Variance (same size)

############################################################
# Multiple testing on real data
.mt <- function(Data, grp, B = 5000){
  n <- nrow(Data)
  mtTest <- mt.maxT(Data, classlabel = grp, B = B)
  mtTest <- mtTest[order(mtTest$index),]
  best <- mtTest$index[mtTest$adjp<0.001]
  nRandom <- sum(grepl('random', rownames(Data)))
  foundRandom <- sum(grepl('random', rownames(Data)[best]))
  nSignif <- sum(grepl('signif', rownames(Data)))
  foundSignif <- sum(grepl('signif', rownames(Data)[best]))
  cat('\n')
  cat('#Signature:', length(best), 'of', n, '\tprop:', round(length(best)/n, 4)*100, '\n')
  cat('#Random:', foundRandom, 'of', nRandom, '\tprop:', round(foundRandom/nRandom, 4)*100, '\n')
  cat('#Signif:', foundSignif, 'of', nSignif, '\tprop:', round(foundSignif/nSignif, 4)*100, '\n')
  return(mtTest)
}
grp <- samples$TumourType
fullTest <- .mt(eset, grp)

select <- pcaSelect(score, 0.05)
filtTest <- .mt(eset[select,], grp)

q <- quantile(finalS, probs = (nrow(eset)-length(select))/nrow(eset))
highVar <- which(finalS >= q)
varTest <- .mt(eset[highVar,], grp)


# Built a Venn diagram here with the 3 signatures.
thresh = 1e-3
  # Global overlap
vdGlobal <- list(full = rownames(fullTest)[which(fullTest$adjp<thresh)],
           byPCA = rownames(filtTest)[which(filtTest$adjp<thresh)],
           byVar = rownames(varTest)[which(varTest$adjp<thresh)])
venn.diagram(vdGlobal, fill = c("red", "green", "blue"),
             alpha = c(0.25, 0.25, 0.25), cex = 2, cat.fontface = 4, lty =2, fontfamily =3,
             filename = 'Lee_vennDiagram_global.png')
  # Added probes overlap
vdAdded <- list(full = rownames(fullTest)[which(fullTest$adjp<thresh & grepl('signif', rownames(fullTest)))],
                 byPCA = rownames(filtTest)[which(filtTest$adjp<thresh & grepl('signif', rownames(filtTest)))],
                 byVar = rownames(varTest)[which(varTest$adjp<thresh & grepl('signif', rownames(varTest)))])
venn.diagram(vdAdded, fill = c("red", "green", "blue"),
             alpha = c(0.25, 0.25, 0.25), cex = 2, cat.fontface = 4, lty =2, fontfamily =3,
             filename = 'Lee_vennDiagram_added.png')

# Show the differences: added signif probes found in bestFilt but not in bestVar
full <- (fullTest$adjp<thresh)*1
byPCA <- is.element(rownames(fullTest), rownames(filtTest))*1
byPCA[byPCA == 1] <- (filtTest$adjp<thresh)*1
byVar <- is.element(rownames(fullTest), rownames(varTest))*1
byVar[byVar == 1] <- (varTest$adjp<thresh)*1

.modelGlm <- function(idx){
  x <- as.numeric(eset[idx,])
  y <- c(0,1)[samples$TumourType]
  model <- glm(y ~ x, family = binomial)
  pvalue <- summary(model)$coefficients[2,4]
  return(pvalue)
}
fullOnlyIdx <- which(full==1 & byPCA+byVar == 0); length(fullOnlyIdx)
fullOnly <- cbind.data.frame(probeId = rownames(eset)[fullOnlyIdx],
                             pvalue = do.call(c, lapply(fullOnlyIdx, function(idx){.modelGlm(idx)})))
pcaOnlyIdx <- which(byPCA==1 & full+byVar == 0); length(pcaOnlyIdx)
pcaOnly <- cbind.data.frame(probeId = rownames(eset)[pcaOnlyIdx],
                             pvalue = do.call(c, lapply(pcaOnlyIdx, function(idx){.modelGlm(idx)})))
pcaNotVarIdx <- which(byPCA==1 & byVar == 0); length(pcaNotVarIdx)
pcaNotVar <- cbind.data.frame(probeId = rownames(eset)[pcaNotVarIdx],
                               pvalue = do.call(c, lapply(pcaNotVarIdx, function(idx){.modelGlm(idx)})))
fullAndPcaIdx <- which(byVar==0 & full*byPCA == 1); length(fullAndPcaIdx)
fullAndPca <- cbind.data.frame(probeId = rownames(eset)[fullAndPcaIdx],
                               pvalue = do.call(c, lapply(fullAndPcaIdx, function(idx){.modelGlm(idx)})))

# # Visualize the fullOnly
# y <- c(0,1)[samples$TumourType]
# par(mfrow = c(2, 2))
# for(i in fullOnlyIdx[1:4]){
#   x <- as.numeric(eset[i,])
# #  boxplot(x~y)
#    model <- glm(y ~ x, family = binomial)
#    newx <- seq(min(x), max(x), len = 100)
#   fit <- predict(model, newdata = list(x = newx), type = 'response')
#   plot(x, y, col = Cols)
#   lines(newx, fit, lwd = 3)
# }
# par(op)

# Show the signatures on the MS plot
png('Lee_Signatures_on_MSplot.png', width = 800, height = 600)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              xlab = 'Means', ylab = 'Log10(sDev)', main = 'Full GE signature')
points(finalM[full==1], log10(finalS[full==1]), pch = 19, cex = 0.2, col = 'green3')
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              xlab = 'Means', ylab = 'Log10(sDev)', main = 'PCA filtering signature')
points(finalM[byPCA==1], log10(finalS[byPCA==1]), pch = 19, cex = 0.2, col = 'green3')
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              xlab = 'Means', ylab = 'Log10(sDev)', main = 'Variance filtering signature')
points(finalM[byVar==1], log10(finalS[byVar==1]), pch = 19, cex = 0.2, col = 'green3')
smoothScatter(finalM, log10(finalS), colramp = colorRampPalette(c("white", "blue4")),
              xlab = 'Means', ylab = 'Log10(sDev)', main = 'On PCA-filering only')
points(finalM[pcaOnlyIdx], log10(finalS[pcaOnlyIdx]), pch = 19, cex = 0.2, col = 'green3')
par(op)
dev.off()

# Show the signatures on the pcaProbes plot
png('Lee_Signatures_on_pcaProbes.png', width = 1000, height = 400)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
plot(pcaProbes$x[,2:3], col = rgb(0.5, 0.6, 0.9, 0.2), cex = 0.2, main = 'Full GE signature')
points(pcaProbes$x[grepl('random', rownames(eset)), 2:3], pch = 19, cex = 0.2, col = 'grey75')
points(pcaProbes$x[full==1, 2:3], pch = 19, cex = 0.2, col = 'darkblue')

plot(pcaProbes$x[,2:3], col = rgb(0.5, 0.6, 0.9, 0.2), cex = 0.2, main = 'PCA-filtered signature')
points(pcaProbes$x[grepl('random', rownames(eset)), 2:3], pch = 19, cex = 0.2, col = 'grey75')
points(pcaProbes$x[byPCA==1, 2:3], pch = 19, cex = 0.2, col = 'darkblue')

plot(pcaProbes$x[,2:3], col = rgb(0.5, 0.6, 0.9, 0.2), cex = 0.2, main = 'Var-filtered signature')
points(pcaProbes$x[grepl('random', rownames(eset)), 2:3], pch = 19, cex = 0.2, col = 'grey75')
points(pcaProbes$x[byVar==1, 2:3], pch = 19, cex = 0.2, col = 'darkblue')
par(op)
dev.off()

# Variances densities:
png('Lee_VarDensities.png', width = 800, height = 600)
par(mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
dS <- density(finalS)
plot(dS, lwd = 2, xlim = range(min(dS$x), 3.5), main = 'Variance densities')
polygon(density(randomS), lwd = 2, col = rgb(0.25, 0.25, 0.25, 0.2))
polygon(density(signifS), lwd = 2, col = rgb(1, 0, 0, 0.25))
polygon(density(finalS[pcaSelect(score, 0.05)]), lwd = 2, col = rgb(0, 0, 1, 0.25))
polygon(density(finalS[highVar]), lwd = 2, col = rgb(0, 0.5, 1, 0.25))
par(op)
dev.off()

# Use gaussian mixture model on Var
png('Lee_gaussianMixture.png', width = 1000, height = 600)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
sampS <- finalS
model <- Mclust(sampS)
n = length(sampS); nG = model$G; m <- model$parameters$mean; p <- model$parameters$pro
s <- sqrt(model$parameters$variance$sigmasq); if(length(s)==1) s <- rep(s, nG)
dS <-density(sampS)
plot(dS, lwd = 2, xlab = 'Variance', main = 'Gaussian mixture on Var')
for (i in 1:nG){
  tmp <- rnorm(n*p[i], m[i], s[i])
  tmpD <- density(tmp, na.rm = T)
  polygon(tmpD$x, tmpD$y*p[i], col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
}
# Use gaussian mixture model on Log10(Var)
sampS <- log(finalS)
model <- Mclust(sampS)
n = length(sampS); nG = model$G; m <- model$parameters$mean; p <- model$parameters$pro
s <- sqrt(model$parameters$variance$sigmasq); if(length(s)==1) s <- rep(s, nG)
dS <-density(sampS)
plot(dS, lwd = 2, xlab = 'Log10(Variance)', main = 'Gaussian mixture on Log10(Var)')
for (i in 1:nG){
  tmp <- rnorm(n*p[i], m[i], s[i])
  tmpD <- density(tmp, na.rm = T)
  polygon(tmpD$x, tmpD$y*p[i], col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
}
par(op)
dev.off()

# Visualize the 4 best pcaNotVar: added probes
png('Lee_pcaNotVar_added.png', width = 800, height = 600)
par(mfrow = c(2, 2), las = 1, mar = c(5, 5, 4, 3.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
y <- c(0,1)[samples$TumourType]
k = 1
for(i in pcaNotVarIdx[order(pcaNotVar$pvalue, decreasing = FALSE)]){
  if(grepl('signif', rownames(eset)[i]) & k<5){
    x <- as.numeric(eset[i,])
    probeName <- rownames(eset)[i]
    pTtest <- t.test(x~y)$p.value
    model <- glm(y ~ x, family = binomial)
    pGlm <- summary(model)$coefficients[2,4]
    newx <- seq(min(x), max(x), len = 100)
    fit <- predict(model, newdata = list(x = newx), type = 'response')
    plot(x, y, cex = 1.25, col = Cols, xlim = range(quantile(x, probs = c(0.02, 1))),
         xlab = 'Log10(int)', ylab = 'Probability',
         main = paste0('probe ',probeName, ' (t.test p: ', signif(pTtest, 3),')'))
    lines(newx, fit, lwd = 3)
    legend('right', legend = paste('logist. reg.\np:',signif(pGlm,3)), cex = 1.25, bty = 'n')
    k = k+1
    axis(side = 4, at = c(0,1), labels = c('Adk', 'Sq'))
  }
}
par(op)
dev.off()

# Visualize the 4 best pcaNotVar: real probes
png('Lee_pcaNotVar_real.png', width = 800, height = 600)
par(mfrow = c(2, 2), las = 1, mar = c(5, 5, 4, 3.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
y <- c(0,1)[samples$TumourType]
k = 1
for(i in pcaNotVarIdx[order(pcaNotVar$pvalue, decreasing = FALSE)]){
  if(!grepl('signif', rownames(eset)[i]) & k<5){
    x <- as.numeric(eset[i,])
    probeName <- rownames(eset)[i]
    pTtest <- t.test(x~y)$p.value
    model <- glm(y ~ x, family = binomial)
    pGlm <- summary(model)$coefficients[2,4]
    newx <- seq(min(x), max(x), len = 100)
    fit <- predict(model, newdata = list(x = newx), type = 'response')
    plot(x, y, cex = 1.25, col = Cols, xlim = range(quantile(x, probs = c(0.02, 1))),
         xlab = 'Log10(int)', ylab = 'Probability',
         main = paste0('probe ',probeName, ' (t.test p: ', signif(pTtest, 3),')'))
    lines(newx, fit, lwd = 3)
    legend('right', legend = paste('logist. reg.\np:',signif(pGlm,3)), cex = 1.25, bty = 'n')
    k = k+1
    axis(side = 4, at = c(0,1), labels = c('Adk', 'Sq'))
  }
}
par(op)
dev.off()

# Visualize the 4 best fullOnly
png('Lee_fullNotPCA.png', width = 800, height = 600)
par(mfrow = c(2, 2), las = 1, mar = c(5, 5, 4, 3.5), cex.main = 2, cex.lab = 1.75, cex.axis = 1.5)
y <- c(0,1)[samples$TumourType]
for(i in fullOnlyIdx[1:min(length(fullOnlyIdx), 4)]){
  x <- as.numeric(eset[i,])
  probeName <- rownames(eset)[i]
  pTtest <- t.test(x~y)$p.value
  model <- glm(y ~ x, family = binomial)
  pGlm <- summary(model)$coefficients[2,4]
  newx <- seq(min(x), max(x), len = 100)
  fit <- predict(model, newdata = list(x = newx), type = 'response')
  plot(x, y, cex = 1.25, col = Cols, xlim = range(quantile(x, probs = c(0.02, 1))),
       xlab = 'Log10(int)', ylab = 'Probability',
       main = paste0('probe ',probeName, ' (t.test p: ', signif(pTtest, 3),')'))
  lines(newx, fit, lwd = 3)
  legend('right', legend = paste('logist. reg.\np:',signif(pGlm,3)), cex = 1.25, bty = 'n')
  axis(side = 4, at = c(0,1), labels = c('Adk', 'Sq'))
}
par(op)
dev.off()

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
