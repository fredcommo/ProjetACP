# PCA filtering on Chemores miR
# v0, 2013_05_17
# FC

###################################
# Build RDS file
# miR T & N, paired
# Path <- '/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_Chemores_miR_simul/'
# data.address <- paste0(Path, "dataCombined_NQ_noCtrl_noFlag_redMir_NA15_KNN_copie_v4 (hsa).txt")
# info.address <- paste0(Path, "Chemores_miR_Samples.txt")
# eset <- read.table(data.address, header = T, sep = "\t")
# infos <- read.table(info.address, header = T, sep = "\t")
# infos <- infos[order(infos$FileName),]
# miRs <- eset[,1:2]
# eset <- eset[,-c(1:2)]
# eset <- eset[,order(colnames(eset))]
# all(colnames(eset) == infos$FileName)
# 
# colnames(eset) <- infos$BareCode
# 
# chemores_miR <- list(eset = eset, annot = miRs, samples = infos)
# saveRDS(chemores_miR, paste0(Path, 'chemores_miR.rds'))
###################################

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

setwd('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_Chemores_miR_simul/')
Chem <- readRDS('chemores_miR.rds')
eset <- Chem$eset
samples <- Chem$samples

# Filter on tumor cell prop
idx <- which(samples$TumCellPercent>=50)
samples <- samples[idx,]
eset <- eset[,idx]

# Consider tumors only
idx <- which(samples$Status == 'Tumor')
samples <- samples[idx,]
eset <- eset[,idx]

# Add nonspecific probes
Random <- .generateRandom1(eset, nrow(eset)*.25)
colnames(Random) <- colnames(eset)
eset <- rbind(eset, Random)

# Add specific probes
M <- apply(eset, 1, mean, na.rm = TRUE)
S <- apply(eset, 1, sd, na.rm = TRUE)

grps <- samples$Status
addTN <- .generateGrps(M, S, ncol(eset), nrow(eset)*.01, grps = grps, minP = 0.5, maxP = 0.8)
tnProbes <- addTN$Data
colnames(tnProbes) <- colnames(eset)
rownames(tnProbes) <- paste0('status', seq(1, nrow(tnProbes)))

grps <- ifelse(grepl('AC|SCC', samples$Disease), as.character(samples$Disease), NA)
addDis <- .generateGrps(M, S, ncol(eset), nrow(eset)*.01, grps = grps, minP = 0.5, maxP = 0.8)
disProbes <- addDis$Data
colnames(disProbes) <- colnames(eset)
rownames(disProbes) <- paste0('disease', seq(1, nrow(disProbes)))

eset <- rbind(eset, tnProbes, disProbes)

# Visualize the construction
  # S vs. M distribution
trueDist <- .computeBounds(M, S)
par(mar = c(5, 6, 4, 2), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
boxplot(log10(trueDist$S) ~ factor(trueDist$M), names = round(unique(trueDist$M), 2), 
        outpch = 19, outcex = 0.25, col = 'steelblue2', outcol = 'steelblue4',
        xlab = 'Means', ylab = 'Log10(Sdev)', main = 'Sdev Vs. mean distribution')
par(op)

  # Visualize the added probes
randomM <- apply(Random, 1, mean, na.rm = TRUE)
randomS <- apply(Random, 1, sd, na.rm = TRUE)
tnProbesM <- apply(tnProbes, 1, mean, na.rm = TRUE)
tnProbesS <- apply(tnProbes, 1, sd, na.rm = TRUE)
disProbesM <- apply(disProbes, 1, mean, na.rm = TRUE)
disProbesS <- apply(disProbes, 1, sd, na.rm = TRUE)

par(mar = c(5, 6, 4, 2), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
smoothScatter(M, log10(S), ylim = range(-1.5, 1),
              colramp = colorRampPalette(c("white", "steelblue")), main = 'Added probes')
points(randomM, log10(randomS), pch = 19, cex = 1, col = 'grey75')
points(randomM, log10(randomS))
points(tnProbesM, log10(tnProbesS), pch = 19, cex = 1, col = 'red3')
points(disProbesM, log10(disProbesS), pch = 19, cex = 1, col = 'blue3')
legend('bottomright', legend = c('original', 'random added', 'tissue added', 'tumType added'),
       pch = c(19,1,19,19), col = c('steelblue', 'black', 'red3', 'blue3'), cex = 1.25, bty = 'n')
par(op)

# Plot original PCA
Cols = ifelse(samples$Disease == 'AC', 'orangered',
              ifelse(samples$Disease == 'SCC', 'darkblue', 'grey'))
Pch <- ifelse(samples$Status == 'Tumor', 19, 1)
pcaSamples <- prcomp(t(eset))
par(mar = c(5, 6, 4, 2), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
plotPCA(pcaSamples, pch = Pch, col = Cols, main = 'PCA on original data set')
legend('bottomright', legend = c('Normal', 'Adk', 'Sq'), pch = c(1, 19, 19),
       col = c('black', 'orangered', 'darkblue'), cex = 1.25, bty = 'n')
par(op)

# PCA Filtering
pcaProbes <- prcomp(eset)
par(mar = c(5, 6, 4, 2), cex.main = 2.5, cex.lab = 1.75, cex.axis = 1.25)
score <- pcaTrace1.1(eset, pcaProbes, main = 'Information curve', lwd = 5)
par(op)
Info <- pcaInfo(score); Info

probeClass <- rep('lightblue', nrow(eset))
probeClass[grep('random',rownames(eset))] <- 'black'
probeClass[grep('status',rownames(eset))] <- 'red3'
probeClass[grep('disease',rownames(eset))] <- 'blue3'
pairs(pcaProbes$x[,1:3], col = probeClass, pch = c(1, 19, 19, 19)[factor(probeClass)], cex = 1)


# Visualize filtered PCA at different cutoffs
Cols = ifelse(samples$Disease == 'AC', 'orangered',
              ifelse(samples$Disease == 'SCC', 'darkblue', 'grey'))
Pch <- ifelse(samples$Status == 'Tumor', 19, 1)
par(mfrow = c(3, 2))
redDot = 10^(1/score$lModel$x.intercept)
for(p in c(0.05, 0.1, 0.25)){
  select <- pcaSelect(score, p)
  n <- length(select)
  pcaS <- prcomp(t(eset[select,]))
  pcaR <- prcomp(t(eset[-select,]))
  plotPCA(pcaS, pch = Pch, col = Cols, main = paste(n,'informative probes'))
  plotPCA(pcaR, pch = Pch, col = Cols, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
          main = paste(nrow(eset)-n,'rejected probes'))  
}
par(op)

# Compare Tumour Vs Normal on the original data set
eset <- Chem$eset
samples <- Chem$samples
idx <- which(samples$TumCellPercent>=50)
samples <- samples[idx,]
eset <- eset[,idx]

TN <- factor(as.character(samples$Status))
testTN <- mt.maxT(eset, classlabel = TN, B = 10000)
bestTN <- testTN$index[testTN$adjp<0.001]
length(bestTN); nrow(eset); length(bestTN)/nrow(eset)
length(intersect(grep('random', rownames(eset)), bestTN))
length(intersect(grep('signif', rownames(eset)), bestTN))

# Compare in Tumors only, Sq Vs AC
idx <- which(samples$Status == 'Tumor' & samples$Disease %in% c('AC', 'SCC'))
disease <- factor(as.character(samples$Disease[idx]))
tumors <- eset[,idx]
testDis <- mt.maxT(tumors, classlabel = disease, B = 10000)
bestDis <- testDis$index[testDis$adjp<0.001]
length(bestDis); nrow(tumors); length(bestDis)/nrow(tumors)
length(intersect(grep('random', rownames(tumors)), bestDis))
length(intersect(grep('signif', rownames(tumors)), bestDis))

# Plot each signif probes on a pairs(pcaProbes)
pcaProbes <- prcomp(eset)
probeClass <- rep('lightblue', nrow(eset))
probeClass[bestTN] <- 'red3'
probeClass[bestDis] <- 'blue3'
pairs(pcaProbes$x[,1:3], col = probeClass, pch = c(19,1,19)[factor(probeClass)], cex = 1)
