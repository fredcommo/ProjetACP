# PCA-based filtering on Platinium spike data: GSE21344
# 2013_05_28
# FC

#########################
# Store in synapse
require(synapseClient)
file <- File(path = 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE21344&format=file',
             name = 'GSE21344', parentId = 'syn1897344', synapseStore = FALSE)
file <- synStore(file)
#########################

#########################
# Download data from synapse
require(synapseClient)
require(affy)
require(multtest)
gse21344 <- synGet('syn1897345')
celPath <- paste0(gsub('\\?.+', '', gse21344@filePath), 'CEL')
untar(gse21344@filePath, exdir = celPath)
setwd(celPath)
rawData <- ReadAffy()

# Normalization
eset.mas5 <- mas5(rawData)
eset <- log2(exprs(eset.mas5))
boxplot(eset)

esetCall <- exprs(mas5calls(rawData))
freqCall <- lapply(1:nrow(esetCall), function(i){
  A <- sum(esetCall[i,]=='A')
  M <- sum(esetCall[i,]=='M')
  P <- sum(esetCall[i,]=='P')
  c(A = A/ncol(esetCall), M = M/ncol(esetCall), P = P/ncol(esetCall))
})
freqCall <- do.call(rbind, freqCall)
freqCall <- apply(freqCall, 1, function(x){names(which.max(x))})
rownames(freqCall) <- rownames(esetCall)

# Get the samples annotations
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaFilters.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
op <- par(no.readonly = TRUE)

setwd('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_platinium_spike')
samples <- read.csv('PlatiniumSpike_samples.txt', header = TRUE, sep = '\t')
trueExp <- read.csv('PlatiniumSpike_probesValues.txt', header = TRUE, sep = '\t')

pcaProbes <- prcomp(eset)
probeCols <- ifelse(freqCall=='A', rgb(0.2, 0,0,0.25),
                    ifelse(freqCall=='M', rgb(0.5, 0,0,0.25), rgb(1, 0,0,0.25)))  
pairs(pcaProbes$x[,1:3], cex = 0.2, col = probeCols)

probeCols <- rep(rgb(.8, .8, .8,0.25), nrow(trueExp))
probeCols[trueExp$expr=='high'] <- rgb(1, 0, 0, 0.25)
probeCols[trueExp$expr=='low'] <- rgb(0, 0, 1, 0.25)
probeCols[trueExp$expr=='normal'] <- rgb(0, 0.5, 0.5, 0.25)
probeCols[grepl('MC|MF', trueExp$expr)] <- rgb(0, 0.5, 0, 0.25)
pairs(pcaProbes$x[,1:3], cex = 0.3, col = probeCols)

score <- pcaTrace1.1(eset, pcaProbes, lwd = 4)
info <- pcaInfo(score)

par(mfrow = c(2, 2))
for(p in c(0.05, 0.1, 0.25, 0.5)){
  select <- pcaSelect(score, p)
  n <- length(select)
  pcaS <- prcomp(t(eset[select,]))
  pcaR <- prcomp(t(eset[-select,]))
#  plotPCA(pcaS, col = c('orangered', 'darkblue')[samples$Condition], pch = 19,
#          main = paste('PCA on', n, 'selected probes'))
  plotPCA(pcaR, col = c('orangered', 'darkblue')[samples$Condition], pch = 19,
          xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
          main = paste('PCA on', nrow(eset)-n, 'rejected probes'))
  tab1 <- table(trueExp$expr[-select])/table(trueExp$expr)
  tab2 <- table(factor(freqCall)[-select])/table(freqCall)
  print(round(tab1, 4)*100)
  print(round(tab2, 4)*100)
}
par(op)

fullTest <- mt.maxT(eset, classlabel=samples$Condition, B = 5000)
fullBest <- fullTest$index[fullTest$adjp<1e-3]
tab1 <- table(trueExp$expr[fullBest])/table(trueExp$expr)
round(tab1, 4)*100

select <- pcaSelect(score, 0.0)
filtTest <- mt.maxT(eset[select,], classlabel=samples$Condition, B = 5000)
filtBest <- filtTest$index[filtTest$adjp<1e-3]
tab2 <- table(trueExp$expr[select][filtBest])/table(trueExp$expr)
round(tab2, 4)*100

select <- pcaSelect(score, 0.1)
filtTest <- mt.maxT(eset[select,], classlabel=samples$Condition, B = 5000)
filtBest <- filtTest$index[filtTest$adjp<1e-3]
tab3 <- table(trueExp$expr[select][filtBest])/table(trueExp$expr)
round(tab3, 4)*100

select <- pcaSelect(score, 0.25)
filtTest <- mt.maxT(eset[select,], classlabel=samples$Condition, B = 5000)
filtBest <- filtTest$index[filtTest$adjp<1e-3]
tab4 <- table(trueExp$expr[select][filtBest])/table(trueExp$expr)
round(tab4, 4)*100

select <- pcaSelect(score, 0.5)
filtTest <- mt.maxT(eset[select,], classlabel=samples$Condition, B = 5000)
filtBest <- filtTest$index[filtTest$adjp<1e-3]
tab5 <- table(trueExp$expr[select][filtBest])/table(trueExp$expr)
round(tab5, 4)*100

filterPerf <- lapply(seq(0, 1, by = 0.05), function(p){
  select <- pcaSelect(score, p)
  cat(p, '\tselect:', length(select),'\n')
  Test <- mt.maxT(eset[select,], classlabel=samples$Condition, B = 5000)
  Best <- Test$index[Test$adjp<1e-2]
  tab <- table(trueExp$expr[select][Best])/table(trueExp$expr)
  cat('\nBest', length(Best),'\n\n')
  return(c(p = p, size = length(select), signSize = length(Best), round(tab, 4)*100))
})
fPerf <- do.call(rbind, filterPerf)
fPerf <- as.data.frame(fPerf)
barplot(fPerf$low, col = rgb(0,0.5,0.8,0.5))
barplot(fPerf$high, col = rgb(0.8,0.5,0,0.5), add = TRUE)
