require(synapseClient)
require(glmnet)
require(venneuler)
require(multtest)
require(corpcor)
require(foreach)

source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaSelect.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaInfo.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace1.1.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/heatmap.3.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/GeneRequest.v7.R')

op <- par(no.readonly = TRUE)

synapseLogin('frederic.commo@sagebase.org', 'Se@ttle7')

ccle <- loadEntity('1671195')

setwd('/Users/fredcommo/Documents/MyProjects/CellLines/')

ccleGE <- ccle$objects$ccle_data$ccle_exp
ccleCNV <- ccle$objects$ccle_data$ccle_cnv
ccleMut <- ccle$objects$ccle_data$ccle_mut

ccleIC <- read.table('CCLE_IC50Norm.txt', sep = '\t')
ccleA <- read.table('CCLE_ActAreaNorm.txt', sep = '\t')
ccleDrug <- read.table('CCLE_DrugList.txt', sep = '\t')
ccleCells <- read.table('CCLE_CellsList.txt', sep = '\t')

levels(ccleCells$CCLE.name)[200:201] <- c('G401_KIDNEY', 'G402_KIDNEY')
levels(ccleCells$CCLE.name)[766] <- c('RKN_OVARY')
colnames(ccleGE)[200] <- colnames(ccleCNV)[196] <- 'G402_KIDNEY'

ccleCells <- ccleCells[order(ccleCells$CCLE.name),]
ccleGE <- ccleGE[,order(colnames(ccleGE))]
ccleCNV <- ccleCNV[,order(colnames(ccleCNV))]
ccleMut <- ccleMut[,order(colnames(ccleMut))]
  
# Match colnames GE & CNV with ccleCells$CCLE.name
# Match ccleIC & ccleA with rownames(ccleCells)

# Get cells for which GE exists
I1 <- intersect(colnames(ccleGE), colnames(ccleCNV))
I2 <- intersect(colnames(ccleMut), I1)
commonCells <- intersect(ccleCells$CCLE.name, I2)
ccleGE <- ccleGE[,colnames(ccleGE) %in% commonCells]
ccleCNV <- ccleCNV[,colnames(ccleCNV) %in% commonCells]
ccleMut <- ccleMut[,colnames(ccleMut) %in% commonCells]
ccleCells <- ccleCells[ccleCells$CCLE.name %in% commonCells,]

geCells <- which(colnames(ccleGE) %in% ccleCells$CCLE.name)
ccleInfo <- which(ccleCells$CCLE.name %in% colnames(ccleGE))
ccleCells <- ccleCells[ccleInfo,]

# Get cells for which IC50 exists and GE exists
cclePrimNames <- toupper(gsub('-| |\\.|\\(|\\)', '', ccleCells$Cell.line.primary.name))
haveIC <- which(rownames(ccleIC) %in% cclePrimNames)
ccleIC <- ccleIC[haveIC, ]
ccleA <- ccleA[haveIC, ]
ccleCells <- ccleCells[which(cclePrimNames %in% rownames(ccleIC)), ]
ccleGE <- ccleGE[,which(colnames(ccleGE) %in% ccleCells$CCLE.name)]
ccleCNV <- ccleCNV[,which(colnames(ccleCNV) %in% ccleCells$CCLE.name)]
ccleMut <- ccleMut[,which(colnames(ccleMut) %in% ccleCells$CCLE.name)]

# Tranform ccleMut into binary values
mutNames <- rownames(ccleMut)
ccleMut <- lapply(1:nrow(ccleMut), function(i){ifelse(ccleMut[i,]=='0', rnorm(1, 0, 0.01), rnorm(1, 1, 0.01))})
ccleMut <- do.call(rbind, ccleMut)
rownames(ccleMut) <- mutNames

# Check the cell names
all(rownames(ccleIC) == toupper(gsub('-| |\\.|\\(|\\)', '', ccleCells$Cell.line.primary.name)))
all(rownames(ccleA) == toupper(gsub('-| |\\.|\\(|\\)', '', ccleCells$Cell.line.primary.name)))
all(ccleCells$CCLE.name == colnames(ccleGE))
all(ccleCells$CCLE.name == colnames(ccleCNV))
all(ccleCells$CCLE.name == colnames(ccleMut))

# Rename each rows in order to avoid duplicated names
rownames(ccleGE) <- paste0(rownames(ccleGE), '_expr')
rownames(ccleCNV) <- paste0(rownames(ccleCNV), '_cnv')
rownames(ccleMut) <- paste0(rownames(ccleMut), '_mut')

setwd('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_on_CCLE/')

# Use CNV only
target <- 'MEK'
drug <- ccleDrug[grep(target, ccleDrug$Target.s.),]
drugName <- as.character(drug$Compound..code.or.generic.name.)

# Which cells
cells <- rownames(ccleA)
GEids <- as.character(ccleCells$CCLE.name[(match(rownames(ccleCells), rownames(ccleA)))])
cbind.data.frame(rownames(ccleA), ccleCells$CCLE.name, GEids)
mekResp <- ccleA[,colnames(ccleA)==drugName[1]]
Data <- ccleGE[ ,match(colnames(ccleGE), GEids)]

  # Restrict to lung
idx <- grep('LUNG', GEids)
Data <- apply(Data[,idx], 1, scale)
Data <- as.data.frame(t(Data))
mekResp <- mekResp[idx]

  # Perform initial PCA on samples
pcaSamples <- prcomp(t(Data))
plotPCA(pcaSamples)

  # Perform initial PCA on probes
pcaProbes <- prcomp(Data)
plotPCA(pcaProbes)

png(file = paste0('CCLE_', target, 'inhib', '_GE_Trace.png'),
    width = 800, height = 600)
par(mar = c(5, 5, 5, 2), cex.main = 1.75, cex.lab = 1.5, cex.axis = 1.25)
score <- pcaTrace1.1(Data, pcaProbes, lwd = 5,
                     main = paste('CCLE_', target, 'inhib', '\ninformation curve'))
par(op)
dev.off()

info <- pcaInfo(score)

par(mfrow = c(3, 2))
for(p in c(0.05, 0.1, 0.25)){
  select <- pcaSelect(score, p)
  n <- length(select)
  pcaS <- prcomp(t(Data[select,]))
  plotPCA(pcaS, pch = 19,
          main = paste('PCA with', n, 'selected probes'))
  pcaR <- prcomp(t(Data[-select,]))
  plotPCA(pcaR, pch = 19, xlim = range(pcaS$x[,1]), ylim = range(pcaS$x[,2]),
          main = paste('PCA with', nrow(Data) - n, 'rejected probes'))
}
par(op)


fullModel <- .netModel1(Data, mekResp, B = 100, alpha = 0.01)
topListFull <- as.numeric(fullModel$boostGenes)
yFitFull <- predict(fullModel$model, newx = t(Data[topListFull,]), s = fullModel$bestL)
pCorFull <- as.numeric(cor(mekResp, yFitFull, use = 'pairwise.complete.obs'))

png(file = paste0(paste0('CCLE_', target, 'inhib', '_GE_glmnet.png')),
    width = 1000, height = 1200)
par(mfrow = c(3, 2), mar = c(5, 5, 6, 2), cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.3)

  plot(mekResp, yFitFull, main = 'model on full GE matrix')
  abline(lm(yFit ~ mekResp), lwd = 3, col = 'blue')
  legend('topleft', legend = paste('pearson cor:', round(pCor, 3)), cex = 1.5)

  for(p in c(0.05, 0.1, 0.25, 0.5, 0.75)){
    select <- pcaSelect(score, p)
    n <- length(select)
    filteredModel <- .netModel1(Data[select,], mekResp, B = 100, alpha = 0.01)
    topList <- as.numeric(filteredModel$boostGenes)
    yFit <- predict(filteredModel$model, newx = t(Data[select[topList],]), s = filteredModel$bestL)
    pCor <- as.numeric(cor(mekResp, yFit, use = 'pairwise.complete.obs'))
    plot(mekResp, yFit, pch = 19, col = 'grey30', cex = 1.25,
       main = paste0('MEK inhib: glmnet on ', n, ' selected genes (p: ', p, '%)\n',
                     nrow(Data) - n, ' rejected'))
    abline(lm(yFit ~ mekResp), lwd = 3, col = 'blue')
    legend('topleft', legend = paste('pearson cor:', round(pCor, 3)), cex = 1.5)
  }
par(op)
dev.off()

# Compute the perf according to cutoff
K = 20

pCorFull <- lapply(1:K, function(k){fullModel <- .netModel1(Data, mekResp, B = 50, alpha = 0.01)
                                    topListFull <- as.numeric(fullModel$boostGenes)
                                    yFitFull <- predict(fullModel$model, newx = t(Data[topListFull,]), s = fullModel$bestL)
                                    as.numeric(cor(mekResp, yFitFull, use = 'pairwise.complete.obs'))
}
)
pCorFull <- do.call(c, pCorFull)
pCorFull <- data.frame(p = rep(0, K), n = rep(nrow(Data), K), pCor = pCorFull)

perf <- lapply(seq(0.1, 1, by = 0.1),function(p){
  select <- pcaSelect(score, p)
  n <- length(select)
  pCor <- rep(NA, K)
  if(n>2){
    tmp <- lapply(1:K, function(k){
      filteredModel <- .netModel1(Data[select,], mekResp, B = 50, Nfold = 10, alpha = 0.01)
      topList <- as.numeric(filteredModel$boostGenes)
      yFit <- predict(filteredModel$model, newx = t(Data[select[topList],]), s = filteredModel$bestL)
      as.numeric(cor(mekResp, yFit, use = 'pairwise.complete.obs'))
    }
    )
    pCor <- do.call(c, tmp)
  }
  cbind(p = rep(p, K), n = rep(n, K), pCor)
}
)
perf <- as.data.frame(do.call(rbind, perf))
perf

Res <- rbind.data.frame(pCorFull, perf)

plot(perf$p, perf$pCor)
abline(h = pCorFull)
