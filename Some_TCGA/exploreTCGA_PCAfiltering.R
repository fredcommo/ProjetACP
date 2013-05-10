# PCA filtering on TCGA pancancer white-lists

source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaSelect.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaInfo.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace1.1.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/pcaTrace2.1.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/plotPCA.R')
source('/Users/fredcommo/Documents/MyProjects/FredScripts/heatmap.3.R')

library(synapseClient)
library(cgdsr)
library(glmnet)
library(parallel)
library(ROCR)
library(impute)
library(corpcor)

synapseLogin()

op <- par(no.readonly = TRUE)

.buildRandom <- function(Data, n){
  rand <- lapply(1:n, function(i){
    x <- Data[sample(1:nrow(Data), 1),]
    return(sample(x))
  })
  rand <- do.call(rbind, rand)
  colnames(rand) <- colnames(Data)
  rownames(rand) <- paste0('rand', seq(1, nrow(rand)))
  return(as.data.frame(rand))
}
.pcaSqDist <- function(PCA, Dim = 2:3){
  X <- as.data.frame(PCA$x[,Dim])
  X <- as.data.frame(scale(X))
  D <- rowSums(X^2)
  return(D)
}
.plotDist <- function(Dist, Names,...){
  plot(density(log(Dist[grep('rand', Names)]), na.rm = TRUE),
       ylim = range(0, 1), xlim = range(log(Dist), na.rm = TRUE),
       lwd = 3, col = 'grey',...)
  leg = c('random probes'); k = 1
  for(f in seq(0.8, 0, by = -0.2)){
    leg <- c(leg, paste('between', f, 'and', f+0.2))
    lines(density(log(Dist[zeroFreq<f+0.2 & zeroFreq>=f]), na.rm = TRUE),
          lwd = 3, col = k)
    k = k+ 1
  }
  legend('topleft', title = 'Prop of zero', legend = leg, cex = 1.5,
         lwd = 2, col = c('grey', 1:k))
}
  
resTCGA <- synapseQuery('select id, name from entity where entity.benefactorId == "syn300013"')
resTCGA <- resTCGA[order(resTCGA$entity.name),]
rnaseqIdx <- grep('(.*)RNASeqV2.geneExp.whitelist_tumor', resTCGA$entity.name)
rnaSeqList <- resTCGA[rnaseqIdx,]

setwd('/Users/fredcommo/Documents/MyProjects/Projet ACP/Some_TCGA/')

for (i in 1:nrow(rnaSeqList)){
  cat('running:', rnaSeqList$entity.id[i])
  synId = rnaSeqList$entity.id[i]
  plotTitle = paste0(rnaSeqList$entity.name[i], '\n', rnaSeqList$entity.id[i])
  e <- downloadEntity(synId)
  eset <- read.table(file.path(e$cacheDir,e$files), header=TRUE,
                   row.names=1, comment="", quote="", sep="\t")
  rand <- .buildRandom(eset, 1000)
  esetRand <- rbind.data.frame(eset, rand)
  zeroFreq <- apply(esetRand, 1, function(x){sum(x==0)/length(x)})
  esetRand <- log(esetRand + 1)
  pcaProbes <- prcomp(esetRand)
  
  Dist <- .pcaSqDist(pcaProbes)
  png(file = paste0(rnaSeqList$entity.name[i], '_', rnaSeqList$entity.id[i], '_sqDist.png'),
      width = 800, height = 800)
    par(mar = c(5, 5, 5, 2), cex.main = 1.75, cex.lab = 1.5, cex.axis = 1.25)
    .plotDist(Dist, Names = rownames(esetRand), main = plotTitle, sub = 'Square Dist. distribution',
            xlab = 'log(Square Dist.)')
    par(op)
  dev.off()
  
  Col = rgb(1-zeroFreq, 0, 0, 1-zeroFreq)
  png(file = paste0(rnaSeqList$entity.name[i], '_', rnaSeqList$entity.id[i], '_PCAfilt.png'),
      width = 800, height = 800)
    par(mar = c(5, 5, 6, 2), cex.main = 1.75, cex.lab = 1.5, cex.axis = 1.25)
    plotPCA(pcaProbes, Axes = 2:3, Cex = 1, l = 1, pch = 19, col = Col,
            main = plotTitle)
    abline(h = 0, v = 0, lty = 2, lwd = 3, col = 'cyan')
    points(pcaProbes$x[,2:3], pch = 19, cex = 1,
         col = rgb(0, 0, zeroFreq, zeroFreq))
    legend('bottomleft', legend = 'Prop of Zero per row', pch = 19, col = 'blue')
    par(op)
  dev.off()
  cat('\n')
}

################
# Original PCAs

for (i in 1:nrow(rnaSeqList)){
  cat('running:', rnaSeqList$entity.id[i], '\n')
  synId = rnaSeqList$entity.id[i]
  tissue <- unlist(strsplit(rnaSeqList$entity.name[i], '_'))[2]
  e <- downloadEntity(synId)
  eset <- read.table(file.path(e$cacheDir,e$files), header=TRUE,
                     row.names=1, comment="", quote="", sep="\t")
  eset <- log(eset + 1)
  pcaSamples <- prcomp(t(eset))
  pc3R <- abs(pcaSamples$x[,3])
  pCol = ifelse(pc3S/max(pc3S)<0.5, 0.5, pc3S/max(pc3S))
  png(file = paste0(rnaSeqList$entity.name[i], '_', rnaSeqList$entity.id[i], '_Original.png'),
      width = 800, height = 800)
    par(mar = c(5, 5, 5, 2), cex.main = 1.75, cex.lab = 1.5, cex.axis = 1.25)
    plotPCA(pcaSamples, Axes = 1:2, Cex = 3, l = 3, pch = 19, col = rgb(0,0,1,pCol),
          main = paste('PCA on', tissue, '-', synId,
                       '\nall', nrow(eset),'features'))
    par(op)
  dev.off()
}

################
# Filtering

for (i in 1:nrow(rnaSeqList)){
  cat('running:', rnaSeqList$entity.id[i], '\n')
  synId = rnaSeqList$entity.id[i]
  tissue <- unlist(strsplit(rnaSeqList$entity.name[i], '_'))[2]
  e <- downloadEntity(synId)
  eset <- read.table(file.path(e$cacheDir,e$files), header=TRUE,
                     row.names=1, comment="", quote="", sep="\t")
  eset <- log(eset + 1)
  pcaProbes <- prcomp(eset)
  
  # Information curve
  png(file = paste0(rnaSeqList$entity.name[i], '_', rnaSeqList$entity.id[i], '_Trace.png'),
      width = 800, height = 600)
    par(mar = c(5, 5, 5, 2), cex.main = 1.75, cex.lab = 1.5, cex.axis = 1.25)
    score <- pcaTrace1.1(eset, pcaProbes, lwd = 5,
                         main = paste(tissue, ':',synId, '\ninformation curve'))
    par(op)
  dev.off()
  
  infos <- pcaInfo(score)
  #redDot = 10^(1/score$lModel$x.intercept)
  png(file = paste0(rnaSeqList$entity.name[i], '_', rnaSeqList$entity.id[i], '_samplePCAfilt.png'),
        width = 1000, height = 1200)
  par(mfrow = c(3, 2), mar = c(5, 5, 6, 2), cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.3)
  for(p in c(0.05, 0.1, 0.25)){
    select <- pcaSelect(score, p)
    pcaOnSelect <- prcomp(t(eset[select,]))
    pc3S <- abs(as.numeric(pcaOnSelect$x[,3]))
    pcaOnReject <- prcomp(t(eset[-select,]))
    pc3R <- abs(pcaOnReject$x[,3])
    xLim = max(abs(pcaOnSelect$x[,1]))
    yLim = max(abs(pcaOnSelect$x[,2]))
    pCol = ifelse(pc3S/max(pc3S)<0.5, 0.5, pc3S/max(pc3S))
    plotPCA(pcaOnSelect, Axes = 1:2, Cex = 3, l = 3, pch = 19, col = rgb(0,0,1,pCol),
            main = paste('PCA on', tissue, '-', synId,
                          '\nfilt:', p, '-', length(select), 'selected probes'))
  pCol = ifelse(pc3R/max(pc3R), 0.25, pc3S/max(pc3R))
  plotPCA(pcaOnReject, Axes = 1:2, Cex = 3, l = 3, pch = 19, col = rgb(0,0,1,pCol),
            xlim = range(-xLim, xLim), ylim = range(-yLim, yLim),
            main = paste('PCA on', tissue, '-', synId,
                         '\nfilt:', p, '-', nrow(eset) - length(select), 'rejeted probes'))
    }
  cat('\n')
  par(op)
  dev.off()
}


#####################################
# Push to synapse
setwd('/Users/fredcommo/Documents/MyProjects/Projet ACP/Some_TCGA/')

Path <- paste0(getwd(), '/')
projectId = 'syn1834953'
myCode <- Code(list(name = "exploreTCGA_PCAfiltering.R", parentId = projectId))
myCode <- addFile(myCode, paste0(Path, 'exploreTCGA_PCAfiltering.R'))
myCode <- storeEntity(myCode)
synIds <- rep(rnaSeqList$entity.id, each = 2)


#####################################
# Push files & create a wiki page
#myCode <- getEntity('syn1834779')

listFiles = list.files()
for (l in 2:length(listFiles)){
  fileName = listFiles[l]
#  synId = synIds[l-1]
#  synId = rnaSeqList$entity.id[i]
  splitFileName <- unlist(strsplit(fileName, '_'))
  tissue <- splitFileName[2]
  synId <- splitFileName[6]
  cat(tissue, synId, '\n')
  
  folderName <- paste(tissue, synId, sep = '_')
  Folders <- synapseQuery('select id, name from entity where entity.parentId == "syn1834953"')
  if(!folderName %in% Folders$entity.name){
    # Create a new folder
    newFolder <- Folder(name = folderName, parentId = projectId)
    newFolder <- synStore(newFolder)
    parentId <- propertyValue(newFolder, 'id')
    # Create the folder wiki uri.    
    folderWikiUri <- sprintf("/entity/%s/wiki", parentId)
    # Start a wiki, and initialize the markdown as empty.
    folderWiki <- list()
    folderWiki$attachmentFileHandleIds <- list()
    folderWiki$markdown <- c()
    folderWiki <- synRestPOST(folderWikiUri, folderWiki)
    cat('New folder created\n')
  }
  
  # Push the file to the current folder
  file <- File(paste0(Path, fileName), parentId = parentId)
  file <- synStore(file)
  cat(fileName, 'pushed\n')
  
  # Add an activity
  rawData <- getEntity(synId)
  used(file)<-list(list(entity = myCode, wasExecuted = TRUE),
                   list(entity = rawData, wasExecuted = FALSE))
  file <- synStore(file)
  
  # Update the folder markdown with the new image to display.
  folderWiki <- synRestGET(folderWikiUri)
  folderUpdateWikiUri <- sprintf("/entity/%s/wiki/%s", parentId, folderWiki$id)
  folderWiki$attachmentFileHandleIds[length(folderWiki$attachmentFileHandleIds)+1] <- file@fileHandle$id
  newMarkdown <-  paste0(paste0(.mkdText(fileName), '\nImage Location: ', propertyValue(file, 'id')),
                               '\n${image?synapseId=', propertyValue(file, 'id'),
                               '&align=None&scale=80}')
  folderWiki$markdown <- paste0(folderWiki$markdown, newMarkdown, '\n')

  # now 'push' the wiki to Synapse
  folderWiki <- synRestPUT(folderUpdateWikiUri, folderWiki)
  cat('wiki pushed\n')
  cat('\n')
}

.mkdText <- function(fileName){
  if(grepl('Original', fileName)) markdownText <- "###Original PCA:\n*PCA on samples using the full matrix.*"
  if(grepl('PCAfilt', fileName)) markdownText <- "###PCA on features:\n*Using 2nd and 3rd axes, high zero frequencies are located at the center.*"
  if(grepl('samplePCAfilt', fileName)) markdownText <- "###PCA according to filtering:\n*PCA on samples using different levels of filtering (left: using selected probes, right: using rejected probes).*"
  if(grepl('sqDist', fileName)) markdownText <- "###Distances from the center:\n*Distribution of distances (log(d^2)), according to the frequence of zero.*"
  if(grepl('Trace', fileName)) markdownText <- "###Information curve:\n*Samples spread increasing according to the information brought by new probes.*"
  return(markdownText)
  }

# #####################################
# # Push files & create a wiki page
# 
# listFiles = list.files()
# for (l in 2:length(listFiles)){
#   fileName = listFiles[l]
#   #  synId = synIds[l-1]
#   #  synId = rnaSeqList$entity.id[i]
#   splitFileName <- unlist(strsplit(tmpName, '_'))
#   tissue <- splitFileName[2]
#   synId <- splitFileName[6]
#   cat(tissue, synId, '\n')
#   
#   folderName <- paste(tissue, synId, sep = '_')
#   Folders <- synapseQuery('select id, name from entity where entity.parentId == "syn1709808"')
#   if(!folderName %in% Folders$entity.name){
#     newFolder <- Folder(name = folderName, parentId = projectId)
#     parentId <- propertyValue(newFolder, 'id')
#   }
#   else
#     parentId <- Folders$entity.id[Folders$entity.name == folderName]
#   
#   file <- File(paste0(Path, fileName), parentId = parentId)
#   file <- synStore(file)
#   fileHandleId <- file@fileHandle$id
#   fileWikiUri <- sprintf("/entity/%s/wiki", propertyValue(file, "id"))
#   
#   # we start a wiki
#   fileWiki<-list()
#   # we add to our wiki the ID of the previously uploaded file
#   fileWiki$attachmentFileHandleIds<-list(fileHandleId)
#   
#   # in the markdown we say to display the image.
#   fileWiki$markdown <- paste0('${image?synapseId=',
#                               propertyValue(file, 'id'),
#                               '&align=None&scale=80}')
#   
#   # now 'push' the wiki to Synapse
#   fileWiki<-synRestPOST(fileWikiUri, fileWiki)
#   
#   # Add an activity
#   rawData <- getEntity(synId)
#   used(file)<-list(list(entity = myCode, wasExecuted = TRUE),
#                    list(entity = rawData, wasExecuted = FALSE))
#   file <- synStore(file)
#   cat('\n')
# }


# .addWiki <- function(parentId, Path, File){
#   # demo of how to create a file then use the uploaded file in a wiki
#   # note, we recommend using "File" rather than "Data"  
#   file <- File(paste0(myPath, myImage), parentId = parentId)
#   file <- synStore(file)
#   fileHandleId <- file@fileHandle$id
#   fileWikiUri <- sprintf("/entity/%s/wiki", propertyValue(file, "id"))
#   # we start a wiki
#   fileWiki<-list()
#   # we add to our wiki the ID of the previously uploaded file
#   fileWiki$attachmentFileHandleIds<-list(fileHandleId)
#   # in the markdown we say to display the image.  Note, 'fileName' is the URLEncoded version of the file chosen above.
#   graphName <- gsub('.', '%2E', File)
#   graphName <- gsub('_', '%5F', File)
#   graphName <- paste0('${image?fileName=', graphName,'}')
#   fileWiki$markdown <- paste0('${image?fileName=', graphName,'}')
#   # now 'push' the wiki to Synapse
#   fileWiki<-synRestPOST(fileWikiUri, fileWiki)
# }


# listFiles = list.files()
# for (l in 2:length(listFiles)){
#   tmpName = listFiles[l]
#   synId = synIds[l-1]
#   cat(synId)
#   myData <- Data(list(name = tmpName, parentId = parentId))
#   myData <- addFile(myData, tmpName)
#   myData <- storeEntity(myData)
#   # Add an activity
#   rawData <- getEntity(synId)
#   used(myData)<-list(list(entity = myCode, wasExecuted = TRUE),
#                      list(entity = rawData, wasExecuted = FALSE))
#   myData <- storeEntity(myData)
#   cat('\n\n')
# }
