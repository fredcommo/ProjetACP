# Filtering on Kim Lung data

eset <- read.csv('/gluster/home/jcommo/Kim_Lung/Kim_data2.txt', header = TRUE, sep = '\t')
samples <- read.csv('/gluster/home/jcommo/Kim_Lung/Kim_patients.txt', header = TRUE, sep = '\t')

source('~/Fred_Scripts/pcaSelect.R')
source('~/Fred_Scripts/pcaInfo.R')
source('~/Fred_Scripts/pcaTrace.R')

PCA <- prcomp(eset)
score <- pcaTrace(eset, PCA)
Info <- pcaInfo(score); Info
informProp <- Info$Inform/(Info$Inform+Info$nonInform)

select <- pcaSelect(score, 0.1); length(select)
final <- prcomp(t(eset[select,]))
unselected <- prcomp(t(eset[-select,]))

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
pairs(final$x[,1:5], cex = pcex, col = c('red', 'blue')[samples$TumourType],
      pch = c(1, 19)[samples$Recurrence+1])

graph3D.8(final, class1 = samples$TumourType, class2 = samples$Recurrence)

require(survival)
Time <- samples$RecFreeSurv_month
status <- samples$Recurrence
st <- Surv(Time, status)
km <- survfit(st ~ samples$Sex)
plot(km)