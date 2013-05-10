######################
library(synapseClient)
library(affy)
synapseLogin("frederic.commo@sagebase.org", "Se@ttle7")


# Upload Normalized GE-Agilent4502a
ent <- loadEntity('syn313607')

# Access to expression matrix
eset <- exprs(ent$objects$eset)
dim(eset)

# Upload clinical data associated with syn285162 (TGCA_Lung_ADK parent Id)
metadata <- loadEntity('syn673169')
metaDat <- metadata$objects$metadata

setwd("/Users/fredcommo/Documents/Sage/Projects/TCGA_Lung_ADK")

# Find the column containing eset Ids: colnames(eset)
metaDat <- metadata$objects$metadata
columnMatchTotals <- apply(metaDat, 2, function(x){sum(!is.na(match(colnames(eset), x)))})
columnId <- which(columnMatchTotals > 0); length(columnId)
columnId; columnMatchTotals[columnId]
matchDat <- match(colnames(eset), metaDat[, columnId])
matchDat <- matchDat[!is.na(matchDat)]
metaDat <- metaDat[matchDat,]
metaDat <- metaDat[,which(colSums(!is.na(metaDat)) > 0)]
rownames(metaDat) <- metaDat[,columnId]

metaDat <- metaDat[,c(14, 48:52, 54:56, 57, 64, 68, 75, 78, 80:85, 87, 126, 127, 132, 159, 161:164)]

# searchBrak <- apply(metaDat, 2, function(x){sum(substr(x, 1, 1) == "[")})
# emptyCol <- which(searchBrak == nrow(metaDat))
# searchUnc <- which(substr(colnames(metaDat), 1, 3) == "unc")
# emptyCol <- which(searchBrak == nrow(metaDat))

# write.table(metaDat, "TCGA_Lung_ADK_syn673169.xls", sep = "\t", row.names = T)

matchEset <- match(rownames(metaDat), colnames(eset))
eset <- eset[,matchEset]

# Verify order
Err <- ifelse(colnames(eset) == rownames(metaDat), "ok", "err")
length(which(Err == "err"))

##########################################
op <- par(no.readonly = TRUE)
# memory.limit(size = 4000)
library(MASS)
library(KernSmooth)
library(mvtnorm)
library(lattice)
library(limma)
library(multtest)
library(gplots)
library(mclust)
library(imputation)
library(impute)
library(survival)
pathScripts = "/Users/fredcommo/Documents/Stats/Fred_Scripts"
source(paste(pathScripts, "/graph3D.9.R", sep = ""))
source(paste(pathScripts, "/rb.colors.R", sep = ""))
# source(paste(pathScripts, "/plot.express.R", sep = ""))
source(paste(pathScripts, "/heatmap.3.R", sep = ""))
source(paste(pathScripts, "/rot3D.2.R", sep = ""))
source(paste(pathScripts, "/BoxPoints5.R", sep = ""))
source(paste(pathScripts, "/Resamp.T.v3.R", sep = ""))
source(paste(pathScripts, "/myPCAplot.R", sep = ""))

mycol <- c("royalblue2", "indianred1", "seagreen3", "goldenrod2", "blue3", "purple1", "burlywood2", "firebrick1", "green4", "deeppink1", "deepskyblue2")

##############################################################################

	# Kim lung
data.address <- "/Volumes/FREECOM HDD/Databases GE/Data sets Poumon/Kim/Kim.data2.txt"
info.address <- "/Volumes/FREECOM HDD/Databases GE/Data sets Poumon/Kim/Kim.patients.txt"
eset <- read.table(data.address, header = T, sep = "\t")
infos <- read.table(info.address, header = T, sep = "\t")
infos$Recurrence <- factor(infos$Recurrence)

# ! v?rifier !
verif <- ifelse(colnames(eset)==infos$Id, "ok", "error")
length(which(verif=="error"))

Time<- infos$RecFreeSurv_month
status <- c(0, 1)[infos$Recurrence]
st <- Surv(Time, status)

grp <- rep("A", nrow(metaDat)) 	# factor(infos$TumourType)# grp <- factor(infos$IHC_status1)
if(!is.factor(grp)) grp <- as.factor(grp)

##############################################################################
	# ACPs initiales
if(any(is.na(eset))){
	eset <- impute.knn(as.matrix(eset))$data
  }
eset <- as.data.frame(eset)
acp1 <- prcomp(scale(eset))			# genes = obs ; exp = variables
acp2 <- prcomp(scale(t(eset)))			# exp = obs ; genes = variables

par(mfrow=c(2,2))
	plot(acp1, main = "Genes")
	plot(acp2, main = "Experiments")
	myPCAplot(acp1$x[,1:3], Fact = rep("A", nrow(eset)), main = "Genes", Leg = FALSE)
	RGB <- myPCAplot(acp2$x[,1:3], Fact = grp, main = "Experiments")
	# legend("topright", legend = levels(grp), pch = 19, col = mycol[1:nlevels(grp)], bty = "n")
par(op)

################################################################################
# Launch Diag by Trace v3
################################################################################

yTarg = (best.Fmax + best.Fb)* 0.95
xTarg = best.c - best.b*log(((best.Fmax - best.Fb)/(yTarg - best.Fb))^(1/best.d) - 1)

aTarg <- 10^(xTarg)
alpha = 10^(-aTarg)
Pmax = 1 - alpha
Q <- qchisq(p = Pmax, df = ncol(X))
D <- apply(X^2, 1, sum)
inform <- which(D>=Q)
length(inform)

sub1 <- eset[inform,]
acp3 <- prcomp(sub1)			# genes = obs ; exp = variables
acp4 <- prcomp(t(sub1))			# exp = obs ; genes = variables

par(mfrow=c(2,2))
	PC3 <- acp2$x[,3]
	mp = min(PC3)
	Mp = max(PC3)
	pcex = 2*(PC3-mp)/(Mp - mp) + 0.5
	myPCAplot(acp2$x, Fact = grp, Leg = F, Col = RGB, main = "Experiments (all probes)")

	plot(acp3$x[,1:2], cex = 0.2, col = "royalblue4", main = "", asp = 1)
	title(main = paste("Number of informative probes\n", length(inform), "of", nrow(eset), "(alpha =", alpha, ")"))

	plot.new()
	# legend("center", legend = levels(grp), col= mycol[1:nlevels(grp)], pch = 19, cex = 1.5, bty = "n" )
	
	myPCAplot(acp4$x, Fact = grp, Leg = F, Col = RGB, main = "Experiments (filtered probes)")
par(op)


# graph3D.9(acp4, bold = T, class1 = grp)

heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2, RowSideColors = mycol[grp])


# T.Test validation by testing PCA clusters
myClust <- hclust(dist(t(sub1)), "ward")
myCut <- cutree(myClust, k = 3)
table(myCut)
which(myCut == 2)

rownames(sub1)[1:10]


