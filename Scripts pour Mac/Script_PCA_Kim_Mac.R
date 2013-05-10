# Recherche de genes informatifs par ACP
# A partir de la matrice virtuelle initiale

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
source(paste(pathScripts, "/plot.express.R", sep = ""))
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

grp <- factor(infos$TumourType)# grp <- factor(infos$IHC_status1)


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


	# Filtre et ACPs finales
# a = 7; alpha = 10^(-a)

yTarg = (best.Fmax + best.Fb)* 0.85
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
	plot(acp2$x, asp = 1, pch = 19, cex = pcex, col = mycol[grp], main = "Experiments (all probes)")

	plot(acp3$x[,1:2], cex = 0.2, col = "royalblue4", main = "", asp = 1)
	title(main = paste("Number of informative probes\n", length(inform), "of", nrow(eset), "(alpha =", alpha, ")"))

	plot.new()
	legend("center", legend = levels(grp), col= mycol[1:nlevels(grp)], pch = 19, cex = 1.5, bty = "n" )

	PC3 <- acp4$x[,3]
	mp = min(PC3)
	Mp = max(PC3)
	pcex = 2*(PC3-mp)/(Mp - mp) + 0.5
	Cols <- rep(NA, length(pcex))
	Cols[which(myCut == 1)] <- rgb(0.2, 0.4, 1, alpha = pcex[which(myCut == 1)]/max(pcex))
	Cols[which(myCut == 2)] <- rgb(1, 0.4, 0.2, alpha = pcex[which(myCut == 2)]/max(pcex))
	# plot(acp4$x, pch = 19, cex = pcex, col = mycol[grp], main = "Experiments (informative probes)", asp = 1)
	plot(acp4$x, pch = 19,
		cex = pcex,
		col = Cols,
		main = "Experiments (informative probes)", asp = 1)
par(op)


# graph3D.9(acp4, bold = T, class1 = grp)

heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2, RowSideColors = mycol[grp])


# T.Test validation by testing PCA clusters
myClust <- hclust(dist(t(sub1)), "ward")
myCut <- cutree(myClust, k = 2)

heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2, RowSideColors = mycol[as.factor(myCut)])

mt <- mt.maxT(sub1, classlabel = as.factor(myCut))
mt <- mt[order(mt$index), ]

hm <- heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2,
	RowSideColors = mycol[as.factor(myCut)], ColSideColors = ifelse(mt$adjp>0.05, "grey75", "grey25"))


# T.Test validation by testing Histology

mt2 <- mt.maxT(eset, classlabel = as.factor(infos$TumourType))
mt2 <- mt2[order(mt2$index), ]

heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2,
	RowSideColors = mycol[as.factor(infos$TumourType)], ColSideColors = ifelse(mt2$adjp[inform] > 0.05, "grey75", "grey25"))

plot(-log10(mt2$adjp[hm$colInd]))

isInform <- rep("No", nrow(eset))
isInform[inform] <- "Yes"
isSignif <- ifelse(mt2$adjp>0.001, "NS", "S")
table(isInform, isSignif)

BoxPoints5(isInform, mt2$teststat)

heatmap.3(t(eset[which(isSignif == "S"),]), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2,
	RowSideColors = mycol[as.factor(infos$TumourType)])


km1 <- survfit(st ~ myCut)
cox1 <- coxph(st ~ myCut)
km2 <- survfit(st ~ infos$TumourType)
cox2 <- coxph(st ~ infos$TumourType)

summary(cox1); summary(cox2)

par(mfrow = c(1, 2))
plot(km1)
plot(km2)
par(op)

