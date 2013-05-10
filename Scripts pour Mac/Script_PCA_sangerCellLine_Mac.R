######################
op <- par(no.readonly = TRUE)
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


library(synapseClient)
library(affy)
synapseLogin("frederic.commo@sagebase.org", "Se@ttle7")


# Upload the mRNA Expression data set: syn427896
ent <- loadEntity('syn427896')
eset <- exprs(ent$objects$eset)
dim(eset)	#; eset[1:10, 1:10]

# Upload GE metadata : syn424134
metadata <- loadEntity('syn424134')
metaDat <- metadata$objects$pData
dim(metaDat)	#; metaDat[1:10, 1:10]
colnames(eset) <- metaDat$SampleName

# Upload the supervised CNV data set: syn464292
cnvEnt <- loadEntity('syn464292')
cnv <- exprs(cnvEnt$objects$eset)
dim(cnv)	#; cnv[1:10, 1:10]

# Upload CNV metadata : syn464191
metaCNV <- loadEntity('syn464191')
cnvDat <- metaCNV$objects$pData
dim(cnvDat)	#; cnvDat[1:10, ]
colnames(cnv) <- cnvDat$SampleName

# Attention: réplicats dans le GE
checkNames <-  table(metaDat$SampleName)
Duplic <- which(checkNames>1)

R <- c()
for(i in 1:length(Duplic)){
	ij = which(colnames(eset) == names(Duplic)[i])
	metaDat <- metaDat[-ij[2:length(ij)],]
	eset[,ij[1]] <- rowMeans(eset[,ij])
	eset <- eset[,-ij[2:length(ij)]]
	}
 
dim(eset); dim(metaDat)

matchCNV <- which(!is.na(match(cnvDat$Name, metaDat$SampleName))); length(matchCNV) 	# which CNV are present in GE
matchGE <- which(!is.na(match(metaDat$SampleName, cnvDat$Name))); length(matchGE) 	    # which GE are present in CNV

cnv <- cnv[,matchCNV]
cnv <- cnv[,order(colnames(cnv))]
cnvDat <- cnvDat[matchCNV,]
cnvDat <- cnvDat[order(cnvDat$Name),]
eset <- eset[,matchGE]
eset <- eset[,order(colnames(eset))]
metaDat <- metaDat[matchGE,]
metaDat <- metaDat[order(metaDat$SampleName),]

cbind.data.frame(colnames(cnv), cnvDat$Name, colnames(eset), metaDat$SampleName)

boxplot(eset[,sample(1:ncol(eset), 100)])

mycol <- c("royalblue2", "indianred1", "seagreen3", "goldenrod2", "blue3", "purple1", "burlywood2", "firebrick1", "green4", "deeppink1", "deepskyblue2")

grp <- factor(rep(1, nrow(metaDat))) # grp <- factor(infos$IHC_status1)
grp <- metaDat$PrimaryHist	# grp <- metaDat$PrimarySite

mycol <- c("royalblue2", "indianred1", "seagreen3", "goldenrod2", "blue3", "purple1", "burlywood2", "firebrick1", "green4", "deeppink1", "deepskyblue2")

if(nlevels(grp)>10){
	mycol <- colors()[sample(c(8:152, 365:657), nlevels(grp))]
}

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
	plot(acp1$x, asp = 1, cex = 0.25, col = "royalblue4", main = "Genes")
	plot(acp2$x, asp = 1, pch = 19, col = mycol[grp], main = "Experiments")
	legend("topright", legend = levels(grp), pch = 19, cex = 0.2, col = mycol[1:nlevels(grp)], bty = "n")
par(op)



################################################################################
# Launch Diag by Trace v3
################################################################################


	# Filtre et ACPs finales
# a = 7; alpha = 10^(-a)

informTable

yTarg = (best.Fmax + best.Fb)* 0.5
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
	myPCAplot(acp2$x, grp, main = "Experiments (all probes)")

	plot(acp3$x[,2:3], cex = 0.2, col = "royalblue4", main = "", asp = 1)
	title(main = paste("Number of informative probes\n", length(inform), "of", nrow(eset), "(alpha =", signif(alpha, 3), ")"))

	plot.new()
	legend("center", legend = levels(grp), col= mycol[1:nlevels(grp)], pch = 19, cex = 1.5, bty = "n" )

	myPCAplot(acp4$x, grp, main = "Experiments (informative probes)")
par(op)


# graph3D.9(acp4, bold = T, class1 = grp)

heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2, RowSideColors = mycol[grp])

grp <- metaDat$PrimaryHist	# grp <- metaDat$PrimarySite

mycol <- c("royalblue2", "indianred1", "seagreen3", "goldenrod2", "blue3", "purple1", "burlywood2", "firebrick1", "green4", "deeppink1", "deepskyblue2")

if(nlevels(grp)>10){
	mycol <- colors()[sample(c(8:152, 365:657), nlevels(grp))]
}

heatmap.3(t(sub1), scale = "column", Method = "ward", col = greenred(min(100, ncol(sub1))), key = TRUE, trace = "none", density = "none",
	breaks = seq(-2.5, 2.5, len = min(100, ncol(sub1))+1),
	labRow = grp, cexRow = 0.5, labCol = "", cexCol = 0.2, RowSideColors = mycol[grp])

# Links between classes & covariates
output <- c()
for(K in 2:20){
	myClust <- hclust(dist(t(sub1)), "ward")
	myCut <- cutree(myClust, k = K)
	myCut <- as.factor(myCut)
	
	for(V in c(8, 9, 11, 12:16)){
		myTab <- table(myCut, as.factor(metaDat[,V]))
		pVal <- chisq.test(myTab)$p.value
		output <- rbind(output, c(nClass = K, colNum = V, nLevels = nlevels(as.factor(metaDat[,V])), pValue = pVal))
		}
}
as.data.frame(output)

V = 8
