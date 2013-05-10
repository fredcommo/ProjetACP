
	# Diagnostic sur la valeur de alpha
X <- as.data.frame(acp1$x[,2:3])	# X <- as.data.frame(acp2$x[,2:3]) # s?rie Daniel gautheret
X <- scale(X)
X <- as.data.frame(X)
D <- rowSums(X^2)

# a.values <- c(0.025, 0.05, 0.1, 0.25, 0.5, seq(1, 5))
a.values <- c(0.025, 0.5, 5)

par(mfrow = c(3, 2))
Score <- S2 <- sum.trace <- c()
for(a in a.values){
	alpha = 10^(-a)
	Pmax = 1 - alpha
	Q <- qchisq(p = Pmax, df = ncol(X))
	# pt.col <- ifelse(D>=Q, "red", "grey")

	theta <- seq(0, 90, by = 1)
	x1 <- sqrt(Q)*cos(theta)
	y1 <- sqrt(Q)*sin(theta)

	plot(X, cex = 0.2, pch = 19, col = ifelse(D<=Q, "royalblue4", "grey"))
	points(y1~x1, col = "red", cex = 0.1, pch = 8, xlim = range(X))
	inform <- which(D<=Q)
	title(main = paste("Q =", signif(Q, 3), ";", length(inform), "probes"))
	tmp.n.select <- tmp.S2 <- tmp.sum.trace <- NA
	if(length(inform)>0){
		tmp.n.select <- length(inform)
		sub.inf <- t(eset[inform,])
		acpTest <- prcomp(sub.inf)
		tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
		# plot(prcomp(sub.inf)$x[, 1:2], pch = 19, cex = 1.5, col = mycol[grp], 
		# 	xlim = range(-100, 100), ylim = range(-100, 100), main = "Samples segregation")
		myPCAplot(prcomp(sub.inf)$x[, 1:3], Fact = grp, Col = RGB, 
			xlim = range(-50, 50), ylim = range(-100, 100), main = paste("Samples segregation", "; Trace", round(tmpTrace))
		# plot(prcomp(sub.inf)$x[, 1:2], pch = 19, cex = 1.5)	#, col = mycol[grp])
		# tmp.S2 <- var(as.vector(as.matrix(sub.inf)))
		# tmp.sum.trace <- sum(diag(var(sub.inf)))
		}
	# Score <- rbind(Score, c(a, tmp.n.select, tmp.sum.trace, tmp.S2))
	# S2 <- c(S2, tmp.S2)
	#  <- c(sum.trace, tmp.sum.trace)
	}
# plot(sum.trace, a.values)
par(op)


	# calcul de la trace

	# Diagnostic sur la valeur de alpha
X <- as.data.frame(acp1$x[,2:3])	# X <- as.data.frame(acp2$x[,2:3]) # s?rie Daniel gautheret
X <- scale(X)
X <- as.data.frame(X)
D <- rowSums(X^2)

# a.values <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.025, 0.05, 0.1, 0.25, 0.5, seq(1, 22, by = 3))
a.values <- c(1e-4, 1e-3, 1e-2, 0.025, 0.05, 0.1, 0.25, 0.5,  1, 2, 4, 8, 16, 24, 32)

useTrace = T
Score <- c()

for(a in a.values){
	alpha = 10^(-a)
	Pmax = 1 - alpha
	Q <- qchisq(p = Pmax, df = ncol(X))
	inform <- which(D<=Q)
	lInf <- length(inform)

	tmpScore <- NA

	if(lInf>0 & lInf<nrow(eset)){
		tmpNselect <- length(inform)
		subInf <- t(eset[inform,])
		tmpS2 <- var(as.vector(as.matrix(subInf)))
		tmpScore <- c(aValues = a, nProbes = tmpNselect, S2 = tmpS2)

			# calculer les variances sur les axes1 et 2 de la nouvelle acp ?
		if(useTrace){
			acpTest <- prcomp(subInf)
			tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
			tmpScore <- c(tmpScore, Trace = tmpTrace)
			}
		}
		if(lInf == nrow(eset)) tmpScore <- c(aValues = a, Score[nrow(Score), -1])
		Score <- rbind(Score, tmpScore)
	}

rownames(Score) <- seq(1, nrow(Score))
Score <- as.data.frame(Score)
Score

# colnames(Score) <- c("aValues", "nProbes", "S2", "Trace")

plot(log10(Score$aValues), Score$Trace, xlab = "Log[-Log(alpha)]", ylab = "Trace of Cov matrix")

# Fonction logistique 5PL
	Richard <- function(x, Fb, Fmax, b, c, d){
		y <- Fb + (Fmax - Fb)/(1 + exp(-(x-c)/b))^d
		return(y)
		}

# Fonction sce (somme carr? r?sidus) avec pond?rations
	sce.5P <- function(param, xobs, yobs, w) {
		# Fb <- param[1]
		Fb = 0
		Fmax <- param[2]
		b <- param[3]
		c <- param[4]
		d <- param[5]
		ytheo <- Richard(xobs, Fb, Fmax, b, c, d)
		sq.res <- (yobs - ytheo)^2
		weights <- 1/sq.res^w
		return(sum(weights*sq.res))
		}

# Fonction sce (somme carr? r?sidus) avec pond?rations
	sce.5P.diag <- function(yobs, ytheo, w) {
		sq.res <- (yobs - ytheo)^2
		weights <- 1/sq.res^w
		return(weights)
		}



x <- as.numeric(log10(Score$aValues))
y <- as.numeric(Score$Trace)

if(any(is.na(y) | is.na(x))){
  na.index <- which(is.na(y) | is.na(x))
  y <- y[-na.index]
  x <- x[-na.index]
}

# initialisation des parametres
Fb.ini = 0	                                                 # min(y)
Fmax.ini = max(y)	                                           #*1.05
c.ini = 0                                                    #(max(x) + min(x))/2
z <- (y)/(Fmax.ini - y)
if (any(abs(z)==Inf)) z[abs(z)==Inf] <- NA
b.ini = coef(lm(x ~ log(z)))[2]
d.ini = 1
init <- c(Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)

w = 0.25

# Estimation du modele
	best <- nlm(f = sce.5P, p = as.vector(init), xobs = x, yobs = y, w = w)
  best

# R?cup?ration des param?tres
	best.Fb <- best$estimate[1]
	best.Fmax <- best$estimate[2]
	best.b <- best$estimate[3]
	best.c <- best$estimate[4]
	best.d <- best$estimate[5]

# Diagnostic de r?gression
	newx <- seq(min(log(a.values)), max(log(a.values)), len = 100)
	yfit <- Richard(newx, best.Fb, best.Fmax, best.b, best.c, best.d)
  yfitInit <- Richard(x, Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)

	plot(x, y, lwd = 2, cex = 1.25, col = "grey25", xlab = "Log[-Log(alpha)]", ylab = "Trace of Cov matrix")
	lines(yfitInit ~ x, lwd = 3, lty = 3, col = "grey75")
	lines(yfit ~ newx, lwd = 3, col = "dodgerblue2")
	# abline(h = c(best.Fb, best.Fmax))

	informTable <- c()
	for(i in seq(0.05, 1, by = 0.05)){
		yTarg = (best.Fmax + best.Fb)*i
		xTarg = best.c - best.b*log(((best.Fmax - best.Fb)/(yTarg - best.Fb))^(1/best.d) - 1)

		aTarg <- 10^(xTarg)
		alpha = 10^(-aTarg)
		Pmax = 1 - alpha
		Q <- qchisq(p = Pmax, df = ncol(X))

		inform <- which(D>=Q)
		nInform <- length(inform)
		nNonInform <- nrow(eset) - nInform
		informTable <- rbind(informTable, c(Prop = i, Alpha = signif(alpha, 3), nInf = nInform, nonInf = nNonInform))
	}

	informTable
	
	
		# S?lection
	
	yTarg = (best.Fmax + best.Fb)* 0.5
	xTarg = best.c - best.b*log(((best.Fmax - best.Fb)/(yTarg - best.Fb))^(1/best.d) - 1)
	points(xTarg, yTarg, pch = 19, cex = 1.5, col = "red3")

	aTarg <- 10^(xTarg)

	alpha = 10^(-aTarg)
	Pmax = 1 - alpha
	Q <- qchisq(p = Pmax, df = ncol(X))
	# pt.col <- ifelse(D>=Q, "red", "grey")

	# theta <- seq(0, 90, by = 1)
	# x1 <- sqrt(Q)*cos(theta)
	# y1 <- sqrt(Q)*sin(theta)


	# plot(X, cex = 0.2, pch = 19, col = ifelse(D>=Q, "royalblue4", "grey"), asp = 1)
	# points(y1~x1, col = "red", cex = 0.25, pch = 8, xlim = range(X), asp = 1)
	inform <- which(D>=Q)
	nInform <- length(inform); cat("Information in", nInform, "\n")
	nNonInform <- nrow(eset) - nInform; cat("No information in", nNonInform, "\n")

	esetInform <- eset[inform,]
	esetNonInform <- eset[-inform,]

	acpInf <- prcomp(t(esetInform))			# exp = obs ; genes = variables
	acpNonInf <- prcomp(t(esetNonInform))			# exp = obs ; genes = variables

	par(mfrow = c(1, 2))
	PC3 <- acpNonInf$x[,3]
	mp = min(PC3)
	Mp = max(PC3)
	pcex = 2*(PC3-mp)/(Mp - mp) + 0.25
	plot(acpNonInf$x[, 1:2], pch = 1, cex = pcex, col = c(2:3)[grp], xlim = range(acpInf$x[,1]), ylim = range(acpInf$x[,2]), main = paste("PCA on", nNonInform, "non informative probes"), asp = 1) # 

	PC3 <- acpInf$x[,3]
	mp = min(PC3)
	Mp = max(PC3)
	pcex = 2*(PC3-mp)/(Mp - mp) + 0.25
	plot(acpInf$x[, 1:2], pch = 1, cex = pcex, col = c(2:3)[grp], main = paste("PCA on", nInform, "informative probes"), asp = 1)
	par(op)


		# Visualisation 3D de la s?lection
	Information <- rep("No", nrow(eset))
	Information[inform] <- "Yes"
	graph3D.9(acp1, class1 = Information, size = 2)
	graph3D.9(acp1, class1 = grp, size = 2)


		# Visualisation de la s?gr?gation
	graph3D.9(acpInf, size = 1)
	graph3D.9(acpNonInf, size = 1)

	# Visualisation de la s?gr?gation
	# F <- as.factor(rnaIds$X)
	# n <- nlevels(F)
	graph3D.9(acpInf, class1 = grp, size = 0.5)
	graph3D.9(acpNonInf, class1 = grp, size = 0.5)
	legend("bottomleft", legend = levels(F), pch = 19, col = mycol[1:n])

	# Visualisation de la s?gr?gation
	F <- as.factor(rnaIds$X)
	n <- nlevels(F)
	graph3D.9(acp1, class1 = F, size = 1.25)
	legend("bottomleft", legend = levels(F), pch = 19, col = mycol[1:n])

	table(rnaIds$X, mycol[rnaIds$X])

	# Dendrogramme sur les ?chantillons
	Clust <- hclust(dist(acpInf$x[,1:2]), "ward")
	clustRes <- cbind.data.frame(org = colnames(eset), Ord = Clust$order)
	clustRes <- clustRes[order(clustRes$Ord),]

	ncol = 100
	hd <- heatmap.3(as.matrix(esetInform), scale = "row", Method = "ward", col = greenred(ncol), breaks = seq(-3, 3, len = ncol + 1), key = TRUE, trace = "none", density = "none",
		labRow = "", RowSideColors = mycol[grp[inform]], cexRow = 0.5, labCol = "", cexCol = 0.2) # , labRow = grp, RowSideColors = mycol[grp]
	


	##################################################

	# Kohonen classif
	source("E:\\Stats\\Doc R\\Scripts R\\carteKoho.R")
	n = 10
	m = 10
	koho <- carteKoho(esetInform , n, m, pch = 19, cex = 1, col = mycol[grp])
	Class <- koho$classifier
	myTab <- as.matrix(t(table(grp[inform], Class$class)))
	myTab <- cbind(Num = seq(1, n*m), myTab)

	# write.table(myTab, "Kohonen_Freq.xls", sep = "\t", row.names = F)

	selectNeur <- sort(unique(Class$class[grp == "misc_RNA"]))
	selectNeur

	index <- which(Class$class == 23)

		# Heatmap sur la s?lection
	ncol = 100
	heatmap.3(as.matrix(esetInform), scale = "row", Method = "ward", col = greenred(ncol), breaks = seq(-3, 3, len = ncol + 1), key = TRUE, trace = "none", density = "none",
		labRow = "", RowSideColors = mycol[grp[inform]], cexRow = 0.5, labCol = "", cexCol = 0.2) # , labRow = grp, RowSideColors = mycol[grp]
	

	heatmap.3(as.matrix(esetNonInform), scale = "row", Method = "ward", col = greenred(ncol), breaks = seq(-3, 3, len = ncol + 1), key = TRUE, trace = "none", density = "none",
		labRow = "", RowSideColors = mycol[grp[-inform]], cexRow = 0.5, labCol = "", cexCol = 0.2) # , labRow = grp, RowSideColors = mycol[grp]
	
	heatmap.3(as.matrix(esetInform[index,]), scale = "row", Method = "ward", col = greenred(ncol), breaks = seq(-3, 3, len = ncol + 1), key = TRUE, trace = "none", density = "none",
		labRow = "", RowSideColors = mycol[grp[inform][index]], cexRow = 0.5, labCol = "", cexCol = 0.2) # , labRow = grp, RowSideColors = mycol[grp]
	

	##################################################

PlotSelect <- function(X, Eset, Score, Prop = c(0.05, 0.25, 0.5), Range = c(-100, 0)){

	N <- nrow(eset)
	x <- as.numeric(log10(Score$aValues))
	y <- as.numeric(Score$Trace)

	par(mfrow = c(3, 3))
	for(p in Prop){
		yTarg = (best.Fmax + best.Fb)*p
		xTarg = best.c - best.b*log(((best.Fmax - best.Fb)/(yTarg - best.Fb))^(1/best.d) - 1)
		aTarg <- 10^(xTarg)

		alpha = 10^(-aTarg)
		Pmax = 1 - alpha
		Q <- qchisq(p = Pmax, df = ncol(X))
		theta <- seq(0, 90, by = 1)
		x1 <- sqrt(Q)*cos(theta)
		y1 <- sqrt(Q)*sin(theta)

		# trace curve
		plot(x, y, xlab = "Log[-Log(alpha)]", ylab = "Trace of Cov matrix")
		lines(yfit~newx, lwd = 3, col = "grey25")
		points(xTarg, yTarg, pch = 19, cex = 1.5, col = "red3")
		title(main = paste("Prop of trace:", p))

		# probes visualization (selection)
		plot(X, cex = 0.2, pch = 19, col = ifelse(D<=Q, "royalblue4", "grey"))
		points(y1~x1, col = "red", cex = 0.1, pch = 8, xlim = range(X))
		inform <- which(D<=Q)
		n <- length(inform)
		title(main = paste(length(inform), " probes (left ", N - n, " probes)", sep = ""))
		tmp.n.select <- tmp.S2 <- tmp.sum.trace <- NA
		if(length(inform)>0){
			subInf <- t(eset[inform,])
			# subNonInf <- t(eset[-inform,])
			acpI <- prcomp(subInf)
			# acpNI <- prcomp(subNonInf)
			plot(acpI$x[, 1:2], pch = 19, col = "grey25", cex = 1.5, xlim = Range, ylim = Range)
			title(main = paste("Suppressed:", n, ", kept:", N-n))
			}
		}
par(op)
}

PlotSelect(X = X, Eset = eset, Score = Score, Prop = c(0.025, 0.25, 0.75), Range = c(-90, 90))


