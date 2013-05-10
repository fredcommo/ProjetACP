acp.trace <- function(Data, Dims = c(2, 3), target = seq(0.05, 0.95, by = 0.05)){

	# Diagnostic sur la valeur de alpha
X <- prcomp(Data)$x[,Dims]
X <- scale(X)
X <- as.data.frame(X)
D <- rowSums(X^2)

a.values <- c(0.025, 0.05, 0.1, 0.25, 0.5, seq(1, 5))

cat("Performing traces calculation...\n")
n.select <- S2 <- sum.trace <- c()
for(a in a.values){
	alpha = 10^(-a)
	Qmax = 1 - alpha
	Q <- qchisq(p = Qmax, df = 1)
	n.inform <- which(D<Q)

	tmp.n.select <- tmp.S2 <- tmp.sum.trace <- NA
	if(length(n.inform)>0){
		sub.n.inf <- t(Data[n.inform,])
		tmp.n.select <- length(n.inform)
		tmp.S2 <- var(as.vector(as.matrix(sub.n.inf)))
		tmp.sum.trace <- sum(diag(var(sub.n.inf)))
		}
	n.select <- rbind(n.select, c(a, tmp.n.select))
	S2 <- c(S2, tmp.S2)
	sum.trace <- c(sum.trace, tmp.sum.trace)
	}

x <- as.numeric(log(a.values))
y <- as.numeric(sum.trace)

if(any(is.na(y) | is.na(x))){
	na.index <- which(is.na(y) | is.na(x))
	y <- y[-na.index]
	x <- x[-na.index]
	}

# Fonction logistique 5PL
	Richard <- function(x, Fb, Fmax, b, c, d){
		y <- Fb + (Fmax - Fb)/(1 + exp(-(x-c)/b))^d
		return(y)
		}

# Fonction sce (somme carré résidus) avec pondérations
	sce.5P <- function(param, xobs, yobs, w) {
		Fb <- param[1]
		Fmax <- param[2]
		b <- param[3]
		c <- param[4]
		d <- param[5]
		ytheo <- Richard(xobs, Fb, Fmax, b, c, d)
		sq.res <- (yobs - ytheo)^2
		weights <- 1/sq.res^w
		return(sum(weights*sq.res))
		}

# Fonction sce (somme carré résidus) avec pondérations
	sce.5P.diag <- function(yobs, ytheo, w) {
		sq.res <- (yobs - ytheo)^2
		weights <- 1/sq.res^w
		return(weights)
		}

cat("Performing model adjustment...\n")
# initialisation des parametres
	Fb.ini = min(y)
	Fmax.ini = max(y)	#*1.05
	c.ini = (max(x) + min(x))/2
	z <- (y)/(Fmax.ini - y)
	if (any(abs(z)==Inf)) z[abs(z)==Inf] <- NA
	b.ini = coef(lm(x~log(z)))[2]
	d.ini = 1
	init <- c(Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)
	w = 0.25

# Estimation du modele
	best <- nlm(f = sce.5P, p = init, xobs = x, yobs = y, w = w)

# Récupération des paramètres
	best.Fb <- best$estimate[1]
	best.Fmax <- best$estimate[2]
	best.b <- best$estimate[3]
	best.c <- best$estimate[4]
	best.d <- best$estimate[5]

# Diagnostic de régression
	newx <- seq(min(log(a.values)), max(log(a.values)), len = 100)
	yfit <- Richard(newx, best.Fb, best.Fmax, best.b, best.c, best.d)

	plot(x, y)
	lines(yfit~newx)
	abline(h = c(best.Fb, best.Fmax))
	Ytarget = best.Fmax * target
	Xtarget = best.c - best.b*log(((best.Fmax - best.Fb)/(Ytarget - best.Fb))^(1/best.d) - 1)

	return(list(Xfit = newx, Yfit = yfit, Ytarget = Ytarget, Xtarget = Xtarget))
}