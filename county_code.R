rm(list=ls())

source("util_functions.R")

library(MASS)
library(GIGrvg)

mydata <- read.table(file="data/data2.txt", header=TRUE)

Y <- mydata$y
m <- nrow(mydata)
X <- cbind(rep(1, m), mydata$x)
D <- mydata$D
p <- ncol(X)


n_burn <- 50000
n_sim <- 50000
n_thin <- 10
	
	
U_lm <- residuals(lm(Y ~ X - 1))
XD <- X/sqrt(D)
sigma_Beta <- solve(t(XD) %*% XD)
PXD <- XD %*% solve(t(XD) %*% XD) %*% t(XD)
omega_u <- diag(1/sqrt(D)) %*% (diag(m) - PXD) %*% diag(1/sqrt(D))

	
##### Hierarchical Gamma (HG) prior model
set.seed(12345)
Beta <- rnorm(p)
U <- U_lm
Lambda2 <- rep(1, m)
tau2 <- 1
a <- runif(1)
b <- 1
tau2_prior_par1 <- 10^-10
tau2_prior_par2 <- 10^-10
HG_res <- HG_sampler(Beta, U, Lambda2, tau2, a, b, Y, X, D, n_burn, n_sim, n_thin, tau2_prior_par1, tau2_prior_par2)

# LA
set.seed(12345)
LA_results <- NG_sampler(1, X, Y, D, n_burn, n_sim, n_thin)

# HS
set.seed(12345)
HS_results <- TPB_sampler(0.51, 0.5, X, Y, D, n_burn, n_sim, n_thin)

# FH
set.seed(12345)
FH_results <- FH_sampler(X, Y, D, n_burn, n_sim, n_thin)

save(X, Y, m, D, p, HG_res, LA_results, HS_results, FH_results, file="county_data_results.RData")

load("county_data_results.RData")

n_method <- 4
method_names <- c("HG", "FH", "HS", "LA")
DIC_res <- rep(NA, n_method)
names(DIC_res) <- method_names

## HG
a_median <- median(HG_res$a)
b_median <- median(HG_res$b)
round(c(a_median, b_median), 2)

HG_Theta.chain <- X%*%HG_res$Beta + HG_res$U
DIC_res["HG"] <- dic(Y, HG_Theta.chain, D)
Theta_est_HG <- apply(HG_Theta.chain, 1, mean)
Theta_sd_HG <- apply(HG_Theta.chain, 1, sd)

## FH
FH_Theta.chain <- FH_results$Theta.chain
DIC_res["FH"] <- dic(Y, FH_Theta.chain, D)
Theta_est_FH <- apply(FH_Theta.chain, 1, mean)
Theta_sd_FH <- apply(FH_Theta.chain, 1, sd)

## HS
HS_Theta.chain <- X%*%HS_results$Beta.chain + HS_results$U.chain
DIC_res["HS"] <- dic(Y, HS_Theta.chain, D)
Theta_est_HS <- apply(HS_Theta.chain, 1, mean)
Theta_sd_HS <- apply(HS_Theta.chain, 1, sd)# HS : -15751.12

## LA
LA_Theta.chain <- X%*%LA_results$Beta.chain + LA_results$U.chain
DIC_res["LA"] <- dic(Y, LA_Theta.chain, D)
Theta_est_LA <- apply(LA_Theta.chain, 1, mean)
Theta_sd_LA <- apply(LA_Theta.chain, 1, sd)


## DIC results
round(DIC_res, 2)
plot(Theta_est_HG, Theta_est_LA, pch=16, cex=0.5)
abline(a=0, b=1, lty=2, col="blue")


## plot figure
pdf("plot_county_res.pdf", width=10.5, height=3.5)
par(mfrow=c(1,3), mar=c(4.5, 4.5, 1.5, 1.5))
# a
plot(density(HG_res$a), xlab=expression(a), ylab="density", main="")
# b
plot(density(HG_res$b), xlab=expression(b), ylab="density", main="")
# theta sd
plot(Theta_sd_HG, Theta_sd_LA, pch=16, cex=0.5,xlab=expression(paste("sd", group("(",theta[HG],")"))), ylab=expression(paste("sd", group("(",theta[LA],")"))))
abline(a=0, b=1, lty=2, col="blue")
dev.off()


