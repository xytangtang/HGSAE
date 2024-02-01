rm(list=ls())

source("util_functions.R")

library(MASS)
library(GIGrvg)


# load data
mydata <- read.table(file="../data/data1.txt", header=TRUE)
theta_true <- scan(file="../data/data1_theta0.txt")

m <- nrow(mydata)
X <- cbind(rep(1, m), mydata$X1, mydata$X2, mydata$X3)
Y <- mydata$Y
D <- mydata$d
p <- 4

n_burn <- 50000
n_sim <- 50000
n_thin <- 10

U_lm <- residuals(lm(Y ~ X - 1))
XD <- X/sqrt(D)
sigma_Beta <- solve(t(XD) %*% XD)
PXD <- XD %*% solve(t(XD) %*% XD) %*% t(XD)
omega_u <- diag(1/sqrt(D)) %*% (diag(m) - PXD) %*% diag(1/sqrt(D))

#### MCMC chains #####
## Hierarchical Gamma (HG) model
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

# Fay-Herriot (FH)
set.seed(12345)
FH_results <- FH_sampler(X, Y, D, n_burn, n_sim, n_thin)

# Horseshoe (HS)
set.seed(12345)
HS_results <- TPB_sampler(0.5, 0.5, X, Y, D, n_burn, n_sim, n_thin)

# Laplace (LA)
set.seed(12345)
LA_results <- NG_sampler(1, X, Y, D, n_burn, n_sim, n_thin)

save(theta_true, m, X, Y, p, D, HG_res, HS_results, LA_results, FH_results, file="state_data_results.RData")

load("state_data_results.RData")

n_method <- 4
method_names <- c("HG", "FH", "HS", "LA")
DIC_res <- rep(NA, n_method)
names(DIC_res) <- method_names
Dev_res <- matrix(NA, n_method, 4)
dimnames(Dev_res) <- list(method_names, c("AAD", "ASD", "ARB", "ASRB"))

## HG
# posterior median of a and b
a_median <- median(HG_res$a)
b_median <- median(HG_res$b)
round(c(a_median, b_median), 2)


HG_Theta.chain <- X%*%HG_res$Beta + HG_res$U
DIC_res["HG"] <- dic(Y, HG_Theta.chain, D)
Theta_est_HG <- apply(HG_Theta.chain, 1, mean)
Theta_sd_HG <- apply(HG_Theta.chain, 1, sd)
Dev_res["HG", ] <- dev_measure(Theta_est_HG, theta_true)


## FH
FH_Theta.chain <- FH_results$Theta.chain
DIC_res["FH"] <- dic(Y, FH_Theta.chain, D)
Theta_est_FH <- apply(FH_Theta.chain, 1, mean)
Theta_sd_FH <- apply(FH_Theta.chain, 1, sd)
Dev_res["FH", ] <- dev_measure(Theta_est_FH, theta_true)

## HS
HS_Theta.chain <- X%*%HS_results$Beta.chain + HS_results$U.chain
DIC_res["HS"] <- dic(Y, HS_Theta.chain, D)
Theta_est_HS <- apply(HS_Theta.chain, 1, mean)
Theta_sd_HS <- apply(HS_Theta.chain, 1, sd)
Dev_res["HS", ] <- dev_measure(Theta_est_HS, theta_true)

## LA
LA_Theta.chain <- X%*%LA_results$Beta.chain + LA_results$U.chain
DIC_res["LA"] <- dic(Y, LA_Theta.chain, D)
Theta_est_LA <- apply(LA_Theta.chain, 1, mean)
Theta_sd_LA <- apply(LA_Theta.chain, 1, sd)
Dev_res["LA", ] <- dev_measure(Theta_est_LA, theta_true)

## DIC results
round(DIC_res, 2)

## Deviation measures results
round(Dev_res[,c("AAD", "ASD")], 2)

## plot figure
pdf("plot_state_res.pdf", width=10.5, height=3.5)
par(mfrow=c(1,3), mar=c(4.5, 4.5, 1.5, 1.5))
# a
plot(density(HG_res$a), xlab=expression(a), ylab="density", main="")
# b
plot(density(HG_res$b), xlab=expression(b), ylab="density", xlim=c(0, 200), main="")
# theta sd
plot(Theta_sd_HG, Theta_sd_HS, pch=16, cex=0.5,xlab=expression(paste("sd", group("(",theta[HG],")"))), ylab=expression(paste("sd", group("(",theta[HS],")"))))
abline(a=0, b=1, lty=2, col="blue")
dev.off()








