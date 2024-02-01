
log_den_log_a_HG <- function(a, b, lambda2) {
	m <- length(lambda2)
	m * a * log(b) - m * lgamma(a) + (a - 1) * sum(log(lambda2))
}

stop_expand <- function(slice_lower, slice_upper, b, Lambda2, z) {
  if (slice_upper == 1 && slice_lower == 0) return(TRUE)
  if (slice_upper == 1 && slice_lower > 0) {
    if (log_den_log_a_HG(slice_lower, b, Lambda2) < z) return(TRUE)
    return(FALSE)
  }
  if (slice_upper < 1 && slice_lower == 0) {
    if (log_den_log_a_HG(slice_upper, b, Lambda2) < z) return(TRUE)
    return(FALSE)
  }
  if (log_den_log_a_HG(slice_upper, b, Lambda2) < z && log_den_log_a_HG(slice_lower, b, Lambda2) < z) return(TRUE)
  return(FALSE)
}

HG_sampler <- function(Beta, U, Lambda2, tau2, a, b, Y, X, D, n_burn, n_sim, n_thin, tau2_prior_par1, tau2_prior_par2, MAX_LAMBDA2 = 100, slice_width=0.1) {
  
  m <- length(Y)
  beta_chain <- matrix(0, p, n_sim/n_thin)
  u_chain <- matrix(0, m, n_sim/n_thin)
  lambda2_chain <- matrix(0, m, n_sim/n_thin)
  b_chain <- rep(0, n_sim/n_thin)
  a_chain <- rep(0, n_sim/n_thin)
  
  XD <- X/sqrt(D)
  sigma_Beta <- solve(t(XD) %*% XD)
  PXD <- XD %*% solve(t(XD) %*% XD) %*% t(XD)
  omega_u <- diag(1/sqrt(D)) %*% (diag(m) - PXD) %*% diag(1/sqrt(D))
  
  for (index_iter in 1:(n_sim + n_burn)) {
    if (index_iter %% 1000 == 0) cat("Iteration:", index_iter, "\n")
    
    # sample lambda2
    for (i in 1:m) {
    	Lambda2[i] <- rgig(1, lambda = a - 0.5, chi = U[i]*U[i]/tau2, psi = 2*b)
    	# check for infinity 
    	if (is.infinite(Lambda2[i])) Lambda2[i] <- MAX_LAMBDA2
    	if (Lambda2[i] == 0) Lambda2[i] <- 10^-6
    }
    
    # sample a using slice sampling
    log_f_a <- log_den_log_a_HG(a, b, Lambda2)
    a0 <- a
    z <- log_f_a - rexp(1)
    w_prop <- runif(1)
    slice_upper <- min(1, a0 + w_prop * slice_width)
    slice_lower <- max(0, a0 - (1-w_prop) * slice_width)
    
    while(!stop_expand(slice_lower, slice_upper, b, Lambda2, z)) {
      slice_upper <- min(1, slice_upper + slice_width)
      slice_lower <- max(0, slice_lower - slice_width)
    }
    a <- runif(1, min=slice_lower, max=slice_upper)
    log_f_a <- log_den_log_a_HG(a, b, Lambda2)
    while(log_f_a < z) {
      if(a < a0) slice_lower <- a
      else slice_upper <- a
      a <- runif(1, min=slice_lower, max=slice_upper)
      log_f_a <- log_den_log_a_HG(a, b, Lambda2)
    }
    
    # sample b
    b <- rgamma(1, a*m + tau2_prior_par2, sum(Lambda2) + tau2_prior_par2)
    
    # # sample tau2
    # tau2 <- 1/rgamma(1, shape=0.5*m+tau2_prior_par1, rate=sum(U*U/Lambda2/2) + tau2_prior_par2)
    
    # sample U
    var_u <- D*Lambda2*tau2/(D + Lambda2*tau2)
    mean_u <- var_u / D * (Y - X %*% Beta)
    U <- rnorm(m, mean = mean_u, sd = sqrt(var_u))
    
    # sample beta
    mean_Beta <- sigma_Beta %*% t(XD) %*% ((Y - U)/sqrt(D))
    Beta <- mvrnorm(1, mu = mean_Beta, Sigma = sigma_Beta)
    
    
    # save chains
    if (index_iter > n_burn && (index_iter - n_burn) %% n_thin == 0) {
      index_element <- (index_iter - n_burn) / n_thin
      beta_chain[,index_element] <- Beta
      u_chain[,index_element] <- U
      lambda2_chain[,index_element] <- Lambda2
      b_chain[index_element] <- b
      a_chain[index_element] <- a
    }	
  }
  
  list(Beta = beta_chain, U = u_chain, Lambda2 = lambda2_chain, b = b_chain, a = a_chain)
}

FH_sampler <- function(X, Y, D, nburn, nsim, nthin)
{
	# parameters
	m <- nrow(X)
	np <- ncol(X)
	
	# initial values
	Theta <- rep(1,m)
	Beta <- rep(1, np)
	Sigma2 <- 1
	
	# chain containers
	Theta.chain <- array(0, dim=c(m,nsim/nthin))
	Beta.chain <- array(0, dim=c(np,nsim/nthin))
	Sigma2.chain <- rep(0, nsim/nthin)
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		
		# update Sigma2
		Sigma2 <- 1/rgamma(1,shape=m/2-1, rate=sum((Theta-X%*%Beta)^2)/2)
		
		# update Beta
		mean.Beta <- solve(t(X)%*%X)%*%t(X)%*%Theta
		var.Beta <- Sigma2*solve(t(X)%*%X)
		Beta <- mvrnorm(1, mu=mean.Beta, Sigma=var.Beta)
		
		# update Theta
		var.Theta <- 1/(1/Sigma2+1/D)
		mean.Theta <- var.Theta*(Y/D+X%*%Beta/Sigma2)
		Theta <- rnorm(m, mean=mean.Theta, sd=sqrt(var.Theta))
		
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Sigma2.chain[(index-nburn)/nthin] <- Sigma2
			Theta.chain[,(index-nburn)/nthin] <- Theta
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, Theta.chain=Theta.chain, Sigma2.chain=Sigma2.chain)
}

TPB_sampler <- function(p, q, X, Y, D, nburn, nsim, nthin)
{
	# parameters
	m <- nrow(X)
	np <- ncol(X)
	
	a <- 10^-10
	b <- 10^-10
	phi <- 1
	
	# initial values
	Lambda2 = rep(1,m)
	Tau2 = 1
	Z = rep(1,m)
	Beta = rep(1,np)
	U = rep(1,m)
	
	# chain containers
	Lambda2.chain = array(0, dim=c(m,nsim/nthin))
	Tau2.chain = rep(0, nsim/nthin)
	Z.chain = array(0,dim=c(m,nsim/nthin))
	Beta.chain = array(0, dim=c(np,nsim/nthin))
	U.chain = array(0, dim=c(m,nsim/nthin))
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		# update Tau2
		Tau2 <- 1/rgamma(1,shape=(a+m)/2,rate=sum(U^2/Lambda2)/2+b/2)
		
		# update Z
		Z <- rgamma(m, shape=p+q, rate=phi+Lambda2)
		
		# update Lambda2
		for (i in 1:m) Lambda2[i] = rgig(1,p-0.5,U[i]^2/Tau2,2*Z[i])
		
		# update Beta
		Xstd <- X/sqrt(D)
		sigma <- solve(t(Xstd)%*%Xstd)
	    mean <- apply(X*(Y-U)/D,2,sum)
		mean <- sigma%*%mean
		Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
		
		# update U
		sigma2 <- 1/(1/D+1/Lambda2/Tau2)
		mean <- sigma2 * (Y-X%*%Beta)/D
		U <- rnorm(m, mean=mean, sd=sqrt(sigma2))
	
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Tau2.chain[(index-nburn)/nthin] <- Tau2
			Lambda2.chain[,(index-nburn)/nthin] <- Lambda2
			Z.chain[,(index-nburn)/nthin] <- Z
			U.chain[,(index-nburn)/nthin] <- U
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, U.chain=U.chain, Tau2.chain=Tau2.chain, Lambda2.chain=Lambda2.chain)
}

NG_sampler <- function(a, X, Y, D, nburn, nsim, nthin)
{
	m <- nrow(X)
	np <- ncol(X)
	
	b <- 1
	p <- 10^-10
	q <- 10^-10
	
	Lambda2 <- rep(1,m)
	Tau2 <- 1
	Beta <- rep(1,np)
	U <- rep(1,m)
		
	Lambda2.chain <- array(0, dim=c(m,nsim/nthin))
	Tau2.chain <- rep(0, nsim/nthin)
	Beta.chain <- array(0, dim=c(np,nsim/nthin))
	U.chain <- array(0, dim=c(m,nsim/nthin))
	
	for (index in 1:(nsim+nburn))
	{
		if (index%%10000==0) cat(index,"\n")
		# update Tau2
		Tau2 <- 1/rgamma(1,shape=(p+m)/2,rate=sum(U^2/Lambda2)/2+q/2)
			
		# update Lambda2
		for (i in 1:m) Lambda2[i] <- rgig(1,a-0.5,U[i]^2/Tau2,2*b)
		
		# update Beta
		Xstd <- X/sqrt(D)
		sigma <- solve(t(Xstd)%*%Xstd)
	    	mean <- apply(X*(Y-U)/D,2,sum)
		mean <- sigma%*%mean
		Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
		
		# update U
		sigma2 <- 1/(1/D+1/Lambda2/Tau2)
		mean <- sigma2 * (Y-X%*%Beta)/D
		U <- rnorm(m, mean=mean, sd=sqrt(sigma2))
	
		if (index > nburn && (index-nburn)%%nthin==0)
		{
			Tau2.chain[(index-nburn)/nthin] <- Tau2
			Lambda2.chain[,(index-nburn)/nthin] <- Lambda2
			U.chain[,(index-nburn)/nthin] <- U
			Beta.chain[,(index-nburn)/nthin] <- Beta		
		}	
	}
	
	list(Beta.chain=Beta.chain, U.chain=U.chain, Tau2.chain=Tau2.chain, Lambda2.chain=Lambda2.chain)
	
}

dev_measure <- function(thetahat, theta0)
{
	AAD <- mean(abs(thetahat - theta0))
	ASD <- mean((thetahat-theta0)^2)
	ARB <- mean(abs(thetahat-theta0)/abs(theta0))
	ASRB <- mean((thetahat-theta0)^2/theta0^2)
	
	out <- c(AAD, ASD, ARB, ASRB)
	names(out) <- c("AAD", "ASD", "ARB", "ASRB")
	return(out)
}



# DIC = D(\bar \theta) + 2*p_D where p_D = \bar D - D(\bar theta) and D(\theta) = -2*log-likelihood
dic <- function(Y, Theta.chain, D)
{
  n <- dim(Theta.chain)[2]
  Theta.est <- apply(Theta.chain, 1, mean)
  D.chain <- rep(0, n)
  for (i in 1:n) D.chain[i] <- -2*sum(dnorm(Y, mean=Theta.chain[,i], sd=sqrt(D), log=TRUE))
  D.thetabar <- -2*sum(dnorm(Y, mean=Theta.est, sd=sqrt(D), log=TRUE))
  
  2*mean(D.chain) - D.thetabar
}


