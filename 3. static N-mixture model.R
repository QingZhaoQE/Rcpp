#==================================================================================================
# Static N-mixture model
# code for simulating data, defining MCMC algorithm, implementing the model, and check convergence
# written by Qing Zhao, 2024 in Colorado
#==================================================================================================

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

# Basic values
nsite <- 360 # number of sites
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 2   # number of observational covariates

beta <- c(1.8, 1, -0.8)    # intercept and slopes for true abundance
alpha <- c(0.3, -0.6, 0.3) # intercept and slopes for detection probability

# Simulate data
x <- matrix(rnorm(nsite*ncovs, 0, 1), nsite, ncovs) # environmental covariates
lambda <- exp(cbind(1,x) %*% beta) # expectation of abundance
N <- rpois(nsite, lambda) # true abundance

w <- matrix(rnorm(nsite * nreps * npcvs, 0, 1), nsite * nreps, npcvs) # observational covariates

p <- matrix(, nsite, nreps) # detection probability
for (i in 1:nsite) {
  for (j in 1:nreps) {
    p[i,j] <- inv.logit(c(1,w[(i-1)*nreps+j,]) %*% alpha)
  } # j
} # i

y <- matrix(, nsite, nreps) # detection/non-detection data
for (i in 1:nsite) {
  for (j in 1:nreps) {
    y[i,j] <- rbinom(1, N[i], p[i,j])
  } # j
} # i

start_time <- Sys.time() # start time of computing

#===============================
# Define MCMC algorithm in Rcpp
#===============================
library(Rcpp)

cppFunction('

List Nmix_mcmc(IntegerMatrix y, NumericMatrix x, NumericMatrix w, int nmcmc) {
  
  // Setup variables
  int nsite = y.nrow();
  int nreps = y.ncol();
  int ncovs = x.ncol();
  int npcvs = w.ncol();
  IntegerVector ymax(nsite); 
  for (int i=0; i<nsite; ++i) {
    ymax(i) = max(y(i,_));
  } // i

  NumericMatrix beta_save(ncovs+1, nmcmc);
  NumericMatrix alpha_save(npcvs+1, nmcmc);

  // Priors
  NumericVector beta_mean(ncovs+1);
  double beta_sd = 10.0;
  NumericVector alpha_mean(npcvs+1);
  double alpha_sd = 10.0;
  
  // Starting values
  NumericVector beta(ncovs+1);
  NumericVector alpha(npcvs+1);
  IntegerVector N(nsite);
  for (int i=0; i<nsite; ++i) {
    N[i] = ymax[i] * 2 + 2;
  } // i

  // Tuning factors
  double beta_tune = 0.05;
  double alpha_tune = 0.05;
  int N_tune = 1;

  for (int k=0; k<nmcmc; ++k){

    //Sample N
    IntegerVector N_star(nsite);
    for (int i=0; i<nsite; ++i) {
      N_star[i] = R::rpois(N[i] + N_tune); 
      double log_lam = beta[0];
      for (int j=1; j<(ncovs+1); ++j) {
        log_lam += beta[j] * x(i,j-1);
      } // j
      double lam = exp(log_lam);

      double mh1N = 0; 
      double mh2N = 0; 
      mh1N += R::dpois(N[i], N_star[i] + N_tune, true) + R::dpois(N_star[i], lam, true);
      mh2N += R::dpois(N_star[i], N[i] + N_tune, true) + R::dpois(N[i], lam, true);
      for (int j=0; j<nreps; ++j) {
        double logit_p = alpha[0];
        for (int l=1; l<(npcvs+1); ++l) {
          logit_p += alpha[l] * w(i*nreps+j,l-1);
        } // l
        double p = R::plogis(logit_p, 0, 1, true, false);
        mh1N += R::dbinom( y(i,j), N_star[i], p, true ); 
        mh2N += R::dbinom( y(i,j), N[i], p, true ); 
      } // j

      double mhN = exp(mh1N - mh2N);
      if ( mhN > R::runif(0,1) ) {
      if ( N_star[i] >= ymax[i] ) {
        N[i] = N_star[i];
      }
      }
    } // i

    //Sample beta
    NumericVector beta_star(ncovs+1);
    for (int j=0; j<(ncovs+1); ++j) {
      beta_star[j] = R::rnorm(beta[j], beta_tune);
    } // j
    double mh1B = 0; 
    double mh2B = 0; 
    for(int j=0; j<(ncovs+1); ++j) {
      mh1B += R::dnorm(beta_star[j], beta_mean[j], beta_sd, true);
      mh2B += R::dnorm(beta[j], beta_mean[j], beta_sd, true);
    } // j
    for (int i=0; i<nsite; ++i) {
      double log_lam_star = beta_star[0];
      double log_lam = beta[0];
      for (int j=1; j<(ncovs+1); ++j) {
        log_lam_star += beta_star[j] * x(i,j-1);
        log_lam += beta[j] * x(i,j-1);
      } // j
      double lam_star = exp(log_lam_star); 
      double lam = exp(log_lam); 

      mh1B += R::dpois( N(i), lam_star, true );
      mh2B += R::dpois( N(i), lam, true );
    } // i
    double mhB = exp(mh1B - mh2B);
    if ( mhB > R::runif(0,1) ) {
      beta = clone(beta_star);
    }

    //Sample alpha
    NumericVector alpha_star(npcvs+1);
    for (int j=0; j<(npcvs+1); ++j) {
      alpha_star[j] = R::rnorm(alpha[j], alpha_tune);
    } // j
    double mh1A = 0; 
    double mh2A = 0; 
    for(int j=0; j<(npcvs+1); ++j) {
      mh1A += R::dnorm(alpha_star[j], alpha_mean[j], alpha_sd, true);
      mh2A += R::dnorm(alpha[j], alpha_mean[j], alpha_sd, true);
    } // j
    for (int i=0; i<nsite; ++i) {
      for (int j=0; j<nreps; ++j) {
        double logit_p_star = alpha_star[0];
        double logit_p = alpha[0];
        for (int l=1; l<(npcvs+1); ++l) {
          logit_p_star += alpha_star[l] * w(i*nreps+j,l-1);
          logit_p += alpha[l] * w(i*nreps+j,l-1);
        } // j

        double p_star = R::plogis(logit_p_star, 0, 1, true, false);
        double p = R::plogis(logit_p, 0, 1, true, false);
  
        mh1A += R::dbinom( y(i,j), N[i], p_star, true );
        mh2A += R::dbinom( y(i,j), N[i], p, true );
      } // j
    } // i
    double mhA = exp(mh1A - mh2A);
    if ( mhA > R::runif(0,1) ) {
      alpha = clone(alpha_star);
    }

    //Save samples
    beta_save(_,k) = beta;
    alpha_save(_,k) = alpha;

  } //k

  // Write output
  List out = List::create(Named("beta_save") = beta_save, Named("alpha_save") = alpha_save);

  return(out);
}

')

#==========
# Run MCMC
#==========
nmcmc <- 50000 # number of iterations
out <- Nmix_mcmc(y, x, w, nmcmc)

end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

ylab_beta <- c(expression(beta[0]), expression(beta[1]), expression(beta[2]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]))

par(mfrow=c(2,3))
par(mar=c(1,3,3,1))
par(oma=c(4,2.5,0,0))

for (i in 1:(ncovs+1)) {
  tt <- out$beta_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta[i] - yint * 8
  ymax <- beta[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='red')
  abline(h=beta[i], col='grey16', lwd=1.5)
  title(main=ylab_beta[i], cex.main=2, line=1.2)
} # i

for (i in 1:(npcvs+1)) {
  tt <- out$alpha_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i] - yint * 8
  ymax <- alpha[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='red')
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1.2)
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=0.4, outer=T)



