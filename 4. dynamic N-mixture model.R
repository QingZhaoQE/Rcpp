#==================================================================================================
# Dynamic N-mixture model
# code for simulating data, defining MCMC algorithm, implementing the model, and check convergence
# written by Qing Zhao, 2024 in Colorado
#==================================================================================================

#===============
# Simulate data
#===============
library(boot) # for using logit and inv.logit functions

# Basic values
nsite <- 360 # number of sites
nyear <- 20  # number of years
nreps <- 5   # number of within-season replicates
ncovs <- 2   # number of environmental covariates
npcvs <- 3   # number of observational covariates

beta0 <- c(1.8, 1, -0.8)          # intercept and slopes for initial abundance
beta_rho <- c(0, -0.8, 0.6, -0.4) # intercept and slopes for population growth
delta <- 0.5                      # expectation of the number of immigrants
alpha <- c(0.4, -0.6, 0.4, -0.2)  # intercept and slopes for detection probability

# Simulated environmental covariates
x <- matrix(rnorm(nsite * nyear * ncovs, 0, 1), nsite * nyear, ncovs) # environmental covariates

# Simulate data
lambda <- matrix(, nsite, nyear) # expectation of abundance
N <- matrix(, nsite, nyear) # abundance
lambda[,1] <- exp(cbind(1,x[1:nsite,]) %*% beta0) # expectation of initial abundance
N[,1] <- rpois(nsite, lambda[,1]) # initial abundance

rho <- matrix(0, nsite, nyear-1) # population growth rate
for (t in 2:nyear) {
  rho[,t-1] <- exp(cbind(1, (N[,t-1]-lambda[,1])/lambda[,1], x[(t-1)*nsite+c(1:nsite),]) %*% beta_rho)
  lambda[,t] <- N[,t-1] * rho[,t-1] + delta
  N[,t] <- rpois(nsite, lambda[,t])
} # t

w <- matrix(rnorm(nsite * nyear * nreps * npcvs, 0, 1), nsite * nyear, nreps * npcvs) # observational covariates

p <- matrix(, nsite * nyear, nreps) # detection probability
for (t in 1:nyear) {
  for (j in 1:nreps) {
    p[(t-1)*nsite+c(1:nsite),j] <- inv.logit(cbind(1,w[(t-1)*nsite+c(1:nsite),(j-1)*npcvs+c(1:npcvs)]) %*% alpha)
  } # j
} # t

y <- matrix(, nsite * nyear, nreps) # detection/non-detection data
for (t in 1:nyear) {
  for (j in 1:nreps) {
    y[(t-1)*nsite+c(1:nsite),j] <- rbinom(nsite, N[,t], p[(t-1)*nsite+c(1:nsite),j])
  } # j
} # t

start_time <- Sys.time() # start time of computing

#===============================
# Define MCMC algorithm in Rcpp
#===============================
library(Rcpp)

cppFunction('

List dynam_Nmix_mcmc(IntegerMatrix y, NumericMatrix x, NumericMatrix w, int nmcmc, 
                     int nsite, int nyear, int nreps, int ncovs, int npcvs) {
  
  // Setup variables
  NumericMatrix beta0_save(ncovs+1, nmcmc);
  NumericMatrix beta_rho_save(ncovs+2, nmcmc);
  NumericVector delta_save(nmcmc); 
  NumericMatrix alpha_save(npcvs+1, nmcmc);
  IntegerVector ymax(nsite*nyear); 
  for (int i=0; i<(nsite*nyear); ++i) {
    ymax(i) = max(y(i,_));
  } // i

  // Priors
  NumericVector beta0_mean(ncovs+1);
  double beta0_sd = 10.0;
  NumericVector beta_rho_mean(ncovs+2);
  double beta_rho_sd = 10.0;
  double log_delta_mean = 0.0; 
  double log_delta_sd = 10.0; 
  NumericVector alpha_mean(npcvs+1);
  double alpha_sd = 10.0;
  
  // Starting values
  NumericVector beta0(ncovs+1);
  NumericVector beta_rho(ncovs+2);
  double delta = 1; 
  NumericVector alpha(npcvs+1);
  IntegerMatrix N(nsite, nyear);
  for (int i=0; i<nsite; ++i) {
    for (int t=0; t<nyear; ++t) {
      N(i,t) = ymax(t*nsite+i) * 2 + 2;
    } // t
  } // i
  NumericVector lambda0(nsite, 1.0); 
  NumericMatrix rho(nsite, nyear-1);
  for (int i=0; i<nsite; ++i) {
    for (int t=0; t<(nyear-1); ++t) {
      rho(i,t) = 1.0; 
    } // t
  } // i
  NumericMatrix p(nsite*nyear, nreps);
  for (int i=0; i<(nsite*nyear); ++i) {
    for (int j=0; j<nreps; ++j) {
      p(i,j) = 0.5; 
    } // j
  } // i

  // Tuning factors
  double beta0_tune = 0.04;
  double beta_rho_tune = 0.02;
  double delta_tune = 0.05;
  NumericVector alpha_tune = {0.02, 0.01, 0.01, 0.01}; 
  int N_tune = 1;

  for (int k=0; k<nmcmc; ++k){

    //Sample N
    IntegerMatrix N_star(nsite, nyear);
    IntegerMatrix N_keep(nsite, nyear); 

    for (int i=0; i<nsite; ++i) {
      double mh1N = 0; 
      double mh2N = 0; 

      N_star(i,0) = R::rpois(N(i,0) + N_tune); 
      mh1N += R::dpois(N(i,0), N_star(i,0) + N_tune, true) + R::dpois(N_star(i,0), lambda0(i), true);
      mh2N += R::dpois(N_star(i,0), N(i,0) + N_tune, true) + R::dpois(N(i,0), lambda0(i), true);
      for (int j=0; j<nreps; ++j) {
        mh1N += R::dbinom( y(i,j), N_star(i,0), p(i,j), true ); 
        mh2N += R::dbinom( y(i,j), N(i,0), p(i,j), true ); 
      } // j

      double mhN = exp(mh1N - mh2N);
      if ( mhN > R::runif(0,1) ) {
      if ( N_star(i,0) >= ymax(i) ) {
        N_keep(i,0) = 1;
      }
      }

      for (int t=1; t<nyear; ++t) {
        double mh1N = 0; 
        double mh2N = 0; 

        N_star(i,t) = R::rpois(N(i,t) + N_tune); 
        mh1N += R::dpois(N(i,t), N_star(i,t) + N_tune, true) + R::dpois(N_star(i,t), N(i,t-1) * rho(i,t-1) + delta, true);
        mh2N += R::dpois(N_star(i,t), N(i,t) + N_tune, true) + R::dpois(N(i,t), N(i,t-1) * rho(i,t-1) + delta, true);
        for (int j=0; j<nreps; ++j) {
          mh1N += R::dbinom( y(t*nsite+i,j), N_star(i,t), p(t*nsite+i,j), true ); 
          mh2N += R::dbinom( y(t*nsite+i,j), N(i,t), p(t*nsite+i,j), true ); 
        } // j

        double mhN = exp(mh1N - mh2N);
        if ( mhN > R::runif(0,1) ) {
        if ( N_star(i,t) >= ymax(t*nsite+i) ) {
          N_keep(i,t) = 1;
        }
        }
      } // t
    } // i

    for (int i=0; i<nsite; ++i) {
      for (int t=0; t<nyear; ++t) {
        if ( N_keep(i,t) == 1 ) {
          N(i,t) = N_star(i,t); 
        }
      } // t
    } // i

    //Sample beta0
    double mh1B0 = 0; 
    double mh2B0 = 0; 

    NumericVector beta0_star(ncovs+1);
    for (int l=0; l<(ncovs+1); ++l) {
      beta0_star[l] = R::rnorm(beta0[l], beta0_tune);
      mh1B0 += R::dnorm(beta0_star[l], beta0_mean[l], beta0_sd, true);
      mh2B0 += R::dnorm(beta0[l], beta0_mean[l], beta0_sd, true);
    } // l

    NumericVector lambda0_star(nsite);
    for (int i=0; i<nsite; ++i) {
      double log_lambda0_star = beta0_star[0];
      for (int l=1; l<(ncovs+1); ++l) {
        log_lambda0_star += beta0_star[l] * x(i,l-1);
      } // l
      lambda0_star(i) = exp(log_lambda0_star);
      mh1B0 += R::dpois( N(i,0), lambda0_star(i), true );
      mh2B0 += R::dpois( N(i,0), lambda0(i), true );
    } // i

    double mhB0 = exp(mh1B0 - mh2B0);
    if ( mhB0 > R::runif(0,1) ) {
      beta0 = clone(beta0_star);
      lambda0 = clone(lambda0_star); 
    }

    //Sample beta_rho
    double mh1BP = 0; 
    double mh2BP = 0; 

    NumericVector beta_rho_star(ncovs+2);
    for (int l=0; l<(ncovs+2); ++l) {
      beta_rho_star[l] = R::rnorm(beta_rho[l], beta_rho_tune);
      mh1BP += R::dnorm(beta_rho_star[l], beta_rho_mean[l], beta_rho_sd, true);
      mh2BP += R::dnorm(beta_rho[l], beta_rho_mean[l], beta_rho_sd, true);
    } // l

    NumericMatrix rho_star(nsite, nyear-1); 
    for (int i=0; i<nsite; ++i) {
      for (int t=1; t<nyear; ++t) {
        double log_rho_star = beta_rho_star[0] + beta_rho_star[1] * (N(i,t-1) - lambda0(i)) / lambda0(i);
        for (int l=2; l<(ncovs+2); ++l) {
          log_rho_star += beta_rho_star[l] * x(t*nsite+i,l-2);
        } // l
        rho_star(i,t-1) = exp(log_rho_star); 
        mh1BP += R::dpois( N(i,t), N(i,t-1) * rho_star(i,t-1) + delta, true );
        mh2BP += R::dpois( N(i,t), N(i,t-1) * rho(i,t-1) + delta, true );
      } // t
    } // i

    double mhBP = exp(mh1BP - mh2BP);
    if ( mhBP > R::runif(0,1) ) {
      beta_rho = clone(beta_rho_star);
      rho = clone(rho_star);
    }

    //Sample delta
    double log_delta_star = R::rnorm(log(delta), delta_tune); 
    double delta_star = exp(log_delta_star); 
    double mh1D = R::dnorm(log(delta_star), log_delta_mean, log_delta_sd, true);
    double mh2D = R::dnorm(log(delta), log_delta_mean, log_delta_sd, true);
    for (int i=0; i<nsite; ++i) {
      for (int t=1; t<nyear; ++t) {
        mh1D += R::dpois( N(i,t), N(i,t-1) * rho(i,t-1) + delta_star, true );
        mh2D += R::dpois( N(i,t), N(i,t-1) * rho(i,t-1) + delta, true );
      } // t
    } // i
    double mhD = exp(mh1D - mh2D);
    if ( mhD > R::runif(0,1) ) {
      delta = delta_star;
    }

    //Sample alpha
    double mh1A = 0; 
    double mh2A = 0; 

    NumericVector alpha_star(npcvs+1);
    for (int l=0; l<(npcvs+1); ++l) {
      alpha_star[l] = R::rnorm(alpha[l], alpha_tune[l]);
      mh1A += R::dnorm(alpha_star[l], alpha_mean[l], alpha_sd, true);
      mh2A += R::dnorm(alpha[l], alpha_mean[l], alpha_sd, true);
    } // l

    NumericMatrix p_star(nsite*nyear, nreps); 
    for (int i=0; i<nsite; ++i) {
      for (int t=0; t<nyear; ++t) {
        for (int j=0; j<nreps; ++j) {
          double logit_p_star = alpha_star[0];
          for (int l=1; l<(npcvs+1); ++l) {
            logit_p_star += alpha_star[l] * w(t*nsite+i, j*npcvs+l-1);
          } // l
          p_star(t*nsite+i,j) = R::plogis(logit_p_star, 0, 1, true, false);
          mh1A += R::dbinom( y(t*nsite+i,j), N(i,t), p_star(t*nsite+i,j), true );
          mh2A += R::dbinom( y(t*nsite+i,j), N(i,t), p(t*nsite+i,j), true );
        } // j
      } // t
    } // i

    double mhA = exp(mh1A - mh2A);
    if ( mhA > R::runif(0,1) ) {
      alpha = clone(alpha_star);
      p = clone(p_star);
    }

    //Save samples
    beta0_save(_,k) = beta0;
    beta_rho_save(_,k) = beta_rho;
    delta_save(k) = delta; 
    alpha_save(_,k) = alpha;

  } // k

  // Write output
  List out = List::create(Named("beta0_save") = beta0_save, 
                          Named("beta_rho_save") = beta_rho_save, 
                          Named("delta_save") = delta_save, 
                          Named("alpha_save") = alpha_save);

  return(out);
}

')

#==========
# Run MCMC
#==========
nmcmc <- 50000 # number of iterations
out <- dynam_Nmix_mcmc(y, x, w, nmcmc, nsite, nyear, nreps, ncovs, npcvs)

end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta[1]^"[0]"), expression(beta[2]^"[0]"))
ylab_beta_rho <- c(expression(beta[0]^""["["*rho*"]"]), expression(beta[1]^""["["*rho*"]"]), 
                   expression(beta[2]^""["["*rho*"]"]), expression(beta[3]^""["["*rho*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))

par(mfrow=c(4,3))
par(mar=c(1,3,3,1))
par(oma=c(4,2.5,0,0))

for (i in 1:(ncovs+1)) {
  tt <- out$beta0_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0[i] - yint * 1
  ymax <- beta0[i] + yint * 1
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='blue')
  abline(h=beta0[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i], cex.main=2, line=1.2)
} # i

for (i in 1:(ncovs+2)) {
  tt <- out$beta_rho_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_rho[i] - yint * 4
  ymax <- beta_rho[i] + yint * 4
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='blue')
  abline(h=beta_rho[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_rho[i], cex.main=2, line=1.2)
} # i

tt <- out$delta_save
yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
yint <- ifelse (yint < 0.1, 0.1, yint)
ymin <- delta - yint * 4
ymax <- delta + yint * 4
plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
box()
lines(tt, col='blue')
abline(h=delta, col='grey16', lwd=1.5)
title(main=expression(delta), cex.main=2, line=1.2)

for (i in 1:(npcvs+1)) {
  tt <- out$alpha_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- alpha[i] - yint * 1
  ymax <- alpha[i] + yint * 1
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='blue')
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1.2)
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=0.4, outer=T)


