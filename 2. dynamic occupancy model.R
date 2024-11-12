#==================================================================================================
# Dynamic occupancy model
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
npcvs <- 2   # number of observational covariates

beta0 <- c(0.4, 1, -0.4)         # intercept and slopes for initial occupancy
beta_phi   <- c(-0.2, 0.6, -0.3) # intercept and slopes for persistence
beta_gamma <- c(0.2, 0.4, -0.2)  # intercept and slopes for colonization
alpha <- c(0.3, -0.6, 0.3)       # intercept and slopes for detection probability

# Simulated environmental covariates
x <- matrix(rnorm(nsite * nyear * npcvs, 0, 1), nsite * nyear, npcvs) # environmental covariates

# Simulate data
psi <- matrix(, nsite, nyear) # occupancy probability
z <- matrix(, nsite, nyear) # occupancy status
psi[,1] <- inv.logit(cbind(1,x[1:nsite,]) %*% beta0) # initial occupancy probability
z[,1] <- rbinom(nsite, 1, psi[,1]) # initial occupancy status

phi   <- matrix(, nsite, nyear-1) # probability of persistence
gamma <- matrix(, nsite, nyear-1) # probability of colonization
for (t in 2:nyear) {
  phi  [,t-1] <- inv.logit(cbind(1,x[(t-1)*nsite+c(1:nsite),]) %*% beta_phi  )
  gamma[,t-1] <- inv.logit(cbind(1,x[(t-1)*nsite+c(1:nsite),]) %*% beta_gamma)
  psi[,t] <- z[,t-1] * phi[,t-1] + (1 - z[,t-1]) * gamma[,t-1]
  z[,t] <- rbinom(nsite, 1, psi[,t])
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
    y[(t-1)*nsite+c(1:nsite),j] <- rbinom(nsite, 1, z[,t] * p[(t-1)*nsite+c(1:nsite),j])
  } # j
} # t

start_time <- Sys.time() # start time of computing

#===============================
# Define MCMC algorithm in Rcpp
#===============================
library(Rcpp)

cppFunction('

List dynam_occupancy_mcmc(IntegerMatrix y, NumericMatrix x, NumericMatrix w, int nmcmc, 
                          int nsite, int nyear, int nreps, int ncovs, int npcvs) {
  
  // Setup variables
  IntegerVector ysum(nsite*nyear); 
  for (int i=0; i<(nsite*nyear); ++i) {
    for (int j=0; j<nreps; ++j) {
      ysum(i) += y(i,j);
    } // j
  } // i

  NumericMatrix beta0_save(ncovs+1, nmcmc);
  NumericMatrix beta_phi_save(ncovs+1, nmcmc);
  NumericMatrix beta_gamma_save(ncovs+1, nmcmc);
  NumericMatrix alpha_save(npcvs+1, nmcmc);

  // Priors
  NumericVector beta0_mean(ncovs+1);
  double beta0_sd = 10.0;
  NumericVector beta_phi_mean(ncovs+1);
  double beta_phi_sd = 10.0;
  NumericVector beta_gamma_mean(ncovs+1);
  double beta_gamma_sd = 10.0;
  NumericVector alpha_mean(npcvs+1);
  double alpha_sd = 10.0;
  
  // Starting values
  NumericVector beta0(ncovs+1);
  NumericVector beta_phi(ncovs+1);
  NumericVector beta_gamma(ncovs+1);
  NumericVector alpha(npcvs+1);
  IntegerMatrix z(nsite, nyear);
  for (int i=0; i<nsite; ++i) {
    for (int t=0; t<nyear; ++t) {
      if (ysum[t*nsite+i] > 0) {
        z(i,t) = 1;
      }
    } // t
  } // i

  // Tuning factors
  double beta0_tune = 0.2;
  double beta_phi_tune = 0.05;
  double beta_gamma_tune = 0.05;
  double alpha_tune = 0.035;

  for (int k=0; k<nmcmc; ++k){

    //Sample beta0
    NumericVector beta0_star(ncovs+1);
    for (int j=0; j<(ncovs+1); ++j) {
      beta0_star[j] = R::rnorm(beta0[j], beta0_tune);
    } // j
    double mh1B0 = 0; 
    double mh2B0 = 0; 
    for(int j=0; j<(ncovs+1); ++j) {
      mh1B0 += R::dnorm(beta0_star[j], beta0_mean[j], beta0_sd, true);
      mh2B0 += R::dnorm(beta0[j], beta0_mean[j], beta0_sd, true);
    } // j
    for (int i=0; i<nsite; ++i) {
      double logit_psi_star = beta0_star[0];
      double logit_psi = beta0[0];
      for (int j=1; j<(ncovs+1); ++j) {
        logit_psi_star += beta0_star[j] * x(i,j-1);
        logit_psi += beta0[j] * x(i,j-1);
      } // j

      double psi_star = R::plogis(logit_psi_star, 0, 1, true, false);
      mh1B0 += R::dbinom( z(i,0), 1, psi_star, true );

      double psi = R::plogis(logit_psi, 0, 1, true, false);
      mh2B0 += R::dbinom( z(i,0), 1, psi, true );
    } // i
    double mhB0 = exp(mh1B0 - mh2B0);
    if ( mhB0 > R::runif(0,1) ) {
      beta0 = clone(beta0_star);
    }

    //Sample beta_phi
    NumericVector beta_phi_star(ncovs+1);
    for (int j=0; j<(ncovs+1); ++j) {
      beta_phi_star[j] = R::rnorm(beta_phi[j], beta_phi_tune);
    } // j
    double mh1BP = 0; 
    double mh2BP = 0; 
    for(int j=0; j<(ncovs+1); ++j) {
      mh1BP += R::dnorm(beta_phi_star[j], beta_phi_mean[j], beta_phi_sd, true);
      mh2BP += R::dnorm(beta_phi[j], beta_phi_mean[j], beta_phi_sd, true);
    } // j
    for (int i=0; i<nsite; ++i) {
      for (int t=1; t<nyear; ++t) {
        double logit_phi_star = beta_phi_star[0];
        double logit_phi = beta_phi[0];
        double logit_gamma = beta_gamma[0];
        for (int j=1; j<(ncovs+1); ++j) {
          logit_phi_star += beta_phi_star[j] * x(t*nsite+i, j-1);
          logit_phi += beta_phi[j] * x(t*nsite+i, j-1);
          logit_gamma += beta_gamma[j] * x(t*nsite+i, j-1);
        } // j

        double gamma = R::plogis(logit_gamma, 0, 1, true, false);
        double phi_star = R::plogis(logit_phi_star, 0, 1, true, false);
        mh1BP += R::dbinom( z(i,t), 1, z(i,t-1) * phi_star + (1 - z(i,t-1)) * gamma, true );

        double phi = R::plogis(logit_phi, 0, 1, true, false);
        mh2BP += R::dbinom( z(i,t), 1, z(i,t-1) * phi + (1 - z(i,t-1)) * gamma, true );
      } // t
    } // i
    double mhBP = exp(mh1BP - mh2BP);
    if ( mhBP > R::runif(0,1) ) {
      beta_phi = clone(beta_phi_star);
    }

    //Sample beta_gamma
    NumericVector beta_gamma_star(ncovs+1);
    for (int j=0; j<(ncovs+1); ++j) {
      beta_gamma_star[j] = R::rnorm(beta_gamma[j], beta_gamma_tune);
    } // j
    double mh1BG = 0; 
    double mh2BG = 0; 
    for(int j=0; j<(ncovs+1); ++j) {
      mh1BG += R::dnorm(beta_gamma_star[j], beta_gamma_mean[j], beta_gamma_sd, true);
      mh2BG += R::dnorm(beta_gamma[j], beta_gamma_mean[j], beta_gamma_sd, true);
    } // j
    for (int i=0; i<nsite; ++i) {
      for (int t=1; t<nyear; ++t) {
        double logit_gamma_star = beta_gamma_star[0];
        double logit_gamma = beta_gamma[0];
        double logit_phi = beta_phi[0];
        for (int j=1; j<(ncovs+1); ++j) {
          logit_gamma_star += beta_gamma_star[j] * x(t*nsite+i, j-1);
          logit_gamma += beta_gamma[j] * x(t*nsite+i, j-1);
          logit_phi += beta_phi[j] * x(t*nsite+i, j-1);
        } // j

        double phi = R::plogis(logit_phi, 0, 1, true, false);
        double gamma_star = R::plogis(logit_gamma_star, 0, 1, true, false);
        mh1BG += R::dbinom( z(i,t), 1, z(i,t-1) * phi + (1 - z(i,t-1)) * gamma_star, true );

        double gamma = R::plogis(logit_gamma, 0, 1, true, false);
        mh2BG += R::dbinom( z(i,t), 1, z(i,t-1) * phi + (1 - z(i,t-1)) * gamma, true );
      } // t
    } // i
    double mhBG = exp(mh1BG - mh2BG);
    if ( mhBG > R::runif(0,1) ) {
      beta_gamma = clone(beta_gamma_star);
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
      for (int t=0; t<nyear; ++t) {
        for (int j=0; j<nreps; ++j) {
          if (z(i,t) == 1) {
            double logit_p_star = alpha_star[0];
            double logit_p = alpha[0];
            for (int l=1; l<(npcvs+1); ++l) {
              logit_p_star += alpha_star[l] * w(t*nsite+i, j*npcvs+l-1);
              logit_p += alpha[l] * w(t*nsite+i, j*npcvs+l-1);
            } // l

            double p_star = R::plogis(logit_p_star, 0, 1, true, false);
            mh1A += R::dbinom( y(t*nsite+i,j), 1, p_star, true );
    
            double p = R::plogis(logit_p, 0, 1, true, false);
            mh2A += R::dbinom( y(t*nsite+i,j), 1, p, true );
          }
        } // j
      } // t
    } // i
    double mhA = exp(mh1A - mh2A);
    if ( mhA > R::runif(0,1) ) {
      alpha = clone(alpha_star);
    }

    //Sample z
    for (int i=0; i<nsite; ++i) {
      double qprod = 1; 
      for (int j=0; j<nreps; ++j) {
        double logit_p = alpha[0];
        for (int l=1; l<(npcvs+1); ++l) {
          logit_p += alpha[l] * w(i, j*npcvs+l-1);
        } // j
        double p = R::plogis(logit_p, 0, 1, true, false);
        qprod = qprod * ( 1 - p );
      } // j
      double logit_s_temp = beta0[0];
      for (int j=1; j<(ncovs+1); ++j) {
        logit_s_temp += beta0[j] * x(i,j-1);
      } // j
      double s_temp = R::plogis(logit_s_temp, 0, 1, true, false);
      double num_temp = s_temp * qprod;
      double psi_temp = num_temp / (num_temp + (1 - s_temp));
      if (ysum(i) == 0) {
        z(i,0) = R::rbinom(1, psi_temp);
      }
      for (int t=1; t<nyear; ++t) {
        double qprod = 1; 
        for (int j=0; j<nreps; ++j) {
          double logit_p = alpha[0];
          for (int l=1; l<(npcvs+1); ++l) {
            logit_p += alpha[l] * w(t*nsite+i, j*npcvs+l-1);
          } // j
          double p = R::plogis(logit_p, 0, 1, true, false);
          qprod = qprod * ( 1 - p );
        } // j
        double logit_phi = beta_phi[0];
        double logit_gamma = beta_gamma[0];
        for (int j=1; j<(ncovs+1); ++j) {
          logit_phi += beta_phi[j] * x(t*nsite+i, j-1);
          logit_gamma += beta_gamma[j] * x(t*nsite+i, j-1);
        } // j
        double phi = R::plogis(logit_phi, 0, 1, true, false);
        double gamma = R::plogis(logit_gamma, 0, 1, true, false);
        double s_temp = z(i,t-1) * phi + (1 - z(i,t-1)) * gamma; 
        double num_temp = s_temp * qprod;
        double psi_temp = num_temp / (num_temp + (1 - s_temp));
        if (ysum(t*nsite+i) == 0) {
          z(i,t) = R::rbinom(1, psi_temp);
        }
      } // t
    } // i

    //Save samples
    beta0_save(_,k) = beta0;
    beta_phi_save(_,k) = beta_phi;
    beta_gamma_save(_,k) = beta_gamma;
    alpha_save(_,k) = alpha;

  } //k

  // Write output
  List out = List::create(Named("beta0_save") = beta0_save, 
                          Named("beta_phi_save") = beta_phi_save, 
                          Named("beta_gamma_save") = beta_gamma_save, 
                          Named("alpha_save") = alpha_save);

  return(out);
}

')

#==========
# Run MCMC
#==========
nmcmc <- 50000 # number of iterations
out <- dynam_occupancy_mcmc(y, x, w, nmcmc, nsite, nyear, nreps, ncovs, npcvs)

end_time <- Sys.time() # end time of computing
time_taken <- end_time - start_time # computing time
time_taken

ylab_beta0 <- c(expression(beta[0]^"[0]"), expression(beta[1]^"[0]"), expression(beta[2]^"[0]"))
ylab_beta_phi <- c(expression(beta[0]^""["["*phi*"]"]), expression(beta[1]^""["["*phi*"]"]), expression(beta[2]^""["["*phi*"]"]))
ylab_beta_gamma <- c(expression(beta[0]^""["["*gamma*"]"]), expression(beta[1]^""["["*gamma*"]"]), expression(beta[2]^""["["*gamma*"]"]))
ylab_alpha <- c(expression(alpha[0]), expression(alpha[1]), expression(alpha[2]))

par(mfrow=c(4,3))
par(mar=c(1,3,3,1))
par(oma=c(4,2.5,0,0))

for (i in 1:(ncovs+1)) {
  tt <- out$beta0_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta0[i] - yint * 8
  ymax <- beta0[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='blue')
  abline(h=beta0[i], col='grey16', lwd=1.5)
  title(main=ylab_beta0[i], cex.main=2, line=1.2)
} # i

for (i in 1:(ncovs+1)) {
  tt <- out$beta_phi_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_phi[i] - yint * 8
  ymax <- beta_phi[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='blue')
  abline(h=beta_phi[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_phi[i], cex.main=2, line=1.2)
} # i

for (i in 1:(ncovs+1)) {
  tt <- out$beta_gamma_save[i,]
  yint <- ceiling(sd(tt[(nmcmc*0.2+1):nmcmc]) * 10) / 10
  yint <- ifelse (yint < 0.1, 0.1, yint)
  ymin <- beta_gamma[i] - yint * 8
  ymax <- beta_gamma[i] + yint * 8
  plot(1, xlim=c(0,nmcmc), ylim=c(ymin, ymax), type='n', xlab='', ylab='', axes=F)
  rect(0, -100, nmcmc*0.2, 100, col = 'grey88', border = NA)
  axis(1, at=seq(0,nmcmc,length.out=6), labels=rep('',6))
  axis(2, at=round(seq(ymin, ymax, yint*4), digits=1), las=2)
  box()
  lines(tt, col='blue')
  abline(h=beta_gamma[i], col='grey16', lwd=1.5)
  title(main=ylab_beta_gamma[i], cex.main=2, line=1.2)
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
  lines(tt, col='blue')
  abline(h=alpha[i], col='grey16', lwd=1.5)
  title(main=ylab_alpha[i], cex.main=2, line=1.2)
} # i

title(xlab='Iteration', cex.lab=2.5, line=2.6, outer=T)
title(ylab='Posterior Value', cex.lab=2.5, line=0.4, outer=T)


