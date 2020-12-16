# Github Repository - Homework 4

### Functions and Libraries
{
  ### Libraries
  {
    library(dplyr)
    library(RColorBrewer)
    library(mosaic)
    library(reshape2)
    library(readr)
    library(ggplot2)
    library(mgcv)
    library(kableExtra)
    library(tmvtnorm)
    library(AER)
    library(car)
    library(tidyverse)
    library(statmod)
    library(rstanarm)
    library(tidyr)
    library(cowplot)
    library(mvtnorm)
    library(truncnorm)
    library(cowplot)
  }
  
  
  ### GLM 209 Functions
  {
    inv.logit = function(x) {
      exp(x) / (1 +exp(x))
    }
    
    pplc = function(yinew, yi, k) {
      var.sum = sum(apply(yinew, 2, var))
      sse = (yi-apply(yinew, 2, mean))^2
      var.sum + k/(k+1)*sum(sse)
    }
    
    predict.pi = function(dose_seq, model){
      
      pi.hat = matrix(NA,ncol=1,nrow=length(dose_seq))
      
      for(i in 1:length(dose_seq)){
        lin.pred = model$coefficients[1]+model$coefficients[2]*dose_seq[i]
        lin.pred = inv.logit(lin.pred)
        pi.hat[i,1] = lin.pred
      }
      out=list(pi.hat=pi.hat)
    }
    
    predict.pi.alpha = function(dose_seq, beta1, beta2, alpha){
      
      pi.hat = as.data.frame(matrix(NA,ncol=length(alpha),
                                    nrow=length(dose_seq)))
      
      for(j in 1:length(alpha)) {
        
        for(i in 1:length(dose_seq)){
          lin.pred = beta1 + beta2 * dose_seq[i]
          lin.pred = exp(alpha[j]*lin.pred) / (1 + exp(lin.pred)) ^ (alpha[j])
          pi.hat[i,j] = lin.pred
        }
      }
      pi.hat = pi.hat %>%
        mutate(Dose = dose_seq)
      out=list(pi.hat=pi.hat)
    }
    
    predict.pi.no.model = function(dose_seq, beta1, beta2){
      
      pi.hat = matrix(NA,ncol=1,nrow=length(dose_seq))
      
      for(i in 1:length(dose_seq)){
        lin.pred = beta1 + beta2 * dose_seq[i]
        lin.pred = inv.logit(lin.pred)
        pi.hat[i,1] = lin.pred
      }
      out=list(pi.hat=pi.hat)
    }
    
    deviance.resid = function(observed, fitted, size) {
      
      ith_dev = 2 * (observed * log(observed / (fitted * size)) +
                       (size - observed) * log((size - observed) / 
                                                 (size - (fitted*size))))
      
      dev = sign(observed-fitted*mi) * sqrt(ith_dev)
      
      for(i in 1:length(dev)){
        if(is.nan(dev[i])) {
          dev[i] = 0
        }
      }
      
      return(dev)
    }
    
    plot.resids.fitted = function(residuals, fitted, model, link){
      data = data.frame("residuals" = residuals, "fitted" = fitted)
      ggplot(data, aes(x=residuals, y=fitted)) +
        geom_point() +
        geom_hline(yintercept = 0, color="red") +
        theme_bw(base_size = 16) + 
        ggtitle(label = model, subtitle = link)
    }
    
    conf.int = function(mean, var, zscore) {
      lower = mean - zscore * sqrt(var)
      upper = mean + zscore * sqrt(var)
      
      out=list(low = lower, up = upper)
      return(c(lower,upper))
    }
    
    variance_funct = function(b1,b2,data,x0,dispersion=NULL){
      
      if(is.null(dispersion)){
        j.11 = sum(exp(b1 + b2 * data))
        j.12 <- j.21 <- sum(data * exp(b1 + b2 * data))
        j.22 = sum(data ^ 2 * exp(b1 + b2 * data))
        J_mat = matrix(c(j.11, j.12, j.21, j.22),
                       nrow = 2,
                       ncol = 2,
                       byrow = TRUE)
        var_hat = solve(J_mat)
        var_b1 = var_hat[1, 1]
        var_b2 = var_hat[2, 2]
        cov_hat = 2 * x0 * var_hat[1, 2]
        var_hat = var_b1 + (x0 ^ 2) * var_b2 + cov_hat
      }
      
      else if(!is.null(dispersion)){
        j.11 = (1/dispersion) * sum(exp(b1 + b2 * data))
        j.12 <- j.21 <- (1/dispersion) * sum(data * exp(b1 + b2 * data))
        j.22 = (1/dispersion) * sum(data ^ 2 * exp(b1 + b2 * data))
        J_mat = matrix(c(j.11, j.12, j.21, j.22),
                       nrow = 2,
                       ncol = 2,
                       byrow = TRUE)
        var_hat = solve(J_mat)
        var_b1 = var_hat[1, 1]
        var_b2 = var_hat[2, 2]
        cov_hat = 2 * x0 * var_hat[1, 2]
        var_hat = var_b1 + (x0 ^ 2) * var_b2 + cov_hat
      }
      
      return(var_hat)
    }
  }
  
  # Bayesian Stuff
  
  ### Binomial MCMC
  {
    binom.glm.MCMC = function(initial_beta, covar.mat, var.tuning, niter, log.post) {
      accept = 0
      beta = initial_beta
      beta.store <- matrix(NA, niter, length(beta))
      
      for(i in 1:niter){
        prop <- mvtnorm::rmvnorm(1, beta, var.tuning*covar.mat)
        ratio <- log.post(prop) - log.post(beta)
        if(log(runif(1)) < ratio){
          beta = prop
          accept = accept + 1
        }
        else{
          beta = beta;
          accep = accept
        }
        beta.store[i,] = beta
      }
      out=list(beta.store=beta.store, accept=accept)
    }
    
  }
  
  ### Gamma Random Component
  {
    gamma.llh.invlink = function(beta) {
      
      b1=beta[1]; b2=beta[2]; b3=beta[3]
      
      prior.mu = c(b1.gam.mle, b2.gam.mle)
      prior.Sigma = 100*diag(2)
      prior.shape = 0.0001
      prior.rate = 0.0001
      link = b1+b2*xi
      if (any((link)<0)) return(-Inf)
      
      gam.prior = 1/dgamma(exp(b3), prior.shape, prior.rate, log=TRUE)
      norm.prior = mvtnorm::dmvnorm(c(b1, b2), prior.mu, prior.Sigma, log=TRUE)
      likelihood = sum(-lgamma(exp(b3)) + exp(b3)*b3 + exp(b3)*log(link) + 
                         (exp(b3)-1)*log(yi)-yi*exp(b3)*(link))
      jacobian = b3
      density = gam.prior + norm.prior + likelihood + jacobian
      return(density)
    }
    
  }
  
  #### Inverse Gaussian Random Component
  {
    loglike.ig = function(beta){
      
      b1=beta[1];b2=beta[2];b3=beta[3]
      prior.mu = c(b1.ig.mle, b2.ig.mle)
      prior.Sigma = 100*diag(2)
      prior.shape = 0.0001
      prior.rate = 0.0001
      phi = exp(b3)
      link = b1+b2*xi
      if (any(link<0)) return(-Inf)
      gam.prior = dgamma(exp(b3), prior.shape, prior.rate, log=TRUE)
      norm.prior = mvtnorm::dmvnorm(c(b1, b2), prior.mu, prior.Sigma, log=TRUE)
      likelihood = sum(-1/2*log(2*pi*phi*yi^3) - yi/(2*phi)*(link)^2 + 
                         1/phi*(link) - 1/(2*phi*yi))
      jacobian = b3
      density = gam.prior + norm.prior + likelihood + jacobian
      return(density)
    }
    
  }
  
  ### Normal RW MCMC
  
  my.mcmc = function(initial_beta, log.post, cov.mat, niter, var.tune){
    
    current = initial_beta
    beta.store = matrix(NA, niter, length(current))
    accept=0 
    
    for(i in 1:niter){
      
      proposed = mvtnorm::rmvnorm(1, current, cov.mat*var.tune)
      #print(proposed); #browser()
      ratio = log.post(proposed) - log.post(current)
      if(log(runif(1)) < ratio) { 
        current = proposed
        accept = accept+1
      }
      else{
        current = current
        accept = accept
      } 
      beta.store[i,] = current
    } 
    #print(accept/niter)
    out=list(beta.store = beta.store, accept = accept)
  }
  
  ### Three Parameter MCMC (positive truncation)
  mcmc.3param = function(current, covar.mat, var.tuning, niter, log.post) { 
    accept = 0
    beta = c(current[1], current[2])
    alpha = current[3]
    
    param_store = matrix(nrow = niter, ncol = 3)
    
    for (i in 1:niter){
      prop = rtmvnorm(1, c(beta,alpha), covar.mat*var.tuning,
                      lower=c(-Inf,-Inf,0))
      ratio = log.post(prop) - log.post(c(beta,alpha))
      
      if(log(runif(1)) < ratio){
        beta = prop[1:2]
        alpha = prop[3]
        accept = accept + 1
      }
      param_store[i,]=c(beta,alpha)
    }
    
    out=list(beta.store=param_store, accept=accept)
  }
  
  mh.within.gibbs = function(lambda.condit,
                             beta.condit,
                             lambda.prior,
                             niter,
                             mu.init,
                             beta.init,
                             lambda.init,
                             lambda.var,
                             beta.cov.mat) {
    
    beta.store = matrix(NA, nrow = (niter + 1), ncol = 2)
    mu.store = matrix(NA, nrow = (niter + 1), ncol = 32)
    lambda.store = matrix(NA, nrow = (niter + 1), ncol = 1)
    
    mu.store[1,] = mu.init; lambda.store[1,] = lambda.init; beta.store[1,] = beta.init
    
    accept.lambda = 0
    accept.beta = 0 
    
    for(i in 2:(niter + 1)){
      
      mu.current = mu.store[(i - 1), ]
      lambda.current = lambda.store[(i - 1), ]
      beta.current = beta.store[(i - 1), ]
      
      ### sample mui
      b1 = beta.current[1]; b2 = beta.current[2]
      gamma = exp(b1 + b2*xi)
      mu = rgamma(32, yi + lambda.current, 1 + lambda.current/gamma)
      
      ### sample lambda
      lambda.prop = rnorm(1, log(lambda.current), lambda.var)
      lambda.prop = exp(lambda.prop)
      log.ratio <- lambda.condit(lambda=lambda.prop, mu=mu.current, beta=beta.current) + lambda.prior(lambda.prop) - 
        lambda.condit(lambda=lambda.current, mu=mu.current, beta=beta.current) - lambda.prior(lambda.current)
      
      if(log(runif(1)) < log.ratio) {
        lambda = lambda.prop
        accept.lambda =  accept.lambda + 1
      }
      else {
        lambda = lambda.current
      }
      
      ### sample beta1, beta2
      beta.prop = mvrnorm(1, beta.current, beta.cov.mat)
      log.ratio.beta = beta.condit(beta = beta.prop, mu = mu.current, lambda = lambda.current) - 
        beta.condit(beta = beta.current, mu = mu.current, lambda = lambda.current)
      
      if(log(runif(1)) < log.ratio.beta) {
        beta = beta.prop
        accept.beta = accept.beta + 1 
      }
      else beta = beta.current
      
      ### Store Updated Parameters
      mu.store[i,] = mu
      lambda.store[i,] = lambda
      beta.store[i,] = beta
    }
    out = list(mu.store = mu.store, lambda.store = lambda.store, beta.store = beta.store,
               accept.lambda = accept.lambda, accept.beta = accept.beta)
  }
  
}

### Question 2
{
  alligator = data.frame(
    length = c(1.30, 1.32, 1.32, 1.40, 1.42, 1.42, 1.47, 1.47, 1.50, 1.52, 1.63, 1.65, 1.65, 
               1.65, 1.65, 1.68, 1.70, 1.73, 1.78, 1.78, 1.80, 1.85, 1.93, 1.93, 1.98, 2.03, 
               2.03, 2.31, 2.36, 2.46, 3.25, 3.28, 3.33, 3.56, 3.58, 3.66, 3.68, 3.71, 3.89, 
               1.24, 1.30, 1.45, 1.45, 1.55, 1.60, 1.60, 1.65, 1.78, 1.78, 1.80, 1.88, 2.16, 
               2.26, 2.31, 2.36, 2.39, 2.41, 2.44, 2.56, 2.67, 2.72, 2.79, 2.84),
    sex = c(rep("M", 39), rep("F", 24)),
    choice = c("I", "F", "F", "F", "I", "F", "I", "F", "I", "I", "I", "O", "O", "I",
               "F", "F","I","O","F","O","F","F","I","F","I","F","F", "F","F","F","O","O",
               "F","F","F","F","O","F","F","I", "I", "I", "O", "I", "I", "I", "F", "I",
               "O", "I", "I", "F", "F", "F", "F", "F", "F", "F", "O", "F","I", "F", "F")
  )
  
  x = alligator$length
  
  ggplot(data=alligator) + geom_boxplot(aes(x=choice,y=length, color=choice)) +
    xlab("Food Choice") + ylab("Length") + 
    scale_color_discrete(name="Food Choice", labels=c("Fish", "Invertebrate", "Other")) + 
    theme_bw()
  
  vglm.fit = vglm(choice ~ length, data=alligator, family=multinomial())
  
  y.other = as.numeric(alligator$choice == "O");
  y.fish = as.numeric(alligator$choice == "F");
  y.inv = as.numeric(alligator$choice == "I")
  
  multinom.llh = function(theta) {
    
    alpha1 = theta[1];alpha2 = theta[2]; beta1 = theta[3]; beta2 = theta[4]
    alpha = c(alpha1, alpha2); beta = c(beta1, beta2)
    
    sum_probs = exp(alpha1 + beta1*x) + exp(alpha2+beta2*x)
    pi1 = exp(alpha1 + beta1*x)/(1+sum_probs)
    pi2 = exp(alpha2 + beta2*x)/(1+sum_probs)
    prior.mu = theta; prior.Sigma = 100*diag(4)
    
    llh = sum(y.other*log(pi1) + y.inv*log(pi2) + y.fish*log(1-pi1-pi2), na.rm=TRUE)
    beta.prior = mvtnorm::dmvnorm(c(theta), prior.mu, prior.Sigma, log=TRUE)
    multinom.llh = llh + beta.prior
    return(multinom.llh)
  }
  
  init = c(coef(vglm.fit)[1], coef(vglm.fit)[2], coef(vglm.fit)[3], coef(vglm.fit)[4])
  
  cov.mat = optim(fn = multinom.llh, par = init,hessian = TRUE,control=list(fnscale=-1))
  fisher.info = solve(-cov.mat$hessian)
  
  iter = 20000
  
  test.mcmc = my.mcmc(initial_beta = c(1,1,1,1), log.post = multinom.llh, cov.mat = fisher.info,
                      niter=iter, var.tune=1)
  
  theta.store = data.frame(alpha1 = test.mcmc$beta.store[-c(1:1500),1], 
                           alpha2 = test.mcmc$beta.store[-c(1:1500),2],
                           beta1 = test.mcmc$beta.store[-c(1:1500),3], 
                           beta2 = test.mcmc$beta.store[-c(1:1500),4])
  
  
  lag1.alpha1 = acf(theta.store$alpha1,plot = FALSE)[1]; lag1.alpha2 = acf(theta.store$alpha2,plot = FALSE)[1]
  lag1.beta1 = acf(theta.store$beta1,plot = FALSE)[1]; lag1.beta2 = acf(theta.store$beta2,plot = FALSE)[1]
  eff.size = effectiveSize(test.mcmc$beta.store)
  
  par(mfrow=c(2,2))
  acf(theta.store$alpha1, main = expression("ACF for "*alpha[1]));  
  acf(theta.store$alpha2, main = expression("ACF for "*alpha[2]));
  acf(theta.store$beta1, main = expression("ACF for "*beta[1]));
  acf(theta.store$beta2, main = expression("ACF for "*beta[2]));
  
  plot(theta.store$alpha1, type="l", ylab=expression(alpha[1]))
  abline(h=median(test.mcmc$beta.store[-c(1:1500),1]),col="red")
  plot(theta.store$alpha2, type="l", ylab=expression(alpha[2]))
  abline(h=median(test.mcmc$beta.store[-c(1:1500),2]),col="red")
  plot(theta.store$beta1, type="l", ylab=expression(beta[1]))
  abline(h=median(test.mcmc$beta.store[-c(1:1500),3]),col="red")
  plot(theta.store$beta2, type="l", ylab=expression(beta[2]))
  abline(h=median(test.mcmc$beta.store[-c(1:1500),4]),col="red")
  
  xgrid = seq(min(x), max(x), len=100)
  pi1.preds.store = sapply(xgrid, function(xx) exp(theta.store$alpha1 + theta.store$beta1*xx) / 
                             (1 + exp(theta.store$alpha1 + theta.store$beta1*xx) + 
                                exp(theta.store$alpha2 + theta.store$beta2*xx)))
  
  pi2.preds.store = sapply(xgrid, function(xx) exp(theta.store$alpha2 + theta.store$beta2*xx) / 
                             (1 + exp(theta.store$alpha1 + theta.store$beta1*xx) + 
                                exp(theta.store$alpha2 + theta.store$beta2*xx)))
  
  pi3.preds.store = 1 - pi1.preds.store - pi2.preds.store
  
  log.odds1 = log(pi1.preds.store/pi3.preds.store)
  log.odds2 = log(pi2.preds.store/pi3.preds.store)
  
  mu.ests.other = apply(pi1.preds.store, 2, mean)
  mu.ests.invert = apply(pi2.preds.store, 2, mean)
  mu.ests.fish = apply(pi3.preds.store, 2, mean)
  lower.other = apply(pi1.preds.store, 2, quantile, probs=0.025); 
  up.other = apply(pi1.preds.store, 2,quantile, probs=0.975)
  lower.invert = apply(pi2.preds.store, 2, quantile, probs=0.025); 
  up.invert = apply(pi2.preds.store, 2,quantile, probs=0.975)
  lower.fish = apply(pi3.preds.store, 2, quantile, probs=0.025); 
  up.fish = apply(pi3.preds.store, 2,quantile, probs=0.975)
  
  ggplot() + xlab("length") + ylab("Probability of Choice") +
    geom_line(aes(x=xgrid, y=mu.ests.other, color="Other"), lwd=1) +
    geom_line(aes(x=xgrid, y=mu.ests.invert, color="Invertebrate"), lwd=1) + 
    geom_line(aes(x=xgrid, y=mu.ests.fish, color="Fish"), lwd=1) + 
    geom_ribbon(aes(x=xgrid, ymin=lower.other, ymax=up.other), alpha=0.1, fill="green") + 
    geom_ribbon(aes(x=xgrid, ymin=lower.invert, ymax=up.invert), alpha=0.2, fill="pink") + 
    geom_ribbon(aes(x=xgrid, ymin=lower.fish, ymax=up.fish), alpha=0.2, fill="blue") + 
    scale_color_manual(name = "Food", 
                       values = c("Other" = "springgreen2", "Invertebrate" = "deeppink2",
                                  "Fish" = "deepskyblue1")) +
    ggtitle("Food Choice as a Function of Length") + 
    theme_bw()
  
  # Point Estimates
  pi1.pt.store = sapply(x, function(xx) exp(theta.store$alpha1 + theta.store$beta1*xx) / (1 + exp(theta.store$alpha1 + theta.store$beta1*xx) + exp(theta.store$alpha2 + theta.store$beta2*xx)))
  pi2.pt.store = sapply(x, function(xx) exp(theta.store$alpha2 + theta.store$beta2*xx) / (1 + exp(theta.store$alpha1 + theta.store$beta1*xx) + exp(theta.store$alpha2 + theta.store$beta2*xx)))
  pi3.pt.store = 1 - pi1.pt.store - pi2.pt.store
  pi1.pt.low = apply(pi1.pt.store, 2, quantile, probs=0.025)
  pi1.pt.up = apply(pi1.pt.store, 2, quantile, probs=0.975)
  pi2.pt.low = apply(pi2.pt.store, 2, quantile, probs=0.025)
  pi2.pt.up = apply(pi2.pt.store, 2, quantile, probs=0.975)
  pi3.pt.low = apply(pi3.pt.store, 2, quantile, probs=0.025)
  pi3.pt.up = apply(pi3.pt.store, 2, quantile, probs=0.975)
  pi1.pt = apply(pi1.pt.store, 2, mean)
  pi2.pt = apply(pi2.pt.store, 2, mean)
  pi3.pt = 1 - pi1.pt - pi2.pt
  
  pt.est.table = round(cbind(pi1.pt.low, pi1.pt, pi1.pt.up, pi2.pt.low, 
                             pi2.pt, pi2.pt.up, pi3.pt.low, pi3.pt, pi3.pt.up),4)
  colnames = c(rep(c("Lower 2.5%", "Point Estimate", "Upper 97.5%"), 3))
  
  x.male = as.numeric(alligator$sex == "M");
  x.female = as.numeric(alligator$sex == "F");
  vglm.fit.m2 = vglm(choice ~ length + sex, data=alligator, family=multinomial())
  
  multinom.llh.m2 = function(theta) {
    
    alpha1 = theta[1];alpha2 = theta[2]; beta1 = theta[3]; beta2 = theta[4];
    beta3 = theta[5]; beta4 = theta[6]
    alpha = c(alpha1, alpha2); beta = c(beta1, beta2, beta3, beta4)
    
    sum_probs = exp(alpha1 + beta1*x + beta3*x.male) + exp(alpha2 + beta2*x + beta4*x.male)
    pi1 = exp(alpha1 + beta1*x+ beta3*x.male)/(1+sum_probs)
    pi2 = exp(alpha2 + beta2*x + beta4*x.male)/(1+sum_probs)
    prior.mu = theta; prior.Sigma = 100*diag(6)
    
    llh = sum(y.other*log(pi1) + y.inv*log(pi2) + y.fish*log(1-pi1-pi2), na.rm=TRUE)
    beta.prior = mvtnorm::dmvnorm(c(theta), prior.mu, prior.Sigma, log=TRUE)
    multinom.llh = llh + beta.prior
    return(multinom.llh)
  }
  
  init = c(coef(vglm.fit.m2)[1], coef(vglm.fit.m2)[2], coef(vglm.fit.m2)[3], coef(vglm.fit.m2)[4],
           coef(vglm.fit.m2)[5], coef(vglm.fit.m2)[6])
  
  cov.mat = optim(fn = multinom.llh.m2, par = init,hessian = TRUE,control=list(fnscale=-1))
  fisher.info = solve(-cov.mat$hessian)
  
  iter = 40000
  
  test.mcmc.m2 = my.mcmc(initial_beta = init, log.post = multinom.llh.m2, cov.mat = fisher.info,
                         niter=iter, var.tune=1)
  
  theta.store.m2 = data.frame(alpha1 = test.mcmc.m2$beta.store[-c(1:1500),1], 
                              alpha2 = test.mcmc.m2$beta.store[-c(1:1500),2],
                              beta1 = test.mcmc.m2$beta.store[-c(1:1500),3], 
                              beta2 = test.mcmc.m2$beta.store[-c(1:1500),4], 
                              beta3 = test.mcmc.m2$beta.store[-c(1:1500),5], 
                              beta4 = test.mcmc.m2$beta.store[-c(1:1500),6])
  
  par(mfrow=c(3,2))
  acf(theta.store.m2$alpha1, main = expression("ACF for "*alpha[1]));  
  acf(theta.store.m2$alpha2, main = expression("ACF for "*alpha[2]));
  acf(theta.store.m2$beta1, main = expression("ACF for "*beta[1]));
  acf(theta.store.m2$beta2, main = expression("ACF for "*beta[2]));
  acf(theta.store.m2$beta3, main = expression("ACF for "*zeta[1]));
  acf(theta.store.m2$beta4, main = expression("ACF for "*zeta[2]));
  # acf(theta.store.m2$alpha1, plot=FALSE)[1];acf(theta.store.m2$alpha2, plot=FALSE)[1] 
  # acf(theta.store.m2$beta1, plot=FALSE)[1];acf(theta.store.m2$beta2, plot=FALSE)[1] 
  # acf(theta.store.m2$beta3, plot=FALSE)[1];acf(theta.store.m2$beta4, plot=FALSE)[1] 
  # eff.size = effectiveSize(test.mcmc.m2$beta.store)
  # print(eff.size)
  
  plot(theta.store.m2$alpha1, type="l", ylab=expression(alpha[1]))
  abline(h=median(test.mcmc.m2$beta.store[-c(1:1500),1]),col="red")
  plot(theta.store.m2$alpha2, type="l", ylab=expression(alpha[2]))
  abline(h=median(test.mcmc.m2$beta.store[-c(1:1500),2]),col="red")
  plot(theta.store.m2$beta1, type="l", ylab=expression(beta[1]))
  abline(h=median(test.mcmc.m2$beta.store[-c(1:1500),3]),col="red")
  plot(theta.store.m2$beta2, type="l", ylab=expression(beta[2]))
  abline(h=median(test.mcmc.m2$beta.store[-c(1:1500),4]),col="red")
  plot(theta.store.m2$beta3, type="l", ylab=expression(zeta[1]))
  abline(h=median(test.mcmc.m2$beta.store[-c(1:1500),5]),col="red")
  plot(theta.store.m2$beta4, type="l", ylab=expression(zeta[2]))
  abline(h=median(test.mcmc.m2$beta.store[-c(1:1500),6]),col="red")
  
  xgrid = seq(min(x), max(x), len=100)
  
  # Males 
  {
    pi1.preds.male = sapply(xgrid, function(xx) exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx + theta.store.m2$beta3) / 
                              (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx + theta.store.m2$beta3) + 
                                 exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx + theta.store.m2$beta4)))
    
    pi2.preds.male = sapply(xgrid, function(xx) exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx + theta.store.m2$beta4) / 
                              (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx + theta.store.m2$beta3) + 
                                 exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx + theta.store.m2$beta4)))
    
    pi3.preds.male = 1 - pi1.preds.male - pi2.preds.male
    
    mu.male.ests.other = apply(pi1.preds.male, 2, mean)
    mu.male.ests.invert = apply(pi2.preds.male, 2, mean)
    mu.male.ests.fish = apply(pi3.preds.male, 2, mean)
    
    lower.other.male = apply(pi1.preds.male, 2, quantile, probs=0.025); 
    up.other.male = apply(pi1.preds.male, 2,quantile, probs=0.975)
    lower.invert.male = apply(pi2.preds.male, 2, quantile, probs=0.025); 
    up.invert.male = apply(pi2.preds.male, 2,quantile, probs=0.975)
    lower.fish.male = apply(pi3.preds.male, 2, quantile, probs=0.025); 
    up.fish.male = apply(pi3.preds.male, 2,quantile, probs=0.975)
  }
  
  # Females
  {
    pi1.preds.female = sapply(xgrid, function(xx) exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx ) / 
                                (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx) + 
                                   exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx)))
    
    pi2.preds.female = sapply(xgrid, function(xx) exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx) / 
                                (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx) + 
                                   exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx)))
    
    pi3.preds.female = 1 - pi1.preds.female - pi2.preds.female
    
    mu.female.ests.other = apply(pi1.preds.female, 2, mean)
    mu.female.ests.invert = apply(pi2.preds.female, 2, mean)
    mu.female.ests.fish = apply(pi3.preds.female, 2, mean)
    
    lower.other.female = apply(pi1.preds.female, 2, quantile, probs=0.025); 
    up.other.female = apply(pi1.preds.female, 2,quantile, probs=0.975)
    lower.invert.female = apply(pi2.preds.female, 2, quantile, probs=0.025); 
    up.invert.female = apply(pi2.preds.female, 2,quantile, probs=0.975)
    lower.fish.female = apply(pi3.preds.female, 2, quantile, probs=0.025); 
    up.fish.female = apply(pi3.preds.female, 2,quantile, probs=0.975)
  }
  
  male_prbs = ggplot() + xlab("Length") + ylab("Probability of Choice") +
    geom_line(aes(x=xgrid, y=mu.male.ests.other, color="Other"), lwd=1) +
    geom_line(aes(x=xgrid, y=mu.male.ests.invert, color="Invertebrate"), lwd=1) + 
    geom_line(aes(x=xgrid, y=mu.male.ests.fish, color="Fish"), lwd=1) + 
    geom_ribbon(aes(x=xgrid, ymin=lower.other.male, ymax=up.other.male), alpha=0.1, fill="green") + 
    geom_ribbon(aes(x=xgrid, ymin=lower.invert.male, ymax=up.invert.male), alpha=0.2, fill="pink") + 
    geom_ribbon(aes(x=xgrid, ymin=lower.fish.male, ymax=up.fish.male), alpha=0.2, fill="blue") + 
    scale_color_manual(name = "Food", 
                       values = c("Other" = "springgreen2", "Invertebrate" = "deeppink2",
                                  "Fish" = "deepskyblue1")) +
    ggtitle("Alligator Food Choice", "Males") + 
    theme_bw() + 
    theme(legend.position="bottom")
  
  female_prbs = ggplot() + xlab("Length") + ylab("Probability of Choice") +
    geom_line(aes(x=xgrid, y=mu.female.ests.other, color="Other"), lwd=1) +
    geom_line(aes(x=xgrid, y=mu.female.ests.invert, color="Invertebrate"), lwd=1) + 
    geom_line(aes(x=xgrid, y=mu.female.ests.fish, color="Fish"), lwd=1) + 
    geom_ribbon(aes(x=xgrid, ymin=lower.other.female, ymax=up.other.female), alpha=0.1, fill="green") + 
    geom_ribbon(aes(x=xgrid, ymin=lower.invert.female, ymax=up.invert.female), alpha=0.2, fill="pink") + 
    geom_ribbon(aes(x=xgrid, ymin=lower.fish.female, ymax=up.fish.female), alpha=0.2, fill="blue") + 
    scale_color_manual(name = "Food", 
                       values = c("Other" = "springgreen2", "Invertebrate" = "deeppink2",
                                  "Fish" = "deepskyblue1")) +
    ggtitle(" ", "Females") +
    theme_bw() +
    theme(legend.position="bottom")
  
  
  plot_grid(male_prbs, female_prbs)
  
  # Point Estimates
  pi1.pt.store.male = sapply(x, function(xx) exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx + theta.store.m2$beta3) / 
                               (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx + theta.store.m2$beta3) + 
                                  exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx + theta.store.m2$beta4)))
  pi2.pt.store.male = sapply(x, function(xx) exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx + theta.store.m2$beta4) / 
                               (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx + theta.store.m2$beta3) + 
                                  exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx + theta.store.m2$beta4)))
  pi3.pt.store.male = 1 - pi1.pt.store - pi2.pt.store
  
  pi1.pt.low.male = apply(pi1.pt.store.male, 2, quantile, probs=0.025)
  pi1.pt.up.male = apply(pi1.pt.store.male, 2, quantile, probs=0.975)
  pi2.pt.low.male = apply(pi2.pt.store.male, 2, quantile, probs=0.025)
  pi2.pt.up.male = apply(pi2.pt.store.male, 2, quantile, probs=0.975)
  pi3.pt.low.male = apply(pi3.pt.store.male, 2, quantile, probs=0.025)
  pi3.pt.up.male = apply(pi3.pt.store.male, 2, quantile, probs=0.975)
  
  pi1.pt.male = apply(pi1.pt.store.male, 2, mean)
  pi2.pt.male = apply(pi2.pt.store.male, 2, mean)
  pi3.pt.male = 1 - pi1.pt.male - pi2.pt.male
  
  pt.est.table.male = round(cbind(pi1.pt.low.male, pi1.pt.male, pi1.pt.up.male, 
                                  pi2.pt.low.male, pi2.pt.male, pi2.pt.up.male, 
                                  pi3.pt.low.male, pi3.pt.male, pi3.pt.up.male),4)
  colnames = c(rep(c("Lower 2.5%", "Point Estimate", "Upper 97.5%"), 3))
  pt.est.table.male %>% kable(col.names = colnames, caption = "Point Estimates and Confidence Intervals - Males",
                              format = "latex") %>%
    add_header_above(c( "pi 1"=3, "pi 2"=3, "pi 3"=3))
  
  # Point Estimates
  pi1.pt.store.female = sapply(x, function(xx) exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx) / 
                                 (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx) + 
                                    exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx)))
  pi2.pt.store.female = sapply(x, function(xx) exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx) / 
                                 (1 + exp(theta.store.m2$alpha1 + theta.store.m2$beta1*xx) + 
                                    exp(theta.store.m2$alpha2 + theta.store.m2$beta2*xx)))
  pi3.pt.store.female = 1 - pi1.pt.store - pi2.pt.store
  
  pi1.pt.low.female = apply(pi1.pt.store.female, 2, quantile, probs=0.025)
  pi1.pt.up.female = apply(pi1.pt.store.female, 2, quantile, probs=0.975)
  pi2.pt.low.female = apply(pi2.pt.store.female, 2, quantile, probs=0.025)
  pi2.pt.up.female = apply(pi2.pt.store.female, 2, quantile, probs=0.975)
  pi3.pt.low.female = apply(pi3.pt.store.female, 2, quantile, probs=0.025)
  pi3.pt.up.female = apply(pi3.pt.store.female, 2, quantile, probs=0.975)
  
  pi1.pt.female = apply(pi1.pt.store.female, 2, mean)
  pi2.pt.female = apply(pi2.pt.store.female, 2, mean)
  pi3.pt.female = 1 - pi1.pt.female - pi2.pt.female
  
  pt.est.table.female = round(cbind(pi1.pt.low.female, pi1.pt.female, pi1.pt.up.female, 
                                    pi2.pt.low.female, pi2.pt.female, pi2.pt.up.female, 
                                    pi3.pt.low.female, pi3.pt.female, pi3.pt.up.female),4)
  colnames = c(rep(c("Lower 2.5%", "Point Estimate", "Upper 97.5%"), 3))
  pt.est.table.female %>% kable(col.names = colnames, caption = "Point Estimates and Confidence Intervals - females",
                                format = "latex") %>%
    add_header_above(c( "pi 1"=3, "pi 2"=3, "pi 3"=3))
}

### Question 3
{
  toxicity = data.frame(
    concentration = c(0, 62.5, 125, 250, 500),
    dead = c(15,17,22,38,144), malform = c(1, 0, 7, 59, 132),
    normal = c(281, 225, 283, 202, 9), total = c(297, 242, 312, 299, 285)
  )
  yi1 = toxicity$dead; yi2 = toxicity$malform; yi3 = toxicity$normal
  mi = toxicity$total
  xi = toxicity$concentration
  
  toxicity = toxicity %>%
    mutate(pi1 = dead / total) %>%
    mutate(pi2 = malform / total) %>%
    mutate(pi3 = normal / total)
  
  ggplot(toxicity) + ylab("Probability of Outcome") + 
    geom_line( aes(x=concentration, y=pi1, color="pi1")) +
    geom_line( aes(x=concentration, y=pi2, color="pi2")) +
    geom_line( aes(x=concentration, y=pi3, color="pi3")) + 
    scale_color_manual(name = "Outcome", 
                       values = c("pi1" = "orange", "pi2" = "skyblue",
                                  "pi3" = "forestgreen"),
                       label = c("Dead", "Malformation", "Normal")) +
    ggtitle("Observed Probabilities of Outcomes", "Function of Concentration Dose") + 
    theme_bw()
  
  fit.pi1 = glm(cbind(yi1, mi-yi1) ~ concentration, data = toxicity,
                family = binomial(link = "logit"))
  fit.pi2 = glm(cbind(yi2, (mi - yi1) - yi2) ~ concentration, data = toxicity,
                family = binomial(link = "logit"))
  
  alpha1.mle = coef(fit.pi1)[1]; alpha2.mle = coef(fit.pi2)[1]
  beta1.mle = coef(fit.pi1)[2]; beta2.mle = coef(fit.pi2)[2]
  
  xgrid = seq(min(xi), max(xi), len=100)
  t.val = qt(p = 0.95,df = 4)
  pi_1_ests = predict(fit.pi1, newdata = data.frame("concentration"=xgrid), 
                      type="response")
  pi_2_ests = predict(fit.pi2, newdata = data.frame("concentration"=xgrid), 
                      type="response")
  pi_2_ests = (1 - pi_1_ests) * pi_2_ests
  pi_3_ests = 1 - pi_2_ests - pi_1_ests
  
  shapes = c("s1" = 16, "s2" = 3, "s3" = 2)
  ggplot() + ylab("Probability of Choice") +
    geom_point(data=toxicity, aes(x=concentration, y=pi1, shape="s1", color="pi1ests")) +
    geom_point(data=toxicity, aes(x=concentration, y=pi2, shape = "s2", color="pi2ests")) +
    geom_point(data=toxicity, aes(x=concentration, y=pi3, shape = "s3", color="pi3ests")) +
    geom_line(aes(x=xgrid, y=pi_1_ests, color="pi1ests")) + 
    geom_line(aes(x=xgrid, y=pi_2_ests, color="pi2ests")) + 
    geom_line(aes(x=xgrid, y=pi_3_ests, color="pi3ests")) + 
    scale_color_manual(name = "Predicted Fit", 
                       values = c("pi1ests" = "#80B1D3", "pi2ests" = "#FB8072",
                                  "pi3ests" = "#B3DE69"),
                       label = c("Dead", "Malformation", "Normal")) +
    scale_shape_manual(name="Observed", breaks = c("s1", "s2", "s3"), values = shapes,
                       labels=c("Dead", "Malformation", "Normal")) + 
    ggtitle("Classical GLM Fit - Observed Probabilities of Outcomes", "Function of Concentration Dose") + 
    guides(shape = guide_legend(override.aes = list(colour = c("#80B1D3","#FB8072", "#B3DE69")))) +
    theme_bw()
  pi1.fit = fit.pi1$fitted.values; pi2.fit = fit.pi2$fitted.values * (1-pi1.fit)
  pi3.fit = 1 - pi1.fit - pi2.fit
  
  multi.binom.llh = function(theta){
    
    a1 = theta[1]; a2 = theta[2]; b1 = theta[3]; b2 = theta[4]
    
    pi1 = logitlink(a1+b1*xi, inverse = TRUE)
    pi2 = (1 - pi1) * logitlink(a2 + b2*xi, inverse = TRUE)
    
    rho1 = pi1; rho2 = pi2 / (1 - pi1)
    
    prior.mu = theta; prior.Sigma = 100*diag(4)
    beta.prior = mvtnorm::dmvnorm(c(theta), prior.mu, prior.Sigma, log=TRUE)
    
    yi2givyi1.llh = sum(dbinom(yi2, mi-yi1, rho2, log = TRUE))
    yi1.llh = sum(dbinom(yi1, mi, rho1, log = TRUE))
    
    loglike = yi1.llh + yi2givyi1.llh + beta.prior
    return(loglike)
    
  }
  theta = c(fit.pi1$coefficients[1], fit.pi2$coefficients[1], 
            fit.pi1$coefficients[2], fit.pi2$coefficients[2])
  
  hess = optim(fn = multi.binom.llh, par = theta, 
               hessian = TRUE, control=list(fnscale=-1))
  fisher.info = solve(-hess$hessian)
  
  iter = 20000
  
  mcmc.binom.multi = my.mcmc(initial_beta = theta, log.post = multi.binom.llh, 
                             cov.mat = fisher.info, niter=iter, var.tune=1)
  
  theta.store = data.frame(alpha1 = mcmc.binom.multi$beta.store[-c(1:1500),1], 
                           alpha2 = mcmc.binom.multi$beta.store[-c(1:1500),2],
                           beta1 = mcmc.binom.multi$beta.store[-c(1:1500),3], 
                           beta2 = mcmc.binom.multi$beta.store[-c(1:1500),4])
  
  par(mfrow=c(2,2))
  acf(theta.store$alpha1, main = expression("ACF for "*alpha[1]));  
  acf(theta.store$alpha2, main = expression("ACF for "*alpha[2]));
  acf(theta.store$beta1, main = expression("ACF for "*beta[1]));
  acf(theta.store$beta2, main = expression("ACF for "*beta[2]));
  
  plot(theta.store$alpha1, type="l", ylab=expression(alpha[1]))
  abline(h=median(mcmc.binom.multi$beta.store[-c(1:1500),1]),col="red")
  plot(theta.store$alpha2, type="l", ylab=expression(alpha[2]))
  abline(h=median(mcmc.binom.multi$beta.store[-c(1:1500),2]),col="red")
  plot(theta.store$beta1, type="l", ylab=expression(beta[1]))
  abline(h=median(mcmc.binom.multi$beta.store[-c(1:1500),3]),col="red")
  plot(theta.store$beta2, type="l", ylab=expression(beta[2]))
  abline(h=median(mcmc.binom.multi$beta.store[-c(1:1500),4]),col="red")
  
  alpha1.bayes = theta.store$alpha1; alpha2.bayes = theta.store$alpha2
  beta1.bayes = theta.store$beta1; beta2.bayes = theta.store$beta2
  
  xgrid = seq(min(xi), max(xi), len=100)
  
  eta1 = sapply(xgrid, function(xx) alpha1.bayes + beta1.bayes*xx)
  eta2 = sapply(xgrid, function (xx) alpha2.bayes + beta2.bayes*xx)
  
  pi1.bayes.store = logitlink(eta1, inverse=TRUE)
  pi2.bayes.store = logitlink(eta2, inverse = TRUE) *(1 - pi1.bayes.store)
  pi3.bayes.store = 1 - pi1.bayes.store - pi2.bayes.store
  
  pi1.bayes = apply(pi1.bayes.store, 2, mean)
  pi1.bayes.low = apply(pi1.bayes.store, 2, quantile, probs=0.025)
  pi1.bayes.up = apply(pi1.bayes.store, 2, quantile, probs=0.975)
  pi2.bayes = apply(pi2.bayes.store, 2, mean)
  pi2.bayes.low = apply(pi2.bayes.store, 2, quantile, probs=0.025)
  pi2.bayes.up = apply(pi2.bayes.store, 2, quantile, probs=0.975)
  pi3.bayes = apply(pi3.bayes.store, 2, mean)
  pi3.bayes.low = apply(pi3.bayes.store, 2, quantile, probs=0.025)
  pi3.bayes.up = apply(pi3.bayes.store, 2, quantile, probs=0.975)
  
  shapes = c("s1" = 16, "s2" = 3, "s3" = 2)
  ggplot() + ylab("Probability of Choice") +
    geom_point(data=toxicity, aes(x=concentration, y=pi1, shape="s1", color="pi1ests")) +
    geom_point(data=toxicity, aes(x=concentration, y=pi2, shape = "s2", color="pi2ests")) +
    geom_point(data=toxicity, aes(x=concentration, y=pi3, shape = "s3", color="pi3ests")) +
    geom_line(aes(x=xgrid, y=pi1.bayes, color="pi1ests")) + 
    geom_line(aes(x=xgrid, y=pi2.bayes, color="pi2ests")) + 
    geom_line(aes(x=xgrid, y=pi3.bayes, color="pi3ests")) + 
    geom_ribbon(aes(x=xgrid, ymin=pi1.bayes.low, ymax=pi1.bayes.up),
                alpha=0.3,fill="#80B1D3") + 
    geom_ribbon(aes(x=xgrid, ymin=pi2.bayes.low, ymax=pi2.bayes.up),
                alpha=0.3,fill="#FB8072") + 
    geom_ribbon(aes(x=xgrid, ymin=pi3.bayes.low, ymax=pi3.bayes.up),
                alpha=0.3,fill="#B3DE69") + 
    guides(shape = guide_legend(override.aes = list(colour = c("#80B1D3","#FB8072", "#B3DE69")))) +
    scale_color_manual(name = "Bayesian Fit", 
                       values = c("pi1ests" = "#80B1D3", "pi2ests" = "#FB8072",
                                  "pi3ests" = "#B3DE69"),
                       label = c("Dead", "Malformation", "Normal")) +
    scale_shape_manual(name="Observed", breaks = c("s1", "s2", "s3"), values = shapes,
                       labels=c("Dead", "Malformation", "Normal")) + 
    ggtitle("Bayesian GLM - Observed Probabilities of Outcomes", "Function of Concentration Dose") + 
    theme_bw()
  
  ## Point Estimates:
  
  eta1.pt = sapply(xi, function(xx) alpha1.bayes + beta1.bayes*xx)
  eta2.pt = sapply(xi, function (xx) alpha2.bayes + beta2.bayes*xx)
  
  pi1.bayes.store.pt = logitlink(eta1.pt, inverse=TRUE)
  pi2.bayes.store.pt = logitlink(eta2.pt, inverse = TRUE) *(1 - pi1.bayes.store.pt)
  pi3.bayes.store.pt = 1 - pi1.bayes.store.pt - pi2.bayes.store.pt
  pi1.bayes.store.pt = apply(pi1.bayes.store.pt, 2, mean)
  pi2.bayes.store.pt = apply(pi2.bayes.store.pt, 2, mean)
  pi3.bayes.store.pt = apply(pi3.bayes.store.pt, 2, mean)
  
}

### Question 4
{
  fabfail=read.delim(file = "~/Documents/STAT 209 Generalized Linear Models/FabricFaults.txt", 
                     header = TRUE, sep = "\t")
  ggplot(fabfail) + 
    geom_point(aes(x=length, y=faults), size=2, color="blue") +
    xlab("Length") + ylab("Number of Faults") +
    theme_bw()
  
  yi = fabfail$faults; xi = fabfail$length
  
  glm.pois = glm(formula = faults ~ length, data = fabfail, family = poisson())
  s=summary(glm.pois)
  
  pois.llh = function(theta) {
    
    b1 = theta[1]; b2 = theta[2]
    link = exp(b1 + b2*xi)
    
    prior = log(1)
    llh = sum(dpois(x = yi, lambda = link, log = TRUE))
    loglike = llh + prior
    return(loglike)
  }
  
  iter = 50000
  init = c(coef(glm.pois)[1],coef(glm.pois)[2])
  s = summary(glm.pois)
  cov.mat = s$cov.unscaled
  
  pois.mcmc = my.mcmc(initial_beta = init, log.post = pois.llh,
                      cov.mat = cov.mat, niter = iter, var.tune = 2)
  theta.store = data.frame(beta1 = pois.mcmc$beta.store[-c(1:1500),1], 
                           beta2 = pois.mcmc$beta.store[-c(1:1500),2])
  
  par(mfrow=c(2,2))
  acf(theta.store$beta1, main = expression("ACF for "*beta[1]));
  acf(theta.store$beta2, main = expression("ACF for "*beta[2]));
  plot(theta.store$beta1, type="l", ylab=expression(beta[1]))
  abline(h=median(pois.mcmc$beta.store[-c(1:1500),1]),col="red")
  plot(theta.store$beta2, type="l", ylab=expression(beta[2]))
  abline(h=median(pois.mcmc$beta.store[-c(1:1500),2]),col="red")
  lag1.beta1 = acf(theta.store$beta1,plot = FALSE)[1]; lag1.beta2 = acf(theta.store$beta2,plot = FALSE)[1]
 
  b1.model1 = theta.store$beta1; b2.model1 = theta.store$beta2
  xgrid = seq(min(xi), max(xi), len=100)
  
  mu_i_store = sapply(xgrid, function(xx) exp(b1.model1+b2.model1*xx))
  
  mui = apply(mu_i_store, 2, mean)
  mui_low = apply(mu_i_store, 2, quantile, probs=0.025)
  mui_up = apply(mu_i_store, 2, quantile, probs=0.975)
  
  ggplot() + geom_point(data=fabfail, aes(x=length, y=faults)) + 
    geom_line(aes(x=xgrid, y=mui, color="mui_ests"), lty=2) +
    geom_ribbon(aes(x=xgrid, ymin = mui_low, ymax = mui_up), 
                fill = "red", alpha = 0.3) + 
    scale_color_manual(name = "Predicted Fit", 
                       values = c("mui_ests" = "red"),
                       label = c("Poission GLM")) +
    theme_bw() + xlab("Length") + ylab("Number of Faults") +
    ggtitle("Bayesian GLM Fit", "Flat Prior")
  
  posterior.yi = sapply(xi, function(x) rpois(length(b1.model1), 
                                              lambda = exp(b1.model1 + b2.model1*x)))
  post.yi.mean = apply(posterior.yi, 2, mean)
  
  boxplot(posterior.yi, pch=20, cex=0.1, xlab="Observation", ylab="Predicted Fault Number",
          main = "Posterior Predictive Samples")
  legend("topleft", legend="Observed", col="red", pch=16)
  points(x=(1:length(xi)), y=yi, col="red", pch=16)
  
  posterior.resids = data.frame(sweep(posterior.yi, 2 , yi))
  colnames(posterior.resids) = seq(1:32)
  boxplot(posterior.resids, pch=20, cex=0.1, xlab="Observation", ylab="Residual",
          main = "Posterior Predictive Residuals")
  abline(h=0, col="red", lwd=2)
  
  glm.fit.gamma = glm(faults~length,data=fabfail,family=Gamma(link="log"))
  
  pareto.lambda.prior = function(lambda){
    return(-2*log(lambda+1))
  } 
  gamma.lambda.prior = function(lambda){
    prior = dgamma(lambda, shape = 0.001, rate=0.001, log=TRUE)
    return(prior)
  } 
  weibull.lambda.prior = function(lambda){
    prior = dweibull(lambda, shape = 0.001, scale = 0.001, log=TRUE)
    return(prior)
  } 
  
  full.condit.lambda <- function(lambda, mu, beta){
    
    b1 = beta[1]; b2 = beta[2]
    gamma = exp(b1 + b2*xi)
    llh = sum(lambda*(log(lambda) - log(gamma) + log(mu))-
                lgamma(lambda) - (lambda*mu)/gamma)
    post = llh
    return(post)
  } 
  full.condit.beta <- function(beta,mu,lambda) {
    
    b1 = beta[1]; b2 = beta[2]
    gamma = exp(b1 + b2*xi)
    beta.prior = log(1)
    llh = sum(-lambda * log(gamma) - (lambda * mu)/gamma)
    post = llh + beta.prior
    return(post)
  }
  
  b.init = c(0,0); l.init = 1
  beta.covar = vcov(glm.fit.gamma)
  iter = 100000
  
  hier.mcmc.pareto =  mh.within.gibbs(lambda.condit = full.condit.lambda, beta.condit = full.condit.beta, 
                                      lambda.prior = pareto.lambda.prior, beta.init = b.init, lambda.init = l.init, 
                                      mu.init = fabfail$faults, niter = iter, lambda.var = 1, 
                                      beta.cov.mat = 2.5*beta.covar)
  
  hier.mcmc.gamma =  mh.within.gibbs(lambda.condit = full.condit.lambda, beta.condit = full.condit.beta, 
                                     lambda.prior = gamma.lambda.prior, beta.init = b.init, lambda.init = l.init, 
                                     mu.init = fabfail$faults, niter = iter, lambda.var = 1, 
                                     beta.cov.mat = 2.5*beta.covar)
  
  hier.mcmc.weibull =  mh.within.gibbs(lambda.condit = full.condit.lambda, beta.condit = full.condit.beta, 
                                       lambda.prior = weibull.lambda.prior, beta.init = b.init, lambda.init = l.init, 
                                       mu.init = fabfail$faults, niter = iter, lambda.var = 1, 
                                       beta.cov.mat = 2.5*beta.covar)
  
  burnin = 10000
  
  b1.hier = hier.mcmc.pareto$beta.store[-c(1:burnin),1] ;
  b2.hier = hier.mcmc.pareto$beta.store[-c(1:burnin),2];
  lambda.hier = hier.mcmc.pareto$lambda.store[-c(1:burnin), 1]
  mu.hier = hier.mcmc.pareto$mu.store[-c(1:burnin),]
  
  xgrid = seq(min(xi), max(xi), len=100)
  
  mu.preds.hier = sapply(xgrid, function(xx) rgamma(n = length(b1.hier), lambda.hier, lambda.hier/(exp(b1.hier+b2.hier*xx))))
  
  mui.hier = apply(mu.preds.hier, 2, mean)
  mui_low_hier = apply(mu.preds.hier, 2, quantile, probs=0.025)
  mui_up_hier = apply(mu.preds.hier, 2, quantile, probs=0.975)
  
  dev.off()
  
  # Replicate Values
  yi.reps = matrix(NA, nrow=length(b1.hier), ncol=length(xi))
  for(i in 1:ncol(yi.reps)){
    yi.reps[,i] = rpois(length(b1.hier), mu.hier[,i])
  }
  
  # Pt est and CI
  yi.pt.est.hier = apply(yi.reps, 2, mean)
  yi.low.hier = apply(yi.reps, 2, quantile, probs=0.025)
  yi.upper.hier = apply(yi.reps, 2, quantile, probs=0.975)
  posterior.resids = data.frame(sweep(yi.reps, 2 , yi))
  colnames(posterior.resids) = seq(1:32)
  
  curve(pplc(yi.reps, yi, x), 
        xlim=c(0,10),
        ylim=c(400,1000),
        col="blue",
        xlab="k", lty = 3,
        ylab="PPLC (k)", lwd = 2,
        main="Quadratic Loss Measure")
  curve(pplc(posterior.yi, yi, x), add = TRUE,
        col = "purple", lwd = 2, lty= 2)
  legend("topleft", legend = c("Hierarchical", "Non-Hierarchical"), 
         col = c("blue","purple"), lty=c(3,2), lwd=2)
  
  pplc.hier = pplc(yi.reps, yi, k=1)
  pplc.nohier = pplc(posterior.yi, yi, k=1)
  pplc.tab = round(cbind(pplc.hier, pplc.nohier),1)%>%`rownames<-`("PPLC (k = 1)")

}