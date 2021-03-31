library(tidyverse)
library(data.table)
library(parallel)
library(doParallel)

#' Approximate pdf for the Wiener diffusion model 
#' 
#' (from Lee, Fuss & Navarro, 2007)
#' @param rt the reaction time, 
#' @param theta is the threshold, 
#' @param mu is the drift rate, 
#' @param delta delta is the non-decision time (or ter)
ddiffusion <- function(rt,theta,delta,mu) { 
  1/sqrt(2*pi)*(2*theta+mu*(rt-delta))/(2*(rt-delta)^1.5)*exp(-(2*theta-mu*(rt-delta))^2/(2*(rt-delta))) }

#' Logarithm of the DDM PDF
#'
#' optimized compared to using log(ddiffusion) by cancelling the exponential against the logarithm
#' @param rt the reaction time, 
#' @param theta is the threshold, 
#' @param mu is the drift rate, 
#' @param delta delta is the non-decision time (or ter)
log_ddiffusion <- function(rt,theta,delta,mu) {  log((2*theta+mu*(rt-delta))/(2*(rt-delta)^1.5)) - log(sqrt(2*pi)) - (2*theta-mu*(rt-delta))^2/(2*(rt-delta)) }

#' Recinormal PDF
#' 
#' @param rt the reaction time
#' @param mu the reciprocal mu
#' @param sig the reciprocal sig
drecinorm <- function(rt, mu, sig){
  1/(rt*rt*sqrt(2*pi*sig*sig))*exp(-(mu*rt-1)^2/2/sig^2/rt^2)
}

# Read the data and add column names, factor labels, etc. Also scale RTs to seconds rather than ms for compatibility with existing code from previous project.
readData <- function (filename) {
  fullname <- file.path('.','data',filename)
  d = read.table(fullname)
  d<-as.data.table(d)
  
  names(d) <- c("tno","config","config2","tside","tpos","d1pos","d2pos","d3pos","tcol","tori","d1ori","d2ori","d3ori",
                "tposprime","tcolprime","toriprime","rt","error")
  d$sub<-substr(filename,1,2)
  d$session<-substr(filename,3,3)
  
  d$config<-factor(d$config, levels=c(1,2), labels=c("square", "diamond"))
  d$tside<-factor(d$tside, levels=c(0,1,2), labels=c("vertical","left", "right"))
  d$tcol<-factor(d$tcol, levels=c(0,1), labels=c("red","green"))
  d$tori<-factor(d$tori, levels=c(0,1), labels=c("top", "bottom"))
  d$d1ori<-factor(d$d1ori, levels=c(0,1), labels=c("top", "bottom"))
  d$d2ori<-factor(d$d2ori, levels=c(0,1), labels=c("top", "bottom"))
  d$d3ori<-factor(d$d3ori, levels=c(0,1), labels=c("top", "bottom"))
  d$tposprime<-factor(d$tposprime, levels=c(0,1,3,2), labels=c("NA", "TT", "TN", "TD"), ordered = TRUE) 
  d$tcolprime<-factor(d$tcolprime, levels=c(0,1), labels=c("same", "different"))
  d$toriprime<-factor(d$toriprime, levels=c(0,1), labels=c("same", "different"))
  d$error<-factor(d$error, levels=c(0,1,2), labels=c("no", "yes", "timeout"))
  d$sub<-factor(d$sub)
  d$rt<-d$rt/1000 
  
  d
}

#' so the Bayesian updating framework can be separated as follows:
#' 1. Calculate trial-wise prior and posterior updates
#'    two methods: simulation, and analytical approach with hyper parameter integration 
#'    This includes memory component (forgetting): $u_i = (1-m) u_0 + m * u_{i-1}$ (similar to Kalman filter)
#' 2. Estimate starting point $S_0$ based one logPrior
#' 3. Estimate likelihood of all trials for given parameters and models

#' prior updates with parameter integration
#' 
#' Updating priors based on Beta distribution
#' When m = 1, keep all updates (full memory). Otherwise, partial leakage
#' In this model, we assume participants has some knowlege of the mean of the prior (i.e., mu)
#' But it is somehow forgetting or partially integrated. 
#' Participants still use some original prior (or partially integrated). Mathematically:
#' mu_update = (1-m) * mu_0 + m * mu
#' Another parameter - sample size v remains. Thus, Beta distribution parameters alpha and beta
#' can be expressed as : alpha = mu * v, beta = (1-u)v
#' It turns out that two beta parameters are updated as follows:
#' alpha_update = (1-m)* alpha_0 + m * alpha
#' beta_update = (1-m)*beta_0 + m * beta
#' 
#' The above approach is very similar to the simulation approach where two prior distributions are fused. 
#' It is a normal approximation of fusion of two beta distribution. 
#' @param targets An array indicate target or non-target. logical FALSE, TRUE
#' @param a Beta distribution parameter a
#' @param b Beta distribution parameter b
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
updatePriorA <- function(targets, a = 1, b = 1, m = 0.5){
  a0 = a
  b0 = b
  # with memory approach
  mu = rep(0,length(targets))
  mu[1] = a/(a+b) # initial prior: equal (a/(a+b))
  for(idx in 2: (length(targets))) {
    # update parameters after trial n
    if (targets[idx-1]){
      a = a + 1
    } else {
      b = b + 1
    }
    mu_cur = a/(a+b) # current prior
    # partial reset (forgetting)
    mu[idx] = (1-m)*mu[1] + mu_cur*m
    
    # similarly, update a and b 
    a = (1-m)*a0 + m*a
    b = (1-m)*b0 + m*b
  }
  return(mu)
}

#' prior updates - simulation approach
#' 
#' Updating priors based on Beta distribution
#' @param targets An array indicate target or non-target. It must be 0,1 or FALSE, TRUE
#' @param a Beta distribution parameter a
#' @param b Beta distribution parameter b
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
updatePrior <- function(targets, a = 1, b = 1, m = 0.5){
  x<-seq(0,1,0.001)
  p0<-dbeta(x,a,b) # assume a and b are the same initial
  p0 <- ifelse(p0==Inf, 20, p0) # replace Inf with 20 (this happens for Jeffrey prior)
  pd <- p0
  mu = rep(a/(a+b),length(targets)+1)
  idx = 2
  for(target in targets) {
    # Update the distribution based on the Bernoulli likelihood 
    pd <- pd*x*target + (1-x)*pd*(1-target)
    # Implement with 'forgetting'
    pd <- (1-m)*p0 + m * pd
    # Normalize the probability distribution
    pd<-pd/mean(pd) 
    mu[idx] <- mean(pd*x)
    idx  = idx + 1
  }
  # return priors
  return(mu[1:length(targets)])
}

#' @param targets An array indicate target or non-target. It must be 0,1 or FALSE, TRUE
#' @param a Beta distribution parameter a
#' @param b Beta distribution parameter b
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
updatePosPrior <- function(targets, tpos, a = 1, b = 1, m = 0.5, d=0.5, multi=0) {
  x <- seq(0,1,0.001)
  p0 <- dbeta(x,a,b)
  p0 <- ifelse(p0==Inf, 20, p0) # replace Inf with 20 (this happens for Jeffrey prior)  
  p0 <- matrix(rep(p0,each=8),8,length(x))
  pd <- p0
  mu = rep(a/(a+b),length(targets)+1)
  idx = 2
  for(target in targets) {
    # Update the distribution based on the Bernoulli likelihood 
    pd[tpos[idx-1],] <- pd[tpos[idx-1],]*x*target + (1-x)*pd[tpos[idx-1],]*(1-target)
    if (multi==1) {
      pd[(tpos[idx-1]+1-1) %% 8 + 1,]<-pd[(tpos[idx-1]+1-1) %% 8 + 1,]*(d*x+(1-d)/2)*target + 
        pd[(tpos[idx-1]+1-1) %% 8 + 1,]*(1 - (d*x+(1-d)/2))*(1-target)
      pd[(tpos[idx-1]-1-1) %% 8 + 1,]<-pd[(tpos[idx-1]-1-1) %% 8 + 1,]*(d*x+(1-d)/2)*target +
        pd[(tpos[idx-1]-1-1) %% 8 + 1,]*(1 - (d*x+(1-d)/2))*(1-target)
      pd[(tpos[idx-1]+2-1) %% 8 + 1,]<-pd[(tpos[idx-1]+2-1) %% 8 + 1,]*(d^2*x+(1-d^2)/2)*target + 
        pd[(tpos[idx-1]+2-1) %% 8 + 1,]*(1 - (d^2*x+(1-d^2)/2))*(1-target) 
      pd[(tpos[idx-1]-2-1) %% 8 + 1,]<-pd[(tpos[idx-1]-2-1) %% 8 + 1,]*(d^2*x+(1-d^2)/2)*target + 
        pd[(tpos[idx-1]-2-1) %% 8 + 1,]*(1 - (d^2*x+(1-d^2)/2))*(1-target)
      pd[(tpos[idx-1]+3-1) %% 8 + 1,]<-pd[(tpos[idx-1]+3-1) %% 8 + 1,]*(d^3*x+(1-d^3)/2)*target + 
        pd[(tpos[idx-1]+3-1) %% 8 + 1,]*(1 - (d^3*x+(1-d^3)/2))*(1-target) 
      pd[(tpos[idx-1]-3-1) %% 8 + 1,]<-pd[(tpos[idx-1]-3-1) %% 8 + 1,]*(d^3*x+(1-d^3)/2)*target + 
        pd[(tpos[idx-1]-3-1) %% 8 + 1,]*(1 - (d^3*x+(1-d^3)/2))*(1-target)
      pd[(tpos[idx-1]+4-1) %% 8 + 1,]<-pd[(tpos[idx-1]+4-1) %% 8 + 1,]*(d^4*x+(1-d^4)/2)*target + 
        pd[(tpos[idx-1]+4-1) %% 8 + 1,]*(1 - (d^4*x+(1-d^4)/2))*(1-target)
    }
    # Implement 'forgetting'
    pd <- (1-m)*p0 + m * pd
    # Normalize the probability distribution
    pd <- pd/rowMeans(pd)
    mu[idx] <- mean(pd[tpos[idx],]*x)
    idx  = idx + 1
  }
  # return priors
  return(mu[1:length(targets)])
}

#' @param targets An array indicate target or non-target. It must be 0,1 or FALSE, TRUE
#' @param a Beta distribution parameter a
#' @param b Beta distribution parameter b
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
updatePosPriorSpread <- function(targets, tpos, a = 1, b = 1, m = 0.5, d=0.5) {
  x <- seq(0,1,0.001)
  p0 <- dbeta(x,a,b)
  p0 <- ifelse(p0==Inf, 20, p0) # replace Inf with 20 (this happens for Jeffrey prior)  
  p0 <- matrix(rep(p0,each=8),8,length(x))
  pd <- p0
  mu = rep(a/(a+b),length(targets)+1)
  idx = 2
  
  w <- function(x) {
    return(dnorm(x,0,d*4))
  }
  
  for(target in targets) {
    # Update the distribution based on the Bernoulli likelihood 
    pd[tpos[idx-1],] <- pd[tpos[idx-1],]*x*target + (1-x)*pd[tpos[idx-1],]*(1-target)
    
    # Implement 'forgetting'
    pd <- (1-m)*p0 + m * pd
    # Normalize the probability distribution
    pd <- pd/rowMeans(pd)
    mu[idx] <- mean(pd[tpos[idx],]*x)
    
    # Implement 'spreading'
    pd_old <- pd 
    for(ind in 1:8) {
      pd[ind,]<-(pd_old[ind,]*w(0) + pd_old[(ind+1-1) %% 8 + 1,]*w(1)+pd_old[(ind-1-1) %% 8 + 1,]*w(1)+
                   pd_old[(ind+2-1) %% 8 + 1,]*w(2)+pd_old[(ind-2-1) %% 8 + 1,]*w(2)+
                   pd_old[(ind+3-1) %% 8 + 1,]*w(3)+pd_old[(ind-3-1) %% 8 + 1,]*w(3)+
                   pd_old[(ind+4-1) %% 8 + 1,]*w(4))/(w(0)+2*w(1)+2*w(2)+2*w(3)+w(4))
    }
    idx  = idx + 1
  }
  # return priors
  return(mu[1:length(targets)])
}

#' prior update - non-bayesian approach
#' @param targets An array indicate target or non-target. It must be 0,1 or FALSE, TRUE
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
#' @param delta determines the amount that the starting point is changed by on each trial
updateS0 <- function(targets, m = 0.5, delta=0.1){
  targets <- 2*(targets-0.5) 
  s = rep(0, length(targets))
  for (i in 2:length(targets)) {
    s[i] = s[i-1] + targets[i-1]*delta 
    s[i] = s[i]*m
  }
  p=exp(s)/(1+exp(s))
  return(p)
}

#' prior update - non-bayesian approach
#' @param targets An array indicate target or non-target. It must be 0,1 or FALSE, TRUE
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
#' @param delta determines the amount that the starting point is changed by on each trial
updatePosS0 <- function(targets, tpos, m = 0.5, delta=0.1, d=0.5, multi=0){
  targets <- 2*(targets-0.5) 
  s <- rep(0, length(targets))
  ss <- matrix(rep(s,each=8),8,length(s))
  for (i in 2:length(targets)) {
    ss[tpos[i-1],i] <- ss[tpos[i-1],i-1] + targets[i-1]*delta 
    if(multi==1) {
      ss[(tpos[i-1]+1-1) %% 8 + 1,i] <- ss[(tpos[i-1]+1-1) %% 8 + 1,i] + targets[i-1]*d*delta
      ss[(tpos[i-1]-1-1) %% 8 + 1,i] <- ss[(tpos[i-1]-1-1) %% 8 + 1,i] + targets[i-1]*d*delta
      ss[(tpos[i-1]+2-1) %% 8 + 1,i] <- ss[(tpos[i-1]+2-1) %% 8 + 1,i] + targets[i-1]*d^2*delta
      ss[(tpos[i-1]-2-1) %% 8 + 1,i] <- ss[(tpos[i-1]-2-1) %% 8 + 1,i] + targets[i-1]*d^2*delta
      ss[(tpos[i-1]+3-1) %% 8 + 1,i] <- ss[(tpos[i-1]+3-1) %% 8 + 1,i] + targets[i-1]*d^3*delta
      ss[(tpos[i-1]-3-1) %% 8 + 1,i] <- ss[(tpos[i-1]-3-1) %% 8 + 1,i] + targets[i-1]*d^3*delta
      ss[(tpos[i-1]+4-1) %% 8 + 1,i] <- ss[(tpos[i-1]+4-1) %% 8 + 1,i] + targets[i-1]*d^4*delta
    }
    ss[,i] <- ss[,i]*m
    s[i] <- ss[tpos[i],i]
  }
  p=exp(s)/(1+exp(s))
  return(p)
}

#' update of the drift rate
#' 
#' Updating drift rate so that it is reduced on dimension switch trials, this version has a memory of one trial back and doesn't
#' care whether there has been many repeats/switches in a row. 
#' @param intertrial An array that indicate dimension repeat (1), switch (0) or target absent (NA)
#' @param sc Scaling parameter, should be <= 1: the drift rate is scaled down by this amount after a dimension switch
#' When m = 1, keep all updates (full memory). Otherwise, partial leakage
updateRate <- function(intertrial, sc) {
  scale <- intertrial + (1-intertrial)*sc
  scale[is.na(scale)] <- 1
  return(scale)
}

#' update of the drift rate - with memory
#' 
#' Updating drift rate so that it is reduced on dimension switch trials, this version has a memory of more than one trial back
#' the scaling factor is reduced after a switch and increased after a repeat, but with exponential discounting of old trials 
#' @param intertrial An array that indicate dimension repeat (1), switch (0) or target absent (NA)
#' @param mem memory parameter, m=0 means no updating, m=1 means old trials matter as much as recent trials
#' @param delta determines the amount that the scaling factor is changed by on each trial
updateRateMem <- function(intertrial, mem, delta) {
  intertrial <- 2*(intertrial-0.5) # rescale intertrial so that a switch is represented by -1 instead of 0
  intertrial[is.na(intertrial)] <- 0 # Target absent trials are neither a dimension switch nor a repeat
  s = rep(1, length(intertrial))
  for (i in 2:length(intertrial)) {
    s[i] = max(0.001, s[i-1] + intertrial[i]*delta) # make sure it is positive
    s[i] = s[i]*mem + (1-mem)
  }
  s[intertrial==0] <- 1 # Dimension switch costs don't affect target absent trials
  return(s)
}

#' update of the drift rate - dimension weighting
#' 
#' Updating drift rate so that it is reduced on dimension switch trials, this version attempts to do something closer to the spirit
#' of dimension weighting by updating based on the dimension itself rather than whether there has been a repeat or switch: 
#' ds reflects how much "weight" has been shifted from orientation to color, consequently the rate is scaled by s+ds on a color trial
#' and by s-ds on an orientation trial.
#' 
#' @param dim An array that indicate the target dimension color, orientation or absent (NA)
#' @param mem memory parameter, m=0 means no updating, m=1 means old trials matter as much as recent trials
#' @param delta determines the amount that the scaling factor is changed by on each trial
updateRateWeight <- function(dim, mem, delta) {
  dim[dim>0] <- 2*(dim[dim>0]-1.5) # rescale intertrial so that color orientation and absent are -1, 1 and 0
  ds = rep(0, length(dim))
  s = rep(1, length(dim))
  for (i in 2:length(dim)) {
    ds[i] = ds[i-1] + dim[i-1]*delta
    ds[i] = ds[i]*mem 
    s[i]<-s[i]+dim[i]*ds[i]
    s[i]=max(0,s[i])
  }
  return(s)
}

updateRateWeightPosition <- function(tarpos, mem, delta) {
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  s = rep(1, length(tarpos))
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] + delta
    ds[i, -tarpos[i-1]] <- ds[i-1, -tarpos[i-1]] - delta/7
    ds[i, ] <- ds[i, ] * mem
    s[i] <- s[i] + ds[i, tarpos[i]]
  }
  return(s)
}

updateRateWeightPositionDist <- function(tarpos, distpos, mem, delta_tar, delta_dist) {
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  s = rep(1, length(tarpos))
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] + delta_tar
    ds[i, -tarpos[i-1]] <- ds[i-1, -tarpos[i-1]] - delta_tar/7
    ds[i, distpos[,i-1]] <- ds[i, distpos[,i-1]] - delta_dist/3
    ds[i, -distpos[,i-1]] <- ds[i, -distpos[,i-1]] + delta_dist/5
    ds[i, ] <- ds[i, ] * mem
    s[i] <- s[i] + ds[i, tarpos[i]]
  }
  return(s)
}

updateRateWeightPositionMatch <- function(tarpos, distpos, mem, delta) {
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  s = rep(1, length(tarpos))
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] + delta
    ds[i, distpos[,i-1]] <- ds[i, distpos[,i-1]] - delta/3
    ds[i, ] <- ds[i, ] * mem
    s[i] <- s[i] + ds[i, tarpos[i]]
  }
  return(s)
}

updateNDTWeight <- function(dim, mem = 0.5, delta=0.02){
  dim <- 2*(dim-0.5)   
  ds = rep(0, length(dim))
  s = rep(0, length(dim))
  for (i in 2:length(dim)) {
    ds[i] <- ds[i-1] - dim[i-1]*delta
    ds[i] <- ds[i]*mem 
    s[i] <- ds[i]*dim[i]
  }
  return(s)
}

updateNDTPosition <- function(tarpos, mem, delta) {
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  s = rep(0, length(tarpos))
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] - delta
    ds[i, -tarpos[i-1]] <- ds[i-1, -tarpos[i-1]] + delta/7
    ds[i, ] <- ds[i, ] * mem
    s[i] <- ds[i, tarpos[i]]
  }
  return(s)
}

updateNDTPositionDist <- function(tarpos, distpos, mem, delta_tar, delta_dist) {
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  s = rep(0, length(tarpos))
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] - delta_tar
    ds[i, -tarpos[i-1]] <- ds[i-1, -tarpos[i-1]] + delta_tar/7
    ds[i, distpos[,i-1]] <- ds[i, distpos[,i-1]] + delta_dist/3
    ds[i, -distpos[,i-1]] <- ds[i, -distpos[,i-1]] - delta_dist/5
    ds[i, ] <- ds[i, ] * mem
    s[i] <- ds[i, tarpos[i]]
  }
  return(s)
}

updateNDTPositionMatch <- function(tarpos, distpos, mem, delta) {
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  s = rep(0, length(tarpos))
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] - delta
    ds[i, distpos[,i-1]] <- ds[i, distpos[,i-1]] + delta/3
    ds[i, ] <- ds[i, ] * mem
    s[i] <- ds[i, tarpos[i]]
  }
  return(s)
}

updatePosNDT <- function(targets, tpos, m = 0.5, delta = 0.01, d = 0.5, multi = 0){
  targets <- 2*(targets-0.5) 
  s <- rep(0, length(targets))
  ss <- matrix(rep(s,each=8),8,length(s))
  for (i in 2:length(targets)) {
    ss[tpos[i-1],i] <- ss[tpos[i-1],i-1] - targets[i-1]*targets[i]*delta 
    if(multi==1) {
      ss[(tpos[i-1]+1-1) %% 8 + 1,i] <- ss[(tpos[i-1]+1-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d*delta
      ss[(tpos[i-1]-1-1) %% 8 + 1,i] <- ss[(tpos[i-1]-1-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d*delta
      ss[(tpos[i-1]+2-1) %% 8 + 1,i] <- ss[(tpos[i-1]+2-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d^2*delta
      ss[(tpos[i-1]-2-1) %% 8 + 1,i] <- ss[(tpos[i-1]-2-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d^2*delta
      ss[(tpos[i-1]+3-1) %% 8 + 1,i] <- ss[(tpos[i-1]+3-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d^3*delta
      ss[(tpos[i-1]-3-1) %% 8 + 1,i] <- ss[(tpos[i-1]-3-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d^3*delta
      ss[(tpos[i-1]+4-1) %% 8 + 1,i] <- ss[(tpos[i-1]+4-1) %% 8 + 1,i] - targets[i-1]*targets[i]*d^4*delta
    }
    ss[,i] <- ss[,i]*m
    s[i] <- ss[tpos[i],i]
  }
  return(s)
}

updatePosRateWeight <- function(tpos, dim, mem, delta, d = 0.5, multi = 0) {
  dim <- 2*(dim-1.5) # rescale intertrial so that red and green are -1 and 1 
  ds = rep(0, length(dim))
  ds <- matrix(rep(ds,each=8),8,length(ds))
  s = rep(1, length(dim))
  ss <- matrix(rep(s,each=8),8,length(s))
  for (i in 2:length(dim)) {
    ds[tpos[i-1],i] = ds[tpos[i-1],i-1] + dim[i-1]*delta
    ds[,i] = ds[,i]*mem 
    ss[tpos[i-1],i]<-ss[tpos[i-1],i]+dim[i]*ds[tpos[i-1],i]
    if(multi==1) {
      ss[(tpos[i-1]+1-1) %% 8 + 1,i] <- ss[(tpos[i-1]+1-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d
      ss[(tpos[i-1]-1-1) %% 8 + 1,i] <- ss[(tpos[i-1]-1-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d
      ss[(tpos[i-1]+2-1) %% 8 + 1,i] <- ss[(tpos[i-1]+2-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d^2
      ss[(tpos[i-1]-2-1) %% 8 + 1,i] <- ss[(tpos[i-1]-2-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d^2
      ss[(tpos[i-1]+3-1) %% 8 + 1,i] <- ss[(tpos[i-1]+3-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d^3
      ss[(tpos[i-1]-3-1) %% 8 + 1,i] <- ss[(tpos[i-1]-3-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d^3
      ss[(tpos[i-1]+4-1) %% 8 + 1,i] <- ss[(tpos[i-1]+4-1) %% 8 + 1,i] + dim[i]*ds[tpos[i-1],i]*d^4
    }
    ss[tpos[i-1],i]=max(0,ss[tpos[i-1],i])
    s[i] <- ss[tpos[i],i]
  }
  return(s)
}

#' Estimate loglikelihood
#'
#' Estimate log likelihood for a given model (theta, mu, sig)
#' @param data input data, which should have one column rt, one column of updated prior, called beta, and one column target
#' @param theta decision boundary (theta1, theta2)
#' @param mu drift rate / ergodic rate (mu1, mu2)
#' @param sig sigma of the drift rate (sig1, sig2)
#' @param ter nondecision time (assume it is independent from response/target condition)
#' @param later_diffusion a flag for using LATER or Diffusion model as RT distribution
logll <- function(data, theta=4, mu=5,  sig=1, ter = 0, later_diffusion =1) {
  
  # data <-   data %>% ungroup(.) %>% filter(!(error | outlier)) %>%
  #   mutate( #beta = max(0.001, min(0.999, beta)), # constrain beta in 0.01 - 0.99  # min/max not working in mutate!!!!
  #     s0 = log(beta/(1-beta)) * ((targets>0) *2 -1),
  #     mu = mu*scale,  sig = sig,  theta = theta,
  #     rs = 1/(rt - ter), # add non-decision time for LATER model
  #     delta =  theta-s0) 
  # LATER model  
  rt = data$rt
  scale = data$scale
  delta = theta - data$s0
  ter = data$ndt
  if(min(rt - ter)<0) {
    return(1e10)
  }
  if(later_diffusion==1) {
    # use recinormal pdf (i.e., rt pdf, comparible to DDM)
    nll = -sum(log(drecinorm(rt-ter, mu*scale/delta, sig/delta)))
    # nll = -sum(log(dnorm(1/(rt - ter), mu*scale/delta, sig/delta)))
  } else {
    # put back scaling parameter sig back, given that delta is not scaled
    nll = -sum(log_ddiffusion(rt, delta/sig, ter, mu/sig*scale))
  }
  if (nll == Inf | nll == -Inf) # avoid optim error with L-BFGS-B method
    nll = 1e10
  return(nll)
  #  if(later_diffusion==1) {
  #    data %>% mutate(delta = theta - s0, ll =  log(dnorm(1/(rt - ter), mu*scale/delta, sig/delta))) %>%
  #      summarise(logll = -sum(ll))
  #  } else {
  #    data %>% mutate(delta = theta - s0, ll = log_ddiffusion(rt, delta, ter, mu*scale)) %>%
  #      summarise(logll = -sum(ll))    
  #  }
  
}

#' for given parameters calculate negative log likelihood
#' 
#' @param x to-be-optimzed parameter
#' @param data data frame contain RT data
#' @param fixed pass fixed parameters here. Same length as x. NA for to-be-optimzed parameter
#' @param op logical flag for optimization or returning parameters
#' @return negloglikelihood return negative log likelihood
findParameters <- function(par1, data, fixed, op = TRUE) {
  tryCatch({
    ll<-0 # log likelihood
    # recombine to-be-fitted parameterd and fixed parameters
    x = fixed
    x[is.na(fixed)] = par1
    data$beta <- 0.5
    data$scale <- 1
    data$ndt <- x[["ter"]]
    
    data$inttrial_resp <- ifelse(data$toriprime=="same", 1, 0)
    data$inttrial_col <- ifelse(data$tcolprime=="same", 1, 0)
    data$inttrial_pos <- ifelse(data$tposprime=="TT", 1, 0)
    
    distpos <- t(matrix(c(data$d1pos, data$d2pos, data$d3pos),nrow(data),3))
    
    # Response based updating 
    switch(x['resp_update'], # If 0 do nothing
           data$beta <- updateS0(data$tori == "top", delta = x["u_delta_resp"], m = x["u_mem_resp"]), # 1: Direct S0 updating with forgetting
           data$beta <- updatePrior(data$tori == "top", a = x["beta_resp"], b=x["beta_resp"], m = x["mem_resp"]), # 2: Bayesian S0 updating with forgetting
           data$beta <- updatePosPrior(data$tori == "top", tpos=data$tpos, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]], multi=0), # 3: Position linked S0 updating with forgetting
           data$beta <- updatePosPrior(data$tori == "top", tpos=data$tpos, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]], d = x[["d_resp"]], multi=1), # 4: Multiple position linked S0 updating with forgetting
           data$scale <- updateRate(data$inttrial_resp, sc=x["scale_resp"]), # 5: single trial back switch cost on DR
           data$scale <- updateRateMem(data$inttrial_resp, mem=x["u_mem_resp"], delta=x["u_delta_resp"]), # 6: switch cost on DR with longer memory
           data$scale <- updateRateWeight(as.numeric(data$tori == "top")+1, mem=x["u_mem_resp"], delta=x["u_delta_resp"]), # 7 response based DR weighting 
           data$beta <- updatePosS0(data$tori == "top", tpos=data$tpos, delta = x["u_delta_resp"], m = x["u_mem_resp"], multi=0), # 8 Position linked direct S0 updating with forgetting
           data$beta <- updatePosS0(data$tori == "top", tpos=data$tpos, delta = x["u_delta_resp"], m = x["u_mem_resp"], d=x[['d_resp']], multi=1), # 9 Multiple position linked direct S0 updating with forgetting 
           data$ndt[!is.na(data$inttrial_resp)] <- data$ndt[!is.na(data$inttrial_resp)] - 
             data$inttrial_resp[!is.na(data$inttrial_resp)]*x[['delta_ndt_resp']], # 10: NDT binary
           data$ndt <- data$ndt + updateNDTWeight(data$tori == "top", mem = x[["m_ndt_resp"]], delta = x[["delta_ndt_resp"]]), # 11: weighted NDT
           data$ndt <- data$ndt + updatePosNDT(data$tori == "top", tpos=data$tpos, delta = x[["delta_ndt_resp"]], m = x[["m_ndt_resp"]], multi=0), # 12 Position linked NDT updating
           data$ndt <- data$ndt + updatePosNDT(data$tori == "top", tpos=data$tpos, delta = x[["delta_ndt_resp"]], m = x[["m_ndt_resp"]], d=x[['d_resp']], multi=1), # 13 Multiple position linked NDT updating
           data$beta <- updatePosPriorSpread(data$tori == "top", tpos=data$tpos, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]], d = x[["d_resp"]]) # 14: Multiple position linked S0 updating with forgetting and spreading
    )
    # Color based updating
    switch(x['col_update'], # If 0 do nothing
           data$scale <- data$scale*updateRate(data$inttrial_col, sc=x["scale_col"]), # 1: single trial back switch cost on DR
           data$scale <- data$scale*updateRateMem(data$inttrial_col, mem=x["u_mem_col"], delta=x["u_delta_col"]), # 2: switch cost on DR with longer memory
           data$scale <- data$scale*updateRateWeight(as.numeric(data$tcol), mem=x["u_mem_col"], delta=x["u_delta_col"]), # 3: response based DR weighting 
           data$ndt[!is.na(data$inttrial_col)] <- data$ndt[!is.na(data$inttrial_col)] -
             data$inttrial_col[!is.na(data$inttrial_col)]*x[['delta_ndt_col']], # 4: NDT binary
           data$ndt <- data$ndt + updateNDTWeight(as.numeric(data$tcol=="red"), mem = x[["m_ndt_col"]], delta = x[["delta_ndt_col"]]), # 5: weighted NDT 
           data$scale <- data$scale*updatePosRateWeight(data$tpos, as.numeric(data$tcol), mem=x["u_mem_col"], delta=x["u_delta_col"]), # 6: Position-linked weighted rate 
           data$scale <- data$scale*updatePosRateWeight(data$tpos, as.numeric(data$tcol), mem=x["u_mem_col"], delta=x["u_delta_col"], d=x[['d_col']], multi=1) # 7: Position-linked weighted rate 
    )
    # Position based updating
    switch(x['pos_update'], # If 0 do nothing
           data$scale <- data$scale*updateRate(data$inttrial_pos, sc=x["scale_pos"]), # 1: single trial back switch cost on DR
           data$scale <- data$scale*updateRateMem(data$inttrial_pos, mem=x["u_mem_pos"], delta=x["u_delta_pos"]), # 2: switch cost on DR with longer memory
           data$scale <- data$scale*updateRateWeightPosition(as.numeric(data$tpos), mem=x["u_mem_pos"], delta=x["u_delta_pos"]), # 3: position based DR weighting 
           data$scale <- data$scale*updateRateWeightPositionDist(as.numeric(data$tpos), distpos, mem=x["u_mem_pos"], delta_tar=x["u_delta_pos"], delta_dist=x["u_delta_pos_dist"]), # 4: position based DR weighting with distractor inhibition 
           data$ndt[!is.na(data$inttrial_pos)] <- data$ndt[!is.na(data$inttrial_pos)] -
             data$inttrial_pos[!is.na(data$inttrial_pos)]*x[['delta_ndt_pos']], # 5: NDT binary
           data$ndt <- data$ndt + updateNDTPosition(as.numeric(data$tpos), mem = x[["m_ndt_pos"]], delta = x[["delta_ndt_pos"]]), # 6: position based NDT weighting
           data$ndt <- data$ndt + updateNDTPositionDist(as.numeric(data$tpos), distpos, mem = x[["m_ndt_pos"]], delta_tar = x[["delta_ndt_pos"]], delta_dist = x[["delta_ndt_pos_dist"]]), # 7: position based NDT weighting with distractor inhibition
           data$scale <- data$scale*updateRateWeightPositionMatch(as.numeric(data$tpos), distpos, mem=x[["u_mem_pos"]], delta = x[["u_delta_pos"]]), # 8: position based DR weighting with matched distractor inhibition 
           data$ndt <- data$ndt + updateNDTPositionMatch(as.numeric(data$tpos), distpos, mem = x[["m_ndt_pos"]], delta = x[["delta_ndt_pos"]]) # 9: position based NDT weighting with matched distractor inhibition
    )
    
    # later or diffusion
    # calculate log likelihood based how many diffusion /later processes 
    # the targets has implicit indication of process (e.g., targets contains 0,1,2 indicates 3 processes)
    # each find optimal parameters here (inner optimization)
    # Inner parameters: theta=x[1], mu=x[2],  sig=x[3], ter  = x[4] , later_diffusion =x[5]
    pars = x # parameters
    inner_par0 = c(theta = 3, mu = 6, sig = 1, ter = x[['ter']], later_diffusion = fixed[['later_diffusion']])
    inner_fixed = c(NA, NA, NA, x['ter'], fixed['later_diffusion']) 
    
    data <- data %>% ungroup(.) %>% filter(error=="no") %>%
      mutate( s0 = log(beta/(1-beta)) * ((tori == "top") *2 -1)) 
    # add contraints for parameters - theta, mu, sig, ter 
    # not working at the moment: return L-BFGS-B needs finite values of 'fn'
    #i_lowers = c(max(abs(data$s0)), 0.001, 0.001) #theta must be greater than s0
    #i_uppers = c(50, 50, 50) 
    ui1 = rbind(c(1,0,0),c(0,1,0),c(0,0,1))
    ci1 = c(max(abs(data$s0)), 0.001, 0.001)
    ui2 = -ui1
    ci2 = c(-50, -50, -50)
    ui = rbind(ui1, ui2)
    ci = c(ci1, ci2)
    
    if(sum(is.na(data$s0)>0))
    {
      ll <- 1e10
      return(ll)
    } else if(max(abs(data$s0)>inner_par0[1])) {
      ll <- 1e10
      return(ll)
    } else if(min(data$scale)<=0) {
      ll <- 1e10
      return(ll)
    }
    
    inner_par1 <- inner_par0[is.na(inner_fixed)] # find to-be-fixed parameters
    par <- constrOptim(inner_par1,  findInnerParameters, NULL,
                       ui = ui, ci = ci, 
                       data = data, fixed = inner_fixed)
    
    #      par <- optim(inner_par1,  findInnerParameters, data = subdata, fixed = inner_fixed)
    #      par <- optim(inner_par1,  findInnerParameters, data = subdata, fixed = inner_fixed,
    #                   lower = i_lowers, upper = i_uppers, method = 'BFGS')
    ll <- ll + par$value
    pars = c(pars, par$par)
    
    pars = c(pars, nll = ll)
    if (op)
      return(ll)
    else
      return(pars)
  }, 
  error = function(e){
    print(e)
    stop(e)
  })
}

# Generates a sequence of starting points and  drift rates
genSeq <- function(x, data) {
  
  data$beta <- 0.5
  data$scale <- 1
  data$ndt <- x[["ter"]]
  
  data$inttrial_resp <- ifelse(data$toriprime=="same", 1, 0)
  data$inttrial_col <- ifelse(data$tcolprime=="same", 1, 0)
  data$inttrial_pos <- ifelse(data$tposprime=="TT", 1, 0)
  
  distpos <- t(matrix(c(data$d1pos, data$d2pos, data$d3pos),nrow(data),3))
  
  # Response based updating 
  switch(x[['resp_update']], # If 0 do nothing
         data$beta <- updateS0(data$tori == "top", delta = x[["u_delta_resp"]], m = x[["u_mem_resp"]]), # 1: Direct S0 updating with forgetting
         data$beta <- updatePrior(data$tori == "top", a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]]), # 2: Bayesian S0 updating with forgetting
         data$beta <- updatePosPrior(data$tori == "top", tpos=data$tpos, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]], multi=0), # 3: Position linked S0 updating with forgetting
         data$beta <- updatePosPrior(data$tori == "top", tpos=data$tpos, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]], d = x[["d_resp"]], multi=1), # 4: Multiple position linked S0 updating with forgetting
         data$scale <- updateRate(data$inttrial_resp, sc=x[["scale_resp"]]), # 5: single trial back switch cost on DR
         data$scale <- updateRateMem(data$inttrial_resp, mem=x[["u_mem_resp"]], delta=x[["u_delta_resp"]]), # 6: switch cost on DR with longer memory
         data$scale <- updateRateWeight(as.numeric(data$tori == "top")+1, mem=x[["u_mem_resp"]], delta=x[["u_delta_resp"]]), # 7 response based DR weighting 
         data$beta <- updatePosS0(data$tori == "top", tpos=data$tpos, delta = x[["u_delta_resp"]], m = x[["u_mem_resp"]], multi=0), # 8 Position linked direct S0 updating with forgetting
         data$beta <- updatePosS0(data$tori == "top", tpos=data$tpos, delta = x[["u_delta_resp"]], m = x[["u_mem_resp"]], d=x[['d_resp']], multi=1), # 9 Multiple position linked direct S0 updating with forgetting 
         data$ndt[!is.na(data$inttrial_resp)] <- data$ndt[!is.na(data$inttrial_resp)] - 
           data$inttrial_resp[!is.na(data$inttrial_resp)]*x[['delta_ndt_resp']], # 10: NDT binary
         data$ndt <- data$ndt + updateNDTWeight(data$tori == "top", mem = x[["m_ndt_resp"]], delta = x[["delta_ndt_resp"]]), # 11: weighted NDT
         data$ndt <- data$ndt + updatePosNDT(data$tori == "top", tpos=data$tpos, delta = x[["delta_ndt_resp"]], m = x[["m_ndt_resp"]], multi=0), # 12 Position linked NDT updating
         data$ndt <- data$ndt + updatePosNDT(data$tori == "top", tpos=data$tpos, delta = x[["delta_ndt_resp"]], m = x[["m_ndt_resp"]], d=x[['d_resp']], multi=1), # 13 Multiple position linked NDT updating
         data$beta <- updatePosPriorSpread(data$tori == "top", tpos=data$tpos, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]], d = x[["d_resp"]]) # 14: Multiple position linked S0 updating with forgetting and spreading
  )
  # Color based updating
  switch(x[['col_update']], # If 0 do nothing
         data$scale <- data$scale*updateRate(data$inttrial_col, sc=x[["scale_col"]]), # 1: single trial back switch cost on DR
         data$scale <- data$scale*updateRateMem(data$inttrial_col, mem=x[["u_mem_col"]], delta=x[["u_delta_col"]]), # 2: switch cost on DR with longer memory
         data$scale <- data$scale*updateRateWeight(as.numeric(data$tcol), mem=x[["u_mem_col"]], delta=x[["u_delta_col"]]), # 3: response based DR weighting 
         data$ndt[!is.na(data$inttrial_col)] <- data$ndt[!is.na(data$inttrial_col)] - 
           data$inttrial_col[!is.na(data$inttrial_col)]*x[['delta_ndt_col']], # 4: NDT binary
         data$ndt <- data$ndt + updateNDTWeight(as.numeric(data$tcol=="red"), mem = x[["m_ndt_col"]], delta = x[["delta_ndt_col"]]), # 5: weighted NDT 
         data$scale <- data$scale*updatePosRateWeight(data$tpos, as.numeric(data$tcol), mem=x[["u_mem_col"]], delta=x[["u_delta_col"]]), # 6: Position-linked weighted rate 
         data$scale <- data$scale*updatePosRateWeight(data$tpos, as.numeric(data$tcol), mem=x[["u_mem_col"]], delta=x[["u_delta_col"]], d=x[['d_col']], multi=1) # 7: Position-linked weighted rate 
  )
  # Position based updating
  switch(x[['pos_update']], # If 0 do nothing
         data$scale <- data$scale*updateRate(data$inttrial_pos, sc=x[["scale_pos"]]), # 1: single trial back switch cost on DR
         data$scale <- data$scale*updateRateMem(data$inttrial_pos, mem=x[["u_mem_pos"]], delta=x[["u_delta_pos"]]), # 2: switch cost on DR with longer memory
         data$scale <- data$scale*updateRateWeightPosition(as.numeric(data$tpos), mem=x[["u_mem_pos"]], delta=x[["u_delta_pos"]]), # 3: position based DR weighting 
         data$scale <- data$scale*updateRateWeightPositionDist(as.numeric(data$tpos), distpos, mem=x[["u_mem_pos"]], delta_tar=x[["u_delta_pos"]], delta_dist=x[["u_delta_pos_dist"]]), # 4: position based DR weighting with distractor inhibition 
         data$ndt[!is.na(data$inttrial_pos)] <- data$ndt[!is.na(data$inttrial_pos)] -
           data$inttrial_pos[!is.na(data$inttrial_pos)]*x[['delta_ndt_pos']], # 5: NDT binary
         data$ndt <- data$ndt + updateNDTPosition(as.numeric(data$tpos), mem = x[["m_ndt_pos"]], delta = x[["delta_ndt_pos"]]), # 6: position based NDT weighting
         data$ndt <- data$ndt + updateNDTPositionDist(as.numeric(data$tpos), distpos, mem = x[["m_ndt_pos"]], delta_tar = x[["delta_ndt_pos"]], delta_dist = x[["delta_ndt_pos_dist"]]), # 7: position based NDT weighting with distractor inhibition
         data$scale <- data$scale*updateRateWeightPositionMatch(as.numeric(data$tpos), distpos, mem=x[["u_mem_pos"]], delta=x[["u_delta_pos"]]), # 8: position based DR weighting with matched distractor inhibition 
         data$ndt <- data$ndt + updateNDTPositionMatch(as.numeric(data$tpos), distpos, mem = x[["m_ndt_pos"]], delta = x[["delta_ndt_pos"]]) # 9: position based NDT weighting with matched distractor inhibition
  )
  data$rate <- x$mu*data$scale
  data$theta <- x$theta
  data$sigma <- x$sig
  data <- mutate(data, s0 = log(beta/(1-beta)) * ((tori == "top") *2 -1)) 
  
  t <- seq(0, 3, 0.01)
  delta <- data$theta - data$s0
  sigma <- data$sigma
  mu <- data$rate 
  ter <- data$ndt
  predRT <- rep(0, nrow(data))
  for(n in 1:nrow(data)) {
    if(x$later_diffusion=="LATER") {
      dist <- drecinorm(t-ter[n], mu[n]/delta[n], sigma[n]/delta[n])
    } else {
      dist <- ddiffusion(t, delta[n]/sigma[n], ter[n], mu[n]/sigma[n])
    }
    predRT[n] <- mean(t*dist, na.rm=TRUE)/mean(dist, na.rm=TRUE)
  }
  data$predRT <- predRT
  
  return(data)
}

# return inner LATER/DIFFUSION parameters
findInnerParameters <- function(par, data, fixed){
  x = fixed
  x[is.na(fixed)] = par
  logll(data, theta=x[1], mu=x[2],  sig=x[3], ter  = x[4] , later_diffusion =x[5]) 
}


#' fit parameters with given models
#' 
#' This function will find optimal parameters for a given model and give subject data
#' @param sub subject rt data frame
#' @param modelname model names. 
#' model name rules: 1 ("L" or "D"): First letter indicates which model is used L = LATER model, D = DDM
#' 2-3 ("RX"): indicates which type of updating is done based on the response history 
#' 4-5 ("CX"): indicates which type of updating is done based on the color history
#' 6-7 ("PX"): indicates which type of updating is done based on the position history
fitPar <- function(sub, modelname = 'LR0C0P0') {
  # parameters for outer fit: mem, beta = parameters of prior updating (memory and beta prior),  
  # scale,  u_mem, u_delta = parameters of rate updating,  ter = non-decision time,
  # flags: later_diffusion, resp_update (response related updating, 0: no update, 1-2: prior updating, 3-5 rate updating), 
  # dim_update (dimension related updating, 0: no update, 1-2: prior updating, 3-5 rate updating)
  par0 <- c(mem_resp = 0.5, beta_resp = 2, d_resp=0.5, mem_col = 0.5, beta_col = 2, mem_pos = 0.5, beta_pos = 0.5,
            scale_resp = 1, u_mem_resp =0.5, u_delta_resp = 0.1,
            scale_col = 1, d_col=0.5, u_mem_col = 0.5, u_delta_col = 0.1,
            scale_pos = 1, u_mem_pos = 0.5, u_delta_pos = 0.1, u_delta_pos_dist = 0.1,
            ter = 0.01, m_ndt_resp = 0.5, delta_ndt_resp = 0.01, m_ndt_col = 0.5, delta_ndt_col = 0.01, 
            m_ndt_pos = 0.5, delta_ndt_pos = 0.01, delta_ndt_pos_dist = 0.01, later_diffusion = 1, 
            resp_update = 2, col_update = 1, pos_update = 1)
  fixed = c(mem_resp=NA, beta_resp=NA, d_resp=NA, mem_col=NA, beta_col=NA, mem_pos=NA, beta_pos=NA,
            scale_resp=NA, u_mem_resp=NA,  u_delta_resp=NA, 
            scale_col=NA, d_col=NA, u_mem_col=NA,  u_delta_col=NA, 
            scale_pos=NA, u_mem_pos=NA,  u_delta_pos=NA, u_delta_pos_dist = NA,
            ter=NA, m_ndt_resp = NA, delta_ndt_resp = NA, m_ndt_col = NA, delta_ndt_col = NA, 
            m_ndt_pos = NA, delta_ndt_pos = NA, delta_ndt_pos_dist = NA, later_diffusion=1, 
            resp_update=0, col_update=0, pos_update=0)
  
  # Later model or DDM
  switch(substr(modelname, 1, 1),
         'L' = {# LATER 
           fixed['later_diffusion'] = 1
         },
         'D' = {# Diffusion 
           fixed['later_diffusion'] = 0
         }
  )
  # Response history based updating: 0 = no updating, 1 = Bayesian S0 updating with full memory, 2 = Bayesian S0 updating with forgetting,
  # 3 = single trial back switch cost on DR, 4 = switch cost on DR with longer memory, 5 = response based DR weighting 
  switch(substr(modelname, 2, 3),
         'R0' = { # No updating based on response history
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, 0, 0, 0, 0, 0)
         }, 
         'R1' = { # Direct S0 updating with forgetting
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 1)
         },
         'R2' = { # Bayesian S0 updating with forgetting
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(NA, NA, 0, 1, 0, 0, 0, 0, 2)
         },
         'R3' = { # Bayesian S0 updating with forgetting
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(NA, NA, 0, 1, 0, 0, 0, 0, 3)
         },
         'R4' = { # Bayesian S0 updating with forgetting
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(NA, NA, NA, 1, 0, 0, 0, 0, 4)
         }, 
         'R5' = { # single trial back switch cost on DR
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, NA, 0, 0, 0, 0, 5)
         },
         'R6' = { # switch cost on DR with longer memory
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 6)
         },
         'R7' = { # response based DR weighting 
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 7)
         },
         'R8' = { # response based direct S0 updating
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 8)
         },
         'R9' = { # response based direct S0 updating
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, NA, 1, NA, NA, 0, 0, 9)
         },
         'RA' = { # response based binary NDT
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, 0, 0, 0, NA, 10)
         },
         'RB' = { # response based NDT weighting 
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, 0, 0, NA, NA, 11)
         },
         'RC' = { # response based binary NDT
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, 0, 1, 0, 0, NA, NA, 12)
         },
         'RD' = { # response based NDT weighting 
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(0, 1, NA, 1, 0, 0, NA, NA, 13)
         },
         'RE' = { # Bayesian S0 updating with forgetting
           fixed[c('mem_resp','beta_resp', 'd_resp', 'scale_resp','u_mem_resp','u_delta_resp',
                   'm_ndt_resp','delta_ndt_resp','resp_update')] <- c(NA, NA, NA, 1, 0, 0, 0, 0, 14)
         }
  )
  
  # Target color based updating: 0 = no updating, 1 = Bayesian S0 updating with full memory, 2 = Bayesian S0 updating with forgetting,
  # 3 = single trial back switch cost on DR, 4 = switch cost on DR with longer memory, 5 = dimension based DR weighting 
  switch(substr(modelname, 4, 5),
         'C0' = { #  no updating
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, 1, 0, 0, 0, 0, 0)
         },
         'C1' = { # single trial back switch cost on DR
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, NA, 0, 0, 0, 0, 1)
         }, 
         'C2' = { # switch cost on DR with longer memory
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 2)
         }, 
         'C3' = { # dimension based DR weighting 
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 3)
         },
         'C4' = { # dimension based DR weighting 
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, 1, 0, 0, 0, NA, 4)
         },
         'C5' = { # dimension based DR weighting 
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, 1, 0, 0, NA, NA, 5)
         },
         'C6' = { # dimension based DR weighting 
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, 0, 1, NA, NA, 0, 0, 6)
         },
         'C7' = { # dimension based DR weighting 
           fixed[c('mem_col','beta_col','d_col','scale_col','u_mem_col','u_delta_col', 'm_ndt_col','delta_ndt_col',
                   'col_update')] <- c(0, 1, NA, 1, NA, NA, 0, 0, 7)
         }       
  )
  
  # Target position based updating: 0 = no updating, 1 = Bayesian S0 updating with full memory, 2 = Bayesian S0 updating with forgetting,
  # 3 = single trial back switch cost on DR, 4 = switch cost on DR with longer memory, 5 = dimension based DR weighting 
  switch(substr(modelname, 6, 7),
         'P0' = { #  no updating
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0)
         },
         'P1' = { # single trial back switch cost on DR
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <- 
             c(0, 1, NA, 0, 0, 0, 0, 0, 0, 1)
         }, 
         'P2' = { # switch cost on DR with longer memory
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, NA, NA, 0, 0, 0, 0, 2)
         }, 
         'P3' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, NA, NA, 0, 0, 0, 0, 3)
         }, 
         'P4' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, NA, NA, NA, 0, 0, 0, 4)
         }, 
         'P5' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, 0, 0, 0, 0, NA, 0, 5)
         }, 
         'P6' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <- 
             c(0, 1, 1, 0, 0, 0, NA, NA, 0, 6)
         }, 
         'P7' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <- 
             c(0, 1, 1, 0, 0, 0, NA, NA, NA, 7)
         }, 
         'P8' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, NA, NA, 0, 0, 0, 0, 8)
         }, 
         'P9' = { # dimension based DR weighting 
           fixed[c('mem_pos','beta_pos','scale_pos','u_mem_pos','u_delta_pos','u_delta_pos_dist',
                   'm_ndt_pos','delta_ndt_pos','delta_ndt_pos_dist','pos_update')] <-
             c(0, 1, 1, 0, 0, 0, NA, NA, 0, 9)
         }
  )
  
  par1 <- par0[is.na(fixed)] # find to-be-fixed parameters
  #  t0 = proc.time()
  
  # add constraints on parameters (mem, beta, u_scale, u_mem, u_delta, non_decision)
  lowers = c(mem_resp=0, beta_resp=0.5, d_resp=0, mem_col=0, beta_col=0.5, mem_pos=0, beta_pos=0.5, scale_resp=0, 
             u_mem_resp=0,  u_delta_resp=0, scale_col=0, d_col=0, u_mem_col=0,  u_delta_col=0,
             scale_pos=0, u_mem_pos=0,  u_delta_pos=0, u_delta_pos_dist=0, ter=0, 
             m_ndt_resp = 0, delta_ndt_resp = 0, m_ndt_col = 0, delta_ndt_col = 0, 
             m_ndt_pos = 0, delta_ndt_pos = 0, delta_ndt_pos_dist = 0)
  lowers = lowers[is.na(fixed[1:26])] # contraints for to-be-fixed params
  uppers = c(mem_resp=1, beta_resp=100, d_resp=1, mem_col=1, beta_col=100,  mem_pos=1, beta_pos=100, scale_resp=1, 
             u_mem_resp=1, u_delta_resp=1, scale_col=1, d_col=1, u_mem_col=1, u_delta_col=1, 
             scale_pos=1, u_mem_pos=1, u_delta_pos=1, u_delta_pos_dist=1, ter=min(sub$rt-0.01),
             m_ndt_resp = 1, delta_ndt_resp = 0.1, m_ndt_col = 1, delta_ndt_col = 0.1, 
             m_ndt_pos = 1, delta_ndt_pos = 0.1, delta_ndt_pos_dist = 0.1)
  uppers = uppers[is.na(fixed[1:26])]
  
  par <- optim(par1,  findParameters, data = sub, fixed = fixed,
               lower = lowers, upper = uppers, method = 'L-BFGS-B')
  #  proc.time() - t0
  # get all parameters
  fixed_par = fixed
  fixed_par[is.na(fixed_par)] = par$par # all fixed
  parameters <- findParameters(fixed_par, sub, fixed_par, op = FALSE)
  N <- nrow(sub)
  npar <- length(par1) + 3*length(unique(sub$targets))
  nll <- parameters[length(parameters)]
  parameters <- c(parameters, AIC=2*nll + 2*npar)
  parameters <- c(parameters, BIC=2*nll + log(N)*npar)
  parameters <- c(parameters, sub = sub$sub[1], model = modelname)
  return(parameters)
}

#' Parallel compute parameters
#' 
#' Parallel compute parameters for each subject each model
#' @param lsubs list of subject data, which should include 'targets' and 'inttrial'
#' @param modelnames a list of modelnames (see optimization function)
#' @return paras return a list of parameters
parEstimate <- function(lsubs, modelnames){
  no_cores <- min(detectCores() - 1, 10)
  print(no_cores)
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(dplyr) })
  clusterExport(cl, c('findParameters','findInnerParameters', 'updateRate','updateRateMem', 'updateRateWeight',
                      'fitPar','drecinorm','updatePriorA','updatePrior','logll','log_ddiffusion','updateS0',
                      'updatePosPrior',  'updateRateWeightPosition',  'updateRateWeightPositionDist', 
                      "updatePosS0", "updateNDTWeight", "updateNDTPosition", "updateNDTPositionDist", 
                      "updatePosNDT","updatePosRateWeight", "updateRateWeightPositionMatch",
                      "updateNDTPositionMatch", "updatePosPriorSpread"))
  
  t0 = proc.time()
  paras <- clusterMap(cl, fitPar, lsubs, modelnames)
  stopCluster(cl)
  print(proc.time()-t0)
  return(paras)
}
