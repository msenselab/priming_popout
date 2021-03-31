library(tidyverse)
library(data.table)
library(cowplot) 

updatePosDemo <- function(data, param) {
  
  tarpos <- data$tpos
  distpos <- t(matrix(c(data$d1pos, data$d2pos, data$d3pos),nrow(data),3))
  delta_tar <- param$u_delta_pos
  delta_dist <- param$u_delta_pos_dist
  mem <- param$u_mem_pos
  
  ds = matrix(rep(0, length(tarpos)*8), length(tarpos), 8)
  for (i in 2:length(tarpos)) {
    ds[i, tarpos[i-1]] <- ds[i-1, tarpos[i-1]] + delta_tar
    ds[i, -tarpos[i-1]] <- ds[i-1, -tarpos[i-1]] - delta_tar/7
    ds[i, distpos[,i-1]] <- ds[i, distpos[,i-1]] - delta_dist/3
    ds[i, -distpos[,i-1]] <- ds[i, -distpos[,i-1]] + delta_dist/5
    ds[i, ] <- ds[i, ] * mem
  }
  
  ds <- data.frame(ds)
  ds$tno <- 1:nrow(ds)
  ds$tpos <- data$tpos
  ds$d1pos <- data$d1pos
  ds$d2pos <- data$d2pos
  ds$d3pos <- data$d3pos
  
  ds <- gather(ds,pos,weight,-tno,-tpos,-d1pos,-d2pos,-d3pos)
  
  ds$weight <- ds$weight + 1 
  
  ds$pos <- factor(ds$pos, labels=c("1", "2", "3", "4", "5", "6", "7", "8"))
  
  ds %>% ggplot(aes(x=tno, y=pos, fill=weight)) + geom_tile() + geom_point(aes(y=tpos), shape=84, size=3) + 
    geom_point(aes(y=d1pos), shape=68, size=3) + geom_point(aes(y=d2pos), shape=68, size=3) + 
    geom_point(aes(y=d3pos), shape=68, size=3) + theme_classic() + labs(x="Trial", y="Position", 
                                                                        fill="Weight") +  
    scale_fill_gradient(low = "white", high = "gray25")  + theme(text = element_text(size=15)) -> fig_hmap
  
  fig_line <- filter(ds, tpos==pos) %>% ggplot(aes(x=tno, y=weight*param$mu, group=1)) + 
    geom_point() + geom_line()+ labs(x="Trial", y="Rate")  + theme_classic() +
    geom_hline(aes(yintercept=param$mu), linetype=2) + theme(text = element_text(size=15)) 
  
  fig_tot <- plot_grid(fig_hmap, fig_line, nrow=2, labels=c("A", "B"), rel_heights = c(2,1))
  
  return(list(ds, fig_tot))
}

updateColDemo <- function(data, param) {
  
  dim <- as.numeric(data$tcol)
  mem <- param$u_mem_col
  delta <- param$u_delta_col
  
  dim[dim>0] <- 2*(dim[dim>0]-1.5) # rescale intertrial so that color orientation and absent are -1, 1 and 0
  ds = rep(0, length(dim))
  for (i in 2:length(dim)) {
    ds[i] = ds[i-1] + dim[i-1]*delta
    ds[i] = ds[i]*mem 
  }
  
  N <- length(ds)
  ds <- matrix(1 + c(ds, -ds), N, 2)
  ds <- data.frame(ds)
  ds$tno <- 1:N
  ds$tcol <- data$tcol
  ds$dcol <- factor(data$tcol, levels=c("green", "red"), labels=c("red", "green"))
  ds <- gather(ds,col,weight,-tno,-tcol,-dcol)
  ds$col <- factor(ds$col, labels=c("green", "red"))
  
  
  ds %>% ggplot(aes(x=tno, y=col, fill=weight)) + geom_tile() + geom_point(aes(y=tcol), shape=84, size=3) + 
    geom_point(aes(y=dcol), shape=68, size=3) + theme_classic() + labs(x="Trial", y="Color", fill="Weight") +  
    scale_fill_gradient(low = "white", high = "gray25") + theme(text = element_text(size=15)) -> fig_hmap
  
  fig_line <- filter(ds, tcol==col) %>% ggplot(aes(x=tno, y=weight*param$mu, group=1)) + 
    geom_point() + geom_line() + labs(x="Trial", y="Rate") + theme_classic()  +
    geom_hline(aes(yintercept=param$mu), linetype=2) + theme(text = element_text(size=15))
  
  fig_tot <- plot_grid(fig_hmap, fig_line, nrow=2, labels=c("A", "B"), rel_heights = c(1,1))
  
  return(list(ds, fig_tot))
}

updateRespDemo <- function(data, param) {
  
  targets <- as.numeric(data$tori=="top") 
  tpos <- data$tpos
  m <- param$mem_resp
  a <- param$beta_resp
  b <- param$beta_resp
  d <- param$d_resp
  
  x <- seq(0,1,0.001)
  p0 <- dbeta(x,a,b)
  p0 <- matrix(rep(p0,each=8),8,length(x))
  pd <- p0
  mu = matrix(rep(a/(a+b), (length(targets)+1)*8), length(targets)+1, 8)
  idx = 2
  
  for(target in targets) {
    # Update the distribution based on the Bernoulli likelihood 
    pd[tpos[idx-1],] <- pd[tpos[idx-1],]*x*target + (1-x)*pd[tpos[idx-1],]*(1-target)
    
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
    
    # Implement 'forgetting'
    pd <- (1-m)*p0 + m * pd
    # Normalize the probability distribution
    pd <- pd/rowMeans(pd)
    for(i in 1:8) {
      mu[idx,i] <- mean(pd[i,]*x)
    }
    idx  = idx + 1
  }
  mu <- mu[1:(idx-2),]
  N <- dim(mu)[1]
  mu <- data.frame(mu)
  names(mu) <- 1:8
  mu$tno <- 1:N
  mu$tpos <- tpos
  mu$tori <- factor(targets, levels=c(0,1), labels=c("bottom", "top"))
  mu <- gather(mu,pos,S0,-tno,-tpos,-tori) %>% mutate(S0=log(S0/(1-S0)))
  
  fig_hmap <- mu %>% ggplot(aes(x=tno, y=pos, fill=S0)) + geom_tile() + 
    geom_point(aes(y=tpos, shape=tori)) + theme_classic() +
    labs(x="Trial", y="Position", fill="S0",shape="RCF") +  
    scale_fill_gradient(low = "white", high = "gray25") + 
    theme(text = element_text(size=15)) + scale_shape_manual(values=c(6,2)) 
  
  fig_line <- filter(mu, tpos==pos) %>% ggplot(aes(x=tno, y=S0, group=1)) + 
    geom_point() + geom_line() + labs(x="Trial", y="Starting point") +
    theme(text = element_text(size=15)) + theme_classic() +
    geom_hline(aes(yintercept=0), linetype=2)
  
  fig_tot <- plot_grid(fig_hmap, fig_line, nrow=2, labels=c("A", "B"), rel_heights = c(2,1))
  
  return(list(mu, fig_tot))
}


