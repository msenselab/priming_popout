library(tidyverse)
library(data.table)
library(cowplot) 

source('bayesian_updates.R')
source('figure_functions.R')


exp_pars_col <- readRDS('./data/fitpars_opt_r2_col.rds')
exp_pars_col$type <- "c"
exp_pars_resp <- readRDS('./data/fitpars_opt_r2_resp.rds')
exp_pars_resp$type <- "r"
exp_pars_resp$d_col <- 0
exp_pars_pos <- readRDS('./data/fitpars_opt_pos.rds')
exp_pars_pos$type <- "p"
exp_pars_pos$d_col <- 0
exp_pars <- rbind(exp_pars_col, exp_pars_resp, exp_pars_pos)

exp_pars$later_diffusion <- factor(exp_pars$later_diffusion, levels=0:1, labels=c("DDM", "LATER"))
exp_pars$resp_update <- factor(exp_pars$resp_update, levels=0:14, labels=c("No update", "PI approx. Bayesian S0", "PI Bayesian S0", "PD Bayesian S0", 
                                                                           "PG Bayesian S0", "PI binary rate", "PI step rate", 
                                                                           "PI weighted rate", "PD approx. Bayesian S0", "PG approx. Bayesian S0",
                                                                           "PI binary NDT", "PI weighted NDT", "PD step NDT", "PG step NDT", 
                                                                           "PS Bayesian S0"))
exp_pars$col_update <- factor(exp_pars$col_update, levels=0:7, labels=c("No update", "PI binary rate", "PI step rate", "PI weighted rate", "PI binary NDT", 
                                                                        "PI weighted NDT", "PD weighted rate", "PG weighted rate"))
exp_pars$pos_update <- factor(exp_pars$pos_update, levels=0:9, labels=c("No update", "Binary rate", "Step rate", "Weighted rate", "Weighted rate with DI",
                                                                        "Binary NDT", "Weighted NDT", "Weighted NDT with DI", "Matched weighted rate",
                                                                        "Matched weighted NDT"))

exp_pars <- unite(exp_pars, sub_session, sub, session)

# Create lists of models, sorted by total AIC
group_by(exp_pars, model) %>% summarize(mAIC=mean(as.numeric(AIC)), mBIC=mean(as.numeric(BIC))) %>% 
  arrange(mAIC) -> model_list

# Create lists of best and worst models based on total AIC
model_list[1:5,] -> top5_models
model_list[(nrow(model_list)-4):nrow(model_list),] -> bottom5_models
top_model <- top5_models$model[1]
filter(exp_pars, model==top_model, type=="c") %>% select(sub_session, AIC) -> best_models
names(best_models)[2] <- "topAIC"
exp_pars <- full_join(exp_pars, best_models) %>% mutate(rAIC=AIC-topAIC)

model_params <- dplyr::filter(exp_pars, model==top_model, type=="c")
model_params$resp_update <- as.numeric(model_params$resp_update)-1
model_params$col_update <- as.numeric(model_params$col_update)-1
model_params$pos_update <- as.numeric(model_params$pos_update)-1

# Reorder factor levels
exp_pars$resp_update <- factor(exp_pars$resp_update, levels(exp_pars$resp_update)[c(1:3, 6:8, 11:12, 9:10, 4:5, 15, 13:14)])
exp_pars$pos_update <- factor(exp_pars$pos_update, levels(exp_pars$pos_update)[c(1:5, 9, 6:8, 10)])

# Top lists based on number of participants for which a model is the best fitting one in terms of AIC
topmod = list()

subs <- unique(exp_pars$sub_session)

N <- length(unique(subs))

for (i in 1:28) { dplyr::filter(exp_pars, sub_session==subs[i]) %>% arrange(AIC) %>% dplyr::filter(AIC==min(AIC)) %>% 
    dplyr::select(model) -> topmod[i] }

# Find "second best" models for position-based updating
topmod = list()
stopmod = list()
best_pos_mods <- data.frame(sub=numeric(0), topmod=numeric(0), second=numeric(0))
for (i in 1:28) { dplyr::filter(exp_pars, sub_session==subs[i]) %>% arrange(AIC) -> topmods
  topmod[i] <- topmods %>% filter(AIC==min(AIC)) %>% select(pos_update)
  topmod[i] <- as.numeric(topmod[[i]][1]) - 1
  j <- 1 
  while(as.numeric(topmods[j,]$pos_update)-1 == topmod[i]) {
    j <- j + 1
  }
  stopmod[i] <- as.numeric(topmods[j,]$pos_update)-1
  best_pos_mods <- rbind(best_pos_mods, data.frame(sub=i, topmod=topmod[[i]], second=stopmod[[i]]))
}

# Look at which proportion of second best models are of the same type 
prop_rate_rate <- filter(best_pos_mods, topmod<=5) %>% summarize(mean(second<=5)) 
prop_NDT_NDT <- filter(best_pos_mods, topmod>5) %>% summarize(mean(second >5))

# Approximate Bayesian rules
approx_bayes <- c('PI approx. Bayesian S0', 'PD approx. Bayesian S0', 'PG approx. Bayesian S0')

cols <- c("S0" = "skyblue2", "Rate" = "orange", "NDT" = "seagreen3")
bgcol <- data.frame(xmin=c(1.5,2.5,5.5,7.5,10.5),xmax=c(2.5,5.5,7.5,10.5,12.5), 
                    Type=c("S0", "Rate", "NDT", "S0", "NDT"), ymin=-5, ymax=70)
fig_resp_aic <- filter(exp_pars, type=="r", !(resp_update %in% approx_bayes)) %>%
  separate(sub_session, c("sub", "session")) %>%
  group_by(resp_update, later_diffusion, sub) %>% summarize(rAIC=mean(rAIC)) %>%
  summarize(mAIC = mean(rAIC), se_AIC=sd(rAIC)/sqrt(N-1)) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=later_diffusion, color=later_diffusion, shape=later_diffusion)) + theme_bw() +
  geom_errorbar(aes(x = resp_update, ymin=mAIC-se_AIC, ymax=mAIC+se_AIC, color=later_diffusion), width=0.2) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
  labs(x="RCF based updating", y="average relative AIC", color="EA model", shape="EA model", fill="Updating variable") + 
  geom_vline(aes(xintercept=7.5), linetype=2) +
  geom_rect(data=bgcol, inherit.aes=FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Type), color="transparent",  alpha=0.3) +
  coord_cartesian(ylim = c(-1, 66)) + scale_fill_manual(values=cols)

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/resp_aic.png', fig_resp_aic, width = 7.5, height = 3.5)
}

bgcol <- data.frame(xmin=c(1.5,4.5,6.5),xmax=c(4.5,6.5,8.5), 
                    Type=c("Rate", "NDT", "Rate"), ymin=-5, ymax=70)
fig_col_aic <- filter(exp_pars, type=="c") %>% separate(sub_session, c("sub", "session")) %>%
  group_by(col_update, later_diffusion, sub) %>% summarize(rAIC=mean(rAIC)) %>%
  summarize(mAIC = mean(rAIC), se_AIC=sd(rAIC)/sqrt(N-1)) %>% 
  ggplot(aes(x = col_update, y=mAIC, group=later_diffusion, color=later_diffusion, shape=later_diffusion)) + theme_bw() +
  geom_errorbar(aes(x = col_update, ymin=mAIC-se_AIC, ymax=mAIC+se_AIC, color=later_diffusion), width=0.2) + 
  geom_line() + geom_vline(aes(xintercept=6.5), linetype=2) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
  labs(x="Color based updating", y="average relative AIC", color="EA model", shape="EA model", fill="Updating variable") +
  geom_rect(data=bgcol, inherit.aes=FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Type), color="transparent",  alpha=0.3) +
  coord_cartesian(ylim = c(-1, 65)) + scale_fill_manual(values=cols)

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/col_aic.png', fig_col_aic, width = 7, height = 3.5)
}

bgcol <- data.frame(xmin=c(1.5,6.5),xmax=c(6.5,10.5), 
                    Type=c("Rate", "NDT"), ymin=-5, ymax=55)
fig_pos_aic <- filter(exp_pars, type=="p") %>% separate(sub_session, c("sub", "session")) %>%
  group_by(pos_update, later_diffusion, sub) %>% summarize(rAIC=mean(rAIC)) %>%
  summarize(mAIC = mean(rAIC), se_AIC=sd(rAIC)/sqrt(N-1)) %>% 
  ggplot(aes(x = pos_update, y=mAIC, group=later_diffusion, color=later_diffusion, shape=later_diffusion)) + theme_bw() +
  geom_errorbar(aes(x = pos_update, ymin=mAIC-se_AIC, ymax=mAIC+se_AIC, color=later_diffusion), width=0.2) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
  labs(x="Position based updating", y="average relative AIC", color="EA model", shape="EA model", fill="Updating variable") +
  geom_rect(data=bgcol, inherit.aes=FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Type), color="transparent",  alpha=0.3) +
  coord_cartesian(ylim = c(-1, 50)) + scale_fill_manual(values=cols)


if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/pos_aic.png', fig_pos_aic, width = 7.5, height = 4)
}

subfiles=dir('data', '*.dat')
data <- data.table(do.call('rbind',lapply(subfiles, readData)))
data <- unite(data, "sub_session", sub, session)

data_pred <- data.frame(sub_session=character(0))
for(s in subs) {
  dp <- genSeq(dplyr::filter(model_params, sub_session==s), dplyr::filter(data, sub_session==s))
  dp$tno <- 1:nrow(dp)
  data_pred <- full_join(data_pred, dp)
}

# Plot examples of starting point and rate updating

sample_sub_session <- "01_A"
sample_sub <- substr(sample_sub_session, 1,2)
sample_session <- substr(sample_sub_session, 4,4)

sample <- filter(data_pred, sub_session==sample_sub_session)[11:110,] %>% mutate(s0=s0 * ((tori=="top") * 2 - 1)) %>% 
  mutate(tar_start = tno, tar_end = tno+1)

p_s0 <- sample %>% ggplot(aes(x=tno, y=s0))  + geom_line() + theme_bw() +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x = NULL, y = NULL, ymin = -0.65, ymax = 0.65, fill = tori), alpha=0.5) + 
  xlab('Trial no.') + theme(legend.position="bottom")

rate <- dplyr::filter(exp_pars, sub_session==sample_sub_session, model==top_model)$mu 

p_rate <- sample %>% ggplot(aes(x=tno, y=rate)) + geom_line() + theme_bw()  +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x =NULL, y = NULL, ymin = 16, ymax = 21, fill = tcol), alpha=0.5) + 
  xlab('Trial no.') + theme(legend.position="bottom")


col_priming <- function(data) {
  col_prime <- data.frame()
  for(i in 1:8) {
    cp <- group_by(data, sub_session) %>% mutate(cn=lag(tcol,i), en=lag(error,i), rn=if_else(tcol==cn, "Repeat", "Switch")) %>% 
      filter(error=="no", !is.na(rn)) %>% group_by(rn) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    cp$d=-i
    col_prime <- rbind(col_prime, cp)
  }
  col_prime
}

ss <- unique(data_pred$sub_session)
N <- length(ss)/2

col_prime <- data.frame()
for(s in ss) {
  cp <- col_priming(filter(data_pred, sub_session==s))
  cp$sub_session <- s
  col_prime <- rbind(col_prime, cp)
}

isub_cprime_plot <- ggplot(col_prime, aes(x=d, y=mRT, color=rn)) +
  geom_point(size=3) + geom_line(aes(y=mpRT, group=rn)) + theme_bw() + facet_wrap(~ sub_session)

sdata_col <- col_prime%>% separate(sub_session, c("sub", "session"))%>%  group_by(sub, session) %>%
  mutate(nRT=(mRT-mean(mRT))*1000, npRT=(mpRT-mean(mRT))*1000) %>%
  group_by(d, rn, sub) %>% summarize(nRT=mean(nRT), npRT=mean(npRT)) %>%
  summarize(mmRT=mean(nRT), mmpRT=mean(npRT), seRT=sd(nRT)/sqrt(N-1))

cprime_plot <- ggplot(sdata_col, aes(x=d, y=mmRT, color=rn)) + geom_point(size=3) + geom_line(aes(y=mmpRT, group=rn)) +
  geom_errorbar(aes(ymax=mmRT + seRT*1.96, ymin=mmRT-seRT*1.96), width=0.2) + theme_bw() +
  labs(x="Lag", y="Normalized RT (ms)", color="Color") + 
  theme(strip.background = element_blank(), strip.text=element_blank(), text = element_text(size=16)) 

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/cprime.png', cprime_plot, width = 8, height = 6)
}

pos_priming <- function(data) {
  pos_prime <- data.frame()
  for(i in 1:8) {
    pp <- group_by(data, sub_session) %>% mutate(pn=lag(tpos,i), d1n=lag(d1pos,i), d2n=lag(d2pos,i), d3n=lag(d3pos,i), 
                                                 en=lag(error,i), rn=if_else(tpos==pn, "TT", 
                                                                             if_else(tpos==d1n, 
                                                                                     "TD", if_else(tpos==d2n, "TD", 
                                                                                                   if_else(tpos==d3n, 
                                                                                                           "TD", "TN"))))) %>% 
      filter(error=="no", !is.na(rn)) %>% group_by(rn) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    pp$d=-i
    pos_prime <- rbind(pos_prime, pp)
  }
  pos_prime
}

pos_prime <- data.frame()
for(s in ss) {
  pp <- pos_priming(filter(data_pred, sub_session==s))
  pp$sub_session <- s
  pos_prime <- rbind(pos_prime, pp)
}

isub_pprime_plot <- ggplot(pos_prime, aes(x=d, y=mRT, color=rn))  +
  geom_point(size=3) + geom_line(aes(y=mpRT, group=rn)) + theme_bw() +
  facet_wrap(~ sub_session)

sdata_pos <- pos_prime %>% separate(sub_session, c("sub", "session")) %>% group_by(sub, session) %>% 
  mutate(nRT=(mRT-mean(mRT))*1000, npRT=(mpRT-mean(mRT))*1000) %>%
  group_by(d, rn, sub) %>% summarize(nRT=mean(nRT), npRT=mean(npRT)) %>%
  summarize(mmRT=mean(nRT), mmpRT=mean(npRT), seRT=sd(nRT)/sqrt(N-1))

pprime_plot <- ggplot(sdata_pos, aes(x=d, y=mmRT, color=rn)) + geom_point(size=3) + geom_line(aes(y=mmpRT, group=rn)) +
  geom_errorbar(aes(ymax=mmRT + seRT*1.96, ymin=mmRT - seRT*1.96), width=0.2) + theme_bw() +
  labs(x="Lag", y="Normalized RT (ms)", color="Positional inter-trial condition") + 
  theme(strip.background = element_blank(), strip.text=element_blank(), 
        legend.position="bottom", text = element_text(size=16)) 

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/pprime.png', pprime_plot, width = 8, height = 6.5)
}

resp_priming <- function(data) {
  resp_prime <- data.frame()
  for(i in 1:8) {
    rp <- group_by(data, sub_session) %>% mutate(on=lag(tori,i), en=lag(error,i), rn=if_else(tori==on, "Repeat", "Switch")) %>% 
      filter(error=="no", !is.na(rn)) %>% group_by(rn) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    rp$d=-i
    resp_prime <- rbind(resp_prime, rp)
  }
  resp_prime
}

pos_resp_priming <- function(data) {
  pos_resp_prime <- data.frame()
  for(i in 1:8) {
    prp <- group_by(data, sub_session) %>% mutate(on=lag(tori,i), pn=lag(tpos,i), d1n=lag(d1pos,i), d2n=lag(d2pos,i), 
                                                  d3n=lag(d3pos,i), en=lag(error,i), rpn=if_else(tpos==pn, "TT", 
                                                                                                 if_else(tpos==d1n, 
                                                                                                         "TD", if_else(tpos==d2n, "TD", 
                                                                                                                       if_else(tpos==d3n, 
                                                                                                                               "TD", "TN")))),
                                                  ron=if_else(tori==on, "Repeat", "Switch")) %>% 
      filter(error=="no", !is.na(rpn)) %>% group_by(rpn, ron) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    prp$d=-i
    pos_resp_prime <- rbind(pos_resp_prime, data.frame(prp))
  }
  pos_resp_prime
}

resp_prime <- data.frame()
for(s in ss) {
  rp <- resp_priming(filter(data_pred, sub_session==s))
  rp$sub_session <- s
  resp_prime <- rbind(resp_prime, rp)
}

isub_rprime_plot <- ggplot(resp_prime, aes(x=d, y=mRT, color=rn)) + geom_point(size=3) + 
  geom_line(aes(y=mpRT, group=rn)) + theme_bw() + facet_wrap(~ sub_session)

sdata_resp <- resp_prime %>% separate(sub_session, c("sub", "session"))  %>% group_by(sub, session) %>%
  mutate(nRT=(mRT-mean(mRT))*1000, npRT=(mpRT-mean(mRT))*1000) %>%
  group_by(d, rn, sub) %>% summarize(nRT=mean(nRT), npRT=mean(npRT)) %>%  
  summarize(mmRT=mean(nRT), mmpRT=mean(npRT), seRT=sd(nRT)/sqrt(N-1))

rprime_plot <- ggplot(sdata_resp, aes(x=d, y=mmRT, color=rn)) + geom_point(size=3) + geom_line(aes(y=mmpRT, group=rn)) +
  geom_errorbar(aes(ymax=mmRT + seRT*1.96, ymin=mmRT-seRT*1.96), width=0.4) + theme_bw() +
  labs(x="Lag", y="Normalized RT (ms)", color="Target notch position (RCF)") +
  theme(strip.background = element_blank(), strip.text=element_blank(), 
        legend.position="bottom", text = element_text(size=16)) 

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/rprime.png', rprime_plot, width = 8, height = 6.5)
}


pos_resp_prime <- data.frame()
for(s in ss) {
  prp <- pos_resp_priming(filter(data_pred, sub_session==s))
  prp$sub_session <- s
  pos_resp_prime <- rbind(pos_resp_prime, prp)
}

sdata_pos_resp <- pos_resp_prime %>% separate(sub_session, c("sub", "session")) %>% group_by(sub, session) %>%
  mutate(nRT=(mRT-mean(mRT))*1000, npRT=(mpRT-mean(mRT))*1000) %>%
  group_by(d, rpn, ron, sub) %>% summarize(nRT=mean(nRT), npRT=mean(npRT)) %>%  
  summarize(mmRT=mean(nRT), mmpRT=mean(npRT), seRT=sd(nRT)/sqrt(N-1))

prprime_plot <- ggplot(sdata_pos_resp, aes(x=d, y=mmRT, color=ron)) + geom_point(size=3) + 
  geom_line(aes(y=mmpRT, group=ron)) + geom_errorbar(aes(ymax=mmRT + seRT*1.96, 
                                                         ymin=mmRT-seRT*1.96), width=0.4) + 
  theme_bw() + labs(x="Lag", y="Normalized RT (ms)", color="Target notch position (RCF)") +
  theme(strip.background = element_blank(), legend.position="bottom", 
        text = element_text(size=16)) + facet_wrap(~ rpn)

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/prprime.png', prprime_plot, width = 8, height = 4)
}


# ---- Updating_examples  ----

# Example of color-based updating
ex_sub <- "11_A"
l <- updateColDemo(filter(data_pred, sub_session==ex_sub, tno %in% 1:8), 
                   filter(exp_pars, sub_session==ex_sub, model==top_model, type=="p"))

fig_col_update <- l[[2]]

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/col_update.png', fig_col_update, width = 8, height = 5)
}

# Example of position-based updating

l <- updatePosDemo(filter(data_pred, sub_session==ex_sub, tno %in% 1:8), 
                   filter(exp_pars, sub_session==ex_sub, model==top_model, type=="p"))

fig_pos_update <- l[[2]]

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/pos_update.png', fig_pos_update, width = 8, height = 8)
}

# Example of response-based updating

l <- updateRespDemo(filter(data_pred, sub_session==ex_sub, tno %in% 1:8), 
                   filter(exp_pars, sub_session==ex_sub, model==top_model, type=="p"))

fig_resp_update <- l[[2]]

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/resp_update.png', fig_resp_update, width = 8, height = 8)
}

