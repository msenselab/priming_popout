source('bayesian_updates.R')

# ---- Functions ----

# Read the data and add column names, factor labels, etc.
# Also scale RTs to seconds rather than ms for compatibility with existing code from previous project.
# Finally duplicate data and create folds, each with a different "test" block
readData_crossval <- function (filename) {
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
  d$error<-factor(d$error, levels=c(0,1,2,3), labels=c("no", "yes", "timeout", "test"))
  d$sub<-factor(d$sub)
  d$rt<-d$rt/1000 
  
  d <- d %>% mutate(block=rep(1:(n()/112), each=112))
  d <- do.call("rbind", replicate(max(d$block), d, simplify = FALSE))
  d <- d %>% mutate(fold=rep(1:max(block), each=n()/max(block)), 
                    test = ifelse(block==fold, 1, 0))

  d
}

# Calculate color priming effects at different lags
col_priming <- function(data) {
  col_prime <- data.frame()
  for(i in 1:8) {
    cp <- data %>% mutate(cn=lag(tcol,i), en=lag(error,i), rn=if_else(tcol==cn, "Repeat", "Switch")) %>% 
      filter(error == "test", !is.na(rn)) %>% group_by(rn) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    cp$d=-i
    col_prime <- rbind(col_prime, cp)
  }
  col_prime
}

pos_priming <- function(data) {
  pos_prime <- data.frame()
  for(i in 1:8) {
    pp <- data %>% mutate(pn=lag(tpos,i), d1n=lag(d1pos,i), d2n=lag(d2pos,i), d3n=lag(d3pos,i), 
                                                 en=lag(error,i), rn=if_else(tpos==pn, "TT", 
                                                                             if_else(tpos==d1n, 
                                                                                     "TD", if_else(tpos==d2n, "TD", 
                                                                                                   if_else(tpos==d3n, 
                                                                                                           "TD", "TN"))))) %>% 
      filter(error == "test", !is.na(rn)) %>% group_by(rn) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    pp$d=-i
    pos_prime <- rbind(pos_prime, pp)
  }
  pos_prime
}

pos_resp_priming <- function(data) {
  pos_resp_prime <- data.frame()
  for(i in 1:8) {
    prp <- data %>% mutate(on=lag(tori,i), pn=lag(tpos,i), d1n=lag(d1pos,i), d2n=lag(d2pos,i), 
                                                  d3n=lag(d3pos,i), en=lag(error,i), rpn=if_else(tpos==pn, "TT", 
                                                                                                 if_else(tpos==d1n, 
                                                                                                         "TD", if_else(tpos==d2n, "TD", 
                                                                                                                       if_else(tpos==d3n, 
                                                                                                                               "TD", "TN")))),
                                                  ron=if_else(tori==on, "Repeat", "Switch")) %>% 
      filter(error == "test", !is.na(rpn)) %>% group_by(rpn, ron) %>% summarize(mRT=mean(rt), mpRT=mean(predRT))
    prp$d=-i
    pos_resp_prime <- rbind(pos_resp_prime, data.frame(prp))
  }
  pos_resp_prime
}


# ---- Read data ----

modelnames <- c("LR0C0P0", "LR4C3P4")

subfiles=dir('data', '....dat')
data <- data.table(do.call('rbind',lapply(subfiles, readData_crossval)))
data <- data %>% mutate(sub = paste0(sub, '_', session, '_', fold))

# Mark test blocks for exclusion
data$error[data$test==1 & data$error=="no"] = "test"

N <- length(unique(data$sub))

lsubs <- split(data, data$sub) # split data for each subject and session

# ---- Fit models ----

# g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
# 
# lpars <- parEstimate(g$dat, g$mod)
# 
# # now convert lists to data.frame
# pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
# names(pars.df)[35:36] <- c("AIC", "BIC")
# pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
# pars.df <- separate(pars.df, sub, c("sub", "session", "fold"))
# pars.df <- pars.df %>% mutate(sub_session=paste0(sub, session))
# 
# saveRDS(pars.df, file = './data/cross_val.rds')

# Read results of fitting models from a file
pars.df <- readRDS('./data/cross_val.rds') 

# ---- Read model parameters from previous non-crossvalidation model fits ----

# Read baseline no updating model params
exp_pars_ctrl <- readRDS('./data/fitpars_col.rds') %>% filter(model==modelnames[1])

# Read parameters of best model
exp_pars_top <- readRDS('./data/fitpars_opt_r2_col.rds')  %>% filter(model==modelnames[2])

# ---- Generate model predictions ----

ss <- unique(pars.df$sub_session)

cv_par <- data.frame(sub_session=character(0), base_nll = numeric(0), base_AIC = numeric(0),
                     test_nll = numeric(0), test_AIC = numeric(0))
data_pred_full <- data.frame(sub=character(0))
for(s in ss) {
  # Control model
  
  modpar <- filter(pars.df, sub_session==s, model==modelnames[1])
  modpar$later_diffusion <- factor(modpar$later_diffusion, levels=0:1, labels=c("DDM", "LATER"))
  
  data_pred_ctrl <- data.frame(sub=character(0))
  data_pred_top <- data.frame(sub=character(0))
  
  folds <- modpar$fold
  for(f in folds) {
    sub_session_fold <- paste0(substr(s,1,2), '_', substr(s,3,3), '_', f)
    dp <- genSeq(filter(modpar, fold==f), filter(data, sub==sub_session_fold))
    dp$tno <- 1:nrow(dp)
    data_pred_ctrl <- full_join(data_pred_ctrl, dp)
  } 
  data_pred_ctrl <- mutate(data_pred_ctrl, delta = theta - s0, include=ifelse(error=="test",1,0))
  data_pred_ctrl$model = "Control"
  dat <- filter(data_pred_ctrl,  include==1)
  train_dat <- filter(data_pred_ctrl, error=="no")
  train_nll <- mutate(train_dat, nll = -log(drecinorm(rt-ndt, rate/delta, sigma/delta))) %>%
    group_by(block, fold) %>% summarize(nll=sum(nll)) %>% summarize(nll=mean(nll)) %>% 
    summarize(nll=sum(nll))
  cv_par <- full_join(cv_par,
                      data.frame(sub_session = s, model="Control", 
                                 base_nll=filter(exp_pars_ctrl, sub==substr(s, 1,2), 
                                                 session==substr(s,3,3))$nll,
                                 base_AIC = filter(exp_pars_ctrl, sub==substr(s, 1,2),
                                                   session==substr(s,3,3))$AIC,
                                 train_nll = train_nll$nll,
                                 test_nll = -sum(log(drecinorm(dat$rt-dat$ndt, dat$rate/dat$delta, dat$sig/dat$delta))),
                                 test_AIC = -sum(log(drecinorm(dat$rt-dat$ndt, dat$rate/dat$delta, dat$sig/dat$delta)))*2 + 2*4))
  
  # Best model
  
  modpar=filter(pars.df, sub_session==s, model==modelnames[2])
  modpar$later_diffusion <- factor(modpar$later_diffusion, levels=0:1, labels=c("DDM", "LATER"))
  
  folds <- modpar$fold
  for(f in folds) {
    sub_session_fold <- paste0(substr(s,1,2), '_', substr(s,3,3), '_', f)
    dp <- genSeq(filter(modpar, fold==f), filter(data, sub==sub_session_fold))
    dp$tno <- 1:nrow(dp)
    data_pred_top <- full_join(data_pred_top, dp)
  } 
  data_pred_top <- mutate(data_pred_top, delta = theta - s0, include=ifelse(error=="test",1,0))
  data_pred_top$model = "Top"
  dat <- filter(data_pred_top, include==1)
  train_dat <- filter(data_pred_top, error=="no")
  train_nll <- mutate(train_dat, nll = -log(drecinorm(rt-ndt, rate/delta, sigma/delta))) %>%
    group_by(block, fold) %>% summarize(nll=sum(nll)) %>% summarize(nll=mean(nll)) %>% 
    summarize(nll=sum(nll))
  cv_par <- full_join(cv_par, 
                      data.frame(sub_session = s,  model="Top", 
                                 base_nll = filter(exp_pars_top, sub==substr(s, 1,2),
                                                   session==substr(s,3,3))$nll,
                                 base_AIC = filter(exp_pars_top, sub==substr(s, 1,2),
                                                   session==substr(s,3,3))$AIC,
                                 train_nll = train_nll$nll,
                                 test_nll = -sum(log(drecinorm(dat$rt-dat$ndt, dat$rate/dat$delta, dat$sig/dat$delta))),
                                 test_AIC = -sum(log(drecinorm(dat$rt-dat$ndt, dat$rate/dat$delta, dat$sig/dat$delta)))*2 + 2*12))
  
  data_pred_full <- full_join(data_pred_full, data_pred_ctrl)
  data_pred_full <- full_join(data_pred_full, data_pred_top)
  
}

s_cv_par <- cv_par %>% group_by(model) %>% summarize(m_base_ll = mean(-base_nll),
                                                     m_train_ll = mean(-train_nll),
                                                     m_test_ll = mean(-test_nll)) %>% 
  mutate(diff=m_test_ll-m_train_ll)
s_cv_par 

# ---- Plot cross-validated model predictions ----

col_prime <- data.frame()
pos_prime <- data.frame()
pos_resp_prime <- data.frame()
for(s in ss) {
  dat <- filter(data_pred_full, block == fold, substr(sub, 1,2)==substr(s,1,2), 
                session==substr(s, 3,3), model=="Top")
  
  cp <- col_priming(dat)
  cp$sub_session <- paste0(substr(s,1,2), '_', substr(s,3,3))
  col_prime <- rbind(col_prime, cp)
  
  pp <- pos_priming(dat)
  pp$sub_session <- paste0(substr(s,1,2), '_', substr(s,3,3))
  pos_prime <- rbind(pos_prime, pp)
  
  prp <- pos_resp_priming(dat)
  prp$sub_session <- paste0(substr(s,1,2), '_', substr(s,3,3))
  pos_resp_prime <- rbind(pos_resp_prime, prp)
}

# Divide by two to get number of participants because two sessions/participant
N <- length(ss)/2

sdata_col <- col_prime %>% separate(sub_session, c("sub", "session"))%>%  group_by(sub, session) %>%
  mutate(nRT=(mRT-mean(mRT))*1000, npRT=(mpRT-mean(mRT))*1000) %>%
  group_by(d, rn, sub) %>% summarize(nRT=mean(nRT), npRT=mean(npRT)) %>%
  summarize(mmRT=mean(nRT), mmpRT=mean(npRT), seRT=sd(nRT)/sqrt(N-1))

cprime_plot <- ggplot(sdata_col, aes(x=d, y=mmRT, color=rn)) + geom_point(size=3) + geom_line(aes(y=mmpRT, group=rn)) +
  geom_errorbar(aes(ymax=mmRT + seRT*1.96, ymin=mmRT-seRT*1.96), width=0.2) + theme_bw() +
  labs(x="Lag", y="Normalized RT (ms)", color="Color") + 
  theme(strip.background = element_blank(), strip.text=element_blank(), text = element_text(size=16)) 

sdata_pos <- pos_prime %>% separate(sub_session, c("sub", "session")) %>% group_by(sub, session) %>% 
  mutate(nRT=(mRT-mean(mRT))*1000, npRT=(mpRT-mean(mRT))*1000) %>%
  group_by(d, rn, sub) %>% summarize(nRT=mean(nRT), npRT=mean(npRT)) %>%
  summarize(mmRT=mean(nRT), mmpRT=mean(npRT), seRT=sd(nRT)/sqrt(N-1))

pprime_plot <- ggplot(sdata_pos, aes(x=d, y=mmRT, color=rn)) + geom_point(size=3) + geom_line(aes(y=mmpRT, group=rn)) +
  geom_errorbar(aes(ymax=mmRT + seRT*1.96, ymin=mmRT - seRT*1.96), width=0.2) + theme_bw() +
  labs(x="Lag", y="Normalized RT (ms)", color="Positional inter-trial condition") + 
  theme(strip.background = element_blank(), strip.text=element_blank(), 
        legend.position="bottom", text = element_text(size=16)) 

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

ggsave('suppfigures/cprime_cv.png', cprime_plot, width = 8, height = 6)
ggsave('suppfigures/pprime_cv.png', pprime_plot,width = 8, height = 6.5)
ggsave('suppfigures/prprime_cv.png', prprime_plot, width = 8, height = 4)

