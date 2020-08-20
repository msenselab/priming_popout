library(tidyverse)
library(data.table)
library(ez)

source('bayesian_updates.R')

# Read data
subfiles=dir('data', '*.dat')
data <- data.table(do.call('rbind',lapply(subfiles, readData)))
data <- unite(data, "sub_session", sub, session)

# Summarize over all different intertrial conditions from one trial back to explore inter-trial effects
inttrial <- separate(data, sub_session, c("sub", "session")) %>% 
  group_by(tcolprime, toriprime, tposprime, sub, session) %>% 
  filter(error=="no", tposprime!="NA") %>% summarize(rt=mean(rt)*1000) %>% summarize(rt=mean(rt)) 

N <- length(unique(inttrial$sub))

# Plot summary of inter-trial effects
fig_inttrial <- inttrial %>% summarize(mRT=mean(rt), seRT=sd(rt)/sqrt(N-1)) %>% 
  ggplot(aes(x=toriprime, y=mRT, color=tcolprime, group=tcolprime)) + 
  geom_errorbar(aes(ymin=mRT-seRT*1.96, ymax=mRT+seRT*1.96), width=0.2) +
  geom_point(size=3) + geom_line() + theme_bw() + facet_wrap(~ tposprime) +
  labs(x="Target notch position (RCF)", y="Response time (ms)", color="Color") + 
  theme(text = element_text(size=14)) 

if(exists("savefigures") && savefigures==TRUE) {
  ggsave('figures/inttrial.png', fig_inttrial, width = 8, height = 3)
}
