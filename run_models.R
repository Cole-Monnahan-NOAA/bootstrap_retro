### Run and compare models
library(tidyverse)
library(r4ss)
library(snowfall)
library(ggplot2)
library(R2admb)
## devtools::install_github("afsc-assessments/GOApollock", ref='fix_dat_fns')
library(GOApollock)
theme_set(theme_bw())
packageVersion('r4ss') #  '1.42.0'
source('code/functions.R')
## GOA pollock is a bespoke model and setup differently
pkdatlist <- readRDS("models/GOA_pollock/datfile.RDS")
pkreplist <- readRDS("models/GOA_pollock/repfile.RDS")
source('code/functions_pollock.R')

## For each boostrap data set, run a retrospective analysis
Nreps <- 500
reps <- 0:Nreps # 0 is special code for original data
Npeels <- 14
peels <- 0:-Npeels

## Setup to run parallel, saving a single core free.
cpus <- parallel::detectCores()-2
sfStop()
sfInit( parallel=TRUE, cpus=cpus)
sfExportAll()

### Run full in parallel for all models. This assumes that the
### starter file was modifed to produce 500 replicates (502
### total) and was run. The "blank" files should also have
### started file modified to reduce output to speed things
### up. Finally, the data input must be named data.ss in the main
### files, to facilitate the Miller sampling scheme.

## ## Run one in serial as a test
## test <- run_SS_boot_iteration(1, 'EBS_Pcod', TRUE)
run_model(reps, model.name='GOA_NRS',clean.files=TRUE)
run_model(reps, model.name='GOA_NRS',miller=TRUE, clean.files=TRUE)



run_model(reps, model.name='EBS_Pcod')
run_model(reps, model.name='EBS_Pcod', miller=TRUE)
## run_model(reps, model.name='GOA_Pcod_prior')
## run_model(reps, model.name='GOA_Pcod_prior', miller=TRUE)
run_model(reps, model.name='GOA_Pcod_noprior')
run_model(reps, model.name='GOA_Pcod_noprior', miller=TRUE)

## make sure to delete runs and result files before
## rerunning
## unlink('runs/GOA_pollock/', recursive=TRUE)
## run_pollock_boot_iteration(boot=0, datlist=pkdatlist, replist=pkreplist)
run_pollock_model(reps,datlist=pkdatlist, replist=pkreplist,
                  model.name='GOA_pollock', miller=TRUE,
                  clean.files=TRUE)
run_pollock_model(reps,datlist=pkdatlist, replist=pkreplist,
                  model.name='GOA_pollock', miller=FALSE,
                  clean.files=TRUE)

source('code/process_results.R')
## The models have different end years so the baseyears are
## different, but check for full replicates
results_afsc %>% group_by(model,miller, baseyear) %>%
  filter(metric=='SSB') %>% summarize(count=n()) %>%
  pivot_wider(c(model, miller),  names_from='baseyear', values_from='count')

filter(results_afsc, boot==1 & grepl('GOA_Pcod', x=model) &
                     metric=='SSB' & baseyear=='2015')


mean_CI_byYr=results_afsc %>% group_by(model,miller, baseyear, metric) %>% summarise( lci_rho = CI(rho)[3], mean_rho = mean(rho),uci_rho = CI(rho)[1])
write.csv(mean_CI_byYr,"rho_meanCIbyYear.csv",row.names=FALSE)

min_max_byYr=results_afsc %>% group_by(model,miller, baseyear, metric) %>% summarise(min_rho = min(rho),max_rho = max(rho))

min_max_ssb_byYr=results_afsc %>% filter(metric=='SSB' & miller==FALSE) %>% group_by(model,miller, baseyear, metric) %>%
  summarise(min_rho = min(rho),max_rho = max(rho),width=max_rho-min_rho)

mean_CI=results_afsc %>% group_by(model,miller, metric) %>% summarise(lci_rho = CI(rho)[3], mean_rho = mean(rho),uci_rho = CI(rho)[1])
write.csv(mean_CI,"rho_meanCI.csv",row.names=FALSE)

## Quick plot of Miller vs SS bootstrap
g <- results_afsc %>%
 # filter(metric!='Rec') %>%
  ggplot(aes(baseyear, y=rho, fill=miller)) + geom_violin() +
  facet_grid(metric~model, scales='free') + geom_hline(yintercept=0, col='red')+
  geom_point(data=rho_obs, col='red',pch='-', size=10) +
  coord_cartesian(ylim=c(-.5,.5))
ggsave('plots/results_miller2.png', g, width=12, height=7)

# New facet label names for model variable
g2 <- results_afsc %>%
  filter(metric %in% c('Rec','F','SSB'))%>%
  ggplot(aes(baseyear, y=rho, fill=miller)) + geom_violin() +
  stat_summary(fun = "mean", geom = "pointrange",position = position_dodge(width = 0.9),color = "black",size=0.15) +
  facet_grid(metric~model, scales='free') + geom_hline(yintercept=0, col='red')+
  geom_point(data=rho_obs[rho_obs$metric %in% c('Rec','F','SSB') ,], col='red',pch='-', size=10) +
  coord_cartesian(ylim=c(-.6,.6))+ scale_fill_discrete(name='Method',labels=c('model',"data")) + labs(x="Base year",y='rho') +
  geom_hline(yintercept=c(-0.15,0.2),linetype='dashed')+
  theme(legend.text=element_text(size=12),axis.title=element_text(size=16),axis.text=element_text(size=16),strip.text=element_text(size=12))
ggsave('plots/results_miller_RecFSSB.png', g2, width=18, height=10)

g3 <- results_afsc %>%
  filter(metric %in% c('SSB') & model!='EBS_Pcod')%>%
  ggplot(aes(baseyear, y=rho, fill=miller)) + geom_violin() +
  facet_grid(metric~model, scales='free') + geom_hline(yintercept=0, col='red')+
  geom_point(data=rho_obs[rho_obs$metric %in% c('SSB') & rho_obs$model!='EBS_Pcod',], col='red',pch='-', size=10) +
  coord_cartesian(ylim=c(-.5,.5))+ scale_fill_discrete(name='Method',labels=c('model',"data")) + labs(x="Base year",y='rho') +
  geom_hline(yintercept=c(-0.15,0.2),linetype='dashed')+
  theme(legend.text=element_text(size=16),axis.title=element_text(size=16),axis.text=element_text(size=12),strip.text=element_text(size=16))
ggsave('plots/results_miller_SSB_noEBSPcod.png', g3, width=18, height=8)

results_afsc$metric <-factor(results_afsc$metric, levels=c("F","SSB","Bratio"),labels=c("F","SSB","Bratio"))
results_afsc$model <-factor(results_afsc$model, levels=c("EBS_Pcod", "GOA_NRS", "GOA_Pcod"),labels=c("EBS Pacific cod", "GOA northern rock sole", "GOA Pacific cod"))
rho_obs$metric <-factor(rho_obs$metric, levels=c("F","SSB","Bratio"),labels=c("F","SSB","Bratio"))
rho_obs$model <-factor(rho_obs$model, levels=c("EBS_Pcod", "GOA_NRS", "GOA_Pcod"),labels=c("EBS Pacific cod", "GOA northern rock sole", "GOA Pacific cod"))

g4 <- results_afsc %>%
  filter(metric %in% c('F','SSB'))%>%
  ggplot(aes(baseyear, y=rho, fill=miller)) + geom_violin() +
  stat_summary(fun = "mean", geom = "pointrange",position = position_dodge(width = 0.9),color = "black",size=0.15) +
  facet_grid(metric~model, scales='free') + geom_hline(yintercept=0, col='red')+
  geom_point(data=rho_obs[rho_obs$metric %in% c('F','SSB') ,], col='red',pch='-', size=10) +
  coord_cartesian(ylim=c(-.6,.6))+ scale_fill_discrete(name='Approach',labels=c('Model',"Data")) + labs(x="Base year",y='rho') +
  geom_hline(yintercept=c(-0.15,0.2),linetype='dashed')+
  th
ggsave('plots/results_miller_FSSB.png', g4, width=18, height=10)

#terminal year plot only
rho_obs$baseyear=as.numeric(as.character(rho_obs$baseyear))
rho_obs2=rho_obs %>%
  group_by(model) %>% filter(metric %in% c('F','SSB') & baseyear==max(baseyear))

results_afsc$baseyear=as.numeric(as.character(results_afsc$baseyear))
g5 <- results_afsc %>%
  group_by(model) %>% filter(metric %in% c('F','SSB') & baseyear==max(baseyear))%>%
  ggplot(aes(as.factor(baseyear), y=rho, fill=miller)) + geom_violin() +
  facet_grid(metric~model, scales='free') + geom_hline(yintercept=0, col='red')+
  geom_point(data=rho_obs2, col='red',pch='-', size=10) +
  coord_cartesian(ylim=c(-.5,.5))+ scale_fill_discrete(name='Method',labels=c('model',"data")) + labs(x="Base year",y='rho') +
   geom_hline(yintercept=c(-0.15,0.2),linetype='dashed')+
  theme(legend.text=element_text(size=16),axis.title=element_text(size=16),axis.text=element_text(size=16),strip.text=element_text(size=16))
ggsave('plots/results_miller_FSSB_termYr.png', g5, width=12, height=7)



## ## Are the variances the same?
## results_afsc %>% mutate(miller=replace_na(miller, FALSE),
##                         model=gsub('flathead', 'FHS', model)) %>%
##   group_by(model, metric, miller) %>%
##   summarize(stdev=round(sd(rho),3)) %>% pivot_wider(c(model, metric),
##   names_from=miller, names_prefix='miller=', values_from=stdev)

## Needs updating since adding Miller stuff broke it:
## source('code/make_plots.R')
