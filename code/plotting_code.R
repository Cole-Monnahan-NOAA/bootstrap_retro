#if(!exists('results')) stop('run code/process_results.R first')

message("Set plotting theme...")

old_par <- par("mar")
par(mar=c(1.5,1.5,1,1))
par("mar")
#Set plotting theme
old_theme=theme_set(theme_bw())
theme_set(theme_bw()+theme(text = element_text(size=14),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
model_colors <- c('navyblue', 'gold')
obs_rho_col <- 'mediumvioletred'#'darkorange2'

## The models have different end years so the baseyears are
## different, but check for full replicates
results_afsc %>% dplyr::group_by(model,miller, baseyear) %>%
  dplyr::filter(metric=='SSB') %>% dplyr::summarize(count=n()) %>%
  pivot_wider(c(model, miller),  names_from='baseyear', values_from='count')

results_afsc$baseyear=as.numeric(as.character(results_afsc$baseyear))

## ## Are the variances the same?
#tmp <- results_afsc %>% dplyr::mutate(miller=replace_na(miller, FALSE),
#                       model=gsub('flathead', 'FHS', model)) %>%
#  dplyr::group_by(model, metric, miller,baseyear) %>%
#  dplyr::summarize(mn=mean(rho), med=median(rho), stdev=round(sd(rho),3),
#                   lci = mn-1.96*(stdev/sqrt(length(rho))), uci = mn+1.96*(stdev/sqrt(length(rho)))) %>%
#  dplyr::mutate(miller = case_when(miller ==FALSE ~ 'model',miller == TRUE ~ 'data') ) %>% dplyr::filter(metric=='SSB')

#ggplot(data = tmp,
#       aes(x= as.factor(baseyear), y=med, ymin=lci, ymax=uci, color=miller, fill=miller)) +
#  geom_linerange(size=1, position = position_dodge(width = 0.5)) +
#  geom_point(size=2, position = position_dodge(width = 0.5)) +
#  coord_cartesian(ylim=c(-0.5,0.5))+
#  geom_hline(yintercept=0,color='black', linetype='solid') +
#  geom_hline(aes(yintercept = 0.2, linetype = 'Rule of thumb'), color='black') +
#  geom_hline(aes(yintercept = -0.15, linetype = 'Rule of thumb'), color='black') +
#  geom_point(data =rho_obs_o %>% dplyr::filter(metric== x),
#             aes(x = as.factor(baseyear), y = rho, ymin=rho, ymax=rho), pch='-', size=10)+
#  scale_color_manual(name="Approach", values=c(model_colors,obs_rho_col)) +
#  scale_fill_manual(name="Approach", values=c(model_colors,obs_rho_col)) +
#  scale_linetype_manual(values=2) + labs(linetype = NULL) +
#  facet_wrap(.~model, scales='free_x') +
#  labs(x = 'Year', y = expression(rho))
## Needs updating since adding Miller stuff broke it:
## source('code/make_plots.R')

#Get median and quantile range for terminal year
rho_maxyr <- results_afsc %>% dplyr::group_by(model, miller, metric) %>%
                  dplyr::filter(metric!='Bratio') %>%
                  dplyr::mutate(baseyear = as.numeric(as.character(baseyear))) %>%
                  dplyr::group_by(model,metric,miller) %>%
                  dplyr::filter( baseyear == max(baseyear)) %>%
                  dplyr::group_by(baseyear,model,metric,miller) %>%
                  dplyr::summarise(x = quantile(rho, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), q=c(0.025,0.5,0.975)) %>%
                  ungroup
rho_maxyr2 <- tidyr::spread(rho_maxyr,q,x) %>%
                  dplyr::mutate(miller = case_when(miller ==FALSE ~ 'model',
                                                   miller == TRUE ~ 'data') )
names(rho_maxyr2)[5:ncol(rho_maxyr2)] <- c('lci','med','uci')

#Get observed rho value
rho_obs_o <- rho_obs %>% dplyr::mutate(baseyear = as.numeric(as.character(baseyear))) %>%
                          dplyr::group_by(model,miller) %>%
                          dplyr::filter( baseyear == max(baseyear)) %>%
                          dplyr::filter(miller == FALSE) %>%
                          dplyr::mutate(miller = 'rho obs')

#Plot terminal year quantile range of rho
#Include Hurtado-Ferro rule of thumb range, observed rho, and zero line
lapply(c('SSB','F','Rec'),function(x){
  ggplot(data = rho_maxyr2 %>% dplyr::filter(metric== x),
         aes(x= as.factor(baseyear), y=med, ymin=lci, ymax=uci, color=miller, fill=miller)) +
    geom_linerange(size=1, position = position_dodge(width = 0.5)) +
    geom_point(size=2, position = position_dodge(width = 0.5)) +
    coord_cartesian(ylim=c(-0.5,0.5))+
    geom_hline(yintercept=0,color='black', linetype='solid') +
    geom_hline(aes(yintercept = 0.2, linetype = 'Rule of thumb'), color='black') +
    geom_hline(aes(yintercept = -0.15, linetype = 'Rule of thumb'), color='black') +
    geom_point(data =rho_obs_o %>% dplyr::filter(metric== x),
               aes(x = as.factor(baseyear), y = rho, ymin=rho, ymax=rho), pch='-', size=10)+
    scale_color_manual(name="Approach", values=c(model_colors,obs_rho_col)) +
    scale_fill_manual(name="Approach", values=c(model_colors,obs_rho_col)) +
    scale_linetype_manual(values=2) + labs(linetype = NULL) +
    facet_wrap(.~model, scales='free_x') +
    labs(x = 'Year', y = expression(rho))
  ggsave(here::here('plots',paste0(x,'TerminalYearMohnsRho.png')), width = 6, height = 6, units = 'in')
})

  #Get median and quantile range for terminal year
  rho_yr <- results_afsc %>% dplyr::group_by(model, miller, metric) %>%
    dplyr::filter(metric!='Bratio') %>%
    dplyr::mutate(baseyear = as.numeric(as.character(baseyear))) %>%
    dplyr::group_by(baseyear, model, metric, miller) %>%
    dplyr::summarise(x = quantile(rho, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), q=c(0.025,0.5,0.975)) %>%
    ungroup
  rho_yr2 <- tidyr::spread(rho_yr,q,x) %>%
    dplyr::mutate(miller = case_when(miller ==FALSE ~ 'model',
                                     miller == TRUE ~ 'data') )
  names(rho_yr2)[5:ncol(rho_maxyr2)] <- c('lci','med','uci')

  #Get observed rho value
  rho_obs_yr <- rho_obs %>% dplyr::mutate(baseyear = as.numeric(as.character(baseyear))) %>%
    dplyr::group_by(baseyear,model,miller) %>%
    dplyr::filter(miller == FALSE) %>%
    dplyr::mutate(miller = 'rho obs')

lapply(c('SSB', 'Rec', 'F'), function(x){
  ggplot(data = rho_yr2 %>% dplyr::filter(metric==x),
         aes(x= as.factor(baseyear), y=med, ymin=lci, ymax=uci, color=miller, fill=miller)) +
    geom_linerange(size=1, position = position_dodge(width = 1)) +
    geom_point(size=2, position = position_dodge(width = 1)) +
    coord_cartesian(ylim=c(-0.75,0.75))+
    geom_hline(yintercept=0,color='black', linetype='solid') +
    geom_hline(aes(yintercept = 0.2, linetype = 'Rule of thumb'), color='black') +
    geom_hline(aes(yintercept = -0.15, linetype = 'Rule of thumb'), color='black') +
    geom_point(data =rho_obs_yr %>% dplyr::filter(metric== x),
               aes(x = as.factor(baseyear), y = rho, ymin=rho, ymax=rho), pch='-', size=10)+
    scale_color_manual(name="Approach", values=c(model_colors, obs_rho_col)) +
    scale_fill_manual(name="Approach", values=c(model_colors, obs_rho_col)) +
    scale_linetype_manual(values=2) +
    labs(linetype = NULL) +
    facet_wrap(.~model, scales='free_x') +
    labs(x = 'Year', y = expression(rho))
  ggsave(here::here('plots',paste0(x,'AnnualMohnsRho.png')), width = 10, height = 10, units = 'in')}
)

ts <- readRDS('results/ts.results.RDS')

s <- sample(1:100, 10)
lapply(s, function(s) {
  ip<- ts %>% dplyr::filter(boot==s & name == 'ssb' & year >=1990) %>%
    dplyr::mutate(peel_yr = assess_yr - peel)
  ip$peel_yr <- factor(ip$peel_yr, levels = rev(sort(unique(ip$peel_yr))))

  term_pd <- ip %>% dplyr::group_by(model, miller, peel) %>%
    dplyr::filter(year == max(year))

  ggplot(data = ip, aes(x=year, y=value/1e4, color = peel_yr)) +
    geom_line() +
    geom_point(data=term_pd, aes(x=year, y=value/1e4, color=peel_yr)) +
    expand_limits(y=0) +
    facet_grid(model~miller, scales = 'free') +
    scale_color_viridis_d(name="Base year", direction=1, option='inferno') +
    labs(y='SSB (10,000 mt)', x='Year')
  ggsave(here::here('plots',paste0('retro_timeseries_boot_',s,'.png')),width = 8, height = 8, units='in')
})


ts4 <- ts %>%
  dplyr::filter(boot>0) %>%
  dplyr::group_by(model,year,peel,name,miller) %>%
  dplyr::summarise(x = quantile(re, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975))
ts5 <- tidyr::spread(ts4,q,x)
names(ts5)=c(names(ts5[1:5]),'lci','med','uci')

met=sort(unique(ts5$name))
for(j in 1:length(met))
{
  g=ts5 %>% filter(name==met[j] & peel==0) %>% mutate(peel=as.factor(peel))

  ggplot(data=g,aes(x=year, y=med, color=peel)) + geom_line(aes(x=year, y=med, color=peel))+
    geom_ribbon(aes(ymin=lci,ymax=uci, fill=peel),alpha=0.3)+ facet_grid(miller~model)+
    coord_cartesian(ylim=c(-1,1))+geom_hline(yintercept=0,color="red")+
    ylab(paste("Relative error in ", met[j],sep=" ")) + xlab("Year")
    ggsave(here::here('plots',paste0('RelError_',met[j],'.png')),units='in',width=8,height=8)
}

