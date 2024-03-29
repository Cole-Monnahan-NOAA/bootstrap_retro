## test that the bootstrap code for the Miller approach is
## working

## test with different input files
datlist <- SS_readdat('models/BSAI_FHS/data.ss')
datlist <- SS_readdat('models/GOA_NRS/data.ss')
datlist <- SS_readdat('models/GOA_SRS/data.ss')

## Bootstrap simulations and combine together to get distributions
out <- lapply(1:50, function(boot) sample_miller_boot(boot, datlist, test=TRUE))
indices <- lapply(out, function(x) x$cpue) %>% bind_rows()
lencomps <- lapply(out, function(x) x$lencomp) %>% bind_rows()
## Split age comps into marginal and CAAL
tmp <- lapply(out, function(x) x$agecomp) %>% bind_rows()
agecomps <- filter(tmp, Lbin_lo<0)
caalcomps <- filter(tmp, Lbin_lo>0)

## Now get the observed ones to plot on top
lencomp.long <- datlist$lencomp %>%
  pivot_longer(-(1:6), values_to='proportion') %>%
  mutate(rep=-1, sex=substr(name,0,1), len=as.numeric(gsub('f|m','', name)))
## Split the age comps into marginal and conditional age at
## length after normalizing each row.
tmp <- datlist$agecomp
tmp[,-(1:9)] <- tmp[,-(1:9)]/rowSums(tmp[,-(1:9)])
agecomp.long <- tmp  %>% filter(Lbin_lo<0) %>%
  pivot_longer(-(1:9), values_to='proportion') %>%
  mutate(rep=-1, sex=substr(name,0,1),
         age=as.numeric(gsub('f|m','', name)))
caalcomp.long <- tmp  %>% filter(Lbin_lo>0) %>%
  pivot_longer(-(1:9), values_to='proportion') %>%
  mutate(rep=-1, sex=substr(name,0,1),
         age=as.numeric(gsub('f|m','', name)))
## Do some fancy filtering to account for the different Gender
## options whihc otherwise won't make sense and not match
agecomps <- filter(agecomps, !(Gender==1 & sex=='m') &
                               !(Gender==2 & sex=='f'))
agecomp.long <- filter(agecomp.long, !(Gender==1 & sex=='m') &
                               !(Gender==2 & sex=='f'))
lencomps <- filter(lencomps, !(Gender==1 & sex=='m') &
                               !(Gender==2 & sex=='f'))
lencomp.long <- filter(lencomp.long, !(Gender==1 & sex=='m') &
                               !(Gender==2 & sex=='f'))
caalcomps <- filter(caalcomps, !(Gender==1 & sex=='m') &
                               !(Gender==2 & sex=='f'))
caalcomp.long <- filter(caalcomp.long, !(Gender==1 & sex=='m') &
                               !(Gender==2 & sex=='f'))
## Normalize the raw data since not necessary sum to 1 but boot
## output is
lencomp.long  <- lencomp.long %>%
  group_by(Yr, FltSvy, Gender) %>%
  mutate(proportion=proportion/sum(proportion)) %>% ungroup
agecomp.long  <- agecomp.long %>%
  group_by(Yr, FltSvy, Gender, Lbin_lo) %>%
  mutate(proportion=proportion/sum(proportion)) %>% ungroup
caalcomp.long  <- caalcomp.long %>%
  group_by(Yr, FltSvy, Gender, Lbin_lo) %>%
  mutate(proportion=proportion/sum(proportion)) %>% ungroup

g.index <- ggplot(indices, aes(year, log(obs))) + geom_point(alpha=.25) +
  geom_point(data=datlist$CPUE, col='red', size=2.5) +
  facet_wrap('index') + labs(y='index')
if(nrow(lencomp.long)>0){
  g.lencomps <- ggplot(filter(lencomps, Yr==max(Yr)),
                       aes(len, proportion)) +
  ##  facet_wrap('Yr', scales='free') +
  facet_grid(Yr+FltSvy~sex) +
  geom_jitter(alpha=.25, width=.5, height=0) +
  geom_point(data=filter(lencomp.long, Yr == max(Yr)),
             col='red', size=2.5) +
    labs(y='length comp')
  } else {g.lencomps <- NULL}
if(nrow(agecomp.long)>0){
g.agecomps <- ggplot(filter(agecomps, Yr==max(Yr)), aes(age, proportion)) +
  ##  facet_wrap('Yr', scales='free') +
  facet_grid(Yr+FltSvy~sex) +
  geom_jitter(alpha=.25, width=.15, height=0)+
  geom_point(data=filter(agecomp.long, Yr==max(Yr)),
             col='red', size=2.5) +
  labs(y='age comp')
} else { g.agecomps <- NULL}
if(nrow(caalcomp.long)>0){
  tmp1 <- caalcomps %>%
    filter(Yr==max(Yr)) %>%
    filter(Lbin_lo > 30 & Lbin_lo <40)
    ##filter(Lbin_lo==min(Lbin_lo) | Lbin_lo==max(Lbin_lo))
  tmp2 <- caalcomp.long %>%
    filter(Yr==max(Yr)) %>%
    filter(Lbin_lo > 30 & Lbin_lo <40)
  ##filter(Lbin_lo==min(Lbin_lo) | Lbin_lo==max(Lbin_lo))
  g.caalcomps <-
    ggplot(tmp1, aes(age, proportion)) +
    facet_grid(sex+FltSvy~Lbin_lo, scales='free') +
    ##facet_grid(Yr+FltSvy~sex+Lbin_lo) +
    geom_jitter(alpha=.25, width=.15, height=0) +
    geom_point(data=tmp2, col='red', size=2.5) +
    labs(y='age comp')
} else { g.caalcomps <- NULL}

g.index
g.lencomps
g.agecomps
g.caalcomps

## numerical checks
indices.E <- indices %>% group_by(year) %>%
  summarize(avg=mean(log(obs))) %>% pull(avg)
(log(datlist$CPUE$obs) - indices.E) # appears unbiased

## try a single serial iteration
run_SS_boot_iteration(25, 'BSAI_FHS', clean.files=FALSE, miller=TRUE)

## Try a few iterations in parallel. Load all the global stuff in
## run_models.R before
results.list <- sfLapply(1:15, function(i)
  run_SS_boot_iteration(boot=i, 'fhs', clean.files=FALSE, miller=TRUE))

