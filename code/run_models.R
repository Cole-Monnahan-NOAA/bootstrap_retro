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

