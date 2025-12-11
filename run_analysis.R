### Explore uncertainty for retrospective patterns in assessments

## This script runs the analysis for the in prep manuscript at
## https://github.com/Cole-Monnahan-NOAA/bootstrap_retro

## Step 1 is to prepare the workspace
library(tidyverse)
library(r4ss)
library(snowfall)
library(ggplot2)
library(Rmisc)
theme_set(theme_bw())
packageVersion('r4ss') #  '1.52.1'
R.version              # 4.5.1
source('code/functions.R')           # this has the key functions

## Note that the models are not on the repo and must be obtained
## offline.
list.dirs('models')

## Step 2 is to run retrospectives using two bootstrap
## approaches, tracked via a 'miller' variable. See MS for more
## details. This will run everything in parallel using (available
## cores - 2) cores. Most output is too large to save on the repo.
#source('code/run_models.R')

## Step 3 is to process the outputted results and save them as
## condensed files
source('code/process_results.R')
str(results_afsc)                       # the main data to use


## Step 4 is to make tables/figures
source("code/plotting_code.R")
#source("code/make_plots.R")


## old code to check if updating EBS_Pcod had a big effect or
m1 <- SS_output('models/EBS_Pcod', verbose=FALSE)
m2 <- SS_output('models/EBS_Pcod2', verbose=FALSE, covar=FALSE)
m3 <- SS_output('models/EBS_Pcod3', verbose=FALSE, covar=FALSE)
m4 <- SS_output('models/EBS_Pcod4', verbose=FALSE, covar=TRUE)
m5 <- SS_output('runs/EBS_Pcod4/boot_167/retros/retro0', verbose=FALSE, covar=FALSE)
SS_plots(m4)
x <- SSsummarize(biglist=list(m1=m1,m2=m2,m3=m3,m4,m5))
SSplotComparisons(x, png=TRUE,
                  plotdir='models/ebs_comparison',
                  legendlabels=c('orig', 'reweight', 'DM off',
                                 'bias ramp', 'arb boot'))

