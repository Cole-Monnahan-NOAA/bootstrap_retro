
#FUNCTION to get results from Report.sso
get_res<-function(wd, boot, model, miller)
{
  library(dplyr)
  get_ts <- function(rep, peel){
    ##Get timeseries results
    ##ssb and apical F
    ssb_f_rec_ts<-rep$sprseries[rep$sprseries$Yr<=rep$endyr,]

    ##Create data frames
    ssb<-data.frame(assess_yr=rep$endyr, name='ssb', peel=peel, year=ssb_f_rec_ts$Yr, value=ssb_f_rec_ts$SSB) #this was set to SSBfished, which was incorrect.
    ap_f<-data.frame(assess_yr=rep$endyr, name='f', peel=peel, year=ssb_f_rec_ts$Yr, value=ssb_f_rec_ts$sum_Apical_F)
    recr<-data.frame(assess_yr=rep$endyr, name='recr', peel=peel,year=ssb_f_rec_ts$Yr,value=ssb_f_rec_ts$Recruits)

    ##Combine all timeseries results
    ts_res<-bind_rows(ssb,ap_f,recr)
    return(ts_res)
  }

  get_mles <- function(rep, peel){
    ##Get parameters
    pars<-rep$parameters
    pars_res <- data.frame(assess_yr=rep$endyr, par=pars$Label,
                         peels=peel, value=pars$Value)
    return(pars_res)
  }

  dir <- list.files(file.path(wd, 'retros'), pattern='retro',
                    full.names=TRUE)
  peels <- gsub('retro','', x=list.files(file.path(wd, 'retros'),
                                         pattern='retro')) %>% as.numeric %>% abs
  stopifnot(length(peels)>0)
  Models <- SSgetoutput(dirvec=dir, getcovar=FALSE)
  pars_res <- lapply(1:length(Models), function(ii) {
    get_mles(Models[[ii]], peels[ii])
  }) %>% bind_rows %>% cbind(boot=boot, model=model, miller=miller)
  write.csv(pars_res, file=file.path(wd,"parameter_res.csv"), row.names=FALSE)
  ts_res <- lapply(1:length(Models), function(ii) {
    get_ts(Models[[ii]], peels[ii])
  }) %>% bind_rows %>% cbind(boot=boot, model=model, miller=miller)
  write.csv(ts_res, file=file.path(wd,"ts_res.csv"), row.names=FALSE)
}



#' Generate a Miller bootstrap data list
#'
#' @param boot The bootstrap number (used as seed)
#' @param datlist A data file list as returned by SS_readdat
#' @param test Whether to return output useful for testing
#'   purposes rather than the modified datlist
#' @details Function to do a parametric bootstrap given an SS
#'   input file. This resamples the data given the assumed
#'   distribution of observed data. It does not use the
#'   expectations like the SS bootstrap does. Inspired by Miller
#'   and Legault (2017).
sample_miller_boot <- function(boot, datlist, test=FALSE){
  set.seed(boot)
  if(test) library(tidyverse)
  xnew <- x <- datlist
  ind <- -(1:6) # len comp cols to drop
  ind2 <- -(1:9) # age comps has extra cols

  ncpue <- nrow(x$CPUE)
  if(ncpue>0){
    xnew$CPUE$obs <- rlnorm(nrow(x$CPUE), log(x$CPUE$obs), x$CPUE$se_log)
  }

  ## Carefully resample from a comps, dealing with dummy data and
  ## unscaled inputs
  resample_comps <- function(bins, prob, Nsamp, sex, i, type){
    prob0 <- prob ## save in case want ot inspect during browser
    ## convert bins to an index
    bins <- seq_along(bins)
    if(any(!is.finite(prob))){
      print(prob)
      message('Bad ', type, ' comp probs for row  ', i)
      browser()
    }
    if(sex!=0 & length(prob)!=2*length(bins))
        stop("length of probabilities in ", type, " ", i, " does not equal length of bins")
    Nsamp <- max(Nsamp, 1) # sometimes < 1 which breaks sampler
    if (sex==1 | sex==0){
      ## combined or females only so male columns are dummy, zero them out
      prob[-bins] <- 0
    } else if(sex==2) {
      ## female columns dummy
      prob[bins] <- 0
    } else if(sex!=3){
      ## sex=3 is females then males so leave as is
      stop("sex=", sex, "not setup to work")
    }
    ## Normalize since SS does this internally people enter
    ## whatever. Needed to fix dummy before doing it.
    prob <- prob/sum(prob)
    if(all(prob==0)) browser()
    if(any(!is.finite(prob))) browser()
    tmp <- rmultinom(1, size=Nsamp, prob=prob)[,1]
    if(any(!is.finite(prob)) | any(!is.finite(tmp))){
      print(prob); print(tmp)
      stop("NAs in simualted data")
    }
    tmp <- tmp/sum(tmp)
    return(tmp)
  }
  lc <- x$lencomp
  if(!is.null(lc)){
    nlencomp <- nrow(x$lencomp)
    stopifnot(nlencomp>0)
    lbins <- x$lbin_vector
    for(i in 1:nlencomp){
      ## Skip resampling ghost data
      if(lc$Yr[i]<0 | lc$FltSvy[i]<0) next
      prob <- as.numeric(lc[i,ind])
      xnew$lencomp[i,ind] <-
        resample_comps(lbins, prob, lc$Nsamp[i], lc$Gender[i], i, 'length')
    }
    ## Make long version for plotting
    if(test)
      lencomp.long <- pivot_longer(xnew$lencomp, -(1:6), values_to='proportion') %>%
        mutate(boot=boot, sex=substr(name,0,1), len=as.numeric(gsub('f|m','', name)))
  }

  ## Repeat with ages and CAAL
  ac <- x$agecomp
  if(!is.null(ac)){
    nagecomp <- nrow(ac)
    stopifnot(nagecomp>0)
    abins <- x$agebin_vector
    for(i in 1:nagecomp){
      ## Skip resampling ghost data
      if(ac$Yr[i]<0 | ac$FltSvy[i]<0) next
      prob <- as.numeric(ac[i,ind2])
      xnew$agecomp[i,ind2] <-
        resample_comps(abins, prob, ac$Nsamp[i], ac$Gender[i], i, 'age')
    }
    ## Make long version for plotting
    if(test)
      agecomp.long <- pivot_longer(xnew$agecomp, -(1:9), values_to='proportion') %>%
        mutate(boot=boot, sex=substr(name,0,1), age=as.numeric(gsub('f|m','', name)))
  }
  stopifnot(x$N_meanbodywt==0)
  stopifnot(x$N_MeanSize_at_Age_obs==0)
  ## If testing return the simulated bits in tidy format for
  ## plotting later
  if(test){
    out <- list(cpue=cbind(boot=boot,xnew$CPUE), lencomp=lencomp.long,
                agecomp=agecomp.long)
    return(out)
  }
  ## otherwise return modified datfile to write
  return(xnew)
}



#' Calculation the retrospective metrics for a single
#' bootstrapped data set. Designed to work with parallel
#' execution.
run_SS_boot_iteration <- function(boot, model.name,
                                  clean.files=TRUE, miller=FALSE){
  ## Some of these are global variables
  library(r4ss)
  if(!miller){
    wd <- file.path('runs', model.name, paste0("boot_", boot))
  } else {
    wd <- file.path('runs', model.name, paste0("millerboot_", boot))
  }
  rhofile <- ifelse(miller, 'results_miller_rho.csv', 'results_rho.csv')
  if(file.exists(file.path(wd, rhofile))) {
    message("Skipping ", boot, " of ", model.name, "for miller=", miller,
            " because it exists already. Delete manually if needed.")
    return(NULL)
  }
  ## boot==0 is the original data, so skip resampling the data
  if(boot==0){
    ## The original data, whether miller or not
    dat <- SS_readdat(file.path('models', model.name,'data.ss'),
                      verbose=TRUE, section=1)
  } else {
    ## resample using SS or miller approach
    if(!file.exists(file.path('models',model.name)))
      stop("model not found")
    ## The two types of bootstraps are implemented here
    if(!miller){
      dat <- SS_readdat(file.path('models', model.name,'data.ss_new'),
                        verbose=TRUE, section=2+boot)
    } else {
      ## Original data
      dat0 <- SS_readdat(file.path('models', model.name,'data.ss'),
                         verbose=TRUE, section=1)
      ## This resamples the observed data
      dat <- sample_miller_boot(boot=boot, datlist=dat0, test=FALSE)
    }
  }
  ## Prepare folder to run this iteration
  dir.create(wd, showWarnings=TRUE, recursive=TRUE)
  blank.files <- list.files(file.path('models', model.name,'blank'), full.names=TRUE)
  test <- file.copy(from=blank.files, to=wd, overwrite=TRUE)
  if(!all(test)){
    message(paste(blank.files[!test], collapse= '\n'))
    stop("Some blank files failed to copy for iteration ", boot)
  }
  ## Write new data
  SS_writedat(dat, outfile=paste0(wd, '/data.ss'), verbose=FALSE,
              overwrite=TRUE)
  newfiles <- list.files(wd, full.names=TRUE) # mark for deletion at end
  ## Run retro for this bootstrap one
  retro(masterdir=getwd(), oldsubdir=wd,
             newsubdir=paste0(wd, '/retros'),
             years=peels, extras='-nohess -nox -iprint 1000')
  dirvec <- file.path(paste0(wd, '/retros'), paste0("retro", peels))
  if(length(dirvec)!=Npeels+1)
    stop("Some retro runs missing in", model.name, " iteration ",
         boot)
  ## Read in all peels
  retroModels <- SSgetoutput(dirvec=dirvec, getcovar=FALSE)
  saveRDS(retroModels, file=paste0(file.path(wd,'retroModels.RDS')))
  ## Pluck out a subset of 7 peels and calculate rho
  pdf(file=file.path(wd, 'retro_plots.pdf'), onefile=TRUE,
      width=7, height=9)
  rhos <- list(); k <- 1
  ## Read in and save timeseries and parameter estimates for all
  ## peels
  get_res(wd, boot, model.name, miller)
  for(peelyr in 0:7){
    peels.tmp <- peels[(1+peelyr):(peelyr+8)] # peel #
    peels.ind <- which(peels %in% peels.tmp)  # peel index
    retroSummary <- SSsummarize(retroModels[peels.ind], verbose=FALSE)
    endyrvec <- retroSummary$endyrs + peels.tmp
    if(model.name %in% c('BSAI_GT', 'BSAI_GT2')) {
      ## What is going on here?? hack to fix error
      retroSummary$startyrs <- rep(1961, length(peels))
    }
    ## We have 14 peels so with a 7-peel year we can do 7
    par(mfrow=c(3,1))
    SSplotComparisons(retroSummary, subplots=c(1,3,9), endyrvec=endyrvec, png=FALSE, plot=TRUE,
                      plotdir=wd,  new=FALSE, add=FALSE, uncertainty=FALSE,
                      legendlabels=paste("Data",peels,"years"))
    mtext(paste0('Base year=', endyrvec[1]), side=3, line=-3,
          outer=TRUE)
    ## Calculate Mohn's rho. Not sure why I need to increase
    ## startyr here but errors out for some models if I don't. It
    ## looks like this is only relevant for WHOI rho calcs so seems
    ## fine to ignore.
    rho <- SSmohnsrho(retroSummary, endyrvec=endyrvec, startyr=retroSummary$startyr[1]+1)
    rhos[[k]] <- data.frame(model=model.name, baseyear=endyrvec[1], miller=miller, boot=boot, rho)
    k <- k+1
  }
  dev.off()
  rhos <- do.call(rbind, rhos)
  write.csv(x=rhos, file=file.path(wd, rhofile), row.names=FALSE)
  if(clean.files){
    unlink(file.path(wd, 'retros'), recursive=TRUE)
    ##file.remove(file.path(wd, 'retroSummary.RDS'))
    file.remove(file.path(wd, 'retroModels.RDS'))
    ## trash <-
    ##   file.remove(list.files(file.path(wd), pattern='.exe', full.names=TRUE))
    ## tmp <- file.path(wd, 'wtatage.ss')
    ## if(file.exists(tmp)) file.remove(tmp)
    ## ## this can be large so it'd be nice to keep but deleting for now
    ## file.remove(file.path(wd, 'data.ss'))
    file.remove(newfiles)
  }
  return(rhos)
}

#' Wrapper to run and save a single model
run_model <- function(reps, model.name, miller=FALSE, clean.files=TRUE){
  ## Run all bootstrap results. The clean.files argument is
  ## helpful b/c it's Nreplicates*Npeels SS runs which gets huge
  ## fast.
  trash <- sfLapply(reps, function(i)
    run_SS_boot_iteration(boot=i, model.name=model.name, clean.files=clean.files,
                          miller=miller))

  ## It fails on some. Not sure why this is happening. But hack is
  ## to loop through and figure out which failed and simply rerun
  ## them. Try 5 loops and break out if they all worked.
  tmp <- ifelse(miller, 'results_miller_rho', 'results_rho')

  for(i in 1:5){
    ff <- list.files(path=file.path('runs', model.name),
                     pattern=tmp, recursive=TRUE, full.names=TRUE)
    results <- lapply(ff, read.csv) %>% bind_rows()
    ind <- which(!reps %in% results$boot)
    if(length(ind)>0){
      warning("Rerunning failed ", model.name, " models= ", paste(ind, collapse=','))
      trash <- sfLapply(reps, function(i)
        run_SS_boot_iteration(boot=i, model.name=model.name, clean.files=clean.files,
                              miller=miller))
    }
  }
  ## Read in all final results, including those not necessarily
  ## in reps since they could have been run ea
  ff <- list.files(path=file.path('runs', model.name),
                   pattern=tmp, recursive=TRUE, full.names=TRUE)
  results <- lapply(ff, read.csv) %>% bind_rows()
  f <- paste0(model.name,ifelse(miller, '_millerboot_retros.RDS',
                                '_boot_retros.RDS'))
  saveRDS(results, file=file.path('results', f))
  message("All replicates finished for model=", model.name,
    " and miller=", miller)
}



