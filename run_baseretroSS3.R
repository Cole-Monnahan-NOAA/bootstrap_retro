#Run retrospective analysis and plot retrospective time series with CI
library(here)
#if(!require("remotes"))install.packages("remotes");
#remotes::install_github("r4ss/r4ss")
library(r4ss)

if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

old_par <- par("mar")
par(mar=c(1.5,1.5,1,1))
par("mar")
#Set plotting theme
old_theme=theme_set(theme_bw())
theme_set(theme_bw()+theme(text = element_text(size=14),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

# Directory with full assessment model results
dirs=list.dirs(here('models'),recursive=FALSE)
endyrvec <- 0:-10	   	      #retrospective peels

for(i in 1:length(dirs))
{
  list.of.files <- list.files(here(dirs[i])) # files to copy to Original folder

  # Retrospective sub dir, created within master
  retros="retrospectives"
  #Create main retrospective folder
  if(!dir.exists(here(dirs[i], retros))){
    dir.create(here(dirs[i], retros), showWarnings = TRUE)
  }

  main_dir<-here(dirs[i], retros)

  #create "Original" subfolder and copy all original assessment files into this subfolder
  #needed for SS_doRetro and SS_doRetro_surveRed
  if(!dir.exists(here(dirs[i], "Original"))){
    dir.create(here(dirs[i],"Original"), showWarnings = TRUE)
  }

  file.copy(from = here(dirs[i],list.of.files), to = here(dirs[i],"Original"),overwrite=TRUE) # copy the files to the new folder

  #Typical retrospective
  retro(dir=dirs[i], oldsubdir="Original", newsubdir = retros,
             subdirstart = 'retro', years = endyrvec , overwrite = TRUE,
             exe = "ss", extras = "-nox")

}

#Plot time series of SSB with 95% confidence envelopes
  #Read report files and get derived quant and sd
  #Calc 95% confidence envelope
if (file.exists(here('models',"retro_summarySS3.csv"))) {
  file.remove(here('models',"retro_summarySS3.csv"))
}


  for(i in 1: length(dirs))
  {
    folders<-list.files(here(dirs[i],'retrospectives'),pattern="retro")
    for(j in 1:length(folders))
    {
      Models<-SSgetoutput(dirvec=here(dirs[i],'retrospectives',folders[j]))

      Assess_endyr= Models$replist1$Retro_year
      yrs=Models$replist1$timeseries$Yr[Models$replist1$timeseries$Yr<=Models$replist1$Retro_year]
      ssb_mat = Models$replist1$derived_quants[grep(pattern = "SSB_",Models$replist1$derived_quants$Label),]
      ssb_stdev = as.numeric(ssb_mat[paste0("SSB_",yrs),"StdDev"])

      mod_op <- data.frame(mod = unlist(strsplit(dirs[i],split="/"))[8],
                           peel = unlist(strsplit(folders[j],split="retro"))[2],
                           year = Models$replist1$timeseries$Yr[Models$replist1$timeseries$Yr<=Models$replist1$Retro_year],
                           SSB = as.numeric(Models$replist1$timeseries$SpawnBio[Models$replist1$timeseries$Yr<=Models$replist1$Retro_year]),
                           SSB_stdev = as.numeric(ssb_stdev))

      write.table(mod_op, file=here('models','retro_summarySS3.csv'), append=TRUE, row.names=FALSE, sep=",",
                  col.names = !file.exists(here('models',"retro_summarySS3.csv")))
    }

  }
#Plot
colors <- palette.Martin()
fn <- "retro_summarySS3.csv"
ssb <- read.csv(here('models',fn),header=TRUE)

  ssb <- ssb %>% dplyr::filter(year >= 1990) %>%
                 dplyr::filter(year<=2019) %>%
                 dplyr::mutate(CV = SSB_stdev/SSB,
                               lci = SSB/exp(2.*sqrt(log(1+(SSB_stdev^2)/(SSB^2)))),
                               uci = SSB*exp(2.*sqrt(log(1+(SSB_stdev^2)/(SSB^2))))) %>%
                 dplyr::mutate(peel = as.numeric(as.character(peel)))

  tmp <- ssb %>% dplyr::group_by(mod,peel) %>%
                 dplyr::summarise(base_year = (max(year)))
  ssb <-full_join(ssb,tmp, by=c('mod','peel'))
  rm(tmp)

  ssb <- ssb %>% dplyr::arrange(desc(peel))
  ssb$base_year <- factor(ssb$base_year, levels = 2019:2007)

  ggplot(data = ssb, aes(x=year, y=SSB, color = factor(base_year))) +
      geom_line() +
      geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=factor(base_year)), alpha=0.15,colour=NA) +
      expand_limits(y=0) +
      facet_wrap(~mod, scales = 'free',ncol=1) +
      scale_color_viridis_d(name="Base year", direction=-1) +
      scale_fill_viridis_d(name = "Base year", direction=-1) +
      labs(y='SSB (10,000 mt)', x='Year')
  ggsave(here('retro_timeseries.png'),width = 8, height = 8, units='in')








