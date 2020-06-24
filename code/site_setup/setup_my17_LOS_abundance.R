## MY 17 LOSTINE 
## ANALYSIS OPTIONS FOR ABUNDANCE CALCULATION 


# Current analysis run
  my.cur <- 2017
  my.txt <- "my17"
  site.cur <- "LOS"

# Analysis options
  include.INT.catch <- TRUE #Default = TRUE. Std ELH policy is to include catch numbers from INT days in simple abundance estimates.
  use.prv.strata <- FALSE   #Default = FALSE. Use standard, algorithm-based stratification (FALSE); or use BG's, previous stratification (TRUE).
  sub.sample <- TRUE        #Was there any subsampling at this site/year?
  abund.est <- "Petersen"   #Petersen: Nhat = catch * (mark/rec)  Simplest abundance estimate, as described in annual report and Thedinga etal 1994.
                            #Chapman:  Nhat = catch * (mark+1)/(rec+1)  Apparently corrects for bias in small sample size, esp of rec. 
  med.date.est <- "exp.catch"  #Determines which method for calculating median date to report. {catch, exp.catch, gam} 
    
# Custom mods
  # Usable INT days
  INT.mods <- c("2016-10-24")  #Vector of INT dates determined, based on subjective assessment, to be usable for mark-recap anl.

  # Adjustments to stratification
  strat.mods <- list(chs.c.wks = c(39, 40, 41),  #J.wks manually assigned to stratum "c" for chs
                     chs.e.wks = c(42, 43, 44) )
  
  # Adjustments to CAL day assignment
    # none for this data set
  
# Run standard Setup and functions
  source("code//0_Setup.R")
  source("code//func_CustomFunctions.R")

  
## END ##