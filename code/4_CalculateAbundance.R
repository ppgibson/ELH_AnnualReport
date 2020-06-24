## CALCULATE ABUNDANCE ##
# And variance

# Requires sourced <func_BootstrapVariance.R> ({original.boot}) and <func_CustomFunctions.R> ({apply.boot}). 

# Read in RDS files of data?
# Add a set.seed for bootstrap
# add a cv column - maybe in the exploratory file?
# User can set a "diagnostics" variable to suppress printing the diagnostics plot 

## ---- 0. Read in data and functions ---- 
# Nonsub data
  nonsub.wk <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "NonSub_weekly.RData", sep="_"), sep=""))
  nonsub.day <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "NonSub_daily.RData", sep="_"), sep=""))

# Sub data    
  if (sub.sample == TRUE) {
    sub.wk <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "Sub_weekly.RData", sep="_"), sep=""))
      sub.wk$season <- "spring"
    sub.day <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "Sub_daily.RData", sep="_"), sep=""))
      sub.day$season <- "spring"
  }

# Functions
  source(paste(dir_code, "func_BootstrapVariance.R", sep=""))
  source(paste(dir_code, "func_GAMInterpolation.R", sep=""))
  
  
## ---- 1. Simple estimate ----
# Combine sub and nonsub data (if applicable)
  if (sub.sample==FALSE) {
    full.bywk <- nonsub.wk  #If no subsampling, then nonsub.data is ready to go
    full.bywk <- rename(full.bywk, abund.full = abundance, catch.full=catch)
  }
  
  if (sub.sample==TRUE) {  
    sub.wk <- rename(sub.wk, abund.sub = abundance)
    full.bywk <- merge(nonsub.wk, select(sub.wk, species, j.wk, catch.exp, abund.sub),
                       by = c("species", "j.wk"),
                       all = TRUE)
    
    full.bywk$catch.exp[is.na(full.bywk$catch.exp)] <- 0
    full.bywk$abund.sub[is.na(full.bywk$abund.sub)] <- 0
    
    # Total catch is nonsub catch + expanded sub catch
    full.bywk$catch.full <- round((full.bywk$catch + full.bywk$catch.exp), 0)
    
    # Total abundance is nonsub abund + sub abund
    full.bywk$abund.full <- round((full.bywk$abundance+ full.bywk$abund.sub), 0)
  }

# Estimate variance using coggins etal formula (though this does not account for subsmp variance)
  simp.bywk <- mutate(full.bywk, 
                        var.simp = ((mark+1)*(catch.full+rec+1)*(mark-rec)*catch.full) / (((rec+1)^2)*(rec+2)) )
  
# Sum up by species/season
  simp.byseas <- simp.bywk %>%
    group_by(myear, site, species, season) %>%
    summarize(est.simp = sum(abund.full),
              var.simp = sum(var.simp),
              se.simp = sqrt(var.simp))
  
# Clean up 
  rm(full.bywk)


## ---- 2. Bootstrap variance ----
# Set.seed currently set in original.boot function
  
# Apply the bootstrap routine to the non-sub data
  nonsub.boot <- apply.boot(nonsub.wk, sub.smp=FALSE)

# Apply the bootstrap routine the subsample data, if applicable
  if (sub.sample==FALSE) {
    full.boot <- nonsub.boot
    full.boot <- rename(full.boot, 
                        est.full = estimate, var.full=var, se.full=sd)
  }
  
  if (sub.sample==TRUE) {
    # Now run boot routine
    sub.boot <- apply.boot(sub.wk, sub.smp=TRUE)
    
    # Combine sub and nonsub values, still assuming sub is relevant
    sub.boot <- rename(sub.boot, est.sub = estimate, var.sub=var)
    full.boot <- merge(nonsub.boot, select(sub.boot, species, season, j.wk, est.sub, var.sub),
                       by = c("species", "season", "j.wk"),
                       all = TRUE)
    full.boot[c("est.sub", "var.sub")][is.na(full.boot[c("est.sub", "var.sub")])] <- 0
    
    full.boot <- mutate(full.boot,
                        est.full = estimate + est.sub,
                        var.full = var + var.sub,
                        se.full = sqrt(var.full))
  }
  
  boot.bywk <- full.boot
  rm(full.boot)

# Now sum by season
  boot.byseas <- boot.bywk[boot.bywk$j.wk!="all",] %>%  #boot.bywk includes totals calculated by boot routine
    group_by(species, season) %>%
    summarize(est.boot = sum(est.full),
              var.boot = sum(var.full),
              se.boot = sqrt(var.boot))
  # Note, could alternatively use CIs for total values as calculated within original.boot
  # For reporting full bootstrapped CI, this wd probably be the preferred approach. But, no sub totals.
  boot.totals <- arrange(filter(boot.bywk, j.wk=="all"), species, season) %>% #cf routine-calculated totals for data checking
    select(species, season, est.full, var.full, se.full) %>%
    mutate(myear = my.cur, site=site.cur)
   # select(boot.totals, species, season, est.full, var.full, se.full)
   # boot.byseas
  
  # Add id variables
  boot.byseas$site <- site.cur
  boot.byseas$myear <- my.cur
  boot.byseas <- select(boot.byseas, myear, site, everything())
  
  # Clean up
  rm(nonsub.boot, sub.boot)
  
  
## ---- 3. GAM interpolation ----
# For the GAM, make sure INT days are NA, we want to interpolate for those days
  nonsub.day$catch[nonsub.day$fish.data==FALSE] <- NA
  
# Test whether a value is already set for "diagnostics"; 
# if not then create it. 
  if (!(exists("diagnostics"))) {  #default behavior is to print diagnostic plots, as in original gam routine.
    diagnostics <- TRUE
  } 
  
# Run the procedure for each species/season combination
  if (sub.sample==FALSE) {
    # print("false")
    gam.out <- apply.gam(data.nonsub=nonsub.day, sub.smp = FALSE, diagnostics=diagnostics)
  } 
  
  if (sub.sample==TRUE) {
    gam.out <- apply.gam(data.nonsub=nonsub.day, sub.smp = TRUE, data.sub=sub.day, diagnostics=diagnostics)
    # print("true")
  }
  
  # Add id variables
  gam.out$site <- site.cur
  gam.out$myear <- my.cur

# Pull out totals  
  gam.byseas <- filter(gam.out, grouping=="season") %>%
    rename(est.gam=estimate, se.gam=se) %>%
    select(myear, site, species, season, everything(), -grouping, -j.wk, -j.day)

# Merge daily estimates with full calendar
  gam.byday <- filter(gam.out, grouping=="daily") %>%
    rename(est.gam=estimate, se.gam=se)
  
  calendar.gam <- merge(nonsub.day, 
                        select(gam.byday, -grouping, -j.wk, -ci.low.log, -ci.high.log),
                        by=c("myear", "site", "species", "season", "j.day"), 
                        all=TRUE)
  
## ---- 4. Calculate median migration dates ---- 
# Combine all the byseas estimates in one table  
  # Bring in expanded SUB day catch, if any
  if (sub.sample==TRUE) {
    sub.day <- mutate(sub.day, catch.exp = round(catch / (pss/p24), 0))
    
    calendar.rev <- merge(calendar.gam, 
                          select(sub.day, species, date, catch.exp),
                          by=c("species", "date"), all=TRUE)
    calendar.rev$catch[!is.na(calendar.rev$catch.exp)] <- calendar.rev$catch.exp[!is.na(calendar.rev$catch.exp)] 
    
  } 
  if (sub.sample==FALSE) {
    calendar.rev <- calendar.gam
  }
  
# Daily simple estimated abundance
  if (abund.est == "Petersen") {
    calendar.rev <- mutate(calendar.rev, 
                           est.simp = round(catch/((rec)/(mark)), 0))
  }
  if (abund.est == "Chapman") {
    calendar.rev <- mutate(calendar.rev, 
                           est.simp = round(catch/((rec+1)/(mark+1)), 0))
  }

# Now calculate median migration date, according to simple catch
  plyr::ddply(calendar.rev, c("species", "season"), 
             calc.med, date.col="date", count.col="catch")
  med.catch <- plyr::ddply(calendar.rev, c("species", "season"), 
           calc.med, date.col="date", count.col="catch") %>%
  rename(catch = V1)
  
# Expanded catch: daily catch expanded by TE estimate
  med.exp.catch <- plyr::ddply(calendar.rev, c("species", "season"), 
           calc.med, date.col="date", count.col="est.simp") %>%
  rename(exp.catch = V1)

# GAM: daily abundance as estimated by the GAM
  med.gam <- plyr::ddply(calendar.rev, c("species", "season"), 
           calc.med, date.col="date", count.col="est.gam") %>%
  rename(gam = V1)

# Combine all the calculated medians
  median.dates <- merge(med.catch, med.exp.catch,
                        by=c("species", "season"), all=TRUE)
  median.dates <- merge(median.dates, med.gam, 
                        by=c("species", "season"), all=TRUE)
  rm(med.catch, med.exp.catch, calendar.rev)

  
## ---- 5. Prepare output ---- 
# Combine all the byseas estimates in one table  
  all.byseas <- merge(select(simp.byseas, -var.simp), 
                      select(boot.byseas, -var.boot),
                      by=c("myear", "site", "species", "season"),
                      all=TRUE)
  all.byseas <- merge(all.byseas, gam.byseas,
                      by=c("myear", "site", "species", "season"),
                      all=TRUE)
  # Annotate with estimator used in the analysis
  all.byseas$abund.method <- abund.est
  
# Bring in selected median date value
  all.byseas <- merge(all.byseas, 
                      median.dates[, c("species", "season", med.date.est)],
                      by=c("species", "season"), all=TRUE)
  colnames(all.byseas)[which(colnames(all.byseas)==med.date.est)] <- "med.mig.date"
  all.byseas$med.method <- med.date.est
  
# Rearrange columns
  all.byseas <- select(all.byseas, myear, site, species, season, everything())

# Write output of seasonal totals
  write.csv(all.byseas, paste(dir_out, paste(my.txt, site.cur, "AbundanceEstimates.csv", sep="_"), sep=""), 
            row.names=FALSE)
  
# Write output of daily gam abundance estimates
  calendar.gam <- arrange(calendar.gam, myear, site, species, j.day)
  write.csv(calendar.gam, paste(dir_out, paste(my.txt, site.cur, "AbundanceCalendar.csv", sep="_"), sep=""), 
            row.names=FALSE)

## END ##  