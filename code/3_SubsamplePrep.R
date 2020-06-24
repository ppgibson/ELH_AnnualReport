## PROCESS SUBSAMPLE DATA FOR CALCULATING ABUNDANCE ##


# Read and clean data
  sub.raw <- read.csv(paste(dir_data, paste(my.txt, site.cur, "subsmp.csv", sep="_"), sep=""))  #Counts as prepared by BG, as I have no independent source to use for the hr-by-hr counts. 
  
  colnames(sub.raw) <- tolower(colnames(sub.raw))
  sub.raw$date <- ymd(sub.raw$date)

# Which CAL days should be used to estimate proportions for each SUB day?
# Apply standard routine, then review results.
  sub.cal <- sub.raw %>%
    group_by(myear, site, species) %>%
    do(assign.cal.days(.))
  
# Apply any CAL mods here
  # None for this data set

# Calculate PSS/P24 for each j.week, using assigned CAL days
  # Long form of by-hr ct data
  by.hr <- melt(select(sub.raw, -temp, -gage),
                id.vars = c("myear", "site", "species", "j.wk", "date", "mode"), 
                variable.name="period", value.name="count")
  
  # Counts by week
  # Add a grouping variable for sub smp ("pss") vs pre/post only ("cal")
  by.hr$sub.grp <- NA
    by.hr$sub.grp[by.hr$period=="pre" | by.hr$period =="post"] <- "cal"
    by.hr$sub.grp[by.hr$period %in% c("h1", "h2", "h3", "h4")] <- "pss"
  
  # Sum PSS/P24 catch counts by day
  by.day <- by.hr %>%
    group_by(myear, site, species, j.wk, date, mode, sub.grp) %>%
    summarize(catch = sum(count, na.rm=TRUE))
  
    # Cast back out: get PSS and P24 for CAL days, aka catch.ss and catch.24
    sub.cts.day <- dcast(by.day,  myear+site+species+j.wk+date+mode ~ sub.grp, value.var="catch")
    sub.cts.day <- mutate(sub.cts.day, 
                          p24 = cal + pss)

  # Total SUB catch, by week
  sub <- sub.cts.day[sub.cts.day$mode=="SUB", ] %>%
    group_by(myear, site, species, j.wk) %>%
    summarize(catch = sum(p24))
  
  # Initialize fields to be filled by aggregating CAL data
  sub$p24 <- NA
  sub$pss <- NA
  
  # Bring in associated CAL days, as identified by func above
  # Currently this routine works only for a single site/my at a time
  for (i in 1:nrow(sub)) {  #loop through each week with SUB days
    # Identify associated cal column in sub.cal, and other current variables
    wk.cur <- sub$j.wk[i]
    sp.cur <- sub$species[i]
    col.cur <- which(colnames(sub.cal) == paste("cal.", wk.cur, sep=""))  # name of current TRUE/FALSE column
    dates.cur <- unique(sub.cal$date[(sub.cal[, col.cur])==TRUE]) # CAL dates to apply to current j.wk of SUB
    
    # Extract relevant records from the daily cal/sub counts --> weekly counts
    cal.cur <- filter(sub.cts.day, date %in% dates.cur, species==sp.cur)
    cal.sum.cur <- cal.cur %>%
      group_by(myear, site, species) %>%
      summarize(p24 = sum(p24),
                pss = sum(pss) )
    
    # Store these values with their appropriate SUB j.wk
    sub$pss[i] <- cal.sum.cur$pss
    sub$p24[i] <- cal.sum.cur$p24
  }

# Expand SUB numbers using CAL proportion
  sub <- mutate(sub,
                ss.prop = pss/p24,
                catch.exp = catch/ss.prop)  #Expanded catch (note, not an abundance estimate yet).

# Bring in stratified mark/recap numbers from trap log/non sub data
  nonsub.data <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "NonSub_weekly.RData", sep="_"), sep=""))
  sub.data <- merge(sub, select(nonsub.data, -strat.id, -season, -catch),
               by = c("myear", "site", "species", "j.wk"),
               all.x = TRUE, all.y=FALSE)
  sub.data <- mutate(sub.data, abundance = round(catch.exp/te, 0)) #This estimated abundance is based on expanded catch.

# Prep sub data for GAM: p24/pss numbers attached to each sub date/jday
  gam.data.sub <- sub.cts.day %>%
    filter(mode=="SUB") %>%
    select(myear, site, species, j.wk, date, catch=p24)  #don't actually need any catch data at this point,but going to keep it anyway.
  # Add weekly pss/p24 cts for prop
  gam.data.sub <- merge(gam.data.sub, select(sub.data, species, j.wk, p24, pss),
                        by=c("species", "j.wk"), all=TRUE)
  gam.data.sub$j.day <- yday(gam.data.sub$date)
  gam.data.sub <- gam.data.sub %>%
    select(myear, site, species, date, j.wk, j.day, everything()) %>%
    arrange(species, date)

# Write Output files
  # 1. Basic subsampling counts by TE week, formatted for use in bootstrap routine
    saveRDS(sub.data, paste(dir_drvd, paste(my.txt, site.cur, "Sub_weekly.RData", sep="_"), sep=""))
    write.csv(sub.data, paste(dir_drvd, paste(my.txt, site.cur, "Sub_weekly.csv", sep="_"), sep=""), row.names=FALSE)
    
  # 2. Simple table for GAM (SUB dates plus associated CAL proportions)
    saveRDS(gam.data.sub, paste(dir_drvd, paste(my.txt, site.cur, "Sub_daily.RData", sep="_"), sep=""))

## END ##  