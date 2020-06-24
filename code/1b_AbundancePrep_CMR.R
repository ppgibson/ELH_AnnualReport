## PROCESS TRAP LOG DATA FOR CALCULATING ABUNDANCE: PART 1B ##
#
# Directly follows 1a_AbundancePrep_INT.R
#
# Apply INT custom mods, then calculate cmr counts by TE week.
# Apply standard stratification, by species and season. Read
# in previous stratification and merge it for comparison. 
# Set official stratification, based on option from setup. 

## ---- 2. Dates with usable data, cont. ----
# Apply any custom mods to INT days
for (i in INT.mods) {  #Loop through list of dates, as formally set via "INT.mods" at top of script.
  # print(i)
  row.cur <- which(log.te$date==i)
  log.te$rec.data.good[row.cur] <- TRUE     #This date's recap numbers are usable...
  log.te$mark.data.good[row.cur-1] <- TRUE  #...as are release numbers from the prv date.
}
rm(row.cur, i)


## ---- 3. CMR counts ----
# Now calculate CMR counts by day/week
# Separate/renamed data frame    
  log.cmr <- mutate(log.te, 
                   chs.catch  = chs.new,
                   sts.catch  = sts.new,
                   chs.mark = chs.te.pit + chs.te.clip,  #This method avoids addition errors in trap log's "te.all" fields (these errors seem reasonably common)
                   sts.mark = sts.te.pit + sts.te.clip,
                   chs.rec  = chs.re.pit + chs.re.clip,  #This method avoids addition errors, and preemptive OTE removals
                   sts.rec  = sts.re.pit + sts.re.clip)
  log.cmr$chs.catch[log.cmr$catch.data.good==FALSE] <- NA
  log.cmr$sts.catch[log.cmr$catch.data.good==FALSE] <- NA
  log.cmr$chs.mark[log.cmr$mark.data.good==FALSE] <- NA
  log.cmr$sts.mark[log.cmr$mark.data.good==FALSE] <- NA
  log.cmr$chs.rec [log.cmr$rec.data.good==FALSE] <- NA
  log.cmr$sts.rec [log.cmr$rec.data.good==FALSE] <- NA
  
# Export a data set for later use
  # saveRDS(log.cmr, paste(dir_drvd, my.txt, site.cur, "logCMR.RData", sep="."))
  
# Reformat to long form, split out variables  
  cmr.daily <- select(log.cmr, 
                      myear, site, season, date, jday, j.wk, te.re.wk, trapmode,
                      chs.catch, chs.rec, chs.mark, sts.catch, sts.mark, sts.rec)
  cmr.long <- melt(cmr.daily, id.vars=c("site", "myear", "season", "date", "jday", "j.wk", "te.re.wk","trapmode"),
                   variable.name="type", value.name="count")
  cmr.long$species <- substr(cmr.long$type, start=1, stop=3)
  cmr.long$event   <- substr(cmr.long$type, start=5, stop=nchar(as.character(cmr.long$type)))
  
# Calculate weekly counts
# *This assumes trap log j.wks are correct; might want to revisit this assmption.  
  cm.weekly <- cmr.long[cmr.long$event!="rec", ] %>%  
    group_by(site, myear, season, species, j.wk, event) %>%
    summarize(tot.wk = sum(count, na.rm=TRUE),
              n.data = sum(!is.na(count)))
  # Recast out
  cmr.counts <- dcast(data=cm.weekly, myear+site+species+season+j.wk ~ event, value.var="tot.wk")

  # Calculate recap totals separately, by te.re.wk 
  recs <- filter(cmr.long, event=="rec")
  r.weekly <- recs %>%  
    group_by(myear, site, species, season, te.re.wk, event) %>%
    summarize(tot.wk = sum(count, na.rm=TRUE),
              n.data = sum(!is.na(count)))
  r.weekly <- rename(r.weekly, rec=tot.wk)

  # Recombine recaps with cm data, by j.wk
  cmr.counts <- merge(cmr.counts, select(r.weekly, -event, -n.data), 
                      by.x=c("myear", "site", "species", "season", "j.wk"),
                      by.y=c("myear", "site", "species", "season", "te.re.wk"), 
                      all=TRUE)
  cmr.counts <- arrange(cmr.counts, species, season, j.wk)

  # Clean up
  rm(cmr.long); rm(recs); rm(r.weekly); rm(cm.weekly)


## ---- 4. Stratification ----
# Calculate standard stratification, using basic rule of >= 10 rec per stratum.
  fish.strata <- cmr.counts %>%
    group_by(species, season) %>%
    do(set.strata(.))
  fish.strata <- rename(fish.strata, stratum.std=stratum)

# Compare stratification used in previous analysis
  # Read in file
  strat.prv <- read.csv(paste(dir_data, paste(my.txt, site.cur, "stratification.csv", sep="_"), sep=""))
  
  # Merge both versions
  fish.strata <- merge(fish.strata, strat.prv,
                       by=c("myear", "site", "species", "season", "j.wk"),
                       all=TRUE)
  rm(strat.prv)
  
# Select a stratification to use
  if (use.prv.strata==TRUE) {
    # Use BG's previous stratification
    strat.col.cur <- which(colnames(fish.strata)=="stratum.prv")
  } else {
    strat.col.cur <- which(colnames(fish.strata)=="stratum.std")
  } 
  fish.strata$stratum <- fish.strata[, strat.col.cur]
  rm(strat.col.cur)
  
# Next, apply any custom mods for stratification
  