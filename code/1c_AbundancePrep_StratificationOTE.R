## PROCESS TRAP LOG DATA FOR CALCULATING ABUNDANCE: PART 1C ##
#
# Directly follows 1b_AbundancePrep_CMR.R
#
# Apply custom mods to stratification, then calculate CMR by
# stratum. Bring in data on OTE fish and remove them from counts
# where appropriate. 


## ---- 4. Stratification, cont. ----
# Apply any custom mods to stratification, as given in setup script
for (i in names(strat.mods)) {
  #print(i)
  sp.cur <- substr(i, start=1, stop=3)
  stratum.cur <- substr(i, start=5, stop=5)
  #print(strat.mods[i])
  strat.recs.cur <- which(fish.strata$species==sp.cur & fish.strata$j.wk %in% strat.mods[[i]])
  fish.strata$stratum[strat.recs.cur] <- stratum.cur 
}
  if(exists("strat.recs.cur")) {
    rm(i, sp.cur, stratum.cur)
  }

# After custom mods (if any), use fish.strata.revised  
  fish.strata.rev <- fish.strata

# Unique stratum ids
  fish.strata.rev <- mutate(fish.strata.rev, strat.id = paste(species, season, stratum, sep="."))

# Calculate mark/recap totals and TE by stratum
  strata.sums <- fish.strata.rev %>%
    group_by(myear, site, species, season, strat.id) %>%
    summarize(catch = sum(catch),
              mark = sum(mark),
              rec = sum(rec),
              first.wk = min(j.wk),
              last.wk  = max(j.wk))
  # Two alternative ways to estimate TE; option is set in RMD script
  if (abund.est == "Petersen") {
    strata.sums <- mutate(strata.sums, te = (rec)/(mark))
  }
  if (abund.est == "Chapman") {
    strata.sums <- mutate(strata.sums, te = (rec+1)/(mark+1))
  }

## ---- 5. OTE fish ----
# Read in manually compiled table of OTE fish
  ote <- read.csv(paste(dir_data, paste(my.txt, site.cur, "ote.csv", sep="_"), sep=""))
  ote <- ote[!is.na(ote$myear), ]

# # Which OTE fish were recaptured in the first week of a new stratum?
# # This presumably indicateS that they were released during the previous stratum.
# # Fish that are released and recaptured within the same stratum can be counted as normal.
#   ote.str <- merge(ote, select(strata.sums, myear, site, species, first.wk, strat.id),
#                    by.x = c("myear", "site", "species", "j.wk"),
#                    by.y = c("myear", "site", "species", "first.wk"),
#                    all.x=TRUE, all.y=FALSE)
#   rm(ote)
  
# Find the stratum for each release and recap week. 
# Fish that are released and recaptured within the same stratum can be counted as normal,
# while fish that are released and recaptured in different strata are excluded from 
# both counts as OTE.
  
  # Stratum for the OTE fish's recap week (ie, [j.wk])
  ote.str <- merge(ote, select(fish.strata.rev, myear, site, species, j.wk, strat.id),
                   by = c("myear", "site", "species", "j.wk"),
                   all.x=TRUE, all.y=FALSE) %>%
    rename(strat.rec = strat.id)
  
  # Stratum for the OTE fish's presumed release week (ie, [est.rel.wk])
  ote.str <- merge(ote.str, select(fish.strata.rev, myear, site, species, j.wk, strat.id),
                   by.x = c("myear", "site", "species", "est.rel.wk"), 
                   by.y = c("myear", "site", "species", "j.wk"), 
                   all.x=TRUE, all.y=FALSE) %>%
    rename(strat.rel = strat.id)
  
  # Do the release and recap strata match?
  ote.str$ote <- !(ote.str$strat.rel == ote.str$strat.rec)  #TRUE=OTE=exclude; FALSE=not OTE=don't exclude.
  
  # Sum confirmed OTE *recaps* by stratum
  ote.rec <- ote.str[ote.str$ote==TRUE, ] %>%
    group_by(strat.rec) %>%
    summarize(count = sum(count))
  # Sum confirmed OTE *releases* by stratum
  ote.rel <- ote.str[ote.str$ote==TRUE, ] %>%
    group_by(strat.rel) %>%
    summarize(count = sum(count))

  rm(ote)

# # Take these confirmed OTE fish and remove them as appropriate from 
# # mark and recap counts in strata.sums.
#   strata.sums.rev <- merge(strata.sums, ote.str[ote.str$ote==TRUE, c("strat.rec", "count")], #only include the fish identified as crossing a stratum boundary
#                 by.x="strat.rec", by.y="strat.id", all=TRUE)
#   strata.sums.rev$count[is.na(strata.sums.rev$count)] <- 0
#   
#   # Revise recap numbers (subtract the OTE count for that stratum)  
#   strata.sums.rev$rec <- strata.sums.rev$rec - strata.sums.rev$count
#   
#   # Revise mark numbers (subtract OTE count from mark # in prv stratum)
#   for (i in 1:(nrow(strata.sums.rev)-1)) {
#     strata.sums.rev$mark[i] <- strata.sums.rev$mark[i] - strata.sums.rev$count[i+1]
#   }

# Take these confirmed OTE fish and remove them as appropriate from
# mark and recap counts in strata.sums.
  # Mark (release) counts
  strata.sums.rev <- merge(strata.sums, ote.rel,  
                by.x="strat.id", by.y="strat.rel", all=TRUE)
  strata.sums.rev$count[is.na(strata.sums.rev$count)] <- 0
  strata.sums.rev$mark <- strata.sums.rev$mark - strata.sums.rev$count
  strata.sums.rev <- select(strata.sums.rev, -count)

  # Recap counts
  strata.sums.rev <- merge(strata.sums.rev, ote.rec,  
                by.x="strat.id", by.y="strat.rec", all=TRUE)
  strata.sums.rev$count[is.na(strata.sums.rev$count)] <- 0
  strata.sums.rev$rec <- strata.sums.rev$rec - strata.sums.rev$count
  strata.sums.rev <- select(strata.sums.rev, -count)

# Recalculate seasonal counts
  season.cts.rev <- strata.sums.rev %>%
    group_by(myear, site, species, season) %>%
    summarize(catch = sum(catch),
              mark = sum(mark),
              rec = sum(rec))
  season.cts.rev <- mutate(season.cts.rev, te.all = round(rec/mark,3))
  

# Make a weekly table, formatted for the bootstrap routine
  nonsub.data <- merge(fish.strata.rev[, c("myear", "site", "species", "season", "j.wk", "catch", "strat.id")],
                       strata.sums.rev[, c("myear", "site", "strat.id", "mark", "rec", "te")],
                       by=c("myear", "site", "strat.id"), all=TRUE)
  nonsub.data <- mutate(nonsub.data, 
                        abundance=catch/te )  #Peterson vs Chapman method has already been applied, through calculation of TE.
  
## ---- 6. OUtput files ----  
# Write output files
  # 1. Complete trap log, with annotations
  log.out <- select(log.cmr, 
                    myear, site, season, date, j.wk, jday, trapmode, startstop, 
                    te.re.wk, revolutions, time, gage, temp, 
                    chs.new, chs.catch, chs.mark, chs.rec,
                    sts.new, sts.catch, sts.mark, sts.rec, 
                    catch.data.good, rec.data.good, mark.data.good,
                    comments) #everything(), to include all remaining columns after
  saveRDS(log.out, paste(dir_drvd, paste(my.txt, site.cur, "LogProcessed.RData", sep="_"), sep=""))
  write.csv(log.out, paste(dir_drvd, paste(my.txt, site.cur, "LogProcessed.csv", sep="_"), sep=""), row.names=FALSE)

  # 2. Weekly CMR counts, with strata
  saveRDS(nonsub.data, paste(dir_drvd, paste(my.txt, site.cur, "NonSub_weekly.RData", sep="_"), sep=""))
  write.csv(nonsub.data, paste(dir_drvd, paste(my.txt, site.cur, "NonSub_weekly.csv", sep="_"), sep=""), row.names=FALSE)

## END ##    