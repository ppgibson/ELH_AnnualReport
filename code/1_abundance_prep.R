## DATA PROCESSING OF TRAP LOG DATA

# Set libraries, directories, and options
  source("code//0_Setup.R")
  source("code//func_CustomFunctions.R")

# Variables/options for current model run
  include.INT.catch <- TRUE  #looks like it is std ELH policy to always include straight catch data from INT days.
  use.prv.strata <- FALSE 
  # TRUE = use stratification from BG's previous run, as imported via <stratification> .csv file.
  # FALSE = use standadized stratification protocol, of combine until >= 10 rec per stratum.
  
  # *Outstanding issues with LOS my17: 
  #  - unk status on March 16 - UPDATE: status should be INT, trap was started as SUB[?], then shut 
  #    down after hour 3[?] due to being overwhelmed. 
  #  - account for any TE fish
  
# There are no OTE fish noted in trap log for Chinook.
# Three 

# Read in data
log.raw <- read.csv(paste(dir_data, "LOS_my17_traplog.csv", sep=""), stringsAsFactors=FALSE)

# Adjust dates
log <- log.raw
  log$date <- ymd(log$date)
  log$jday <- yday(log$date)  
  
# Add season
# 28 Jan is end of fall season
  yday("2017-01-28")
  yday("2016-09-01") 
log$season[log$jday<=28 | log$jday>=245] <- "fall"
log$season[log$jday>28 & log$jday<245] <- "spring"
  
# # Generate a calendar table
#   range(log$date)
#   calendar <- data.frame(date = seq(from=min(log$date), to=max(log$date), by=1)) #Sequence of dates covering the range of the current trap log data
#   calendar <- mutate(calendar, 
#                      dow = wday(date, label=TRUE, abbr=TRUE),
#                      j.day = yday(date),
#                      j.week = ceiling(j.day/7))
#   days.in.log <- nrow(calendar)  #Usually but not necessarily 365. Might be longer due to leap year or shorter due to trap log entries. 

# Create a version of log.te that doesn't have all the fish columns in the way 
  fish.cols <- colnames(log)[17:57]
  head(select(log, -fish.cols))
  # trap <- log[, c(1:16, 58:60)]
  
  # Examine INT days; are any of them usable?
  table(log$trapmode)
  trap.int <- filter(filter(select(log, season, date, jday, j.wk, gage, trapmode, revolutions, chs.new, sts.new, comments), trapmode=="INT"))
  usable.days <- c("2016-10-24")  #2017-02-16 also possible, but to be conservative I left it out; gage went way up between 2/15 and 2/16, making it harder to gage # revs.
  log$trapmode[log$date==usable.days] <- "CON"
  table(log$trapmode)
  
# Sort trap log according to which dates have usable data (for fish counts) 
# Option to exclude vs include data from days when trapmode=INT  
  # [xx.data.good] indicates usable day of trap log data: 
  log.te <- mutate(log, 
                   catch.data.good = FALSE,
                   rec.data.good = FALSE)  #Default value in all cases is FALSE
  # catch data
    if (include.INT.catch==TRUE){
      log.te$catch.data.good[log.te$trapmode %in% c("CON", "INT", "CAL")] <- TRUE  
    } else {
      log.te$catch.data.good[log.te$trapmode %in% c("CON", "CAL")] <- TRUE  
    }
  
  # recap data
    log.te$rec.data.good[log.te$trapmode %in% c("CON", "CAL")] <- TRUE  #rec data is never valid on SUB and INT days
    
  # release data
    rec.days <- log.te$trapmode %in% c("CON", "CAL")
    mark.data.good <- c(rec.days[2:length(rec.days)], NA)  #Take the vector of eligible recap days, move it up one, and append an NA to the end.
    log.te <- cbind(log.te, mark.data.good)
    
    rm(rec.days, mark.data.good)

  # Are there any instances of [mark.data.good]=TRUE when [startstop]=STP? (should be none)
  table(log.te$startstop, log.te$mark.data.good)
  filter(select(log.te, -fish.cols), mark.data.good==TRUE, startstop=="STP") #this one is okay because next day is CAL.
  
# For now, lets ignore OTE and matching fin clips. 
# These could be removed from the trap log counts - ?
# Create a df with capture/mark/recapture numbers for each species
  # Remove int days (if specified)
  log.te <- mutate(log.te, 
                   chs.catch  = chs.new,
                   sts.catch  = sts.new,
                   chs.mark = chs.te.pit + chs.te.clip,  #this method avoids addition errors - seems reasonably common
                   sts.mark = sts.te.pit + sts.te.clip,
                   chs.rec  = chs.re.pit + chs.re.clip,  #this method avoids addition errors, and preemptive OTE removals
                   sts.rec  = sts.re.pit + sts.re.clip)
  log.te$chs.catch[log.te$catch.data.good==FALSE] <- NA
  log.te$sts.catch[log.te$catch.data.good==FALSE] <- NA
  log.te$chs.mark[log.te$mark.data.good==FALSE] <- NA
  log.te$sts.mark[log.te$mark.data.good==FALSE] <- NA
  log.te$chs.rec [log.te$rec.data.good==FALSE] <- NA
  log.te$sts.rec [log.te$rec.data.good==FALSE] <- NA
  
# Pull out some relevant columns  
  cmr.daily <- select(log.te, 
                      myear, site, season, date, jday, j.wk, te.re.wk, trapmode,
                      chs.catch, chs.rec, chs.mark, sts.catch, sts.mark, sts.rec)
  cmr.long <- melt(cmr.daily, id.vars=c("site", "myear", "season", "date", "jday", "j.wk", "te.re.wk","trapmode"),
                   variable.name="type", value.name="count")
  cmr.long$species <- substr(cmr.long$type, start=1, stop=3)
  cmr.long$event   <- substr(cmr.long$type, start=5, stop=nchar(as.character(cmr.long$type)))
  
# Calculate weekly counts
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


  rm(cmr.long)
  rm(recs)
  rm(r.weekly)
  
  # Total number of CMR, by species/season
  # season.cts <- cmr.counts %>%
  #   group_by(myear, site, species, season) %>%
  #   summarize(catch = sum(catch),
  #             mark = sum(mark),
  #             rec = sum(rec))
  # season.cts <- mutate(season.cts, te.all = round(rec/mark,3))
  
  # log.te[log.te$te.rel.good==FALSE, c("date", "season", "trapmode", "chs.te.all")]->chs.mark
  # log.te[log.te$te.rel.good==FALSE, c("date", "season", "trapmode", "sts.te.all")]->sts.mark
  # chs.mark <- chs.mark[!is.na(chs.mark$chs.te.all),]
  # sts.mark <- sts.mark[!is.na(sts.mark$sts.te.all),]
  # 
  # log.te[log.te$te.re.wk==39, c("j.wk", "te.re.wk", "date", "trapmode", "sts.new", "sts.te.all", "sts.re.all", "sts.rec")]
  
#### STRATIFICATION ####
# Bring standardized stratification function, which sets strata based 
# only on minimum of 10 recaps per stratum
  # source(paste(dir_code, "func_SetStrata.R", sep=""))
  # source(paste(dir_code, "func_setCalGroups.R", sep=""))
  
# *Initial accounting for TE here?
# *add initial variable controlling which strat to use, or modifications to custom strat
  
# Stratification used in previous analysis
  strat.prv <- read.csv(paste(dir_data, "LOS_my17_stratification.csv", sep=""))

  # First generate std stratification
  fish.strata <- cmr.counts %>%
    group_by(species, season) %>%
    do(set.strata(.))
  fish.strata <- rename(fish.strata, stratum.std=stratum)

  # Merge in prv stratification for comparison
  fish.strata <- merge(fish.strata, strat.prv,
                       by=c("myear", "site", "species", "season", "j.wk"),
                       all=TRUE)
  
  # Select a stratification to use, and modify if appropriate
  if (use.prv.strata==TRUE) {
    # Use BG's previous stratification
    strat.col.cur <- which(colnames(fish.strata)=="stratum.prv")
  } else {
    strat.col.cur <- which(colnames(fish.strata)=="stratum.std")
  } 
  fish.strata$stratum <- fish.strata[, strat.col.cur] 
  ##**Custom mods here
    fish.strata$stratum[fish.strata$j.wk %in% c(40,41) & fish.strata$species=="chs"] <- "c" #*custom mod
    fish.strata$stratum[fish.strata$j.wk %in% c(43,44) & fish.strata$species=="chs"] <- "e" #*custom mod
  ##**

  # Unique stratum ids
  fish.strata <- mutate(fish.strata, strat.id = paste(species, season, stratum, sep="."))

# Calculate mark/recap totals and TE by stratum
  strata.sums <- fish.strata %>%
    group_by(myear, site, species, strat.id) %>%
    summarize(catch = sum(catch),
              mark = sum(mark),
              rec = sum(rec),
              te = rec/mark, 
              first.wk = min(j.wk),
              last.wk  = max(j.wk))

  first.te.wk = min(te.week),
              last.te.wk = max(te.week))
  
# Deal with OTE fish
  strata.rev <- strata.sums #or do this above?
  
  ote <- read.csv(paste(dir_data, "LOS_my17_OTE.csv", sep=""))
  
  # Which OTE fish were recaptured in the first week of a new stratum?
  # (This presumably indicateS that they were released during the previous stratum.)
  ote.str <- merge(ote, select(strata.sums, myear, site, species, first.wk, strat.id),
                   by.x = c("myear", "site", "species", "j.wk"),
                   by.y = c("myear", "site", "species", "first.wk"),
                   all.x=TRUE, all.y=FALSE)
  
  # Take these fish and append the relevant fields to stratasums, for removing from 
  # mark and rec counts
  test <- merge(strata.rev, ote.str[!is.na(ote.str$strat.id), c("strat.id", "count")],
                by="strat.id", all=TRUE)
  test$count[is.na(test$count)] <- 0
  # Revise recap numbers, subtract the OTE count for that stratum
  test$mark.rev <- NA
  test$rec.rev <- NA
  for (i in 2:nrow(test)) {
    test$rec[i] <- test$rec[i] - test$count[i]
    test$mark[i] <- test$mark[i] - test$count[i+1]
  }
  test$mark.rev[1] <- test$mark[1]
  test$rec.rev[1] <- test$rec[1]
  test$mark.rev[nrow(test)] <- test$mark[nrow(test)]
  test <- mutate(test, 
                 rec.rev = rec-count,
                 mark.rev = mark-())
  
  
  
# Make a weekly table for the boot routine
  nonsub.data <- merge(fish.strata[, c("myear", "site", "species", "season", "j.wk", "catch", "strat.id")],
                       strata.sums[, c("myear", "site", "strat.id", "mark", "rec", "te")],
                       by=c("myear", "site", "strat.id"), all=TRUE)
  nonsub.data <- mutate(nonsub.data, abundance=catch/te)
  
#### SUBSAMPLING ####
# How to choose appropriate weeks of CAL for each sub?
  sub.raw <- read.csv(paste(dir_data, "LOS_my17_subsmp.csv", sep=""))
  
  colnames(sub.raw) <- tolower(colnames(sub.raw))
  sub.raw$date <- ymd(sub.raw$date)
  
  sub.cal <- sub.raw %>%
    group_by(myear, site, species) %>%
    do(assign.cal.groups(.))
  

# Using assigned CAL groups, sum PSS/P24 for each j.week
  # long form of by-hr ct data
  by.hr <- melt(select(sub.cal, -temp, -gage, -group),
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
                          p24 = cal + pss,
                          prop = pss/p24)
    sub.cts.day$p24 <- sub.cts.day$cal + sub.cts.day$pss

  # Total SUB counts, by week
  sub <- sub.cts.day[sub.cts.day$mode=="SUB", ] %>%
    group_by(myear, site, species, j.wk) %>%
    summarize(catch = sum(p24))

  # IN preparation, add empty cols to sub, to be filled
  sub$pss <- NA
  sub$p24 <- NA
  
  # Bring in associated CAL days, as identified by func above
    # Currently only for a single site/my/species at a time
  for (i in 1:nrow(sub)) {  #loop through each week with SUB days
    # identify associated cal column in sub.cal
    wk.cur <- sub$j.wk[i]
    sp.cur <- sub$species[i]
    col.cur <- which(colnames(sub.cal) == paste("cal.", wk.cur, sep=""))  # name of current TRUE/FALSE column
    dates.cur <- unique(sub.cal$date[(sub.cal[, col.cur])==TRUE]) # CAL dates to apply to current j.wk of SUB
    cal.cur <- filter(sub.cts.day, date %in% dates.cur, species==sp.cur)
    cal.sum.cur <- cal.cur %>%
      group_by(myear, site, species) %>%
      summarize(pss = sum(pss),
                p24 = sum(p24))
    
    sub$pss[i] <- cal.sum.cur$pss
    sub$p24[i] <- cal.sum.cur$p24
  }
  
  # Now expand sub numbers
  sub <- mutate(sub,
                ss.prop = pss/p24,
                catch.exp = catch/ss.prop)
  
  # Bring in stratified mark/recap numbers
  sub.data <- merge(sub, select(nonsub.data, -strat.id, -season, -catch),
               by = c("myear", "site", "species", "j.wk"),
               all.x = TRUE, all.y=FALSE)
  sub.data <- mutate(sub.data, abundance = round(catch.exp/te, 0))
  
  
  
# test boot
# Function to take two input data tables (nonsub and sub data), apply orig.boot 
# to the appropriate subsets, and report results with appropriate indexing factors
run.boot.abund <- function(std.data, sub.data=NULL, sub.sample=F)  
  
  source
  source("code//func_BootstrapVariance.R")
  
  boot.chs <- original.boot(data = sub[sub$species=="chs", ],
                            time ="week", catch="catch", mark="mark",recap="rec",
                         sub.sample = T, P24 = "p24", PSS = "pss")
  boot.sts.nosub <- original.boot(data=weekly.data[weekly.data$species=="chs", ], 
                                  time = "j.wk", catch="cat", mark="mark", recap="rec", sub.sample=F)
  
  test <- weekly.data %>%
    group_by(myear, site, species, season) %>%
    do(original.boot(data=., time = "j.wk", catch="cat", mark="mark", recap="rec", sub.sample=F))

  test <- plyr::ddply(.data=weekly.data, .variables=c("species", "season"), .fun=original.boot, time = "j.wk", catch="cat", mark="mark", recap="rec", sub.sample=F)
#   # Cast back out: get PSS and P24 for CAL days
#   pss <- dcast(pss.wk[pss.wk$mode=="CAL" & pss.wk,],  site+my+species+week+mode ~ sub.grp)
#   pss <- mutate(pss, 
#                 p24 = cal + pss,
#                 prop = pss/p24)
#   pss$p24 <- pss$cal + pss$pss
#   # bring in sub catch numbers
#   pss <- merge(pss, pss.wk[pss.wk$mode=="SUB" & pss.wk$sub.grp=="pss", ],
#                by = c("site", "my", "species", "week"), all.x=TRUE, all.y=FALSE)
#   
#   # Pull out a clean table toward use in bootstrap
#   sub.table <- select(pss, site, my, species, week, catch, pss, p24, prop)
#   sub.table <- filter(sub.table, !is.na(catch))
#   
#   # merge with mark/recap data from cmr
#   test <- merge(sub.table, fish.strata[, c("site", "myear", "species", "j.wk", "mar", "rec")],
#                 by.x = c("site", "my", "species", "week"),
#                 by.y = c("site", "myear", "species", "j.wk"),
#                 all.x = TRUE, all.y=FALSE)
# # Goal output
#   
  
