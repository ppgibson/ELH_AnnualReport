## 2. CREATE A FULLY DAILY CALENDAR WITH CMR ##
#
# Creates a complete daily calendar (no gaps in the date) 
# sequence and merges it with daily trap log information (trap mode, catch) 
# and with stratified mark and recap numbers. Formatted for use in the 
# GAM routine.
# Inputs: processed trap log and stratification/mark/recap by week and by stratum; all from part 1 scripts.
# Outputs: <NonSub_daily> .R and .csv files.

# Read in pre-processed data
  # Full log
  log <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "LogProcessed.RData", sep="_"), sep=""))
  log <- rename(log, j.wk.log = j.wk)

  # Stratification and associated mark/recap/te values
  fish.strata <- readRDS(paste(dir_drvd, paste(my.txt, site.cur, "NonSub_weekly.RData", sep="_"), sep=""))
  

# Generate a continuous calendar table
  calendar <- data.frame(date = seq(from=min(log$date), to=max(log$date), by=1)) #Sequence of dates covering the range of the current trap log data
  calendar <- mutate(calendar,
                     dow = wday(date, label=TRUE, abbr=TRUE),
                     j.day = yday(date),
                     j.wk = ceiling(j.day/7))

# Merge in relevant data from trap log
  cal2 <- merge(calendar,
                select(log, myear, site, season, date, j.wk.log, trapmode, startstop, 
                       chs.new, sts.new,  #Could use .catch instead of .new here, to exclude SUB counts and be more official.
                       catch.data.good, rec.data.good, mark.data.good), 
                by="date", all=TRUE)
  
# Cut out records between fall and spring season
  last.day.fall <- max(cal2$date[cal2$season=="fall" & !is.na(cal2$season)])  
  first.day.spr <- min(cal2$date[cal2$season=="spring" & !is.na(cal2$season)])
  cal2 <- filter(cal2, date<=last.day.fall | date>=first.day.spr)
  
# Simple indicator for dates when trap was collecting fish data (ie, GAM input vs NAs)
# For a more complicated routine to deal with multi-day samples, see LS scripts
  cal2$fish.data <- FALSE
    cal2$fish.data[cal2$rec.data.good=="TRUE"] <- TRUE  #Any day with "good" recap data means the trap was officially fishing for that date.
    cal2$fish.data[cal2$trapmode=="SUB"] <- TRUE  #SUB days also count as fish data, for GAM purposes.

# Fill out indicator variables
  cal2$myear <- my.cur
  cal2$site <- site.cur
    
  fall.st.day <- yday(paste(my.cur, "-07-01", sep=""))
  cal2$season[cal2$j.day<=28 | cal2$j.day>=fall.st.day] <- "fall"
  cal2$season[cal2$j.day>28 & cal2$j.day<fall.st.day] <- "spring"
  rm(fall.st.day)

# Check for and cut out any non-fish days at the beg and end of the sequence (per season)
# These extra days cause problems in the GAM.
  cal3 <- data.frame()
  for (i in c("fall", "spring")) {
    seas.cur <- i
    first.fish.day <- min((cal2$date)[cal2$fish.data==TRUE & cal2$season==seas.cur])  #First day when trap was fishing
    last.fish.day  <- max((cal2$date)[cal2$fish.data==TRUE & cal2$season==seas.cur])  #Last day when trap was fishing and checked.
    cal.cur <- cal2[cal2$date>=first.fish.day & cal2$date<=last.fish.day & cal2$season==seas.cur, ]  #Limit calendar to only dates within first and last fish day range.
    cal3 <- rbind(cal3, cal.cur)
  }
  rm(first.fish.day, last.fish.day, seas.cur, cal.cur)
  
# Melt into two calendars, one for each species
  cal4 <- melt(cal3, 
               id.vars=c("myear", "site", "season", "date", "dow", "j.day", "j.wk", "trapmode", "fish.data"),
               measure.vars=c("chs.new", "sts.new"),
               variable.name="species", value.name="catch")
  cal4$species <- substr(cal4$species, start=1, stop=3)
  
# Merge in stratum assignments for each week
  cal5 <- merge(cal4, select(fish.strata, -catch, -te, -abundance),
                by=c("myear", "site", "species", "season", "j.wk"),
                all=TRUE) %>%
    arrange(species, season, date) %>%
    select(myear, site, species, season, j.wk, strat.id, date, everything())

# Write output files
  saveRDS(cal5, paste(dir_drvd, paste(my.txt, site.cur, "NonSub_daily.RData", sep="_"), sep=""))
  write.csv(cal5, paste(dir_drvd, paste(my.txt, site.cur, "NonSub_daily.csv", sep="_"), sep=""), row.names=FALSE)
  
## END ##     