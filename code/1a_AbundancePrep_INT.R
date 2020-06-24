## PROCESS TRAP LOG DATA FOR CALCULATING ABUNDANCE: PART 1A ##
#
# Run site setup file first.
#
# Take raw trap log, format dates, apply season, set days of 
# usable data based on standard options.


## ---- 1. Basic data processing ----
# Read in raw log file
log.raw <- read.csv(paste(dir_data, paste(my.txt, site.cur, "traplog.csv", sep="_"), sep=""), stringsAsFactors=FALSE)

# Remove any blank rows or columns
log <- log.raw[!is.na(log.raw$date), ]
log <- log[log$date!="",]
log <- select(log, -starts_with("X."))

# Adjust dates
  log$date <- ymd(log$date)
  log$jday <- yday(log$date) #Precise jday values will be one off during leap years.
rm(log.raw)
  
# Add season
  # Fall season begins July 1, ends Jan 28
fall.st.day <- yday(paste(my.cur, "-07-01", sep=""))
log$season[log$jday<=28 | log$jday>=fall.st.day] <- "fall"
log$season[log$jday>28 & log$jday<fall.st.day] <- "spring"
rm(fall.st.day)

## ---- 2. Dates with usable data ----
# Sort out which days have usable fish data for CMR
# According to std rules and options.
# (Customization applied later.)
  log.te <- mutate(log, 
                   catch.data.good = FALSE,
                   rec.data.good = FALSE)  #Default value in all cases is FALSE
  # 1. catch data [C]
    if (include.INT.catch==TRUE){
      log.te$catch.data.good[log.te$trapmode %in% c("CON", "INT", "CAL")] <- TRUE  #SUB data excluded, this will be processed separately.
    } else {
      log.te$catch.data.good[log.te$trapmode %in% c("CON", "CAL")] <- TRUE  
    }
  
  # 2. recap data [R]
    log.te$rec.data.good[log.te$trapmode %in% c("CON", "CAL")] <- TRUE  #Recap data are never valid from SUB and INT days
    
  # 3. release (mark) data [M]
    rec.days <- log.te$trapmode %in% c("CON", "CAL") #Eligible recap days
    mark.data.good <- c(rec.days[2:length(rec.days)], NA)  #Take the vector of eligible recap days, move it up one, and append an NA to the end.
    log.te <- cbind(log.te, mark.data.good)
    
    rm(rec.days, mark.data.good)
    
  # Next, apply any custom mods for INT days