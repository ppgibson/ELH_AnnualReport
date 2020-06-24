## CUSTOM FUCTIONS TO ASSIST WITH PROCESSING SMOLT TRAP DATA ##
# (1) assign.groups(vector) 
# (2) assign.cal.days(sub.df)
# (3) set.strata(cts.df)
# (4) apply.boot(data, sub.smp=FALSE)
# (5) apply.gam(data.nonsub, sub.smp=FALSE, data.sub="none")


## ---- (1) assign.groups ----
# Function assign.groups: simple, generic function to classify a vector of values 
#     based on repeating values. As soon as the value in the input vector changes,
#     assign a new group. function Function to assign groups based on continuous values
# Input: vector of values to be grouped
# Output: vector of group ids (numeric), of same length as input vector
  assign.groups <- function(vector) {
    groups <- rep(NA, times=length(vector))
    groups[1] <- 1
    for (i in 2:length(vector)) {
      if (vector[i] == vector[i-1]) {  #If value is same a previous record...
        groups[i] <- groups[i-1]   #...then assign same group as previous record
      } else {
        groups[i] <- groups[i-1] + 1
      }
    }
    return(groups)
  }
  

  
## ---- (2) assign.cal.days ----  
# Function assign.cal.days: given a data frame of CAL and SUB counts, with dates,
#     decide which CAL days should be applied to which SUB day for estimating the 
#     subsmp prop. Calculation is based on j.weeks; all SUB days within a j.wk are
#     grouped, and any CAL days within the same TE wk automatically apply to those 
#     SUB days.  Then the function checks for any contiguous CAL days before the first
#     SUB day of the week/sub group, and any after, and includes up to three pre and 
#     three post days even if they are in a different TE week. 
#     Function results should be examined manually; for example, if there is a date 
#     gap between the SUB and following CAL, and the CAL occur in a new j.wk, then
#     those CAL days will not be assigned.
# Input: a dataframe with [date], [j.wk], and [mode] fields (colnames must match).
#        Note this function does not do any subsetting, so data must be grouped/subset
#        before running function.
# Output: same dataframe, now with a [group] column plus one cal.xx field for each
#         j.wk value with TRUE for each CAL day that applies to the prop calculation
#         for that j.wk
assign.cal.days <- function(sub.df) { 
  group <- assign.groups(vector=sub.df$mode)
  sub.df <- cbind(sub.df, group)

  # Loop through each TE week with any SUB days
  for (i in unique(sub.df$j.wk[sub.df$mode=="SUB"])) {
    sub.df <- mutate(sub.df, cal.cur=FALSE)
    sub.df$cal.cur[sub.df$j.wk==i & sub.df$mode=="CAL"] <- TRUE  #all CAL days in the current jwk apply to SUB days in the jwk
    n.cal.cur <- sum(sub.df$cal.cur==TRUE)
    if (n.cal.cur < 6) {
      # Additional post-sub CAL days
      last.sub.group <- max(sub.df$group[sub.df$mode=="SUB" & sub.df$j.wk==i])  #group of latest SUB day in current J wk
      last.sub.day <- max(sub.df$date[sub.df$group==last.sub.group]) #latest SUB day of current group
      post.cal.days <- (filter(sub.df, date>last.sub.day, mode=="CAL"))[1:3,] #3 cal days following latest SUB day
      post.cal.days$t.diff <- post.cal.days$date - last.sub.day
      for (j in 1:sum(!is.na(post.cal.days$site))) {
        if (post.cal.days$t.diff[j]==j) {post.cal.days$cal.cur[j] <- TRUE}
      }
      # Additional pre-sub CAL days
      first.sub.group <- min(sub.df$group[sub.df$mode=="SUB" & sub.df$j.wk==i]) #group of earliest SUB day in current J wk
      first.sub.day <- min(sub.df$date[sub.df$group==first.sub.group])  #earliest SUB day in current group
      pre.cal.days <- arrange(filter(sub.df, date<first.sub.day, mode=="CAL"), desc(date))[1:3,] #3 cal days following latest SUB day
      if (sum(!is.na(pre.cal.days$site)) > 0) {  #in case there are no pre-sub cal days...
        pre.cal.days$t.diff <- first.sub.day - pre.cal.days$date 
        for (k in 1:sum(!is.na(pre.cal.days$site))) {
        if (pre.cal.days$t.diff[k]==k) {pre.cal.days$cal.cur[k] <- TRUE}
        }
      }
      # Get the set of dates where cal.cur should be TRUE for the current j.wk
      dates.cur <- pre.cal.days$date[!is.na(pre.cal.days$date) & pre.cal.days$cal.cur==TRUE]
      dates.cur <- c(dates.cur, post.cal.days$date[!is.na(post.cal.days$date) & post.cal.days$cal.cur==TRUE])
      sub.df$cal.cur[sub.df$date %in% dates.cur] <- TRUE
      which(colnames(sub.df)=="cal.cur")
      colnames(sub.df)[which(colnames(sub.df)=="cal.cur")] <- paste("cal", i, sep=".")  #rename cal.cur column
    }
  }
  
  return(sub.df)
}



## ---- (3) set.strata ----
# Function for automated stratification of TE weeks based on when R.j (recaps) cum sum is at least 10.
#
# This {set.strata} function takes an input data frame of weekly (by TE wk) fish counts with a [rec] 
# field and assigns each week to a stratum (labelled "a", "b", etc.) strictly according to the  
# rule that weeks will be combined until the cumulative rec sum (total recaps) is at least 10. 
# If the final stratum has <10 recaps, it will be combined with the previous stratum.
# 
# The function returns the same df with additional populated fields:
#   [cumsum] gives the running cumulative sum of weekly rec counts within the stratum
#   [stratum] reports the assigned stratum
#
set.strata <- function(cts.df) {   #Input is a data frame with a by-week (or day) rec column
  # Create new, empty fields to be filled
  cts.df$cumsum  <- NA  #for the resetting cum sum
  cts.df$stratum <- NA  #for the assigned stratum code
  
  # Set initial values (easier to do this here than within the loop)
  j <- 1  #This j counter moves through the set of letters to use for naming strata
  cts.df$cumsum[1]  <- cts.df$rec[1]  #For the first row, cum sum is same as current rec value.
  cts.df$stratum[1] <- "a"              #First row will always be assigned to first "a" stratum.
  
  # Loop through each row of the data frame to set [stratum] according to the >= 10 rule.
  for (i in 2:nrow(cts.df)) {  #Start with row 2, row 1 is already set above
    if (cts.df$cumsum[i-1] < 10) {  #Is the previous cumulative sum value < 10? If so then...
      cts.df$cumsum[i] <- cts.df$cumsum[i-1] + cts.df$rec[i]  #...cumsum value for current record is previous cumusm value plus Rj for current record.
      # (j value stays the same)
    } else {  #Alternatively, if the previous record's cumsum value is >= 10 then...
      cts.df$cumsum[i] <- cts.df$rec[i]   #...the current rec count starts a new running cum sum
      j <- j+1   #...and move to the next letter in sequence for setting the stratum id
    }
    cts.df$stratum[i] <- letters[j]  #Fill in the stratum value. The if/then statement controls the value of j.
  }
  
  # If final stratum does not have >= 10 recaps, it should be combined with the previous one.
  if (cts.df$cumsum[nrow(cts.df)] < 10) {   #If final cumsum value is less than 10, then...
    str.cur <- cts.df$stratum[nrow(cts.df)] #...extract the final stratum code as 'str.cur'
    if (j>1) {
      cts.df$stratum[cts.df$stratum==str.cur] <- letters[j-1]  #...and update all instances of this stratum code to the previous one.
    } else {
      cts.df$stratum[cts.df$stratum==str.cur] <- letters[j]  #Special case where cumsum never gets to 10, just use stratum 'a' again. 
    }
  }
  
  # Return the completed dataframe
  return(cts.df)
}

## ---- (4) apply.boot ----
# Function to run original.boot routine for each species/season combination.
# 
# This {apply.boot} function takes an input data frame of non-subsampling or 
# of subsampling data (as prepared in scripts 1 and 3, respectively), and runs
# the original.boot function for estimating variance associated with weekly and
# by-season estimates of smolt abundance. 
#
# The function returns a single data frame combining all original.boot output, 
# indexed by species, season, and j.wk. The "combined population estimate" values
# are labeled with [j.wk] = "all". 
#
# Note that <func_BootstrapVariance.R> must be already sourced in order to run this function.
#
apply.boot <- function(data, sub.smp=FALSE) {
  # Create a table of unique combinations of species and season
  groupings <- unique(data[, c("species", "season")])
  
  # Empty DF to be filled
  boot.results <- data.frame()
  
  for (i in 1:(nrow(groupings))) {
    sp.cur <- groupings[i, "species"]
    seas.cur <- groupings[i, "season"]
    print(paste(sp.cur, seas.cur))
    
    # Subset data for boot run
    dat.cur <- filter(data, species==sp.cur, season==seas.cur)
    
    if (sub.smp==FALSE) {
    # Run boot routine: version for non-sub data
    out.cur <- original.boot(data = dat.cur, time = "j.wk", 
                             catch="catch", mark="mark", recap="rec", sub.sample=FALSE)
    }
    if (sub.smp==TRUE) {
    # Run boot routine: version for sub data
    out.cur <- original.boot(data = dat.cur, time = "j.wk", 
                             catch="catch", mark="mark", recap="rec", 
                             sub.sample=TRUE, P24 = "p24", PSS="pss")
    }
   
    # Annotate results
    bywk.cur <- out.cur[[1]]
    bywk.cur <- mutate(bywk.cur, species=sp.cur, season=seas.cur)
    byseas.cur <- out.cur[[2]]
    byseas.cur <- mutate(byseas.cur, 
                         species=sp.cur, season=seas.cur,
                         J.WK = "all", Est_cap_eff = NA)
    if (sub.smp==TRUE) {
      byseas.cur$Est_prop <- NA
    }
    
    vals.cur <- rbind(bywk.cur, byseas.cur)
    
    # Append to shared DF
    boot.results <- rbind(boot.results, vals.cur)
  }
  
  # Clean the output
  colnames(boot.results) <- tolower(colnames(boot.results))
  boot.results <- select(boot.results, species, season, everything())
  
  # Add variance 
  boot.results$var <- boot.results$sd^2
  
  # Return compiled data frame
  return(boot.results)
}

## ---- (5) apply.gam ----
# Function to run the GAM interpolation routine for each species/season combination.
# 
# This {apply.gam} function takes input data frames of daily CMR data and (if
# applicable) daily subsample proportion data (as prepared in scripts 2 and 3, 
# respectively), and runs the GAM interpolation function for estimating total 
# smolt abundance for each species/season, including interpolation for missed
# trapping days. 
#
# The function returns a single data frame combining all gam function output, 
# indexed by grouping ({daily, weekly, or season}), species, season, and j.wk or j.day. 
# The `diagnostics` argument determines whether diagnostic plots will be printed (TRUE)  
# or not (FALSE). 
#
# Note that <func_GAMInterpolation.R> must be already sourced in order to run this function.
#
apply.gam <- function(data.nonsub, sub.smp = FALSE, data.sub="none", diagnostics=TRUE) {
  # Create a table of unique combinations of species and season
  groupings <- unique(data.nonsub[, c("species", "season")])
  
  # Empty DF to be filled
  gam.results <- data.frame()
  
  for (i in 1:(nrow(groupings))) {
    sp.cur <- groupings[i, "species"]
    seas.cur <- groupings[i, "season"]
    print(paste(sp.cur, seas.cur))
    
    # Subset data for GAM run
    dat.cur <- filter(data.nonsub, species==sp.cur, season==seas.cur)
    if (sub.smp==FALSE) {
      sub.cur <- data.frame()
    } 
    if (sub.smp==TRUE) {
      sub.cur <- filter(data.sub, species==sp.cur, season==seas.cur)
    }
    
    if (nrow(sub.cur)>0) {
      # Run gam routine: version for sub data
      out.cur <- GAM.pSpline.trap.est(data = dat.cur, 
                                      time ="j.day",
                                      groups="j.wk", 
                                      catch="catch", 
                                      mark="mark", 
                                      recap="rec",
                                      sub.sample=TRUE,
                                      sub.sample.data=sub.cur, 
                                      stime="j.day", 
                                      pss="pss", 
                                      p24="p24",
                                      diagnostics=diagnostics)
    }
    if (nrow(sub.cur)==0) {
      # Run gam routine: version for no sub data
      out.cur <- GAM.pSpline.trap.est(data = dat.cur, 
                                      time ="j.day",
                                      groups="j.wk", 
                                      catch="catch", 
                                      mark="mark", 
                                      recap="rec",
                                      sub.sample=FALSE, 
                                      diagnostics=diagnostics)
    }
    
    # Annotate results
    byday.cur <- out.cur[[1]] %>%
      mutate(species=sp.cur, season=seas.cur,
             grouping="daily", J.WK=NA,
             Lower_95_lognorm=NA, Upper_95_lognorm=NA)
    bywk.cur <- out.cur[[2]] %>%
      mutate(species=sp.cur, season=seas.cur,
             grouping="weekly", J.DAY=NA)
    byseas.cur <- out.cur[[3]] %>%
      mutate(species=sp.cur, season=seas.cur,
             grouping="season", J.DAY="all", J.WK="all")

    vals.cur <- rbind(byday.cur, bywk.cur, byseas.cur)
    
    # Append to shared DF
    gam.results <- rbind(gam.results, vals.cur)
    
    # Add some space to help with reviewing diagnostic plots
    print("==========================================================")
  }
  
  # Clean the output
  colnames(gam.results) <- tolower(colnames(gam.results))
  gam.results <- select(gam.results, species, season, grouping, j.wk, j.day, everything()) %>%
    rename(ci.low = lower_95, ci.high=upper_95, 
           ci.low.log=lower_95_lognorm, ci.high.log=upper_95_lognorm)
  gam.results$estimate <- round(gam.results$estimate, 0)
  
  # Return compiled data frame
  return(gam.results)
}


# Function to calculate median date
calc.med <- function(df, date.col, count.col) {
  df <- df[!is.na(df[, count.col]), ]
  seq.cur <- rep(x=df[, date.col], times=df[, count.col])
  return(median(seq.cur))
}