## BASIC SETUP FOR ALL ANALYSES ##
# Source prior to running anything else.


# Libraries and options
  library(dplyr)
  library(reshape2)
  library(lubridate)
  # library(RODBC)
  library(ggplot2)
  options(stringsAsFactors = FALSE)
  
# Make table display NAs
  table = function (..., useNA = 'ifany') base::table(..., useNA = useNA)


# Standard directories within the project
  dir_data <- "data_source//"  
  dir_drvd <- "data_derived//"
  dir_out  <- "output//"
  dir_code <- "code//"

## END ##  