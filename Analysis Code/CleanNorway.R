#Libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(nls.multstart)
library(nlstools)
library(lubridate)
library(segmented)
library(tidylog)
library(readxl)

# Set wd to main R project
source("code/correct_6800_functions.R")
source("code/match_correct.R")
source("code/noneq_correct_full.R")
source("code/curve_fitting_michaletz_2021.R")
source("code/curve_fitting_criteria.R")
source("code/cut.hooks.R")
source("code/Discard.Hooks.R")
source("code/Spname_fixes.R")
source("code/fit_mod_schoolfield.R")
source("code/fit_weibull_breadth.R")

setwd("../data/FasterLicorNorway")

# Load list of file paths for Licor data
at.files <- list.files(full.names = TRUE)
at.ids <- lapply(at.files, function(x) {
  # Check if the filename matches the specified pattern
  if (grepl("\\d{4}-\\d{2}-\\d{2}-\\d+_logdata\\d{4}", x)) {
    # Extract the digits after "_logdata"
    match <- regmatches(x, regexpr("_logdata\\d{4}", x))
    if (length(match) > 0) {
      return(substr(match, 9, 12))
    } else {
      return(NA)
    }
  }
})
# Load Licor files
file.x <- list()
for(i in seq_along(at.files)) {
  x <- at.files[[i]]  # Take a single file name
  
  if(file.info(x)$size > 2000) {  # Weed out files that are empty and throw an error
    y <- read_6800_with_BLC(x)  # Read file
    at.id <- at.ids[[i]]  # Take corresponding individual ID
    
    y$curveID <- rep(at.id, nrow(y))  # Add ID as a column to file
    y$obs <- as.integer(y$obs)
    y$time <- as.double(y$time)
    
    file.x[[i]] <- y  # Add file (now a dataframe) to the list of files/dataframes
  }
}
# Create a list of all columns from each dataframe
all_columns_list <- lapply(file.x, names)

# Merge all unique column names
all_columns <- unique(unlist(all_columns_list))
all_columns
# Fill in missing columns with NA in each dataframe
file.x_filled <- lapply(file.x, function(df) {
  # Get missing columns in the current dataframe
  missing_columns <- setdiff(all_columns, names(df))
  if (length(missing_columns) > 0) {
    # Add missing columns to the dataframe and fill with NA
    for (col in missing_columns) {
      df[[col]] <- NA
    }
  }
  return(df)
})

# Combine all dataframes in the list file.x_filled into one big dataframe
at.df <- do.call(rbind, file.x_filled)
at.df = at.df %>%
  select(curveID, everything())

at.df <- subset(at.df, A >= -5) # remove unreasonable A values

setwd("../../data")
# Add in the key for curveID and barcodes from photosynthesis group
at.meta = read.csv("Norway.Key.csv")
at.meta$curveID <- as.numeric(at.meta$curveID)
at.df$curveID <- as.numeric(at.df$curveID)
at.meta$site <- as.numeric(at.meta$site)

# Join Licor files with the key
at.all <- left_join(at.df, at.meta, by="curveID")

at.all$`Î”Pcham` <- at.all$Î.Pcham


# Replace S (licor variable for area) with correct leaf area and fix species column
at.all.n = at.all %>% 
  mutate(S = area_cutout_cm2)

# Recalculate gas exchange variables using new leaf area
at.corr = calc_licor6800(at.all.n)
# Apply match offset correction
at.corr = match_correct(at.corr) 

# Apply nonequilibrium correction just to the rows with UseDynamic set to False
# The Licor has different settings depnding on the Flow/Flow_r rate. These can be found directly on the machine.
at.corr.noneq.1 = at.corr%>%
  filter(Flow < 275)%>%
  noneq_correct_full(.,dt1_c = 5.53, dt2_c = 3.52, aV_c = 66.44, dt1_h = 3.44, dt2_h = 5.83, aV_h = 94)
at.corr.noneq.2 = at.corr%>%
  filter(Flow > 275)%>%
  noneq_correct_full(.,dt1_c = 4.32, dt2_c = 2.77, aV_c = 66.46, dt1_h = 3.5, dt2_h = 4.86, aV_h = 86.45)

at.corr.noneq.norway <- bind_rows(at.corr.noneq.1, at.corr.noneq.2)
at.corr.noneq.norway <- at.corr.noneq.norway %>%
  rename(Species = taxon)%>%
  select(-species)%>%
  rename(Elevation = Elevation.masl)%>%
  mutate(site=site+5)

write.csv(at.corr.noneq.norway, "raw.at.corr.noneq.norway.csv")
