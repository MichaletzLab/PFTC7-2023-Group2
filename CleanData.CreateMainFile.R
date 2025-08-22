#dependencies: Run CleanNorway.R first

# Code to read LI-6800 files
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
library(rTPC)

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

setwd("../data/FasterLicor") ##Right now these are actually on desktop

# load list of file paths for Licor data
at.files = list.files(full.names=T)

# extract list of IDs from the file names
at.ids <- lapply(at.files, function(x)strsplit(x, split = "[.]")[[1]][3]) #save the individual IDs, extracted from file names

# LOAD LICOR FILES
file.x <- list()
for (i in seq_along(at.files)) {
  x <- at.files[[i]]  # Take a single file name
  
  # Print the file name before processing
  print(paste("Processing file:", x))
  
  # Check if the file size is greater than 2000 bytes
  if (file.info(x)$size > 2000) {
    y <- tryCatch({
      read_6800_with_BLC(x)  # Read file
    }, error = function(e) {
      message(paste("Error reading file:", x))
      message(e)
      return(NULL)
    })
    
    # Only proceed if the file was read successfully
    if (!is.null(y)) {
      at.id <- at.ids[[i]]  # Take corresponding individual ID
      y$curveID <- rep(at.id, nrow(y))  # Add ID as a column to file
      y$obs <- as.integer(y$obs)
      y$time <- as.double(y$time)
      
      file.x[[i]] <- y  # Add file (now a dataframe) to the list of files/dataframes
      print(paste("Successfully processed file:", x))
    } else {
      print(paste("Failed to process file:", x))
    }
  } else {
    print(paste("File size is less than or equal to 2000 bytes, skipping:", x))
  }
}
file.x <- file.x[lengths(file.x) != 0] # remove empty placeholders
at.df <- do.call(bind_rows, file.x) # combine files into 1 dataframe
at.df <- subset(at.df, A >= -5) # remove unreasonable A values

#Change working directory back to project here!

# Add in the key for curveID and barcodes from photosynthesis group
at.meta = read.csv("data/Faster_Key.csv") %>%
  rename(curveID = 1, site = SiteID)
at.meta$curveID <- as.character(at.meta$curveID)
at.meta$site <- as.character(at.meta$site)

# Join Licor files with the key
at.all <- left_join(at.df, at.meta, by="curveID")

at.all$`Î”Pcham` <- at.all$Î.Pcham

# Now join the leaf trait pftc7 file with our file by matching the barcodes
at.meta.full= read.csv("data/PFTC7_SA_cleanish_traits_2023.csv")


####### NEW 5.8.24
#Now do the species name correction
FT <- at.meta.full
tochange <- read_excel("data/Heli_name_changes.xlsx")
FT_step_1 <- correct_leaf_ID(data = FT, changes = tochange, data_ID_column = "id", data_species_column = "species")
heli_naming_system <- read_excel("data/New Helichrysum naming system.xlsx")
#comm <- read_csv("community_data_names_cleaned.csv")

FT_new_name_system <- new_heli_naming_system(data = FT_step_1, naming_system = heli_naming_system, data_species_column = "species")
#comm_new_name_system <- new_heli_naming_system(data = comm, naming_system = heli_naming_system, data_species_column = "species")

write.csv(FT_new_name_system, "data/PFTC7_SA_clean_traits_08May2024.csv") # write a new file with the finalized data
#write.csv(comm_new_name_system, "PFTC7_SA_clean_community_08May2024.csv")

# Create a new dataset that has a new BarcodeLeaf column that matches id
## If you want to include further trait data - modify the select to include other rows (or remove this line to include all of the trait data)
at.meta.full.S = FT_new_name_system%>%
  mutate(BarcodeLeaf = id)%>%
  dplyr::select(id,BarcodeLeaf, leaf_area_cm2, species)

# Join based on BarcodeCutout
merged_df <- at.all %>%
  left_join(at.meta.full.S, by = "BarcodeLeaf")
merged_df <- mutate(merged_df, CutoutArea = ifelse(leaf_area_cm2 < 6, leaf_area_cm2, 6.0))

# Replace S (licor variable for area) with correct leaf area and fix species column
at.all.n = merged_df %>% 
  mutate(S = CutoutArea)%>%
  mutate(S = ifelse(is.na(S), 6, S))%>%
  select(-Species)%>%
  rename(Species = species)


# Recalculate gas exchange variables using new leaf area
at.corr = calc_licor6800(at.all.n)
# Apply match offset correction
at.corr = match_correct(at.corr) 
#write.csv(at.corr, "data/at.corr_deletethis.csv")
#at.corr=read.csv("data/at.corr_deletethis.csv")
# Apply nonequilibrium correction just to the rows with UseDynamic set to False
  # The Licor has different settings depending on the Flow/Flow_r rate. These can be found directly on the machine.
at.corr.noneq.1 = at.corr%>%
   filter(is.na(UseDynamic))%>%
   filter(Flow>= 598 & Flow<= 601)%>%
   filter(Flow_r>= 518 & Flow_r<= 600)%>%
   noneq_correct_full(.,dt1_c = 1.76, dt2_c = 1.24, aV_c = 64.4, dt1_h = 3.44, dt2_h = 5.83, aV_h = 94)
at.corr.noneq.2 = at.corr%>%
  filter(is.na(UseDynamic))%>%
  filter(Flow>= 398 & Flow<= 401)%>%
  filter(Flow_r>= 408 & Flow_r<= 448)%>%
  noneq_correct_full(.,dt1_c = 2.44, dt2_c = 1.56, aV_c = 65.4, dt1_h = 3.5, dt2_h = 4.86, aV_h = 86.45)
at.corr.noneq.3 = at.corr%>%
  filter(is.na(UseDynamic))%>%
  filter(Flow>= 548 & Flow<= 551)%>%
  filter(Flow_r>= 592 & Flow_r<= 627)%>%
  noneq_correct_full(.,dt1_c = 1.76, dt2_c = 1.24, aV_c = 64.4, dt1_h = 3.44, dt2_h = 5.83, aV_h = 94)
at.corr.noneq.4 = at.corr%>%
  filter(!is.na(UseDynamic))

at.corr.noneq.SA <- bind_rows(at.corr.noneq.2, at.corr.noneq.3, at.corr.noneq.4)
write.csv(at.corr.noneq.SA, "data/raw.at.corr.noneq.SA.csv")

at.corr.noneq.SA <- read.csv("data/raw.at.corr.noneq.SA.csv")
at.corr.noneq.norway <- read.csv("data/raw.at.corr.noneq.norway.csv")

cols_to_convert = c("P1_dur","P1_Fmax","P2_dur","P2_dQdt","P3_dur","P3_ΔF","ID.1","Fo.1","F1","T.F1","T.HIR","F2", "T.F2")
at.corr.noneq.norway <- at.corr.noneq.norway %>%
  mutate(across(all_of(cols_to_convert), as.character))

at.corr.noneq.SA.NW <- bind_rows(at.corr.noneq.SA, at.corr.noneq.norway)
write.csv(at.corr.noneq.SA.NW, "data/at.corr.noneq.SA.NW.csv")

#write.csv(at.corr.noneq, "SAforRIA.dat.csv")

#if at.subset1 or 2 fails, run install.packages("grDevices") and library(grDevices) and it fixes...
#Clean data files to only include "nice" curves
at.subset = at.corr.noneq.SA.NW
at.subset1 = cut.multimodal(at.subset)
at.subset2 = cut.topt.out.of.bounds(at.subset1)
at.subset3 = discard.hooks(at.subset2, 0.01) ##47 is stubborn...
#at.subset4 = cut.hooks(at.subset2,"truncate_plots.pdf")
write.csv(at.subset3, 'data/raw.discardHooks_data.csv', row.names = F)
#write.csv(at.subset4, 'raw.cutHooks_data.csv', row.names = F)

#write.csv(at.subset2, 'outputs/rawData.at.subset2.SANW.csv', row.names = F)

species.key = read.csv("data/Faster_Key.csv")%>%
  rename(curveID=Obs)%>%
  rename(site=SiteID)
species.key.n=read.csv("data/Norway.Key.csv")%>%
  mutate(site=site+5)%>%
  rename(Elevation = Elevation.masl)

results.discard.hooks = fit_weibull_breadth(at.subset3) #also if this fails try the grDevices thing
results.discard.hooks$curveID = as.numeric(results.discard.hooks$curveID)
results_meta_discard.hooks = left_join(results.discard.hooks, species.key, by="curveID")
results_meta_discard.hooks = left_join(results_meta_discard.hooks,species.key.n, by = "curveID")
results_meta_discard.hooks = results_meta_discard.hooks%>%
  #mutate(country = case_when(curveID > 1000~"Norway",curveID < 1000~"SAfrica"))%>%
  mutate(site = coalesce(site.x, site.y))%>%
  mutate(Elevation = coalesce(Elevation.x, Elevation.y))
write.csv(results_meta_discard.hooks, 'weibull.discard.hooks.SANW.csv', row.names = F)

Sfield.results.discard.hooks =fit_mod_schoolfield(at.subset3 %>% select(Tleaf, A, curveID), x = "Tleaf", y = "A", T_ref = 25)
Sfield.results.discard.hooks$curveID = as.numeric(Sfield.results.discard.hooks$curveID)
Sfield.results_meta_discard.hooks = left_join(Sfield.results.discard.hooks, species.key, by="curveID")
Sfield.results_meta_discard.hooks = left_join(Sfield.results_meta_discard.hooks,species.key.n, by = "curveID")
Sfield.results_meta_discard.hooks = Sfield.results_meta_discard.hooks%>%
  mutate(country = case_when(curveID > 1000~"Norway",curveID < 1000~"SAfrica"))
write.csv(Sfield.results_meta_discard.hooks, 'schoolfield.discard.hooks.SANW.csv', row.names = F)

#results.cut.hooks = fit_weibull_breadth(at.subset4)
#results.cut.hooks$curveID = as.numeric(results.cut.hooks$curveID)
#results_meta_cut.hooks = merge(results.cut.hooks, species.key, by="curveID")
#write.csv(results_meta_cut.hooks, 'weibull.cut.hooks.csv', row.names = F)

#Sfield.results.cut.hooks =fit_mod_schoolfield(at.subset4 %>% select(Tleaf, A, curveID), x = "Tleaf", y = "A", T_ref = 25)
#Sfield.results.cut.hooks$curveID = as.numeric(Sfield.results.cut.hooks$curveID)
#Sfield.results_meta_cut.hooks = merge(Sfield.results.cut.hooks, species.key, by="curveID")
#write.csv(Sfield.results_meta_cut.hooks, 'schoolfield.cut.hooks.csv', row.names = F)

##### Datafiles are all saved (raw and modeled sets) now. Analysis can proceed. 