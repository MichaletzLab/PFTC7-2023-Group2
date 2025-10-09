# PFTC7-2023
This repository contains:
  A data folder containing all licor and Excel files generated during PFTC7.
  A code folder containing functions based on the Fast Assimilation-Temperature Response (FAsTeR) method (Garen and Michaletz 2024).
  An Analysis Code folder containing all code to open and run.

  To run this code: 
    1. Clone repository to personal device.
    2. Open CleanData.CreateMainFile.R and run this file (set working directory to PFTC7-2023 folder).
    3. Open CleanNorway.R and run this file
    This will generate results2.csv which is the cleaned and complete data file.
    4. Open pftc7_vpd_analysis and run configure.datfile.R
    5. Run the Medlyn.g1.R file and the redblueFigs.R file for generation of Figure 3 and S1
    6. Open Analysis Code and run the DataPrep.R file
    7. Open the rest of the figure files in the Analysis Code folder and run them (starting with PCA)
    8. You should never edit the code folder nor the data folder.
