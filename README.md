# PFTC7-2023-Group2

Analysis code for leaf gas exchange temperature responses across elevation
gradients in Norway and South Africa (PFTC7).

## Repository contents

- `data/` — LI-COR and Excel files generated during PFTC7, and vegetation height
  data from Norway (THREE-D and SeedClim).
- `code/` — functions implementing the Fast Assimilation-Temperature Response
  (FAsTeR) method (Garen and Michaletz 2024).
- `Analysis Code/` — the analysis scripts to open and run.
- `FigureS1and2_VPDanalysis/` — scripts for the VPD analysis figures.
- `outputs/` — created automatically. All figures and generated tables are
  written here.

Do not edit the `code/` or `data/` folders.

## Before running

Open `PFTC7-2023.Rproj`, or otherwise set the working directory to the
repository root. All scripts use paths relative to that root and will fail
elsewhere.

Some raw data files are not tracked in this repository because of their size.
[TODO: list which files, and where to obtain them.]

## Run order

Run each file from the `Analysis Code/` folder unless noted.

1. `CleanNorway.R`
2. `CleanData.CreateMainFile.R`
3. `DataPrep.R`
4. `PCA.R`
5. `Models.R`
6. `Figure1.R`, `Figure2.R`, `Figure3.R`
7. `Model_selection.R` — classifies each temperature response as peaked or
   monotonic and writes the Sharpe-Schoolfield parameter estimates that
   Figure 4 reads. Slow. Skip if `outputs/figure4_ss_parameters.csv` already
   exists.
8. `Figure4.R` — Sharpe-Schoolfield parameters for the peaked responses
   (*A*, iWUE) against PC1.
9. `Figure5.R` — temperature sensitivities and reference values for the
   monotonic responses (*E*, *g*sw, WUE) against PC1.
10. `FigureS4.R` and `FigureS5.R` — the PC2 analogues of Figures 4 and 5.
11. In `FigureS1and2_VPDanalysis/`: `configure.datfile.R`, then `FigureS1.R`,
    then `FigureS2.R`.
12. Back in `Analysis Code/`: `FigureS3_fullEnviron.R`.

`Fig_helpers.R` holds plotting and reporting helpers shared by Figures 4, 5,
S4, and S5. It is sourced automatically and is not run directly.

## Figures 4, 5, S4, and S5

These four figures are produced by two scripts. `Figure4.R` and `Figure5.R`
each take a principal component axis through the variable `PC_AXIS`, which
defaults to `PC1`. `FigureS4.R` and `FigureS5.R` set `PC_AXIS <- "PC2"` and
source the corresponding main script, so the supporting figures are strict
analogues of the main-text figures by construction.

Output names follow the axis automatically, so `Figure4.png` and
`FigureS4.png` cannot be confused.

`PC_AXIS` persists in the R session. After running `FigureS4.R` or
`FigureS5.R`, restart R or set `PC_AXIS <- "PC1"` before running the
main-text figures again. Each script reports the axis it is using when it
starts.

Alongside each figure, these scripts write a reported-statistics table to
`outputs/` giving the panel *p*-values at full precision, the degrees of
freedom, the random-effects structure actually fitted, and, where a site
random effect could not be estimated, a site-level test for comparison.

## References

Garen, J. C. and S. T. Michaletz. 2024. Fast Assimilation-Temperature
Response (FAsTeR) method.
