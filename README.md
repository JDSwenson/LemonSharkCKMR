# LemonSharkCKMR
This is a project to test the sensitivity of close-kin mark-recapture models when applied to simulated and real populations of a long-lived marine species, with a focus on the Lemon Shark and Cownose Ray.

As of 05/20/2021, I've taken a step back and am working with scripts that estimate average abundance over the (approximate) sample period. The main scripts I'm testing are the scripts with shortcuts here: 02_IBS/currently_testing

Specifically:
**CKMR_DoviIBS_Lemon_sharks_AvgN_05.20.2021_Lemon.R**: This script uses Dovi's IBS simulation to simulation a population of Lemon Sharks (no skipped breeding), and samples the population over six years. *This is the primary script I'm working with right now*.
**fishSim_CKMR_sex-specific_and_aggregated_loop_AvgN_6yrs_05.13.2021_Lemon.R**: this script uses fishSim to simulate a population of Lemon Sharks and samples this population over six years. The script returns relatively unbiased abundance estimates, with a median relative bias for males, females, and all adults around 2-4%.
**fishSim_CKMR_sex-specific_and_aggregated_loop_AvgN_6yrs_05.13.2021_CNR.R**: this script uses fishSim to simulate a population of Cownose Rays and samples the population over six years. The script returns abundance estimates that are quite biased, with a median relative bias closer to 20%


*All of these scripts include the kinship probability and likelihood functions as part of the script i.e. they're not sourced from the 00_functions folder*





More broadly, the file structure is organized as follows:
**00_functions:** This folder contains functions that define the prior probability of kinship and specify the likelihood function.
**00_key_resources:**This folder contains key papers outlining the theory and application of CKMR.
**01_Leslie_matrix:** This folder contains code to test CKMR performance when applied to an ideal population of Lemon Sharks simulated from a Leslie matrix. All the scripts run on a loop, storing output iteratively. These scripts examine CKMR performance at the population level, with no individual variation. The model performs well under these ideal situations.
**02_IBS:** This folder contains code to test CKMR performance when applied to a more stochastic population of Lemon Sharks using either a) the individual-based simulation program fishSim (https://github.com/SMBaylis/fishSim), or b) Dovi's IBS script. All the scripts run on a loop, storing output iteratively. This is the step I am currently focused on, as more stochastic demographics lead to a greater bias in CKMR estimates.
**03_Lemon_Shark:** This folder contains code to format data from a real Lemon Shark population and fit a CKMR model.
**04_other:** This folder contains other folders and scripts that have been used for model testing, but are not being actively used right now.