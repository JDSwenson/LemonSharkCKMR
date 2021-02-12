# LemonSharkCKMR
This is a project to test the sensitivity of close-kin mark-recapture models when applied to simulated and real populations of a long-lived marine species, with a focus on the Lemon Shark.
The file structure is organized as follows:
**Lemon_Shark:** This folder contains code to format data from a real Lemon Shark population and fit a CKMR model.
**Leslie_matrix_simulations:** This folder contains code to test CKMR performance when applied to an ideal population of Lemon Sharks simulated from a Leslie matrix. All the scripts run on a loop, storing output iteratively. These scripts examine CKMR performance at the population level, with no individual variation.
**fishSim_simulations:** This folder contains code to test CKMR performance when applied to a more stochastic population of Lemon Sharks using the individual-based simulation program fishSim (https://github.com/SMBaylis/fishSim). All the scripts run on a loop, storing output iteratively.
**functions:** This folder contains functions that define the prior probability of kinship and specify the likelihood function.
**key_resources:** This folder contains key papers outlining the theory and application of CKMR.
**master_working_scripts:** This folder contains the scripts that are actively under development. When I've confirmed that these scripts work, I copy the script to the associated 'simulations' folder and run it on a loop. I also update the files in the reference_scripts folder as necessary.
**reference_scripts:** This folder contains a markdown file with scripts that have produced informative data.
