# LemonSharkCKMR
This is a project to test the sensitivity of close-kin mark-recapture models when applied to simulated and real populations of a long-lived marine species, with a focus on the Lemon Shark and Cownose Ray.

**Here's my process for visualizing results:**
1. Load output files from simulation into the data analysis markdown document
	a) Use the template in 00_MAIN_scripts and save in the model.assessment folder corresponding to the objective.
2. Generate and save HPD summary files and any other intermediate file needed for viz
3. Load intermediate files into markdown file for knitting and knit (esp for manuscript)

## Updates from 2022.06.29
- The model_settings.log file was wiped and updated to reflect the new breakout of objectives. Now, I've split model validation and sample scheme testing into two objectives:
1) Validate the model using informed priors, a uniform prior on lambda, and three fixed values for lambda, all while sampling just 1% of the population and only the YOY portion. This can be Figure S1.
2) Test different sampling schemes. Here, I'll try testing YOY only, sampling all juvenile age classes, and sampling all age classes (plus include POPs in model), as well as different sample sizes. This can be Figure 1.

- I also saved two different data analysis files that I can run to calculate the HPDIs and save the output files. One was made to run on just one trial; and the other was made to run on two trials.
- Going forward, whichever one I use, I should save it in the folder for the corresponding trial, while maintaining the original for replication.

## Updates from 2022.06.20
- The simulation_log file was wiped and updated. It was split into two files: one for holding the population simulation and sampling settings; one for the Leslie matrix and prior settings.

