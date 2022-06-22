# LemonSharkCKMR
This is a project to test the sensitivity of close-kin mark-recapture models when applied to simulated and real populations of a long-lived marine species, with a focus on the Lemon Shark and Cownose Ray.

## Updates from 2022.06.20
- The simulation_log file was wiped and updated. It was split into two files: one for holding the population simulation and sampling settings; one for the Leslie matrix and prior settings.

**Here's my process for visualizing results:**
1. Load output files from simulation into markdown document
	a) Each trial should have a sub-objective and each sub-objective should have its own data analysis script
2. Generate and save HPD summary files and any other intermediate file needed for viz
3. Load intermediate files into markdown file for knitting and knit