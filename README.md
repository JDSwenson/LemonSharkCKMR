# LemonSharkCKMR
This is a project to test the sensitivity of close-kin mark-recapture models when applied to simulated and real populations of a long-lived marine species, with a focus on the Lemon Shark.

My process, beginning with population simulations and proceeding through analysis and vizualisation, starts with the script population_simulation_collated.R in the 01_Data.generating.model folder. This script will generate distinct populations with specified parameters, and output an RDS file containing sampled individuals with associated metadata (e.g. age, parents, year sampled, etc.).

Next, the script model.estimation_collated.R in the 02_Estimation.model folder reads in the RDS file containing samples, creates a pairwise comparison matrix for each group of samples (corresponding to the number of distinct populations that were generated), and fits a CKMR model that is appropriate for the age composition of sampled individuals (half-sibling if only juveniles were sampled; half-sibling + parent-offspring if adults were also sampled). The CKMR model, including all samples from the posterior distribution, is saved as an RDS file (it's a list with each element containing the model and samples for each iteration), and parameter estimates are saved as a csv file.

The folder 03_Lemon_shark_data contains code to filter a longterm genetic dataset from lemon sharks in Bimini, Bahamas, create a pairwise comparison matrix, and fit a CKMR model. The code includes various options for subsetting the dataset, and outputs a csv file with parameter estimates.

Model results (i.e. parameter estimates) are read into the script collate_and_export_results.Rmd in the 04_DataViz folder. This script focuses on collating results from multiple runs and filtering to only include results where the Markov chains converged. These results are combined to give one output file per objective and exported as a csv file.

Finally, the script results_analysis_and_figures_markdown.Rmd - also in the 04_DataViz folder - reads in the collated results for each objective and produces plots and tables.