library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)
library(ggmcmc)


rm(list=ls())
today <- format(Sys.Date(), "%d%b%Y")

purpose.combo <- "skipped_breeding_HS.only_vs_HS.PO" #What am I comparing right now? This will be used to name output files

#----------------Read in files ------------------------------
#------------- MCMC parameters ----------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains
MCMC.settings <- paste0("thin", nt, "_draw", ni, "_burn", nb)

#------------- Simulation parameters and labels  ----------------#
#Initialize values
date.of.simulation1 = purpose1 = purpose1.lab = model.type1 = date.of.simulation2 = purpose2 = purpose2.lab = model.type2 = date.of.simulation3 = purpose3 = purpose3.lab = model.type3 = date.of.simulation4 = purpose4 = purpose4.lab = model.type4 <- NULL

date.of.simulation1 <- "25Apr2022"
purpose1 <- "skipped.breeding_trial.2"
purpose1.lab <- "HS|PO model"
model.type1 <- "HS.PO"

date.of.simulation2 <- "25Apr2022"
purpose2 <- "skipped.breeding_trial.1"
purpose2.lab <- "HS only model"
model.type2 <- "HS.only"

date.of.simulation3 <- "11Apr2022"
purpose3 <- "HS.only_refined.samples_all.comps"
purpose3.lab <- "HS only all comparisons"
model.type3 <- "HS.PO"

date.of.simulation4 <- "11Apr2022"
purpose4 <- "HS.only_refined.samples_downsample"
purpose4.lab <- "HS only downsample"
model.type4 <- "HS.PO"


seeds <- "Seeds2022.04.15"
sim.samples.1 = samples1.lab <- "240.samples"
sim.samples.2 = samples2.lab <- "480.samples"
sim.samples.3 = samples3.lab <- "720.samples"
sim.samples.4 = samples4.lab <- "960.samples"
sim.samples.all <- c(240, 480, 720, 960)

MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.results/"
results_plots_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.results/figures/"
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
parents_prefix <- "parents_breakdown/CKMR_parents.breakdown"
sample.prefix <- "sample_info/CKMR_sample.info"
survival.prefix <- "survival/CKMR_survival"
pop.size.prefix <- "pop_size/CKMR_pop.size"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"

#------------------------------- Results -----------------------------------#
#Set results to NULL and the script will run regardless of how many actual results are being compared
results.1 = results.2 = results.3 = results.4 <- NULL

results.1 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1, ".csv")) %>% 
  mutate(model_type = model.type1, 
         purpose.lab = purpose1.lab)

results.2 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2, ".csv")) %>% 
  mutate(model_type = model.type2,
         purpose.lab = purpose2.lab)

#date.of.simulation <- "05Apr2022"
results.3 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation3, "_", seeds, "_", purpose3, ".csv")) %>% 
  mutate(model_type = model.type3,
         purpose.lab = purpose3.lab)

results.4 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation4, "_", seeds, "_", purpose4, ".csv")) %>% 
    mutate(model_type = model.type4,
           purpose.lab = purpose4.lab)


#------------------------------- MCMC output -----------------------------------#
#date.of.simulation <- "07Apr2022"
s1.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose1))

s1.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose1))

s1.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose1))

s1.4 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose1))

#Trial 2
s2.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose2))

s2.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose2))

s2.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose2))

s2.4 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose2))

#Trial 3
#date.of.simulation <- "05Apr2022"
s3.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose3))

s3.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose3))

s3.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose3))

s3.4 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose3))

#Trial 4
s4.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose4))

s4.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose4))

s4.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose4))

s4.4 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose4))

#------------------------------- Population size details -----------------------------------#
#date.of.simulation <- "07Apr2022"
pop.size.1 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1))

pop.size.2 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2))

pop.size.3 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation3, "_", seeds, "_", purpose3))

pop.size.4 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation4, "_", seeds, "_", purpose4))

#------------------------------ Pairwise comparison matrices -------------------------------#
#Mom
mom.comps.1 <- readRDS(paste0(results_location, mom.comps.prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1))

mom.comps.2 <- readRDS(paste0(results_location, mom.comps.prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2))

mom.comps.3 <- readRDS(paste0(results_location, mom.comps.prefix, "_", date.of.simulation3, "_", seeds, "_", purpose3))

mom.comps.4 <- readRDS(paste0(results_location, mom.comps.prefix, "_", date.of.simulation4, "_", seeds, "_", purpose4))

#Dad
dad.comps.1 <- readRDS(paste0(results_location, dad.comps.prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1))

dad.comps.2 <- readRDS(paste0(results_location, dad.comps.prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2))

dad.comps.3 <- readRDS(paste0(results_location, dad.comps.prefix, "_", date.of.simulation3, "_", seeds, "_", purpose3))

dad.comps.4 <- readRDS(paste0(results_location, dad.comps.prefix, "_", date.of.simulation4, "_", seeds, "_", purpose4))

####------------------------------- Quick analysis of results -----------------------------------####

results.all <- results.1 %>% 
  bind_rows(results.2, results.3, results.4) %>% 
  mutate(cv = (sd/mean)*100) %>% 
  mutate(samples.lab = paste0(total_samples, ".samples"),
         parameter = ifelse(parameter == "Nf", "Nfb", parameter))

results.all$samples.lab <- factor(results.all$samples.lab, levels = c(samples1.lab, samples2.lab, samples3.lab, samples4.lab))
results.all$purpose.lab <- factor(results.all$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

jags_params <- c("Nfa", "Nfb", "Nf", "Nm", "surv", "lam", "psi", "pb") #Specify parameters

head(results.all)
tail(results.all)

#View instances where the model failed to converge
results.all %>% dplyr::filter(Rhat > 1.01) %>% View()

#Remove instances that failed to converge
no.convergence <- results.all %>% dplyr::filter(Rhat > 1.01)


 #-----------Median Relative bias by sample size-------------------------#
results.all %>% group_by(total_samples, parameter, purpose.lab) %>% 
#  dplyr::filter(model_type == "HS.PO") %>%
  dplyr::summarize(median = median(relative_bias), n = n()) %>% 
  arrange(desc(median)) %>% 
  View()

# results2 %>% dplyr::filter(parameter == "lam") %>% 
#   ggplot(aes(relative_bias, colour = which.dup, fill = which.dup)) +
#   geom_density(alpha = 0.5)
# 
# 
# results2 %>%  mutate(cv = (sd/mean)*100) %>% 
#   group_by(parameter, which.dup) %>% 
#   dplyr::summarize(mean.cv = mean(cv), mean.relative.bias = mean(relative_bias)) %>%
# #  dplyr::filter(parameter == "Nf" | parameter == "Nm") %>% 
#   arrange(mean.cv) %>% 
#   write_csv(file = paste0(results_location, "Mean_bias_and_precision_", date.of.simulation, "_", purpose, ".csv"))


#------------------------------- Within HPDI? -----------------------------------#
results.all %>% group_by(purpose.lab, total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100) %>% 
#  filter(parameter != "lam") %>% 
  arrange(desc(percent_in_interval)) %>% 
  View()

#------------------------------- CV -----------------------------------#
results.all %>% group_by(model_type, purpose.lab, total_samples, parameter) %>% 
  dplyr::summarize(mean.cv = mean(cv), 
                   median.bias = median(relative_bias),
                   percent_in_interval = sum(in_interval == "Y")/n() * 100) %>%
#  dplyr::filter(parameter == "Nf" | parameter == "Nm") %>% 
  arrange(mean.cv) %>% 
  write_csv(file = paste0(results_location, "Mean_bias_and_precision_", today, "_", purpose.combo, ".csv"))


####------------------------------- Relative bias density plots -----------------------------------####
#Try looking at distribution of relative bias to get an idea of how biased they are
param.plot.1 <- "Nf"
param.plot.2 <- "Nm"
param.plot.3 <- "survival" #Which parameters are we plotting below?
param.plot.4 <- "lambda"

DensityPlots.file <- paste0(results_plots_location, "Density_plots_", today, "_", purpose.combo, ".pdf") #Name file
pdf(file = DensityPlots.file, width = 14, height = 8) #Open pdf
sim.samples.all <- c(200, 600, 1000) #vector of total sample sizes

#Loop over different sample sizes
#Different plot for each sample size, color distributions by purpose
  for(s in 1:length(sim.samples.all)){
    
    samps <- sim.samples.all[s]

p1 <- results.all %>% dplyr::filter(parameter == "Nf", total_samples == samps) %>% 
  dplyr::select(Q50, mean, truth, relative_bias, purpose.lab) %>% 
  ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("Relative bias: ", param.plot.1, "; ", samps, " samples")) + 
  xlim(-50, 50)
  
p2 <- results.all %>% dplyr::filter(parameter == "Nm", total_samples == samps) %>% 
  dplyr::select(Q50, mean, truth, relative_bias, purpose.lab) %>% 
  ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("Relative bias: ", param.plot.2, "; ", samps, " samples")) + 
  xlim(-50, 50)

p3 <- results.all %>% dplyr::filter(parameter == "surv", total_samples == samps) %>% 
  dplyr::select(Q50, mean, truth, relative_bias, purpose.lab) %>% 
  ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("Relative bias: ", param.plot.3, "; ", samps, " samples")) + 
  xlim(-20, 20)
  
p4 <- results.all %>% dplyr::filter(parameter == "lam", total_samples == samps) %>%
  dplyr::select(Q50, mean, truth, relative_bias, purpose.lab) %>%
  ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
  geom_density(alpha = 0.6) +
  ggtitle(label = paste0("Relative bias: ", param.plot.4, "; ", samps, " samples")) +
  xlim(-20, 20)

print(ggarrange(p1, p2, p3, p4, common.legend = TRUE))

  }

dev.off() #Close pdf

####------------------------------- Relative bias violin plot -----------------------------------####
  # RB_violin.p1 <- results.all %>% dplyr::filter(total_samples == samps & purpose == purpose1) %>% 
  #   ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
  # geom_violin(aes(y=relative_bias), draw_quantiles = 0.5) +
  # #ylim(-50, 160) +
  # geom_hline(yintercept=0, col="black", size=1.25) +
  # #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  # labs(x="Iteration", y="Relative bias", title="Relative bias by iteration") +
  # scale_fill_brewer(palette="Set2") +
  # font("title", size = 10, face = "bold")

####------------------------------- Relative bias box plots -----------------------------------####
#Specify save location for pdf of plots
boxPlots.file <- paste0(results_plots_location, "BoxPlots_", today, "_", purpose.combo, ".pdf") #Name output file
pdf(file = boxPlots.file, width = 12, height = 12) #open pdf
sim.samples.all <- c(200, 600, 1000)

#Loop over different sample sizes and purposes
for(i in 1:length(sim.samples.all)){
  samps <- sim.samples.all[i]
  
  #boxplot plot
  RB_boxplot.s1 <- results.all %>% filter(total_samples == samps, purpose == purpose1.lab) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose1.lab, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 100)
  
  RB_boxplot.s2 <- results.all %>% filter(total_samples == samps, purpose == purpose2.lab) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose2.lab, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 100)

  
  RB_boxplot.s3 <- results.all %>% filter(total_samples == samps, purpose == purpose3.lab) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose3.lab, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 100)

  RB_boxplot.s4 <- results.all %>% filter(total_samples == samps, purpose == purpose4.lab) %>%
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose4.lab, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") +
    theme(legend.title = element_blank()) +
    ylim(-50, 100)
  
print(ggarrange(RB_boxplot.s1, RB_boxplot.s2, RB_boxplot.s3, RB_boxplot.s4, nrow = 2, ncol = 2, common.legend = TRUE))
  }

dev.off() #close pdf

####------------------------------- Calculate different HPDI intervals-------------------------------#### 
#How many estimates fall in different HPD intervals?
#Calculate HPD interval for each iteration
intervals <- c(seq(from = 0.95, to = 0.05, by = -0.05)) #Specify HPDI intervals to calculate
iterations = 200

#Assumes four "purposes" are being compared
#Initialize dataframes for loop below
HPD.240.1 <- NULL
HPD.480.1 <- NULL
HPD.720.1 <- NULL
HPD.960.1 <- NULL
HPD.240.2 <- NULL
HPD.480.2 <- NULL
HPD.720.2 <- NULL
HPD.960.2 <- NULL
HPD.240.3 <- NULL
HPD.480.3 <- NULL
HPD.720.3 <- NULL
HPD.960.3 <- NULL
HPD.240.4 <- NULL
HPD.480.4 <- NULL
HPD.720.4 <- NULL
HPD.960.4 <- NULL

for(i in 1:length(s1.1)){ #Loop over all the posterior samples from each iteration; length should be same as number of iterations
  for(j in 1:length(intervals)){ #Loop over the different intervals for calculating HPDI
    
    #---------------------------Purpose 1--------------------------------#
    #240 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.240.1 <- combine.mcmc(s1.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.240.1 <- post.HPD.240.1[row.names(post.HPD.240.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.240.1 <- rbind(HPD.240.1, post.HPD.240.1)
    
    #480 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.480.1 <- combine.mcmc(s1.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.480.1 <- post.HPD.480.1[row.names(post.HPD.480.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.480.1 <- rbind(HPD.480.1, post.HPD.480.1)
    
    #720 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.720.1 <- combine.mcmc(s1.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.720.1 <- post.HPD.720.1[row.names(post.HPD.720.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.720.1 <- rbind(HPD.720.1, post.HPD.720.1)
    
    #960 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.960.1 <- combine.mcmc(s1.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.960.1 <- post.HPD.960.1[row.names(post.HPD.960.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.960.1 <- rbind(HPD.960.1, post.HPD.960.1)
    
    
    #---------------------------Purpose 2--------------------------------#    
    #240 samples
    #Calcualte HPD interval for iteration i and interval j
    
    post.HPD.240.2 <- combine.mcmc(s2.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.240.2 <- post.HPD.240.2[row.names(post.HPD.240.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.240.2 <- rbind(HPD.240.2, post.HPD.240.2)
    
    #480 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.480.2 <- combine.mcmc(s2.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.480.2 <- post.HPD.480.2[row.names(post.HPD.480.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.480.2 <- rbind(HPD.480.2, post.HPD.480.2)
    
    #720 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.720.2 <- combine.mcmc(s2.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.720.2 <- post.HPD.720.2[row.names(post.HPD.720.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.720.2 <- rbind(HPD.720.2, post.HPD.720.2)
    
    #960 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.960.2 <- combine.mcmc(s2.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.960.2 <- post.HPD.960.2[row.names(post.HPD.960.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.960.2 <- rbind(HPD.960.2, post.HPD.960.2)
    

    # # #---------------------------Purpose 3--------------------------------#
    # #240 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.240.3 <- combine.mcmc(s3.1[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.240.3 <- post.HPD.240.3[row.names(post.HPD.240.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.240.3 <- rbind(HPD.240.3, post.HPD.240.3)
    # 
    # #480 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.480.3 <- combine.mcmc(s3.2[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.480.3 <- post.HPD.480.3[row.names(post.HPD.480.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.480.3 <- rbind(HPD.480.3, post.HPD.480.3)
    # 
    # #720 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.720.3 <- combine.mcmc(s3.2[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.720.3 <- post.HPD.720.3[row.names(post.HPD.720.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.720.3 <- rbind(HPD.720.3, post.HPD.720.3)
    # 
    # #960 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.960.3 <- combine.mcmc(s3.3[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.960.3 <- post.HPD.960.3[row.names(post.HPD.960.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.960.3 <- rbind(HPD.960.3, post.HPD.960.3)
    
    # #---------------------------Purpose 4--------------------------------#
    #240 samples
    #Calcualte HPD interval for iteration i and interval j
    # post.HPD.240.4 <- combine.mcmc(s4.1[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.240.4 <- post.HPD.240.4[row.names(post.HPD.240.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.240.4 <- rbind(HPD.240.4, post.HPD.240.4)
    # 
    # #480 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.480.4 <- combine.mcmc(s4.2[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.480.4 <- post.HPD.480.4[row.names(post.HPD.480.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.480.4 <- rbind(HPD.480.4, post.HPD.480.4)
    # 
    #     #720 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.720.4 <- combine.mcmc(s4.2[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.720.4 <- post.HPD.720.4[row.names(post.HPD.720.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.720.4 <- rbind(HPD.720.4, post.HPD.720.4)
    # 
    # #960 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.960.4 <- combine.mcmc(s4.3[[i]]) %>%
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>%
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.960.4 <- post.HPD.960.4[row.names(post.HPD.960.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.960.4 <- rbind(HPD.960.4, post.HPD.960.4)
    
  }
  print(paste0("Finished with iteration ", i))
}
    


####------------------------------- Create dataframes for HPDI plotting------------------------------#### 
#---------------------------Purpose 1 --------------------------------#    
#Separate out by sample size
HPD.240.1.4viz <- results.1 %>% filter(total_samples == 240) %>% 
  right_join(HPD.240.1, by = c("parameter", "iteration"))

HPD.480.1.4viz <- results.1 %>% filter(total_samples == 480) %>% 
  right_join(HPD.480.1, by = c("parameter", "iteration"))

HPD.720.1.4viz <- results.1 %>% filter(total_samples == 720) %>% 
  right_join(HPD.720.1, by = c("parameter", "iteration"))

HPD.960.1.4viz <- results.1 %>% filter(total_samples == 960) %>% 
  right_join(HPD.960.1, by = c("parameter", "iteration"))

HPD.240.1.summary <- HPD.240.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.240.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.480.1.summary <- HPD.480.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.480.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.720.1.summary <- HPD.720.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.720.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD960.1.summary <- HPD.960.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.960.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes    
HPD.1.summary_all <- HPD.240.1.summary %>% inner_join(HPD.480.1.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD.720.1.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD960.1.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

#Read in previously saved file
#HPD.1.summary_all <- read_csv(file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation1, "_", purpose1, ".csv"))

HPD.1.summary.tidy <- HPD.1.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = model.type1,
         purpose.lab = purpose1.lab)

#Save
#write_csv(HPD.1.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation1, "_", purpose1, ".csv"))


#---------------------------Purpose 2 --------------------------------#    
#Separate out by sample size
HPD.240.2.4viz <- results.2 %>% filter(total_samples == 240) %>% 
  right_join(HPD.240.2, by = c("parameter", "iteration")) %>% 
  mutate(parameter = ifelse(parameter == "Nf", "Nfb", parameter),
         truth = ifelse(parameter == "Nfb", breed.truth, truth))

HPD.480.2.4viz <- results.2 %>% filter(total_samples == 480) %>% 
  right_join(HPD.480.2, by = c("parameter", "iteration")) %>% 
  mutate(parameter = ifelse(parameter == "Nf", "Nfb", parameter),
         truth = ifelse(parameter == "Nfb", breed.truth, truth))

HPD.720.2.4viz <- results.2 %>% filter(total_samples == 720) %>% 
  right_join(HPD.720.2, by = c("parameter", "iteration")) %>% 
  mutate(parameter = ifelse(parameter == "Nf", "Nfb", parameter),
         truth = ifelse(parameter == "Nfb", breed.truth, truth))

HPD.960.2.4viz <- results.2 %>% filter(total_samples == 960) %>% 
  right_join(HPD.960.2, by = c("parameter", "iteration")) %>% 
  mutate(parameter = ifelse(parameter == "Nf", "Nfb", parameter),
         truth = ifelse(parameter == "Nfb", breed.truth, truth))


#Calculate percent estimates in each interval
HPD.240.2.summary <- HPD.240.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.240.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.480.2.summary <- HPD.480.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.480.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.720.2.summary <- HPD.720.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.720.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD960.2.summary <- HPD.960.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.960.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes
HPD.2.summary_all <- HPD.240.2.summary %>% inner_join(HPD.480.2.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD.720.2.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD960.2.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

#Read in previously saved file
#HPD.2.summary_all <- read_csv(file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation2, "_", purpose2, ".csv"))

HPD.2.summary.tidy <- HPD.2.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = model.type2,
         purpose.lab = purpose2.lab)

#Save
#write_csv(HPD.2.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation2, "_", purpose2, ".csv"))


#---------------------------Purpose 3 --------------------------------# 
#Separate out by sample size
HPD.240.3.4viz <- results.3 %>% filter(total_samples == 240) %>% 
  right_join(HPD.240.3, by = c("parameter", "iteration"))

HPD.480.3.4viz <- results.3 %>% filter(total_samples == 480) %>% 
  right_join(HPD.480.3, by = c("parameter", "iteration"))

HPD.720.3.4viz <- results.3 %>% filter(total_samples == 720) %>% 
  right_join(HPD.720.3, by = c("parameter", "iteration"))

HPD.960.3.4viz <- results.3 %>% filter(total_samples == 960) %>% 
  right_join(HPD.960.3, by = c("parameter", "iteration"))


#Calculate percent estimates in each interval
HPD.240.3.summary <- HPD.240.3.4viz %>% group_by(parameter, interval) %>% 
  drop_na() %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.240.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.480.3.summary <- HPD.480.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.480.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.720.3.summary <- HPD.720.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.720.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD960.3.summary <- HPD.960.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.960.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes
HPD.3.summary_all <- HPD.240.3.summary %>% inner_join(HPD.480.3.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD.720.3.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD960.3.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

#Read in previously saved file
#HPD.3.summary_all <- read_csv(file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation3, "_", purpose3, ".csv"))

HPD.3.summary.tidy <- HPD.3.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = model.type3,
         purpose.lab = purpose3.lab)

#Save
#write_csv(HPD.3.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation3, "_", purpose3, ".csv"))


#---------------------------Purpose 4 --------------------------------# 
#Separate out by sample size
HPD.240.4.4viz <- results.4 %>% filter(total_samples == 240) %>%
  right_join(HPD.240.4, by = c("parameter", "iteration"))

HPD.480.4.4viz <- results.4 %>% filter(total_samples == 480) %>%
  right_join(HPD.480.4, by = c("parameter", "iteration"))

HPD.720.4.4viz <- results.4 %>% filter(total_samples == 720) %>%
  right_join(HPD.720.4, by = c("parameter", "iteration"))

HPD.960.4.4viz <- results.4 %>% filter(total_samples == 960) %>%
  right_join(HPD.960.4, by = c("parameter", "iteration"))


#Calculate percent estimates in each interval
HPD.240.4.summary <- HPD.240.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.240.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.480.4.summary <- HPD.480.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.480.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.720.4.summary <- HPD.720.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.720.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD960.4.summary <- HPD.960.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.960.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes
HPD.4.summary_all <- HPD.240.4.summary %>% inner_join(HPD.480.4.summary, by = c("parameter", "interval")) %>%
  inner_join(HPD.720.4.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD960.4.summary, by = c("parameter", "interval")) %>%
  mutate(interval.scale = interval*100)

#Read in previously saved file
#HPD.4.summary_all <- read_csv(file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation4, "_", purpose4, ".csv"))

HPD.4.summary.tidy <- HPD.4.summary_all %>%
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>%
  mutate(model_type = model.type4,
         purpose.lab = purpose4.lab)



#---------------Create dataframes to plot by sample size--------#
#Set eveyrthing to Null so the below scripts will run
HPD.240.1.4viz = HPD.480.1.4viz = HPD.720.1.4viz = HPD.960.1.4viz = HPD.240.2.4viz = HPD.480.2.4viz = HPD.720.2.4viz = HPD.960.2.4viz = HPD.240.3.4viz = HPD.480.3.4viz = HPD.720.3.4viz = HPD.960.3.4viz = HPD.240.4.4viz = HPD.480.4.4viz = HPD.720.4.4viz = HPD.960.4.4viz <- NULL

HPD.240.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "240.samples")
HPD.480.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "480.samples")
HPD.720.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "720.samples")
HPD.960.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "960.samples")

HPD.240.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "240.samples")
HPD.480.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "480.samples")
HPD.720.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "720.samples")
HPD.960.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "960.samples")

HPD.240.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "240.samples")
HPD.480.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "480.samples")
HPD.720.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "720.samples")
HPD.960.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "960.samples")

HPD.240.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "240.samples")
HPD.480.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "480.samples")
HPD.720.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "720.samples")
HPD.960.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "960.samples")


#Combine above dataframes, and re-order purpose as factors
all.240.4viz <- HPD.240.1.4viz %>% bind_rows(HPD.240.2.4viz, HPD.240.3.4viz, HPD.240.4.4viz) 
all.240.4viz$purpose.lab <- factor(all.240.4viz$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

all.480.4viz <- HPD.480.1.4viz %>% bind_rows(HPD.480.2.4viz, HPD.480.3.4viz, HPD.480.4.4viz) 
all.480.4viz$purpose.lab <- factor(all.480.4viz$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

all.720.4viz <- HPD.720.1.4viz %>% bind_rows(HPD.720.2.4viz, HPD.720.3.4viz, HPD.720.4.4viz) 
all.720.4viz$purpose.lab <- factor(all.720.4viz$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

all.960.4viz <- HPD.960.1.4viz %>% bind_rows(HPD.960.2.4viz, HPD.960.3.4viz, HPD.960.4.4viz)
all.960.4viz$purpose.lab <- factor(all.960.4viz$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

#For when we need everything together
all.4viz <- HPD.1.summary.tidy %>% bind_rows(HPD.2.summary.tidy) %>% 
  mutate(parameter = ifelse(parameter == "Nf", "Nfb", parameter))
all.4viz$purpose.lab <- factor(all.4viz$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))



####------------------------------- Make HDPI line graphs ------------------------------#### 
#---------------------Create and store plots----------------------------------------------#
#-----------------------Make figures-----------------------------
#-----------------HPDI plot--------------------#
(p.Nfb <- all.4viz %>% dplyr::filter(parameter == "Nfb") %>% 
   ggplot(aes(x = interval.scale)) +
   geom_point(aes(y = percent_in_interval, color = purpose.lab, fill = purpose.lab, shape = purpose.lab), size = 2.5, alpha = .7) + 
   geom_smooth(aes(y = percent_in_interval, color = purpose.lab, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
   geom_line(aes(y = interval.scale), size = 1.5) + 
   ggtitle("Nfb") + 
   labs(x = "HPDI",
        y = "Percent in HPDI") +
   theme_bw() + 
   theme(legend.title = element_blank()) + 
   scale_colour_brewer(palette = "Set1") +
   facet_wrap(~samples))



(p.Nm <- all.4viz %>% dplyr::filter(parameter == "Nm") %>% 
  ggplot(aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose.lab, fill = purpose.lab, shape = purpose.lab), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose.lab, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("Nm") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~samples))


(p.surv <- all.4viz %>% dplyr::filter(parameter == "surv") %>% 
  ggplot(aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose.lab, fill = purpose.lab, shape = purpose.lab), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose.lab, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~samples))

#May want to add this: purpose.lab != purpose3.lab to the first line to filter out PO from lambda
(p.lam <- all.4viz %>% dplyr::filter(parameter == "lam") %>% 
  ggplot(aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose.lab, fill = purpose.lab, shape = purpose.lab), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose.lab, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~samples))


HPDI.plots <- ggarrange(p.Nfb, p.Nm, p.surv, p.lam, common.legend = TRUE, ncol = 2, nrow = 2, legend = "bottom")
annotate_figure(p = HPDI.plots, fig.lab = "C)", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

#Percent breeders
(p.pb <- all.4viz %>% dplyr::filter(parameter == "pb") %>% 
    ggplot(aes(x = interval.scale)) +
    geom_point(aes(y = percent_in_interval, color = purpose.lab, fill = purpose.lab, shape = purpose.lab), size = 2.5, alpha = .7) + 
    geom_smooth(aes(y = percent_in_interval, color = purpose.lab, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
    geom_line(aes(y = interval.scale), size = 1.5) + 
    ggtitle("Percent breeders") + 
    labs(x = "HPDI",
         y = "Percent in HPDI") +
    theme_bw() + 
    theme(legend.title = element_blank()) + 
    scale_colour_brewer(palette = "Set1") +
    facet_wrap(~samples))

#Total adult females
(p.Nfa <- all.4viz %>% dplyr::filter(parameter == "Nfa") %>% 
    ggplot(aes(x = interval.scale)) +
    geom_point(aes(y = percent_in_interval, color = purpose.lab, fill = purpose.lab, shape = purpose.lab), size = 2.5, alpha = .7) + 
    geom_smooth(aes(y = percent_in_interval, color = purpose.lab, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
    geom_line(aes(y = interval.scale), size = 1.5) + 
    ggtitle("Total adult females") + 
    labs(x = "HPDI",
         y = "Percent in HPDI") +
    theme_bw() + 
    theme(legend.title = element_blank()) + 
    scale_colour_brewer(palette = "Set1") +
    facet_wrap(~samples))

ggarrange(p.pb, p.Nfa, common.legend = TRUE, legend = "bottom")

#Density plots
(p2.Nfb <- results.all %>% dplyr::filter(parameter == "Nfb") %>% 
    ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
    geom_density(alpha = 0.6) + 
    ggtitle("Nfb") + 
    labs(x = "Relative bias of posterior median") +
    geom_vline(xintercept=c(0), linetype="dotted") +
    theme_bw() + 
    theme(legend.title = element_blank()) + 
    scale_fill_brewer(palette = "Set1") +
    xlim(-100, 100) + 
    facet_wrap(~samples.lab))


(p2.Nm <- results.all %>% dplyr::filter(parameter == "Nm") %>% 
  ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
  geom_density(alpha = 0.6) + 
  ggtitle("Nm") + 
  labs(x = "Relative bias of posterior median") +
  geom_vline(xintercept=c(0), linetype="dotted") +
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  scale_fill_brewer(palette = "Set1") +
  xlim(-100, 100) + 
  facet_wrap(~samples.lab))

(p2.surv <- results.all %>% dplyr::filter(parameter == "surv") %>% 
  ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
  geom_density(alpha = 0.6) + 
  ggtitle("Survival") + 
  labs(x = "Relative bias of posterior median") +
  geom_vline(xintercept=c(0), linetype="dotted") +
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  scale_fill_brewer(palette = "Set1") +
  xlim(-50, 50) + 
  facet_wrap(~samples.lab))

(p2.lam <- results.all %>% dplyr::filter(parameter == "lam") %>% 
    ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
    geom_density(alpha = 0.6) + 
    ggtitle("Lambda") + 
    labs(x = "Relative bias of posterior median") +
    geom_vline(xintercept=c(0), linetype="dotted") +
    theme_bw() + 
    theme(legend.title = element_blank()) + 
    scale_fill_brewer(palette = "Set1") +
    xlim(-10, 10) + 
    facet_wrap(~samples.lab))

RelBias.Density.Plots <- ggarrange(p2.Nfb, p2.Nm, p2.surv, p2.lam, common.legend = TRUE, ncol = 2, nrow = 2, legend = "bottom")
annotate_figure(p = RelBias.Density.Plots, fig.lab = "D)", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")


(p2.Nfa <- results.all %>% dplyr::filter(parameter == "Nfa") %>% 
    ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
    geom_density(alpha = 0.6) + 
    ggtitle("Nf - total") + 
    labs(x = "Relative bias of posterior median") +
    geom_vline(xintercept=c(0), linetype="dotted") +
    theme_bw() + 
    theme(legend.title = element_blank()) + 
    scale_fill_brewer(palette = "Set1") +
    xlim(-100, 100) + 
    facet_wrap(~samples.lab))


(p2.pb <- results.all %>% dplyr::filter(parameter == "pb") %>% 
    ggplot(aes(x = relative_bias, color = purpose.lab, fill = purpose.lab)) +
    geom_density(alpha = 0.6) + 
    ggtitle("Percent breeders") + 
    labs(x = "Relative bias of posterior median") +
    geom_vline(xintercept=c(0), linetype="dotted") +
    theme_bw() + 
    theme(legend.title = element_blank()) + 
    scale_fill_brewer(palette = "Set1") +
    xlim(-100, 100) + 
    facet_wrap(~samples.lab))


(p2.psi <- results.all %>% dplyr::filter(parameter == "psi") %>% 
    ggplot(aes(x = factor(parameter), fill = purpose.lab)) +
    geom_boxplot(aes(y = relative_bias)) + 
    ggtitle("Percent biennial breeders") + 
    labs(x = "Relative bias of posterior median") +
    geom_hline(yintercept=c(0), linetype="dotted") +
    theme_bw() + 
    ylim(-10, 10) +
    theme(legend.title = element_blank()) + 
    #scale_fill_brewer(palette = "Set1") +
    facet_wrap(~samples.lab))


ggarrange(p2.psi, p2.Nfa, p2.pb, common.legend = TRUE, legend = "bottom")





#Visualize CV
Bias_and_prec.df <- read_csv(file = paste0(results_location, "Mean_bias_and_precision_", today, "_", purpose.combo, ".csv"))

(CV.1 <- results.all %>% dplyr::filter(purpose.lab == purpose1.lab, samples.lab == "240.samples" | samples.lab == "480.samples", parameter == "Nfa" | parameter == "Nfb" | parameter == "Nm") %>% 
  ggplot(aes(x = HSPs_detected, y = cv, fill = parameter, color = parameter)) +
  geom_point() +
  theme_bw() +
  ylim(0, 100) +
  facet_wrap(~samples.lab) + 
  labs(title = purpose1.lab, x = "HSPs detected", y = "CV"))

(CV.2 <- results.all %>% dplyr::filter(purpose.lab == purpose2.lab, samples.lab == "240.samples" | samples.lab == "480.samples", parameter == "Nfa" | parameter == "Nfb" | parameter == "Nm") %>% 
    ggplot(aes(x = HSPs_detected, y = cv, fill = parameter, color = parameter)) +
    geom_point() +
    theme_bw() +
    ylim(0, 100) +
    facet_wrap(~samples.lab) + 
    labs(title = purpose2.lab, x = "HSPs detected", y = "CV"))

ggarrange(CV.1, CV.2, common.legend = TRUE, legend = "bottom")



####---------------------Density plots of kin detected----------------------------------------------####
#---------------------Prepare data for plotting----------------------------------------------#
#How many POPs detected per sample size?
results.all %>% 
  group_by(purpose.lab, samples.lab) %>% 
  summarize(mean(POPs_detected), mean(HSPs_detected))

#Format PO dataframe
results.PO <- results.all %>% dplyr::select(parameter, total_samples, kin.detected = POPs_detected,  unique_parents_in_sample, iteration, cv, purpose.lab) %>% 
  mutate(type = "PO")

#Format HS dataframe
results.HS <- results.all %>% dplyr::select(parameter, total_samples, kin.detected = HSPs_detected,  unique_parents_in_sample, iteration, cv, purpose.lab) %>% 
  mutate(type = "HS")

#Combine PO and HS dataframe for plotting
results.kin <- rbind(results.PO, results.HS) %>% 
  mutate(total_samples = paste0(total_samples, " samples"))
results.kin$total_samples <- factor(results.kin$total_samples, levels = c("200 samples", "600 samples", "1000 samples"))

#---------------------Create density plots----------------------------------------------#
#Maternal
k1 <- results.kin %>% dplyr::filter(parameter == "Nf", type == "PO") %>%
  ggplot(aes(x = kin.detected, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "MOPs detected") + 
  xlim(0, 150) + 
  facet_wrap(~total_samples, ncol = 1) + 
  theme(legend.title = element_blank())

k2 <- results.kin %>% dplyr::filter(parameter == "Nf", type == "HS") %>%
  ggplot(aes(x = kin.detected, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "MHSPs detected") + 
  xlim(0, 500) + 
  facet_wrap(~total_samples, ncol = 1)

m.kin <- ggarrange(k1, k2, ncol = 2, common.legend = TRUE)
annotate_figure(m.kin, fig.lab = "Mother kin pairs", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

#Paternal
k3 <- results.kin %>% dplyr::filter(parameter == "Nm", type == "PO") %>%
  ggplot(aes(x = kin.detected, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "FOPs detected") + 
  xlim(0, 150) + 
  facet_wrap(~total_samples, ncol = 1) + 
  theme(legend.title = element_blank())

k4 <- results.kin %>% dplyr::filter(parameter == "Nm", type == "HS") %>%
  ggplot(aes(x = kin.detected, color = purpose.lab, fill = purpose.lab)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "FHSPs detected") + 
  xlim(0, 500) + 
  facet_wrap(~total_samples, ncol = 1)

f.kin <- ggarrange(k3, k4, ncol = 2, common.legend = TRUE)
annotate_figure(f.kin, fig.lab = "Father kin pairs", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

#---------------------Create pdf of kin pair density plots----------------------------------------------#
#Specify file
KinPairs.file <- paste0(results_plots_location, "Distribution_of_kin_", today, "_", purpose.combo, ".pdf")
pdf(file = KinPairs.file, width = 11, height = 14) #Open pdf

annotate_figure(m.kin, fig.lab = "Mother kin pairs", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")
annotate_figure(f.kin, fig.lab = "Father kin pairs", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

dev.off() #close pdf



#-----------------------Examine expected vs observed kin pairs--------------------
#Make better labels for plotting
results.all$sample.label <- factor(paste0(results.all$total_samples, " samples"), levels = c("200 samples", "600 samples", "1000 samples"))

#Calcualte % difference in expected vs observed kin
results.all2 <- results.all %>% mutate(HS_exp_obs_diff = (Exp_HSPs - HSPs_detected)/Exp_HSPs,
                                      PO_exp_obs_diff = (Exp_POPs - POPs_detected)/Exp_POPs) %>%
  mutate(PO_exp_obs_diff = replace_na(PO_exp_obs_diff, 0)) %>% 
  mutate(all_exp_obs_diff = HS_exp_obs_diff + PO_exp_obs_diff)


# results.all %>% mutate(HS_exp_obs_diff = Exp_HSPs - HSPs_detected,
#                        PO_exp_obs_diff = Exp_POPs - POPs_detected) %>%
#   mutate(all_exp_obs_diff = HS_exp_obs_diff + PO_exp_obs_diff) %>% 
#   dplyr::select(purpose,
#                 parameter, 
#                 Q50, 
#                 truth, 
#                 Exp_POPs, 
#                 POPs_detected, 
#                 PO_exp_obs_diff, 
#                 Exp_HSPs, 
#                 HSPs_detected, 
#                 HS_exp_obs_diff, 
#                 total_samples, 
#                 relative_bias, 
#                 in_interval) %>% 
#   View()

#----------------Make scatter plot of expected vs observed 
#Purpose 1
EO1 <- results.all2 %>% dplyr::filter(purpose.lab == purpose1.lab) %>% 
  ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
  geom_point() +
  facet_wrap(~sample.label) + 
  labs(title = purpose1.lab, x = "Percent difference", y = "Relative bias") + 
  theme(legend.title = element_blank())

#Purpose 2
EO2 <- results.all2 %>% dplyr::filter(purpose.lab == purpose2.lab) %>% 
  ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
  geom_point() +
  facet_wrap(~sample.label) + 
  labs(title = purpose2.lab, x = "Percent difference", y = "Relative bias")

#Purpose 3
EO3 <- results.all2 %>% dplyr::filter(purpose.lab == purpose3.lab) %>% 
  ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
  geom_point() +
  facet_wrap(~sample.label) + 
  labs(title = purpose3.lab, x = "Percent difference", y = "Relative bias")

#Purpose 4
EO4 <- results.all2 %>% dplyr::filter(purpose.lab == purpose4.lab) %>% 
  ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
  geom_point() +
  facet_wrap(~sample.label) + 
  labs(title = purpose4.lab, x = "Percent difference", y = "Relative bias")


EO.all <- ggarrange(EO1, EO2, EO3, EO4, common.legend = TRUE)
annotate_figure(EO.all, fig.lab = "Relative bias by percent difference in expected vs observed kin pairs", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")


####------------------Posterior predictive distribution------------------
head(s1.1[[1]])

#Add mcmc chains together
s.comb <- combine.mcmc(s1.1[[1]]) %>% 
  as_tibble()

Nf.pre <- s.comb$Nf
Nm.pre <- s.comb$Nm
lam.pre <- s.comb$lam
surv.pre <- s.comb$surv
n.trials <- nrow(s.comb)
mom.n.comps <- mom.comps.1$

  #Put in a loop
rbinom(n = 1, size = mom.n.comps, (surv^mom.mort.yrs[i])/(Nf*(lam^mom.popGrowth.yrs[i])), mom.n.comps[i])





####-------------------Examine comparisons of funny results----------------####
HPD.3.summary.tidy %>% dplyr::filter(parameter == "lam") %>% View()







####---------------------Extra scripts----------------------------------------------####
#---------------------Survival----------------------------------------------#
surv.1 <- readRDS(paste0(results_location, survival.prefix, "_", date.of.simulation, "_", seeds, "_", purpose1)) %>% 
  mutate(iteration = iter) %>% 
  dplyr::filter(year >= 85) %>% 
  distinct(iter, year, .keep_all = TRUE)

#Look at HPDI and relative bias for each of the last five years of the simulation
results.surv <- results.all %>% dplyr::filter(parameter == "surv", total_samples == 1000) %>% 
  dplyr::select(parameter, Q2.5, Q50, Q97.5, HPD2.5, mean, HPD97.5, sd, truth, iteration,
                mean.relative.bias = relative_bias,
                mean.in.interval = in_interval,
                total_samples) %>% 
  inner_join(surv.1, by = "iteration") %>% 
  dplyr::rename(mean.truth = truth, survival.per.yr = survival) %>% 
  mutate(annual.relative.bias = round(((Q50 - survival.per.yr)/survival.per.yr)*100,1)) %>%
  mutate(annual.in.interval = ifelse(HPD2.5 < survival.per.yr & survival.per.yr < HPD97.5, "Y", "N"))


HPD.1000.1.survOnly <- HPD.1000.1 %>% dplyr::filter(parameter == "surv") %>% 
  as_tibble()

HPD.1000.1.4viz.survOnly <- results.surv %>%
  right_join(HPD.1000.1.survOnly, by = c("parameter", "iteration"))



HPD.1000.annual.surv.interval <- HPD.1000.1.4viz.survOnly %>% 
  mutate(annual.in.interval = ifelse(lower < survival.per.yr & survival.per.yr < upper, "Y", "N")) %>% 
  group_by(interval) %>% 
  dplyr::summarize(percent.in.interval = sum(annual.in.interval == "Y")/n() * 100) %>% 
  mutate(type = "annual",
         interval.scale = interval * 100) %>% 
  dplyr::select(interval.scale, percent.in.interval, type)

HPD1000.mean.surv.interval <- all.1000.4viz %>% ungroup() %>% 
  dplyr::filter(parameter == "surv", samples == "1000.samples", purpose == "all.comps") %>% 
  dplyr::select(interval.scale, percent.in.interval = percent_in_interval) %>% 
  mutate(type = "mean")

HPD1000.4viz.survOnly <- rbind(HPD.1000.annual.surv.interval, HPD1000.mean.surv.interval)

ggplot(data = HPD1000.4viz.survOnly,
       aes(x = interval.scale)) +
  geom_point(aes(y = percent.in.interval, color = type, fill = type, shape = type), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent.in.interval, color = type, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("Survival: mean vs annual") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())















#----------Bring in info for parents and samples----------------#
# Breakdown of offspring for each parent
rents.1 <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose1)) %>% 
  mutate(model_type = "HS.PO")

rents.2 <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose2)) %>% 
  mutate(model_type = "HS.only")

rents.3 <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose3)) %>% 
  mutate(model_type = "HS.PO")

rents.4 <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose4)) %>% 
  mutate(model_type = "HS.only")

#Breakdown of samples drawn from simulation
samples1 <- sample.info.1 <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose1)) %>% 
  mutate(model_type = "HS.PO",
         purpose = purpose1)

samples2 <- sample.info.2 <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose2)) %>% 
  mutate(model_type = "HS.only",
         purpose = purpose2)

samples3 <- sample.info.3 <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose3)) %>% 
  mutate(model_type = "HS.PO",
         purpose = purpose3)

samples4 <- sample.info.4 <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose4)) %>% 
  mutate(model_type = "HS.only",
         purpose = purpose4)

head(samples1)

source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO.R")
source("./01_MAIN_scripts/functions/remove_dups.R")

#Make dataframes of pairwise comparisons for each iteration and sample size for further investigation into gap years, pop growth years, and number of detected kin.

mom_comps.samples1 <- dad_comps.samples1 <- NULL

options(dplyr.summarise.inform = FALSE)

for(iter in 1:100){
  print(paste0("working on iteration", iter))

    #50 samples per year
  samples1.temp1 <- samples1 %>% filter(sample.size == 50 & iteration == iter)
  NoDups1.temp1 <- samples1.temp1 %>% split.dups()
  save.last1 <- NoDups1.temp1[[2]]
  filter1.temp1 <- save.last1 %>% filter.samples()
  PO.samples1.list1 <- filter1.temp1[[1]]
  HS.samples1.df1 <- filter1.temp1[[2]]
  pairwise.samples1 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list1, filtered.samples.HS.df = HS.samples1.df1)
  
  mom_comps.temp1 <- pairwise.samples1[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 50)
  dad_comps.temp1 <- pairwise.samples1[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 50)

  #150 samples per year
  samples1.temp2 <- samples1 %>% filter(sample.size == 150 & iteration == iter)
  NoDups1.temp2 <- samples1.temp2 %>% split.dups()
  save.last2 <- NoDups1.temp2[[2]]
  filter1.temp2 <- save.last2 %>% filter.samples()
  PO.samples1.list2 <- filter1.temp2[[1]]
  HS.samples1.df2 <- filter1.temp2[[2]]
  pairwise.samples2 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list2, filtered.samples.HS.df = HS.samples1.df2)
  
  mom_comps.temp2 <- pairwise.samples2[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 150)
  dad_comps.temp2 <- pairwise.samples2[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 150)
  
  #250 samples per year
  samples1.temp3 <- samples1 %>% filter(sample.size == 250 & iteration == iter)
  NoDups1.temp3 <- samples1.temp3 %>% split.dups()
  save.last3 <- NoDups1.temp3[[2]]
  filter1.temp3 <- save.last3 %>% filter.samples()
  PO.samples1.list3 <- filter1.temp3[[1]]
  HS.samples1.df3 <- filter1.temp3[[2]]
  pairwise.samples3 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list3, filtered.samples.HS.df = HS.samples1.df3)
  
  mom_comps.temp3 <- pairwise.samples3[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 250)
  dad_comps.temp3 <- pairwise.samples3[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 250)
  
  
    mom_comps.samples1 <- rbind(mom_comps.samples1, mom_comps.temp1, mom_comps.temp2, mom_comps.temp3) %>% 
      as_tibble()
    dad_comps.samples1 <- rbind(dad_comps.samples1, dad_comps.temp1, dad_comps.temp2, dad_comps.temp3) %>% 
      as_tibble()
  

}
  
head(mom_comps.samples1)
head(dad_comps.samples1)

nrow(mom_comps.samples1)

#Split pairwise comparisons by number of samples so we can plot
mom_comps1.50samps <- mom_comps.samples1 %>% dplyr::filter(samples.per.yr == 50)
mom_comps1.150samps <- mom_comps.samples1 %>% dplyr::filter(samples.per.yr == 150)
mom_comps1.250samps <- mom_comps.samples1 %>% dplyr::filter(samples.per.yr == 250)

dad_comps1.50samps <- dad_comps.samples1 %>% dplyr::filter(samples.per.yr == 50)
dad_comps1.150samps <- dad_comps.samples1 %>% dplyr::filter(samples.per.yr == 150)
dad_comps1.250samps <- dad_comps.samples1 %>% dplyr::filter(samples.per.yr == 250)


#Histograms of mortality years
ggmort.list <- list()
ggmort.list[[1]] <- mom_comps1.50samps %>% gghistogram(x = "mort.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Mom mort years: 50 samples per year") + 
  xlim(0, 50)

ggmort.list[[2]] <- dad_comps1.50samps %>% gghistogram(x = "mort.yrs", fill = "lightblue", add = "mean") + 
  ggtitle("Dad mort years: 50 samples per year") + 
  xlim(0, 50)

ggmort.list[[3]] <- mom_comps1.150samps %>% gghistogram(x = "mort.yrs", fill = "lightblue", add = "mean") + 
  ggtitle("Mom mort years: 150 samples per year")+ 
  xlim(0, 50)

ggmort.list[[4]] <- dad_comps1.150samps %>% gghistogram(x = "mort.yrs", fill = "lightblue", add = "mean") + 
  ggtitle("Dad mort years: 150 samples per year")+ 
  xlim(0, 50)

ggmort.list[[5]] <- mom_comps1.250samps %>% gghistogram(x = "mort.yrs", fill = "lightblue", add = "mean") + 
  ggtitle("Mom mort years: 250 samples per year")+ 
  xlim(0, 50)

ggmort.list[[6]] <- dad_comps1.250samps %>% gghistogram(x = "mort.yrs", fill = "lightblue", add = "mean") + 
  ggtitle("Dad mort years: 250 samples per year")+ 
  xlim(0, 50)

ggarrange(plotlist = ggmort.list, ncol = 2, nrow = 3)



#Histograms of population growth gaps
ggPopGrowth.list <- list()
ggPopGrowth.list[[1]] <- mom_comps1.50samps %>% gghistogram(x = "pop.growth.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Mom pop growth years: 50 samples per year") + 
  xlim(-40, 10)

ggPopGrowth.list[[2]] <- dad_comps1.50samps %>% gghistogram(x = "pop.growth.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Dad pop growth years: 50 samples per year") + 
  xlim(-40, 10)

ggPopGrowth.list[[3]] <- mom_comps1.150samps %>% gghistogram(x = "pop.growth.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Mom pop growth years: 150 samples per year") + 
  xlim(-40, 10)

ggPopGrowth.list[[4]] <- dad_comps1.150samps %>% gghistogram(x = "pop.growth.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Dad pop growth years: 150 samples per year") + 
  xlim(-40, 10)

ggPopGrowth.list[[5]] <- mom_comps1.250samps %>% gghistogram(x = "pop.growth.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Mom pop growth years: 250 samples per year") + 
  xlim(-40, 10)

ggPopGrowth.list[[6]] <- dad_comps1.250samps %>% gghistogram(x = "pop.growth.yrs", fill = "lightblue", add = "mean") +
  ggtitle("Dad pop growth years: 250 samples per year") + 
  xlim(-40, 10)

ggarrange(plotlist = ggPopGrowth.list, ncol = 2, nrow = 3)

#NEXT TO DO:
#What is the distribution of number of sampled offspring per sampled parent per iteration?

head(mom_comps.samples1)
head(dad_comps.samples1)

#Save the total number of comparisons
all.HScomps <- mom_comps.samples1 %>% 
  dplyr::filter(type == "HS") %>% 
  group_by(samples.per.yr, iteration) %>% 
  summarize(all.comps = sum(all), all.pos.comps = sum(yes)) %>% 
  ungroup()

#Examine the number and percentage of comparisons for each sample size with a mortality gap greater than repro.age (which is the point at which it will become difficult to distinguish between HS and grandparent/grandchild)
mom_comps.gaps <- mom_comps.samples1 %>% 
  mutate(big.gap = ifelse(mort.yrs > repro.age, "yes", "no")) %>%
  group_by(samples.per.yr, type, big.gap, iteration) %>%
  summarize(gap.comps = sum(all), gap.pos = sum(yes)) %>% 
  dplyr::filter(type == "HS" & big.gap == "yes") %>% 
  left_join(all.HScomps, by = c("samples.per.yr", "iteration")) %>% 
  mutate(percent.gapcomps = (gap.comps/all.comps)*100) %>% 
  ungroup()

dad_comps.gaps <- dad_comps.samples1 %>% 
  mutate(big.gap = ifelse(mort.yrs > repro.age, "yes", "no")) %>%
  group_by(samples.per.yr, type, big.gap, iteration) %>%
  summarize(gap.comps = sum(all), gap.pos = sum(yes)) %>% 
  dplyr::filter(type == "HS" & big.gap == "yes") %>% 
  left_join(all.HScomps, by = c("samples.per.yr", "iteration")) %>% 
  mutate(percent.gapcomps = (gap.comps/all.comps)*100) %>% 
  ungroup()









#----------------------View prior and posterior----------------------------#
it <- 10 #Simulation iteration to visualize
#make dataframe of draws from prior distribution
tau = 1E-6
norm.prior <- rnorm(n = 10000, mean = 0, sd = sd)
unif.prior <- runif(n = 10000, min = 0, max = 1000)
sd <- sqrt(1/tau)
lam.tau <- 1/(0.02342^2)
lam.sd <- sqrt(1/lam.tau)
#tau <- 1/(sd^2)

N.prior <-  norm.prior %>% 
  as_tibble() %>% 
  mutate(Chain = "prior")

surv.prior <- rbeta(n = 10000, shape1 = 1, shape2 = 2) %>% 
  as_tibble() %>% 
  mutate(Chain = "prior")

lam.prior <- rnorm(n = 10000, mean = 1, sd = lam.sd) %>% 
  as_tibble() %>% 
  mutate(Chain = "prior")

#Prepare posteriors for plotting
Nf.post <- ggs(sim.samples.1_MCMC[[it]]) %>% 
  filter(Parameter == "Nf")

Nm.post <- ggs(sim.samples.1_MCMC[[it]]) %>% 
  filter(Parameter == "Nm")

surv.post <- ggs(sim.samples.1_MCMC[[it]]) %>% 
  filter(Parameter == "surv")

lam.post <- ggs(sim.samples.1_MCMC[[it]]) %>% 
  filter(Parameter == "lam")

Nf.density <- ggs_density(Nf.post) + 
  geom_density(data = N.prior, aes(x = value), alpha = 0.6)

Nm.density <- ggs_density(Nm.post) + 
  geom_density(data = N.prior, aes(x = value), alpha = 0.6)

surv.density <- ggs_density(surv.post) + 
  geom_density(data = surv.prior, aes(x = value), alpha = 0.6)

lam.density <- ggs_density(lam.post) + 
  geom_density(data = lam.prior, aes(x = value), alpha = 0.6)

ggarrange(Nf.density, Nm.density, surv.density, lam.density)








# #-------------------More elaborate/specific line plots-------------------------------#
# #Specify save location for pdf of plots
# coverage.file.400 <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.diagnostics/Plots/Statistical.coverage_", today, "_400Samps.pdf")
# 
# plotlist <- list()
# 
# #Create one page of coverage plots for each iteration
# for(it in 1:iterations){
#   HPD.temp <- HPD.400.4viz %>% filter(iteration == it)
# 
# gg.Fabundance <- ggplot(data = HPD.temp %>% filter(parameter == "Nf"), 
#                         aes(x = interval, y = truth)) + 
#   geom_linerange(aes(ymin = lower, ymax = upper, colour = parameter), size = 2) + 
#   geom_point(aes(fill = parameter), size = 2) + 
#   scale_color_manual(labels = "HPD interval", values = "maroon") + 
#   scale_fill_manual(labels = "Truth", values = "maroon") +
#   ggtitle("Female abundance") + 
#   scale_x_continuous(breaks = intervals) +
#   theme(legend.title = element_blank(),
#         axis.title.y = element_blank(), 
#         axis.title.x = element_text(vjust = -1), 
#         axis.text.x = element_text(angle = 45, vjust = .5))
# 
# gg.Mabundance <- ggplot(data = HPD.temp %>% filter(parameter == "Nm"), 
#                         aes(x = interval, y = truth)) + 
#   geom_linerange(aes(ymin = lower, ymax = upper, colour = parameter), size = 2) + 
#   geom_point(aes(fill = factor(parameter)), size = 2) + 
#   scale_color_manual(labels = "HPD interval", values = "mediumslateblue") + 
#   scale_fill_manual(labels = "Truth", values = "mediumslateblue") +
#   ggtitle("Male abundance") + 
#   scale_x_continuous(breaks = intervals) +
#   theme(legend.title = element_blank(),
#         legend.position = "none",
#         axis.title.y = element_blank(), 
#         axis.title.x = element_text(vjust = -1), 
#         axis.text.x = element_text(angle = 45, vjust = .5))
# 
# gg.survival <- ggplot(data = HPD.temp %>% filter(parameter == "surv"), 
#                       aes(x = interval, y = truth)) + 
#   geom_linerange(aes(ymin = lower, ymax = upper, colour = parameter), size = 2) + 
#   geom_point(aes(fill = factor(parameter), shape = "Truth"), size = 2) + 
#   scale_color_brewer(palette="PRGn") + 
#   ggtitle("Adult survival") + 
#   scale_x_continuous(breaks = intervals) +
#   theme(legend.title = element_blank(), 
#         axis.title.y = element_blank(), 
#         axis.title.x = element_text(vjust = -1), 
#         axis.text.x = element_text(angle = 45, vjust = .5))
# 
# plots <- ggarrange(gg.Fabundance, gg.Mabundance, gg.survival, common.legend = TRUE, legend = "right", vjust = -1)
# plotlist[[it]] <- annotate_figure(plots, top = text_grob(paste0("Statistical coverage for iteration ", it), size = 18, face = "bold"))
# }
# 
# ggexport(plotlist, filename = coverage.file.400)
