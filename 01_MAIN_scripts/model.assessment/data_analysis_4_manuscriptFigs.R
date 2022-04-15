library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)
library(ggmcmc)


rm(list=ls())
today <- format(Sys.Date(), "%d%b%Y")

#----------------Read in files ------------------------------
#------------- MCMC parameters ----------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains
MCMC.settings <- paste0("thin", nt, "_draw", ni, "_burn", nb)

#------------- Simulation parameters and labels  ----------------#
date.of.simulation1 <- "05Apr2022"
purpose1 <- "all.comps"
purpose1.lab <- "random sampling | all comparisons"

date.of.simulation2 <- "07Apr2022"
purpose2 <- "HS.PO_downsample"
purpose2.lab <- "random sampling | downsample"

date.of.simulation3 <- "11Apr2022"
purpose3 <- "HS.PO_refined.samples_all.comps"
purpose3.lab <- "targeted sampling | all comparisons"

date.of.simulation4 <- "11Apr2022"
purpose4 <- "HS.PO_refined.samples_downsample"
purpose4.lab <- "targeted sampling | downsample"

seeds <- "Seeds2022.03.23"
sim.samples.1 <- "200.samples"
sim.samples.2 <- "600.samples"
sim.samples.3 <- "1000.samples"
sim.samples.all <- c(200, 600, 1000)

MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.results/"
results_plots_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.results/figures/"
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
parents_prefix <- "parents_breakdown/CKMR_parents.breakdown"
sample.prefix <- "sample_info/CKMR_sample.info"
survival.prefix <- "survival/CKMR_survival"
pop.size.prefix <- "pop.size/CKMR_pop.size"

#------------------------------- Results -----------------------------------#
results.1 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1, ".csv")) %>% 
  mutate(model_type = "HS.PO", 
         purpose = purpose1,
         purpose.lab = purpose1.lab)

results.2 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2, ".csv")) %>% 
  mutate(model_type = "HS.PO",
         purpose = purpose2,
         purpose.lab = purpose2.lab)

 results.3 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation3, "_", seeds, "_", purpose3, ".csv")) %>% 
   mutate(model_type = "HS.PO",
          purpose = purpose3,
          purpose.lab = purpose3.lab)

 results.4 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation4, "_", seeds, "_", purpose4, ".csv")) %>%
   mutate(model_type = "HS.PO",
          purpose = purpose4,
          purpose.lab = purpose4.lab)


#------------------------------- MCMC output -----------------------------------#
#Trial 1
s1.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose1))

s1.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose1))

s1.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation1, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose1))

#Trial 2
s2.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose2))

s2.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose2))

s2.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation2, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose2))

#Trial 3
 s3.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose3))
 
 s3.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose3))
 
 s3.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation3, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose3))

#Trial 4
  s4.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose4))
  
  s4.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose4))
  
  s4.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation4, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose4))

#------------------------------- Population size details -----------------------------------#
pop.size.1 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1))

pop.size.2 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2))

pop.size.3 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose3))

pop.size.4 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose4))
####------------------------------- Quick analysis of results -----------------------------------####

results.all <- results.1 %>% 
  bind_rows(results.2, results.3, results.4) %>% 
  mutate(relative_bias = round(((Q50 - truth)/truth)*100,1)) %>%
  mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N")) %>% 
  mutate(percent_sampled = round((as.numeric(total_samples)/as.numeric(pop_size_mean)) * 100, 0)) %>% 
  mutate(percent_parents_sampled = as.numeric(unique_parents_in_sample)/as.numeric(mean_unique_parents_in_pop)) %>% 
  mutate(cv = (sd/mean)*100) %>% 
  mutate(purpose2 = ifelse(purpose == purpose1, "Down sampled", "All samples"))

jags_params <- c("Nf", "Nm", "surv", "lam") #Specify parameters

head(results.all)

results.all$purpose <- factor(results.all$purpose.lab, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

 
 #-----------Median Relative bias by sample size-------------------------#
results.all %>% group_by(total_samples, parameter, purpose) %>% 
#  dplyr::filter(model_type == "HS.PO") %>%
  dplyr::summarize(median = median(relative_bias), n = n()) %>% 
  arrange(desc(median))

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
results.all %>% group_by(purpose2, total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100) %>% 
  filter(parameter != "lam") %>% 
  arrange(desc(percent_in_interval))

#------------------------------- CV -----------------------------------#
purpose.combo <- "examine_sampling_combos"

results.all %>% group_by(model_type, purpose2, total_samples, parameter) %>% 
  dplyr::summarize(mean.cv = mean(cv), 
                   median.bias = median(relative_bias)) %>%
#  dplyr::filter(parameter == "Nf" | parameter == "Nm") %>% 
  arrange(mean.cv) #%>% 
  #write_csv(file = paste0(results_location, "Mean_bias_and_precision_", today, "_", purpose.combo, ".csv"))


####------------------------------- Relative bias density plots -----------------------------------####
#Try looking at distribution of relative bias to get an idea of how biased they are
param.plot.1 <- "Nf"
param.plot.2 <- "Nm"
param.plot.3 <- "survival" #Which parameters are we plotting below?
param.plot.4 <- "lambda"

sim.samples.all <- c(200, 600, 1000) #vector of total sample sizes


####------------------------------- Calculate different HPDI intervals-------------------------------#### 
#How many estimates fall in different HPD intervals?
#Calculate HPD interval for each iteration
intervals <- c(seq(from = 0.95, to = 0.05, by = -0.05)) #Specify HPDI intervals to calculate
iterations = 100

#Assumes four "purposes" are being compared
#Initialize dataframes for loop below
HPD.200.1 <- NULL
HPD.600.1 <- NULL
HPD.1000.1 <- NULL
HPD.200.2 <- NULL
HPD.600.2 <- NULL
HPD.1000.2 <- NULL
HPD.200.3 <- NULL
HPD.600.3 <- NULL
HPD.1000.3 <- NULL
HPD.200.4 <- NULL
HPD.600.4 <- NULL
HPD.1000.4 <- NULL

for(i in 1:length(s1.1)){ #Loop over all the posterior samples from each iteration; length should be same as number of iterations
  for(j in 1:length(intervals)){ #Loop over the different intervals for calculating HPDI
    
    #---------------------------Purpose 1--------------------------------#
    #200 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.200.1 <- combine.mcmc(s1.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.1 <- post.HPD.200.1[row.names(post.HPD.200.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200.1 <- rbind(HPD.200.1, post.HPD.200.1)
    
    #600 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.600.1 <- combine.mcmc(s1.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.1 <- post.HPD.600.1[row.names(post.HPD.600.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.1 <- rbind(HPD.600.1, post.HPD.600.1)

    #1000 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.1000.1 <- combine.mcmc(s1.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.1 <- post.HPD.1000.1[row.names(post.HPD.1000.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.1 <- rbind(HPD.1000.1, post.HPD.1000.1)
    
    
    #---------------------------Purpose 2--------------------------------#    
    #200 samples
    #Calcualte HPD interval for iteration i and interval j
    
    post.HPD.200.2 <- combine.mcmc(s2.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.2 <- post.HPD.200.2[row.names(post.HPD.200.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.2 <- rbind(HPD.200.2, post.HPD.200.2)
    
    #600 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.600.2 <- combine.mcmc(s2.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.2 <- post.HPD.600.2[row.names(post.HPD.600.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.2 <- rbind(HPD.600.2, post.HPD.600.2)
    
    #1000 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.1000.2 <- combine.mcmc(s2.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.2 <- post.HPD.1000.2[row.names(post.HPD.1000.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.2 <- rbind(HPD.1000.2, post.HPD.1000.2)
    
    # #---------------------------Purpose 3--------------------------------#
    #200 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.200.3 <- combine.mcmc(s3.1[[i]]) %>%
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>%
      mutate(interval = intervals[j], iteration = i)

    post.HPD.200.3 <- post.HPD.200.3[row.names(post.HPD.200.3) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200.3 <- rbind(HPD.200.3, post.HPD.200.3)

    #600 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.600.3 <- combine.mcmc(s3.2[[i]]) %>%
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>%
      mutate(interval = intervals[j], iteration = i)

    post.HPD.600.3 <- post.HPD.600.3[row.names(post.HPD.600.3) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.600.3 <- rbind(HPD.600.3, post.HPD.600.3)

    #1000 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.1000.3 <- combine.mcmc(s3.3[[i]]) %>%
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>%
      mutate(interval = intervals[j], iteration = i)

    post.HPD.1000.3 <- post.HPD.1000.3[row.names(post.HPD.1000.3) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.1000.3 <- rbind(HPD.1000.3, post.HPD.1000.3)

    # #---------------------------Purpose 4--------------------------------#
    #200 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.200.4 <- combine.mcmc(s4.1[[i]]) %>%
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>%
      mutate(interval = intervals[j], iteration = i)

    post.HPD.200.4 <- post.HPD.200.4[row.names(post.HPD.200.4) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200.4 <- rbind(HPD.200.4, post.HPD.200.4)

    #600 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.600.4 <- combine.mcmc(s4.2[[i]]) %>%
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>%
      mutate(interval = intervals[j], iteration = i)

    post.HPD.600.4 <- post.HPD.600.4[row.names(post.HPD.600.4) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.600.4 <- rbind(HPD.600.4, post.HPD.600.4)

    #1000 samples
    #Calcualte HPD interval for iteration i and interval j
    post.HPD.1000.4 <- combine.mcmc(s4.3[[i]]) %>%
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>%
      mutate(interval = intervals[j], iteration = i)

    post.HPD.1000.4 <- post.HPD.1000.4[row.names(post.HPD.1000.4) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.1000.4 <- rbind(HPD.1000.4, post.HPD.1000.4)
    
  }
  print(paste0("Finished with iteration ", i))
  }
    


####------------------------------- Create dataframes for HPDI plotting------------------------------#### 

#---------------------------Purpose 1 --------------------------------#    
#Separate out by sample size
HPD.200.1.4viz <- results.1 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.1, by = c("parameter", "iteration"))

HPD.600.1.4viz <- results.1 %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.1, by = c("parameter", "iteration"))

HPD.1000.1.4viz <- results.1 %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.1, by = c("parameter", "iteration"))


#Calculate percent estimates in each interval
HPD.200.1.summary <- HPD.200.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.1.summary <- HPD.600.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.1.summary <- HPD.1000.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes    
HPD.1.summary_all <- HPD.200.1.summary %>% inner_join(HPD.600.1.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.1.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.1.summary.tidy <- HPD.1.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
    ) %>% 
  mutate(model_type = "HS.PO",
         purpose = purpose1.lab)

#Save
#write_csv(HPD.1.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation, "_", purpose1, ".csv"))


#---------------------------Purpose 2 --------------------------------#    
#Separate out by sample size
HPD.200.2.4viz <- results.2 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.2, by = c("parameter", "iteration"))

HPD.600.2.4viz <- results.2 %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.2, by = c("parameter", "iteration"))

HPD.1000.2.4viz <- results.2 %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.2, by = c("parameter", "iteration"))


#Calculate percent estimates in each interval
HPD.200.2.summary <- HPD.200.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.2.summary <- HPD.600.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.2.summary <- HPD.1000.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes
HPD.2.summary_all <- HPD.200.2.summary %>% inner_join(HPD.600.2.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.2.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.2.summary.tidy <- HPD.2.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = "HS.only",
                  purpose = purpose2.lab)

#Save
#write_csv(HPD.2.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation, "_", purpose2, ".csv"))


#---------------------------Purpose 3 --------------------------------# 
#Separate out by sample size
HPD.200.3.4viz <- results.3 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.3, by = c("parameter", "iteration"))

HPD.600.3.4viz <- results.3 %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.3, by = c("parameter", "iteration"))

HPD.1000.3.4viz <- results.3 %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.3, by = c("parameter", "iteration"))


#Calculate percent estimates in each interval
HPD.200.3.summary <- HPD.200.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.3.summary <- HPD.600.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.3.summary <- HPD.1000.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes
HPD.3.summary_all <- HPD.200.3.summary %>% inner_join(HPD.600.3.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.3.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.3.summary.tidy <- HPD.3.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = "HS.PO",
         purpose = purpose3.lab)

#Save
#write_csv(HPD.3.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation, "_", purpose3, ".csv"))


#---------------------------Purpose 4 --------------------------------# 
#Separate out by sample size
HPD.200.4.4viz <- results.4 %>% filter(total_samples == 200) %>%
  right_join(HPD.200.4, by = c("parameter", "iteration"))

HPD.600.4.4viz <- results.4 %>% filter(total_samples == 600) %>%
  right_join(HPD.600.4, by = c("parameter", "iteration"))

HPD.1000.4.4viz <- results.4 %>% filter(total_samples == 1000) %>%
  right_join(HPD.1000.4, by = c("parameter", "iteration"))


#Calculate percent estimates in each interval
HPD.200.4.summary <- HPD.200.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.4.summary <- HPD.600.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.4.summary <- HPD.1000.4.4viz %>% group_by(parameter, interval) %>%
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>%
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Combine above dataframes
HPD.4.summary_all <- HPD.200.4.summary %>% inner_join(HPD.600.4.summary, by = c("parameter", "interval")) %>%
  inner_join(HPD1000.4.summary, by = c("parameter", "interval")) %>%
  mutate(interval.scale = interval*100)

HPD.4.summary.tidy <- HPD.4.summary_all %>%
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>%
  mutate(model_type = "HS.only",
         purpose = purpose4.lab)

#Save
#write_csv(HPD.4.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation, "_", purpose4, ".csv"))



#---------------Create dataframes to plot by sample size--------#
HPD.200.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "200.samples")
HPD.600.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "1000.samples")

HPD.200.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "200.samples")
HPD.600.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "1000.samples")

HPD.200.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "200.samples")
HPD.600.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "1000.samples")

HPD.200.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "200.samples")
HPD.600.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "1000.samples")


#Combine above dataframes, and re-order purpose as factors
all.200.4viz <- HPD.200.1.4viz %>% bind_rows(HPD.200.2.4viz, HPD.200.3.4viz, HPD.200.4.4viz) 
all.200.4viz$purpose <- factor(all.200.4viz$purpose, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

all.600.4viz <- HPD.600.1.4viz %>% bind_rows(HPD.600.2.4viz, HPD.600.3.4viz, HPD.600.4.4viz) 
all.600.4viz$purpose <- factor(all.600.4viz$purpose, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))

all.1000.4viz <- HPD.1000.1.4viz %>% bind_rows(HPD.1000.2.4viz, HPD.1000.3.4viz, HPD.1000.4.4viz)

all.1000.4viz$purpose <- factor(all.1000.4viz$purpose, levels = c(purpose1.lab, purpose2.lab, purpose3.lab, purpose4.lab))




#-----------------------Make figures-----------------------------
#-----------------HPDI plot--------------------#
#Nm & surv - 200 samples
#Random (non-targeted) sampling: all comparisons and downsample
p1.200.df <- all.200.4viz %>% filter(parameter == "Nm" | parameter == "surv" | parameter == "Nf") %>% 
  dplyr::filter(purpose == purpose1.lab | purpose == purpose3.lab) %>% 
  mutate(sampling = ifelse(purpose == purpose1.lab, "Indiscriminate sampling", "Targeted sampling"))

p1 <- ggplot(data = p1.200.df, 
             aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = sampling, fill = sampling, shape = sampling), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = sampling, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("A) Percent in HPDI") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank()) + 
  facet_wrap(~parameter)

samps <- 200
p2.200.df <- results.all %>% dplyr::filter(parameter == "Nm" | parameter == "surv" | parameter == "Nf", total_samples == samps, purpose == purpose1.lab | purpose == purpose3.lab) %>% 
  mutate(sampling = ifelse(purpose == purpose1.lab, "Indiscriminate sampling", "Targeted sampling"))

p2 <- p2.200.df %>% ggplot(aes(x = relative_bias, color = sampling, fill = sampling)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("B) Relative bias of posterior median")) + 
  xlim(-75, 75) + 
  facet_wrap(~parameter) +
  theme(legend.title = element_blank())



#Nm & surv - 1000 samples
#random (non-targeted) sampling: all comps and downsample
p3.1000.df <- all.1000.4viz %>% filter(parameter == "Nm" | parameter == "surv" | parameter == "Nf") %>% 
  dplyr::filter(purpose == purpose1.lab | purpose == purpose2.lab) %>% 
  mutate(comparisons = ifelse(purpose == purpose1.lab, "All comparisons", "Downsample"))

p3 <- ggplot(data =  p3.1000.df,
             aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = comparisons, fill = comparisons, shape = comparisons), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = comparisons, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Indiscriminate sampling") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank()) + 
  facet_wrap(~parameter)

#Targeted sampling: all comps and downsample
p4.1000.df <- all.1000.4viz %>% filter(parameter == "Nm" | parameter == "surv" | parameter == "Nf") %>% 
  dplyr::filter(purpose == purpose3.lab | purpose == purpose4.lab) %>% 
  mutate(comparisons = ifelse(purpose == purpose3.lab, "All comparisons", "Downsample"))

p4 <- ggplot(data = p4.1000.df, 
                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = comparisons, fill = comparisons, shape = comparisons), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = comparisons, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Targeted sampling") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank()) + 
  facet_wrap(~parameter)


#1000 samples, indiscriminate sampling
samps <- 1000
p5.1000.df <- results.all %>% dplyr::filter(parameter == "Nm" | parameter == "surv" | parameter == "Nf", total_samples == samps, purpose == purpose1.lab | purpose == purpose2.lab) %>% 
  mutate(comparisons = ifelse(purpose == purpose1.lab, "All comparisons", "Downsample"))

p5 <- p5.1000.df %>% ggplot(aes(x = relative_bias, color = comparisons, fill = comparisons)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("C) Indiscriminate sampling")) + 
  xlim(-50, 50) + 
  facet_wrap(~parameter) +
  theme(legend.title = element_blank())

#1000 samples, targeted sampling
p6.1000.df <- results.all %>% dplyr::filter(parameter == "Nm" | parameter == "surv" | parameter == "Nf", total_samples == samps, purpose == purpose3.lab | purpose == purpose4.lab) %>% 
  mutate(comparisons = ifelse(purpose == purpose3.lab, "All comparisons", "Downsample"))

p6 <- p6.1000.df %>% ggplot(aes(x = relative_bias, color = comparisons, fill = comparisons)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("D) Targeted sampling")) + 
  xlim(-50, 50) + 
  facet_wrap(~parameter) +
  theme(legend.title = element_blank())


ggarrange(p1, p2, common.legend = TRUE, ncol = 1, nrow = 2) #Random sampling across age classes does better when few samples are taken
ggarrange(p3, p4, common.legend = TRUE, ncol = 1, nrow = 2) #Targeted sampling of specific age classes performs better when sampling heavily.
ggarrange(p5, p6, common.legend = TRUE, ncol = 1, nrow = 2, legend = "bottom")

####NEXT: Separate 1000 density plot into indiscriminate and random sampling



####---------------------Density plots of kin detected----------------------------------------------####
#---------------------Prepare data for plotting----------------------------------------------#
#How many POPs detected per sample size?
results.all %>% 
  group_by(purpose2, total_samples) %>% 
  summarize(mean(POPs_detected), mean(HSPs_detected))

#Format PO dataframe
results.PO <- results.all %>% dplyr::select(parameter, total_samples, kin.detected = POPs_detected,  unique_parents_in_sample, iteration, cv, purpose = purpose2) %>% 
  mutate(type = "PO")

#Format HS dataframe
results.HS <- results.all %>% dplyr::select(parameter, total_samples, kin.detected = HSPs_detected,  unique_parents_in_sample, iteration, cv, purpose = purpose2) %>% 
  mutate(type = "HS")

#Combine PO and HS dataframe for plotting
results.kin <- rbind(results.PO, results.HS) %>% 
  mutate(total_samples = paste0(total_samples, " samples"))

results.kin$total_samples <- factor(results.kin$total_samples, levels = c("200 samples", "600 samples", "1000 samples"))
results.kin$purpose <- factor(results.kin$purpose, levels = c("Down sampled", "All samples"))

#---------------------Create density plots----------------------------------------------#
#Maternal
k1 <- results.kin %>% dplyr::filter(parameter == "Nf", type == "PO") %>%
  ggplot(aes(x = kin.detected, color = purpose, fill = purpose)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "MOPs detected") + 
  xlim(0, 60) + 
  facet_wrap(~total_samples, ncol = 1) + 
  theme(legend.title = element_blank())

k2 <- results.kin %>% dplyr::filter(parameter == "Nf", type == "HS") %>%
  ggplot(aes(x = kin.detected, color = purpose, fill = purpose)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "Maternal HSPs detected") + 
  xlim(0, 500) + 
  facet_wrap(~total_samples, ncol = 1) + 
  theme(legend.title = element_blank())

m.kin <- ggarrange(k1, k2, ncol = 2, common.legend = TRUE)
annotate_figure(m.kin, fig.lab = "Mother kin pairs", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")



#Paternal
k3 <- results.kin %>% dplyr::filter(parameter == "Nm", type == "PO") %>%
  ggplot(aes(x = kin.detected, color = purpose, fill = purpose)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = "FOPs detected") + 
  xlim(0, 60) + 
  facet_wrap(~total_samples, ncol = 1) + 
  theme(legend.title = element_blank())

k4 <- results.kin %>% dplyr::filter(parameter == "Nm", type == "HS") %>%
  ggplot(aes(x = kin.detected, color = purpose, fill = purpose)) + 
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
# EO1 <- results.all2 %>% dplyr::filter(purpose == purpose1.lab) %>% 
#   ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
#   geom_point() +
#   facet_wrap(~sample.label) + 
#   labs(title = purpose1.lab, x = "Percent difference", y = "Relative bias") + 
#   theme(legend.title = element_blank())
# 
# #Purpose 2
# EO2 <- results.all2 %>% dplyr::filter(purpose == purpose2) %>% 
#   ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
#   geom_point() +
#   facet_wrap(~sample.label) + 
#   labs(title = purpose2.lab, x = "Percent difference", y = "Relative bias")

#Purpose 3
EO3 <- results.all2 %>% dplyr::filter(purpose == purpose3.lab, parameter != "lam") %>% 
  ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
  geom_point() +
  facet_wrap(~sample.label) + 
  labs(x = "Percent difference in expected vs observed kin pairs", y = "Relative bias")

#Purpose 4
# EO4 <- results.all2 %>% dplyr::filter(purpose == purpose4) %>% 
#   ggplot(aes(x = all_exp_obs_diff, y = relative_bias, fill = parameter, color = parameter)) +
#   geom_point() +
#   facet_wrap(~sample.label) + 
#   labs(title = purpose4.lab, x = "Percent difference", y = "Relative bias")
# 
# 
# 
# EO.all <- ggarrange(EO1, EO2, EO3, EO4, common.legend = TRUE)
# annotate_figure(EO.all, fig.lab = "Relative bias by percent difference in expected vs observed kin pairs", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")
