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
date.of.simulation1 <- "07Apr2022"
purpose1 <- "HS.PO_downsample"

date.of.simulation2 <- "05Apr2022"
purpose2 <- "all.comps"
#purpose3 <- "all.comps"
#purpose4 <- "5yr.comps"

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
         purpose = purpose1)

results.2 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation2, "_", seeds, "_", purpose2, ".csv")) %>% 
  mutate(model_type = "HS.only",
         purpose = purpose2)

# results.3 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose3, ".csv")) %>% 
#   mutate(model_type = "HS.PO",
#          purpose = purpose3)

 # results.4 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose4, ".csv")) %>% 
 #   mutate(model_type = "HS.PO",
 #          purpose = purpose4)


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
# date.of.simulation <- "05Apr2022"
# s3.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose3))
# 
# s3.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose3))
# 
# s3.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose3))

#Trial 4
 # s4.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose4))
 # 
 # s4.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose4))
 # 
 # s4.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose4))

#------------------------------- Population size details -----------------------------------#
pop.size.1 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation1, "_", seeds, "_", purpose1))

pop.size.2 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose2))

pop.size.3 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose3))

pop.size.4 <- readRDS(paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose4))
####------------------------------- Quick analysis of results -----------------------------------####

results.all <- results.1 %>% 
  bind_rows(results.2) %>% 
  mutate(relative_bias = round(((Q50 - truth)/truth)*100,1)) %>%
  mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N")) %>% 
  mutate(percent_sampled = round((as.numeric(total_samples)/as.numeric(pop_size_mean)) * 100, 0)) %>% 
  mutate(percent_parents_sampled = as.numeric(unique_parents_in_sample)/as.numeric(mean_unique_parents_in_pop)) %>% 
  mutate(cv = (sd/mean)*100) %>% 
  mutate(purpose2 = ifelse(purpose == purpose1, "Down sampled", "All samples"))

jags_params <- c("Nf", "Nm", "surv", "lam") #Specify parameters

head(results.all)

results.all$purpose <- factor(results.all$purpose, levels = c(purpose1, purpose2))

 
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
purpose.combo <- "test.downsample"

results.all %>% group_by(model_type, purpose2, total_samples, parameter) %>% 
  dplyr::summarize(mean.cv = mean(cv), 
                   median.bias = median(relative_bias)) %>%
#  dplyr::filter(parameter == "Nf" | parameter == "Nm") %>% 
  arrange(mean.cv) %>% 
  write_csv(file = paste0(results_location, "Mean_bias_and_precision_", today, "_", purpose.combo, ".csv"))


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
    # #200 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.200.3 <- combine.mcmc(s3.1[[i]]) %>% 
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>% 
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.200.3 <- post.HPD.200.3[row.names(post.HPD.200.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.200.3 <- rbind(HPD.200.3, post.HPD.200.3)
    # 
    # #600 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.600.3 <- combine.mcmc(s3.2[[i]]) %>% 
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>% 
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.600.3 <- post.HPD.600.3[row.names(post.HPD.600.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.600.3 <- rbind(HPD.600.3, post.HPD.600.3)
    # 
    # #1000 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.1000.3 <- combine.mcmc(s3.3[[i]]) %>% 
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>% 
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.1000.3 <- post.HPD.1000.3[row.names(post.HPD.1000.3) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.1000.3 <- rbind(HPD.1000.3, post.HPD.1000.3)

    # #---------------------------Purpose 4--------------------------------#    
    # #200 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.200.4 <- combine.mcmc(s4.1[[i]]) %>% 
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>% 
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.200.4 <- post.HPD.200.4[row.names(post.HPD.200.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.200.4 <- rbind(HPD.200.4, post.HPD.200.4)
    # 
    # #600 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.600.4 <- combine.mcmc(s4.2[[i]]) %>% 
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>% 
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.600.4 <- post.HPD.600.4[row.names(post.HPD.600.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.600.4 <- rbind(HPD.600.4, post.HPD.600.4)
    # 
    # #1000 samples
    # #Calcualte HPD interval for iteration i and interval j
    # post.HPD.1000.4 <- combine.mcmc(s4.3[[i]]) %>% 
    #   HPDinterval(prob = intervals[j]) %>%
    #   data.frame() %>% 
    #   mutate(interval = intervals[j], iteration = i)
    # 
    # post.HPD.1000.4 <- post.HPD.1000.4[row.names(post.HPD.1000.4) %in% jags_params,] %>% #Remove deviance
    #   rownames_to_column(var = "parameter")
    # 
    # HPD.1000.4 <- rbind(HPD.1000.4, post.HPD.1000.4)
    
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
         purpose = purpose1)

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
                  purpose = purpose2)

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
         purpose = purpose3)

#Save
write_csv(HPD.3.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation, "_", purpose3, ".csv"))


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
         purpose = purpose4)

#Save
write_csv(HPD.4.summary_all, file = paste0(results_location, "HPD.summaries/HPD.summary_", date.of.simulation, "_", purpose4, ".csv"))



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
all.200.4viz <- HPD.200.1.4viz %>% bind_rows(HPD.200.2.4viz) %>% 
  mutate(purpose = ifelse(purpose == purpose1, "Down sampled", "All samples"))

all.200.4viz$purpose <- factor(all.200.4viz$purpose, levels = c("Down sampled", "All samples"))

all.600.4viz <- HPD.600.1.4viz %>% bind_rows(HPD.600.2.4viz) %>% 
  mutate(purpose = ifelse(purpose == purpose1, "Down sampled", "All samples"))

all.600.4viz$purpose <- factor(all.600.4viz$purpose, levels = c("Down sampled", "All samples"))

all.1000.4viz <- HPD.1000.1.4viz %>% bind_rows(HPD.1000.2.4viz) %>% 
  mutate(purpose = ifelse(purpose == purpose1, "Down sampled", "All samples"))

all.1000.4viz$purpose <- factor(all.1000.4viz$purpose, levels = c("Down sampled", "All samples"))


#-----------------------Make figures-----------------------------
#-----------------HPDI plot--------------------#
#Nf & surv - 200 samples
p1 <- ggplot(data = all.200.4viz %>% filter(parameter == "Nf" | parameter == "surv"), 
                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Percent in HPDI: 200 samples") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank()) + 
  facet_wrap(~parameter)

#Nf & surv - 1000 samples
p2 <- ggplot(data = all.1000.4viz %>% filter(parameter == "Nf" | parameter == "surv"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Percent in HPDI: 1000 samples") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank()) + 
  facet_wrap(~parameter)


#-----------------------relative bias density plots-----------------------------#
samps <- 200

p3 <- results.all %>% dplyr::filter(parameter == "Nf" | parameter == "surv", total_samples == samps) %>% 
  ggplot(aes(x = relative_bias, color = purpose2, fill = purpose2)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("C) Relative bias: ", samps, " samples")) + 
  xlim(-50, 50) + 
  facet_wrap(~parameter) +
  theme(legend.title = element_blank())

samps <- 1000
p4 <- results.all %>% dplyr::filter(parameter == "Nf" | parameter == "surv", total_samples == samps) %>% 
  ggplot(aes(x = relative_bias, color = purpose2, fill = purpose2)) + 
  geom_density(alpha = 0.6) + 
  ggtitle(label = paste0("D) Relative bias: ", samps, " samples")) + 
  xlim(-50, 50) + 
  facet_wrap(~parameter) +
  theme(legend.title = element_blank())



ggarrange(p1, p2, p3, p4, common.legend = TRUE, ncol = 1, nrow = 4)







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
