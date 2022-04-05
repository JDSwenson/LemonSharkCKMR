library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)
library(ggmcmc)


rm(list=ls())
today <- format(Sys.Date(), "%d%b%Y")

#----------------Read in files ------------------------------
#Check results from model diagnostics
# Results
# results <- read_csv("G://My Drive/Personal_Drive/R/CKMR/Model.results/Model.validation/CKMR_results_08Dec2021_longChain.csv")
# 
# head(results)
# 
# results <- results %>% mutate(relative_bias2 = round(((`50` - truth)/truth)*100,1))
# 
# #Median relative bias by sample size
# results %>% group_by(total_samples, parameter) %>% 
#   dplyr::summarize(median = median(relative_bias2), n = n())
repro.age <- 12
max.age <- 50
burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs <- burn.in + Num.years #Total number of simulation years
estimation.year <- n_yrs - 5 # Set year of estimation

####------------- MCMC parameters ----------------####
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains
MCMC.settings <- paste0("thin", nt, "_draw", ni, "_burn", nb)


date.of.simulation <- "04Apr2022"
purpose1 <- "HS.PO_new.compsDF"
purpose2 <- "HS.only_new.compsDF"
purpose3 <- "HS.PO_surv.prior"
#purpose4 <- "HS.only_one.indv.per.parent"

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


#Results
results.1 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose1, ".csv")) %>% 
  mutate(model_type = "HS.PO", 
         purpose = purpose1)

results.2 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose2, ".csv")) %>% 
  mutate(model_type = "HS.only",
         purpose = purpose2)

results.3 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose3, ".csv")) %>% 
  mutate(model_type = "HS.PO",
         purpose = purpose3)

# results.4 <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose4, ".csv")) %>% 
#   mutate(model_type = "HS.only",
#          litter.samples = "one")


#MCMC samples/output
#Trial 1
s1.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose1))

s1.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose1))

s1.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose1))

#Trial 2
s2.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose2))

s2.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose2))

s2.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose2))

#Trial 3
s3.1 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose3))

s3.2 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose3))

s3.3 <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose3))

#Trial 4
# s1.HS.only_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose4))
# 
# s2.HS.only_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose4))
# 
# s3.HS.only_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose4))

#Examine results quickly
results.all <- results.1 %>% 
  bind_rows(results.2, 
            results.3) %>% 
  mutate(relative_bias = round(((Q50 - truth)/truth)*100,1)) %>%
  mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N")) %>% 
  mutate(percent_sampled = round((as.numeric(total_samples)/as.numeric(pop_size_mean)) * 100, 0)) %>% 
  mutate(percent_parents_sampled = as.numeric(unique_parents_in_sample)/as.numeric(mean_unique_parents_in_pop)) %>% 
  mutate(cv = (sd/mean)*100)

# results.all %>% group_by(model_type, total_samples, parameter, in_interval) %>% 
#   summarize(num.in.interval = n()) %>% 
#   arrange(desc(num.in.interval)) %>% 
#   dplyr::filter(parameter != "lam") %>% 
#   View()



#----------------------Quick analysis of results----------------------------
#today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
jags_params <- c("Nf", "Nm", "surv", "lam") #Specify parameters

head(results.all)



#Temporary fix since results saved funny ... 
#Lambda and survival were switched, so switching them back ... 
# HPD2.5_surv.norm <- results.norm %>% filter(parameter == "lam") %>% 
#   pull(HPD2.5)
# 
# HPD2.5_lam.norm <- results.norm %>% filter(parameter == "surv") %>% 
#   pull(HPD2.5)  
# 
# 
# HPD97.5_surv.norm <- results.norm %>% filter(parameter == "lam") %>% 
#   pull(HPD97.5)
# 
# HPD97.5_lam.norm <- results.norm %>% filter(parameter == "surv") %>% 
#   pull(HPD97.5)  
# 
# results.norm[results.norm$parameter == "lam",]$HPD2.5 <- HPD2.5_lam.norm
# results.norm[results.norm$parameter == "surv",]$HPD2.5 <- HPD2.5_surv.norm
# 
# results.norm[results.norm$parameter == "lam",]$HPD97.5 <- HPD97.5_lam.norm
# results.norm[results.norm$parameter == "surv",]$HPD97.5 <- HPD97.5_surv.norm
# 
# results2.norm <- results.norm %>% dplyr::select(-c(in_interval, relative_bias))

#Calculate relative bias
 # results2.norm <- results2.norm %>% 
 #   mutate(relative_bias = round(((Q50 - truth)/truth)*100,1)) %>%
 #   mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N")) 
 #%>% 
  # mutate(percent_sampled = round((total_samples/pop_size_mean) * 100, 0)

 
 #-----------Median Relative bias by sample size-------------------------#
results.all %>% group_by(model_type, total_samples, parameter, purpose) %>% 
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


#-----------Within HPDI?-------------------------#
 #Within HPD interval?
results.all %>% group_by(model_type, purpose, total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100) %>% 
  filter(parameter != "lam") %>% 
  arrange(percent_in_interval)

#---------------CV------------------------#
results.all %>% group_by(model_type, purpose, total_samples, parameter) %>% 
  dplyr::summarize(mean.cv = mean(cv)) %>%
  dplyr::filter(parameter == "Nf" | parameter == "Nm") %>% 
  arrange(mean.cv) %>% 
  View()


#Violin plots
  # RB_violin.p1 <- results.all %>% dplyr::filter(total_samples == samps & purpose == purpose1) %>% 
  #   ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
  # geom_violin(aes(y=relative_bias), draw_quantiles = 0.5) +
  # #ylim(-50, 160) +
  # geom_hline(yintercept=0, col="black", size=1.25) +
  # #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  # labs(x="Iteration", y="Relative bias", title="Relative bias by iteration") +
  # scale_fill_brewer(palette="Set2") +
  # font("title", size = 10, face = "bold")

#---------------------------Box plots----------------------------------
#Specify save location for pdf of plots
boxPlots.file <- paste0(results_plots_location, "BoxPlots_newDF_", today, ".pdf")
pdf(file = boxPlots.file, width = 14, height = 4)
sim.samples.all <- c(200, 600, 1000)
plot.list <- NULL

for(i in 1:length(sim.samples.all)){
  samps <- sim.samples.all[i]
  
  #boxplot plot
  RB_boxplot.s1 <- results.all %>% filter(total_samples == samps, purpose == purpose1) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose1, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 100)
  
  RB_boxplot.s2 <- results.all %>% filter(total_samples == samps, purpose == purpose2) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose2, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 100)

  
  RB_boxplot.s3 <- results.all %>% filter(total_samples == samps, purpose == purpose3) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Parameter", y="Relative bias", title=paste0(purpose3, "_", samps, " samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 100)

# plot.list[[1]] <- RB_boxplot.norm
# plot.list[[j+1]] <- RB_boxplot.uniform
# plot.list[[j+2]] <- RB_boxplot.hier
# print(ggarrange(plotlist = plot.list, common.legend = TRUE, nrow = 1, ncol = 3, widths = c(1, 1, 1))) #Need to wrap

print(ggarrange(RB_boxplot.s1, RB_boxplot.s2, RB_boxplot.s3, nrow = 1, ncol = 3, common.legend = TRUE))
  }

dev.off()

#-------------------Calculate HPD intervals---------------------------
# How many estimates fall in different HPD intervals?
#Calculate HPD interval for each iteration
intervals <- c(seq(from = 0.95, to = 0.05, by = -0.05))
iterations = 100

HPD.200.1 <- NULL
HPD.600.1 <- NULL
HPD.1000.1 <- NULL
HPD.200.2 <- NULL
HPD.600.2 <- NULL
HPD.1000.2 <- NULL
HPD.200.3 <- NULL
HPD.600.3 <- NULL
HPD.1000.3 <- NULL
# HPD.200.4 <- NULL
# HPD.600.4 <- NULL
# HPD.1000.4 <- NULL

for(i in 1:length(s1.1)){
  for(j in 1:length(intervals)){
    #HS.PO_new.compsDF
    #200 samples
    post.HPD.200.1 <- combine.mcmc(s1.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.1 <- post.HPD.200.1[row.names(post.HPD.200.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200.1 <- rbind(HPD.200.1, post.HPD.200.1)
    
    #1
    #600 samples
    post.HPD.600.1 <- combine.mcmc(s1.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.1 <- post.HPD.600.1[row.names(post.HPD.600.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.1 <- rbind(HPD.600.1, post.HPD.600.1)

    #1
    #1000 samples
    post.HPD.1000.1 <- combine.mcmc(s1.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.1 <- post.HPD.1000.1[row.names(post.HPD.1000.1) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.1 <- rbind(HPD.1000.1, post.HPD.1000.1)
    
    #HS.only_new.compsDF
    #200 samples
    post.HPD.200.2 <- combine.mcmc(s2.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.2 <- post.HPD.200.2[row.names(post.HPD.200.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.2 <- rbind(HPD.200.2, post.HPD.200.2)
    
    #2
    #600 samples
    post.HPD.600.2 <- combine.mcmc(s2.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.2 <- post.HPD.600.2[row.names(post.HPD.600.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.2 <- rbind(HPD.600.2, post.HPD.600.2)
    
    #2
    #1000 samples
    post.HPD.1000.2 <- combine.mcmc(s2.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.2 <- post.HPD.1000.2[row.names(post.HPD.1000.2) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.2 <- rbind(HPD.1000.2, post.HPD.1000.2)
    
    #3
    #200 samples
    post.HPD.200.3 <- combine.mcmc(s3.1[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.3 <- post.HPD.200.3[row.names(post.HPD.200.3) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.3 <- rbind(HPD.200.3, post.HPD.200.3)
    
    #3
    #600 samples
    post.HPD.600.3 <- combine.mcmc(s3.2[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.3 <- post.HPD.600.3[row.names(post.HPD.600.3) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.3 <- rbind(HPD.600.3, post.HPD.600.3)
    
    #3
    #1000 samples
    post.HPD.1000.3 <- combine.mcmc(s3.3[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.3 <- post.HPD.1000.3[row.names(post.HPD.1000.3) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.3 <- rbind(HPD.1000.3, post.HPD.1000.3)

  }
  print(paste0("Finished with iteration ", i))
  }
    


#--------------Create dataframes for viz and analysis----------------------

#--------------------------HS.PO_new.compsDF-----------------------------------#
HPD.200.1.4viz <- results.1 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.1, by = c("parameter", "iteration"))

HPD.600.1.4viz <- results.1 %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.1, by = c("parameter", "iteration"))

HPD.1000.1.4viz <- results.1 %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.1, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.1.summary <- HPD.200.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.1.summary <- HPD.600.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.1.summary <- HPD.1000.1.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
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

write_csv(HPD.1.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose1, ".csv"))


#--------------------------HS.only_new.compsDF-----------------------------------#
HPD.200.2.4viz <- results.2 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.2, by = c("parameter", "iteration"))

HPD.600.2.4viz <- results.2 %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.2, by = c("parameter", "iteration"))

HPD.1000.2.4viz <- results.2 %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.2, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.2.summary <- HPD.200.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.2.summary <- HPD.600.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.2.summary <- HPD.1000.2.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
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

write_csv(HPD.2.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose2, ".csv"))


#--------------------------HS.PO_surv.prior-----------------------------------#
HPD.200.3.4viz <- results.3 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.3, by = c("parameter", "iteration"))

HPD.600.3.4viz <- results.3 %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.3, by = c("parameter", "iteration"))

HPD.1000.3.4viz <- results.3 %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.3, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.3.summary <- HPD.200.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.3.summary <- HPD.600.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.3.summary <- HPD.1000.3.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
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

write_csv(HPD.3.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose3, ".csv"))


# #----------------------------HS.only_one.indv---------------------------------#
# HPD.200.4.4viz <- results.4 %>% filter(total_samples == 200) %>% 
#   right_join(HPD.200.4, by = c("parameter", "iteration"))
# 
# HPD.600.4.4viz <- results.4 %>% filter(total_samples == 600) %>% 
#   right_join(HPD.600.4, by = c("parameter", "iteration"))
# 
# HPD.1000.4.4viz <- results.4 %>% filter(total_samples == 1000) %>% 
#   right_join(HPD.1000.4, by = c("parameter", "iteration"))
# 
# 
# #How many iterations fall within the expected HPD interval
# HPD.200.4.summary <- HPD.200.4.4viz %>% group_by(parameter, interval) %>% 
#   mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
#   dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)
# 
# HPD.600.4.summary <- HPD.600.4.4viz %>% group_by(parameter, interval) %>% 
#   mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
#   dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)
# 
# HPD1000.4.summary <- HPD.1000.4.4viz %>% group_by(parameter, interval) %>% 
#   mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
#   dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)
# 
# #Create combined dataframe from above
# HPD.4.summary_all <- HPD.200.4.summary %>% inner_join(HPD.600.4.summary, by = c("parameter", "interval")) %>% 
#   inner_join(HPD1000.4.summary, by = c("parameter", "interval")) %>% 
#   mutate(interval.scale = interval*100)
# 
# HPD.4.summary.tidy <- HPD.4.summary_all %>% 
#   pivot_longer(
#     cols = starts_with("percent"),
#     names_to = "samples",
#     names_prefix = "percent_in_interval.",
#     values_to = "percent_in_interval"
#   ) %>% 
#   mutate(model_type = "HS.only",
#          litter.samples = "one")
# 
# write_csv(HPD.4.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose4, ".csv"))



#---------------create dataframes for plotting by sample size--------
HPD.200.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "200.samples")
HPD.600.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.1.4viz <- HPD.1.summary.tidy %>% filter(samples == "1000.samples")

HPD.200.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "200.samples")
HPD.600.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.2.4viz <- HPD.2.summary.tidy %>% filter(samples == "1000.samples")

HPD.200.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "200.samples")
HPD.600.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "600.samples")
HPD.1000.3.4viz <- HPD.3.summary.tidy %>% filter(samples == "1000.samples")

# HPD.200.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "200.samples")
# HPD.600.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "600.samples")
# HPD.1000.4.4viz <- HPD.4.summary.tidy %>% filter(samples == "1000.samples")



all.200.4viz <- HPD.200.1.4viz %>% bind_rows(HPD.200.2.4viz, HPD.200.3.4viz)

all.600.4viz <- HPD.600.1.4viz %>% bind_rows(HPD.600.2.4viz, HPD.600.3.4viz)

all.1000.4viz <- HPD.1000.1.4viz %>% bind_rows(HPD.1000.2.4viz, HPD.1000.3.4viz)

#-----------------------Viz----------------------------------------------------
#---------------------line plot -----------------------------------------------
#---------------------200 samples----------------------------------------------#
Nf.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "Nf"), 
                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "Nm"), 
       aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "surv"), 
                         aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "lam"), 
                        aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.200 <- ggarrange(Nf.200samps.HPDI.plot, Nm.200samps.HPDI.plot, surv.200samps.HPDI.plot, lam.200samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.200, fig.lab = "200 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")


#---------------------300 samples----------------------------------------------#
Nf.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "Nf"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "Nm"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "surv"), 
                                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "lam"), 
                                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.600 <- ggarrange(Nf.600samps.HPDI.plot, Nm.600samps.HPDI.plot, surv.600samps.HPDI.plot, lam.600samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.600, fig.lab = "600 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

#---------------------1000 samples----------------------------------------------#
Nf.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "Nf"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "Nm"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "surv"), 
                                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "lam"), 
                                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = purpose, fill = purpose, shape = purpose), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = purpose, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.1000 <- ggarrange(Nf.1000samps.HPDI.plot, Nm.1000samps.HPDI.plot, surv.1000samps.HPDI.plot, lam.1000samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.1000, fig.lab = "1000 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")



head(results.all)

#Save to pdf
HPDI.file <- paste0(results_plots_location, "HPDI_model_comparison_", today, ".pdf")
pdf(file = HPDI.file, width = 11, height = 8)

#200 samples
annotate_figure(plots.200, fig.lab = "200 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")
#300 samples
annotate_figure(plots.600, fig.lab = "600 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")
#400 samples
annotate_figure(plots.1000, fig.lab = "1000 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

dev.off()









#----------------Investigate results more closely---------------
#How many POPs detected per sample size?
results.all %>% dplyr::filter(model_type == "HS.PO") %>% 
  group_by(purpose, total_samples) %>% 
  summarize(mean(POPs_detected))

#How many HSPs detected per sample size?
results.all %>%
  group_by(purpose, total_samples) %>% 
  summarize(mean(HSPs_detected))

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









#----------------------View prior and posterior----------------------------
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








#-------------------More elaborate/specific line plots-------------------------------
#Specify save location for pdf of plots
coverage.file.400 <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.diagnostics/Plots/Statistical.coverage_", today, "_400Samps.pdf")

plotlist <- list()

#Create one page of coverage plots for each iteration
for(it in 1:iterations){
  HPD.temp <- HPD.400.4viz %>% filter(iteration == it)

gg.Fabundance <- ggplot(data = HPD.temp %>% filter(parameter == "Nf"), 
                        aes(x = interval, y = truth)) + 
  geom_linerange(aes(ymin = lower, ymax = upper, colour = parameter), size = 2) + 
  geom_point(aes(fill = parameter), size = 2) + 
  scale_color_manual(labels = "HPD interval", values = "maroon") + 
  scale_fill_manual(labels = "Truth", values = "maroon") +
  ggtitle("Female abundance") + 
  scale_x_continuous(breaks = intervals) +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(vjust = -1), 
        axis.text.x = element_text(angle = 45, vjust = .5))

gg.Mabundance <- ggplot(data = HPD.temp %>% filter(parameter == "Nm"), 
                        aes(x = interval, y = truth)) + 
  geom_linerange(aes(ymin = lower, ymax = upper, colour = parameter), size = 2) + 
  geom_point(aes(fill = factor(parameter)), size = 2) + 
  scale_color_manual(labels = "HPD interval", values = "mediumslateblue") + 
  scale_fill_manual(labels = "Truth", values = "mediumslateblue") +
  ggtitle("Male abundance") + 
  scale_x_continuous(breaks = intervals) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.title.x = element_text(vjust = -1), 
        axis.text.x = element_text(angle = 45, vjust = .5))

gg.survival <- ggplot(data = HPD.temp %>% filter(parameter == "surv"), 
                      aes(x = interval, y = truth)) + 
  geom_linerange(aes(ymin = lower, ymax = upper, colour = parameter), size = 2) + 
  geom_point(aes(fill = factor(parameter), shape = "Truth"), size = 2) + 
  scale_color_brewer(palette="PRGn") + 
  ggtitle("Adult survival") + 
  scale_x_continuous(breaks = intervals) +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(vjust = -1), 
        axis.text.x = element_text(angle = 45, vjust = .5))

plots <- ggarrange(gg.Fabundance, gg.Mabundance, gg.survival, common.legend = TRUE, legend = "right", vjust = -1)
plotlist[[it]] <- annotate_figure(plots, top = text_grob(paste0("Statistical coverage for iteration ", it), size = 18, face = "bold"))
}

ggexport(plotlist, filename = coverage.file.400)
