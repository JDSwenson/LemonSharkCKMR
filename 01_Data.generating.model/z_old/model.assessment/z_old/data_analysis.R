library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)
library(ggmcmc)


rm(list=ls())


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
ni <- 30000 # number of post-burn-in samples per chain
nb <- 40000 # number of burn-in samples
nt <- 15     # thinning rate
nc <- 2      # number of chains
MCMC.settings <- paste0("thin", nt, "_draw", ni, "_burn", nb)


date.of.simulation <- "28Mar2022"
purpose1 <- "Test_HS.PO"
purpose2 <- "Test_HS.only"
purpose3 <- "HS.PO_one.indv.per.parent"
purpose4 <- "HS.only_one.indv.per.parent"

seeds <- "Seeds2022.03.23"
sim.samples.1 <- "200.samples"
sim.samples.2 <- "600.samples"
sim.samples.3 <- "1000.samples"

MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.results/"
results_plots_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.construction/Model.results/figures/"
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
parents_prefix <- "parents_breakdown/CKMR_parents.breakdown"
sample.prefix <- "sample_info/CKMR_sample.info"


#Results
results.Test_HS.PO <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose1, ".csv")) %>% 
  mutate(model_type = "HS.PO", 
         litter.samples = "all")

results.Test_HS.only <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose2, ".csv")) %>% 
  mutate(model_type = "HS.only",
         litter.samples = "all")

results.HS.PO_one.indv <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose3, ".csv")) %>% 
  mutate(model_type = "HS.PO",
         litter.samples = "one")

results.HS.only_one.indv <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose4, ".csv")) %>% 
  mutate(model_type = "HS.only",
         litter.samples = "one")


#MCMC samples/output
#Trial 1
s1.Test_HS.PO <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose1))

s2.Test_HS.PO <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose1))

s3.Test_HS.PO <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose1))

#Trial 2
s1.Test_HS.only <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose2))

s2.Test_HS.only <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose2))

s3.Test_HS.only <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose2))

#Trial 3
s1.HS.PO_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose3))

s2.HS.PO_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose3))

s3.HS.PO_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose3))

#Trial 4
s1.HS.only_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose4))

s2.HS.only_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose4))

s3.HS.only_one.indv <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose4))

#Examine results quickly
results.all <- results.Test_HS.PO %>% bind_rows(results.Test_HS.only, results.HS.PO_one.indv, results.HS.only_one.indv) %>% 
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
 results.all %>% group_by(model_type, total_samples, parameter) %>% 
   dplyr::summarize(median = median(relative_bias), n = n()) %>% 
  mutate(abs.bias = abs(median)) %>% 
  arrange(median)
 
 #-----------Within HPDI?-------------------------#
 #Within HPD interval?
results.all %>% group_by(model_type, litter.samples, total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100) %>% 
  filter(parameter != "lam") %>% 
  arrange(desc(percent_in_interval))

#---------------CV------------------------#
results.all %>% group_by(model_type, litter.samples, total_samples, parameter) %>% 
  dplyr::summarize(mean.cv = mean(cv)) %>%
  dplyr::filter(parameter == "Nf" | parameter == "Nm") %>% 
  arrange(mean.cv) %>% 
  View()


#Violin plot
RB_violin <- ggplot(data=results, aes(x=factor(parameter), fill = factor(parameter))) +
  geom_violin(aes(y=relative_bias), draw_quantiles = 0.5) +
  #ylim(-50, 160) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Iteration", y="Relative bias", title="Relative bias by iteration") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")

RB_violin

#---------------------------Box plots----------------------------------
#Specify save location for pdf of plots
boxPlots.file <- paste0(results_plots_location, "BoxPlots_Priorcomparison.pdf")
pdf(file = boxPlots.file, width = 14, height = 4)
sim.samples.all <- c(200, 300, 400)
plot.list <- NULL

for(i in 1:length(sim.samples.all)){
  samps <- sim.samples.all[i]
  
  #boxplot plot
  RB_boxplot.norm <- results2.norm %>% filter(total_samples == samps) %>% 
    ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
    geom_boxplot(aes(y=relative_bias)) +
    #ylim(-50, 160) +
    geom_hline(yintercept=0, col="black", size=1.25) +
    #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
    labs(x="Iteration", y="Relative bias", title=paste0("Normal prior (uninformative) ", samps, " population samples")) +
    scale_fill_brewer(palette="Set2") +
    font("title", size = 10, face = "bold") + 
    theme(legend.title = element_blank()) + 
    ylim(-50, 150)
  
  RB_boxplot.uniform <- results.uniform %>% filter(total_samples == samps) %>% 
  ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
  geom_boxplot(aes(y=relative_bias)) +
  #ylim(-50, 160) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Iteration", y=element_blank(), title=paste0("Uniform prior (uninformative) ",  samps, " population samples")) +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold") + 
  ylim(-50, 150)

  
RB_boxplot.hier <- results.hier %>% filter(total_samples == samps) %>%
  ggplot(aes(x=factor(parameter), fill = factor(parameter))) +
  geom_boxplot(aes(y=relative_bias)) +
  #ylim(-50, 160) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  #annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Iteration", y = element_blank(), title=paste0("Hierarchical model ", samps, " population samples")) +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold") +
  ylim(-50, 150)

# plot.list[[1]] <- RB_boxplot.norm
# plot.list[[j+1]] <- RB_boxplot.uniform
# plot.list[[j+2]] <- RB_boxplot.hier
# print(ggarrange(plotlist = plot.list, common.legend = TRUE, nrow = 1, ncol = 3, widths = c(1, 1, 1))) #Need to wrap

print(ggarrange(RB_boxplot.norm, RB_boxplot.uniform, RB_boxplot.hier, nrow = 1, ncol = 3, common.legend = TRUE))
  }

dev.off()

#-------------------Calculate HPD intervals---------------------------
# How many estimates fall in different HPD intervals?
#Calculate HPD interval for each iteration
intervals <- c(seq(from = 0.95, to = 0.05, by = -0.05))
iterations = 100

#200 samples
HPD.200.Test_HS.PO <- NULL
HPD.600.Test_HS.PO <- NULL
HPD.1000.Test_HS.PO <- NULL
HPD.200.Test_HS.only <- NULL
HPD.600.Test_HS.only <- NULL
HPD.1000.Test_HS.only <- NULL
HPD.200.HS.PO_one.indv <- NULL
HPD.600.HS.PO_one.indv <- NULL
HPD.1000.HS.PO_one.indv <- NULL
HPD.200.HS.only_one.indv <- NULL
HPD.600.HS.only_one.indv <- NULL
HPD.1000.HS.only_one.indv <- NULL

for(i in 1:length(s1.Test_HS.PO)){
  for(j in 1:length(intervals)){
    #Test_HS.PO
    #200 samples
    post.HPD.200.Test_HS.PO <- combine.mcmc(s1.Test_HS.PO[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.Test_HS.PO <- post.HPD.200.Test_HS.PO[row.names(post.HPD.200.Test_HS.PO) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200.Test_HS.PO <- rbind(HPD.200.Test_HS.PO, post.HPD.200.Test_HS.PO)
    
    #Test_HS.PO
    #600 samples
    post.HPD.600.Test_HS.PO <- combine.mcmc(s2.Test_HS.PO[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.Test_HS.PO <- post.HPD.600.Test_HS.PO[row.names(post.HPD.600.Test_HS.PO) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.Test_HS.PO <- rbind(HPD.600.Test_HS.PO, post.HPD.600.Test_HS.PO)

    #Test_HS.PO
    #1000 samples
    post.HPD.1000.Test_HS.PO <- combine.mcmc(s3.Test_HS.PO[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.Test_HS.PO <- post.HPD.1000.Test_HS.PO[row.names(post.HPD.1000.Test_HS.PO) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.Test_HS.PO <- rbind(HPD.1000.Test_HS.PO, post.HPD.1000.Test_HS.PO)
    
    #Test_HS.only
    #200 samples
    post.HPD.200.Test_HS.only <- combine.mcmc(s1.Test_HS.only[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.Test_HS.only <- post.HPD.200.Test_HS.only[row.names(post.HPD.200.Test_HS.only) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.Test_HS.only <- rbind(HPD.200.Test_HS.only, post.HPD.200.Test_HS.only)
    
    #Test_HS.only
    #600 samples
    post.HPD.600.Test_HS.only <- combine.mcmc(s2.Test_HS.only[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.Test_HS.only <- post.HPD.600.Test_HS.only[row.names(post.HPD.600.Test_HS.only) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.Test_HS.only <- rbind(HPD.600.Test_HS.only, post.HPD.600.Test_HS.only)
    
    #Test_HS.only
    #1000 samples
    post.HPD.1000.Test_HS.only <- combine.mcmc(s3.Test_HS.only[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.Test_HS.only <- post.HPD.1000.Test_HS.only[row.names(post.HPD.1000.Test_HS.only) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.Test_HS.only <- rbind(HPD.1000.Test_HS.only, post.HPD.1000.Test_HS.only)

    
    #HS.PO_one.indv
    #200 samples
    post.HPD.200.HS.PO_one.indv <- combine.mcmc(s1.HS.PO_one.indv[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.HS.PO_one.indv <- post.HPD.200.HS.PO_one.indv[row.names(post.HPD.200.HS.PO_one.indv) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.HS.PO_one.indv <- rbind(HPD.200.HS.PO_one.indv, post.HPD.200.HS.PO_one.indv)
    
    #HS.PO_one.indv
    #600 samples
    post.HPD.600.HS.PO_one.indv <- combine.mcmc(s2.HS.PO_one.indv[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.HS.PO_one.indv <- post.HPD.600.HS.PO_one.indv[row.names(post.HPD.600.HS.PO_one.indv) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.HS.PO_one.indv <- rbind(HPD.600.HS.PO_one.indv, post.HPD.600.HS.PO_one.indv)
    
    #HS.PO_one.indv
    #1000 samples
    post.HPD.1000.HS.PO_one.indv <- combine.mcmc(s3.HS.PO_one.indv[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.HS.PO_one.indv <- post.HPD.1000.HS.PO_one.indv[row.names(post.HPD.1000.HS.PO_one.indv) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.HS.PO_one.indv <- rbind(HPD.1000.HS.PO_one.indv, post.HPD.1000.HS.PO_one.indv)
    
    
    #HS.only_one.indv
    #200 samples
    post.HPD.200.HS.only_one.indv <- combine.mcmc(s1.HS.only_one.indv[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.HS.only_one.indv <- post.HPD.200.HS.only_one.indv[row.names(post.HPD.200.HS.only_one.indv) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.HS.only_one.indv <- rbind(HPD.200.HS.only_one.indv, post.HPD.200.HS.only_one.indv)
    
    #HS.only_one.indv
    #600 samples
    post.HPD.600.HS.only_one.indv <- combine.mcmc(s2.HS.only_one.indv[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.600.HS.only_one.indv <- post.HPD.600.HS.only_one.indv[row.names(post.HPD.600.HS.only_one.indv) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.600.HS.only_one.indv <- rbind(HPD.600.HS.only_one.indv, post.HPD.600.HS.only_one.indv)
    
    #HS.only_one.indv
    #1000 samples
    post.HPD.1000.HS.only_one.indv <- combine.mcmc(s3.HS.only_one.indv[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.1000.HS.only_one.indv <- post.HPD.1000.HS.only_one.indv[row.names(post.HPD.1000.HS.only_one.indv) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.1000.HS.only_one.indv <- rbind(HPD.1000.HS.only_one.indv, post.HPD.1000.HS.only_one.indv)
    
  }
}


#--------------Create dataframes for viz and analysis----------------------

#--------------------------Test_HS.PO-----------------------------------#
HPD.200.Test_HS.PO.4viz <- results.Test_HS.PO %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.Test_HS.PO, by = c("parameter", "iteration"))

HPD.600.Test_HS.PO.4viz <- results.Test_HS.PO %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.Test_HS.PO, by = c("parameter", "iteration"))

HPD.1000.Test_HS.PO.4viz <- results.Test_HS.PO %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.Test_HS.PO, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.Test_HS.PO.summary <- HPD.200.Test_HS.PO.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.Test_HS.PO.summary <- HPD.600.Test_HS.PO.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.Test_HS.PO.summary <- HPD.1000.Test_HS.PO.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
HPD.Test_HS.PO.summary_all <- HPD.200.Test_HS.PO.summary %>% inner_join(HPD.600.Test_HS.PO.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.Test_HS.PO.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.Test_HS.PO.summary.tidy <- HPD.Test_HS.PO.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
    ) %>% 
  mutate(model_type = "HS.PO",
         litter.samples = "all")

write_csv(HPD.Test_HS.PO.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose1, ".csv"))


#--------------------------Test_HS.only-----------------------------------#
HPD.200.Test_HS.only.4viz <- results.Test_HS.only %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.Test_HS.only, by = c("parameter", "iteration"))

HPD.600.Test_HS.only.4viz <- results.Test_HS.only %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.Test_HS.only, by = c("parameter", "iteration"))

HPD.1000.Test_HS.only.4viz <- results.Test_HS.only %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.Test_HS.only, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.Test_HS.only.summary <- HPD.200.Test_HS.only.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.Test_HS.only.summary <- HPD.600.Test_HS.only.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.Test_HS.only.summary <- HPD.1000.Test_HS.only.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
HPD.Test_HS.only.summary_all <- HPD.200.Test_HS.only.summary %>% inner_join(HPD.600.Test_HS.only.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.Test_HS.only.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.Test_HS.only.summary.tidy <- HPD.Test_HS.only.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = "HS.only",
                  litter.samples = "all")

write_csv(HPD.Test_HS.only.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose2, ".csv"))


#--------------------------HS.PO_one.indv-----------------------------------#
HPD.200.HS.PO_one.indv.4viz <- results.HS.PO_one.indv %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.HS.PO_one.indv, by = c("parameter", "iteration"))

HPD.600.HS.PO_one.indv.4viz <- results.HS.PO_one.indv %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.HS.PO_one.indv, by = c("parameter", "iteration"))

HPD.1000.HS.PO_one.indv.4viz <- results.HS.PO_one.indv %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.HS.PO_one.indv, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.HS.PO_one.indv.summary <- HPD.200.HS.PO_one.indv.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.HS.PO_one.indv.summary <- HPD.600.HS.PO_one.indv.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.HS.PO_one.indv.summary <- HPD.1000.HS.PO_one.indv.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
HPD.HS.PO_one.indv.summary_all <- HPD.200.HS.PO_one.indv.summary %>% inner_join(HPD.600.HS.PO_one.indv.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.HS.PO_one.indv.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.HS.PO_one.indv.summary.tidy <- HPD.HS.PO_one.indv.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = "HS.PO",
         litter.samples = "one")

write_csv(HPD.HS.PO_one.indv.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose3, ".csv"))


#----------------------------HS.only_one.indv---------------------------------#
HPD.200.HS.only_one.indv.4viz <- results.HS.only_one.indv %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.HS.only_one.indv, by = c("parameter", "iteration"))

HPD.600.HS.only_one.indv.4viz <- results.HS.only_one.indv %>% filter(total_samples == 600) %>% 
  right_join(HPD.600.HS.only_one.indv, by = c("parameter", "iteration"))

HPD.1000.HS.only_one.indv.4viz <- results.HS.only_one.indv %>% filter(total_samples == 1000) %>% 
  right_join(HPD.1000.HS.only_one.indv, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.HS.only_one.indv.summary <- HPD.200.HS.only_one.indv.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.600.HS.only_one.indv.summary <- HPD.600.HS.only_one.indv.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.600.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD1000.HS.only_one.indv.summary <- HPD.1000.HS.only_one.indv.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.1000.samples = sum(in_HPD_interval == "Y")/n() * 100)

#Create combined dataframe from above
HPD.HS.only_one.indv.summary_all <- HPD.200.HS.only_one.indv.summary %>% inner_join(HPD.600.HS.only_one.indv.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD1000.HS.only_one.indv.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.HS.only_one.indv.summary.tidy <- HPD.HS.only_one.indv.summary_all %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(model_type = "HS.only",
         litter.samples = "one")

write_csv(HPD.HS.only_one.indv.summary_all, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose4, ".csv"))



#---------------create dataframes for plotting by sample size--------
Test_HS.PO.200.4viz <- HPD.Test_HS.PO.summary.tidy %>% filter(samples == "200.samples")
Test_HS.PO.600.4viz <- HPD.Test_HS.PO.summary.tidy %>% filter(samples == "600.samples")
Test_HS.PO.1000.4viz <- HPD.Test_HS.PO.summary.tidy %>% filter(samples == "1000.samples")

Test_HS.only.200.4viz <- HPD.Test_HS.only.summary.tidy %>% filter(samples == "200.samples")
Test_HS.only.600.4viz <- HPD.Test_HS.only.summary.tidy %>% filter(samples == "600.samples")
Test_HS.only.1000.4viz <- HPD.Test_HS.only.summary.tidy %>% filter(samples == "1000.samples")

HS.PO_one.indv.200.4viz <- HPD.HS.PO_one.indv.summary.tidy %>% filter(samples == "200.samples")
HS.PO_one.indv.600.4viz <- HPD.HS.PO_one.indv.summary.tidy %>% filter(samples == "600.samples")
HS.PO_one.indv.1000.4viz <- HPD.HS.PO_one.indv.summary.tidy %>% filter(samples == "1000.samples")

HS.only_one.indv.200.4viz <- HPD.HS.only_one.indv.summary.tidy %>% filter(samples == "200.samples")
HS.only_one.indv.600.4viz <- HPD.HS.only_one.indv.summary.tidy %>% filter(samples == "600.samples")
HS.only_one.indv.1000.4viz <- HPD.HS.only_one.indv.summary.tidy %>% filter(samples == "1000.samples")



all.200.4viz <- Test_HS.PO.200.4viz %>% bind_rows(Test_HS.only.200.4viz, HS.PO_one.indv.200.4viz, HS.only_one.indv.200.4viz) %>% 
  mutate(model.label = paste0(model_type, "_", litter.samples))

all.600.4viz <- Test_HS.PO.600.4viz %>% bind_rows(Test_HS.only.600.4viz, HS.PO_one.indv.600.4viz, HS.only_one.indv.600.4viz) %>% 
  mutate(model.label = paste0(model_type, "_", litter.samples))

all.1000.4viz <- Test_HS.PO.1000.4viz %>% bind_rows(Test_HS.only.1000.4viz, HS.PO_one.indv.1000.4viz, HS.only_one.indv.1000.4viz) %>% 
  mutate(model.label = paste0(model_type, "_", litter.samples))

#-----------------------Viz----------------------------------------------------
#---------------------line plot -----------------------------------------------
#---------------------200 samples----------------------------------------------#
Nf.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "Nf"), 
                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "Nm"), 
       aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "surv"), 
                         aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "lam"), 
                        aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
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
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "Nm"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "surv"), 
                                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.600samps.HPDI.plot <- ggplot(data = all.600.4viz %>% filter(parameter == "lam"), 
                                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
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
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "Nm"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "surv"), 
                                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.1000samps.HPDI.plot <- ggplot(data = all.1000.4viz %>% filter(parameter == "lam"), 
                                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = model.label, fill = model.label, shape = model.label), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = model.label, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.1000 <- ggarrange(Nf.1000samps.HPDI.plot, Nm.1000samps.HPDI.plot, surv.1000samps.HPDI.plot, lam.1000samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.1000, fig.lab = "1000 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")



head(results)

#Save to pdf
HPDI.file <- paste0(results_plots_location, "HPDI_model_comparison.pdf")
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
  group_by(litter.samples, total_samples) %>% 
  summarize(mean(POPs_detected))

#How many HSPs detected per sample size?
results.all %>%
  group_by(model_type, litter.samples, total_samples) %>% 
  summarize(mean(HSPs_detected))

#----------Bring in info for parents and samples----------------#
# Breakdown of offspring for each parent
rents.Test_HS.PO <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose1)) %>% 
  mutate(model_type = "HS.PO", 
         litter.samples = "all")

rents.Test_HS.only <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose2)) %>% 
  mutate(model_type = "HS.only", 
         litter.samples = "all")

rents.HS.PO_one.indv <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose3)) %>% 
  mutate(model_type = "HS.PO", 
         litter.samples = "one")

rents.HS.only_one.indv <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose4)) %>% 
  mutate(model_type = "HS.only", 
         litter.samples = "one")

#Breakdown of samples drawn from simulation
samples1 <- sample.info.Test_HS.PO <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose1)) %>% 
  mutate(model_type = "HS.PO", 
         litter.samples = "all")

samples2 <- sample.info.Test_HS.only <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose2)) %>% 
  mutate(model_type = "HS.only", 
         litter.samples = "all")

samples3 <- sample.info.HS.PO_one.indv <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose3)) %>% 
  mutate(model_type = "HS.PO", 
         litter.samples = "one")

samples4 <- sample.info.HS.only_one.indv <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose4)) %>% 
  mutate(model_type = "HS.only", 
         litter.samples = "one")

head(samples1)

source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO.R")


#Make dataframes of pairwise comparisons for each iteration and sample size for further investigation into gap years, pop growth years, and number of detected kin.

mom_comps.samples1 <- dad_comps.samples1 <- NULL
for(iter in 1:100){
  print(paste0("working on iteration", iter))

    #50 samples per year
  samples1.temp1 <- samples1 %>% filter(sample.size == 50 & iteration == iter)
  filter1.temp1 <- samples1.temp1 %>% filter.samples()
  PO.samples1.list1 <- filter1.temp1[[1]]
  HS.samples1.df1 <- filter1.temp1[[2]]
  pairwise.samples1 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list1, filtered.samples.HS.df = HS.samples1.df1)
  
  mom_comps.temp1 <- pairwise.samples1[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 50)
  dad_comps.temp1 <- pairwise.samples1[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 50)

  #150 samples per year
  samples1.temp2 <- samples1 %>% filter(sample.size == 150 & iteration == iter)
  filter1.temp2 <- samples1.temp2 %>% filter.samples()
  PO.samples1.list2 <- filter1.temp2[[1]]
  HS.samples1.df2 <- filter1.temp2[[2]]
  pairwise.samples2 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list2, filtered.samples.HS.df = HS.samples1.df2)
  
  mom_comps.temp2 <- pairwise.samples2[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 150)
  dad_comps.temp2 <- pairwise.samples2[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 150)
  
  #250 samples per year
  samples1.temp3 <- samples1 %>% filter(sample.size == 250 & iteration == iter)
  filter1.temp3 <- samples1.temp3 %>% filter.samples()
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
