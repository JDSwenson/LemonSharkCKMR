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
results <- read_csv("G://My Drive/Personal_Drive/R/CKMR/Model.results/Model.validation/CKMR_results_08Dec2021_longChain.csv")

head(results)

results <- results %>% mutate(relative_bias2 = round(((`50` - truth)/truth)*100,1))

#Median relative bias by sample size
results %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(median = median(relative_bias2), n = n())
####------------- MCMC parameters ----------------####
ni <- 30000 # number of post-burn-in samples per chain
nb <- 40000 # number of burn-in samples
nt <- 15     # thinning rate
nc <- 2      # number of chains
MCMC.settings <- paste0("thin", nt, "_draw", ni, "_burn", nb)


date.of.simulation <- "01Jan2022"
purpose <- "estSurvLam"
seeds <- "Seeds12.27"
sim.samples.1 <- "200.samples"
sim.samples.2 <- "300.samples"
sim.samples.3 <- "400.samples"

MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.results/"
results_plots_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.results/figures/"
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
parents_prefix <- "parents_breakdown/CKMR_parents.breakdown"
sample.prefix <- "sample_info/CKMR_sample.info"


#Results - Normal prior (uninformative)
results.norm <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, ".csv")) %>% 
  mutate(N.prior_max = NA, prior = "normal")

#MCMC samples/output
s1.norm <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))

s2.norm <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))

s3.norm <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))

# Breakdown of offspring for each parent
rents.norm <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

#Breakdown of samples drawn from simulation
sample.info.norm <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))



#Results - uniform prior (uninformative)
date.of.simulation <- "03Jan2022"
purpose <- "testUniformPriors"

results.uniform <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, ".csv")) %>% 
  mutate(prior = "uniform")

#MCMC samples/output
s1.uniform <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))

s2.uniform <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))

s3.uniform <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))

# Breakdown of offspring for each parent
rents.uniform <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

#Breakdown of samples drawn from simulation
sample.info.uniform <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))



#Results - hierarchical approach
date.of.simulation <- "09Jan2022"
purpose <- "testHierarchical2"

results.hier <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, ".csv")) %>% 
  mutate(N.prior_max = NA, prior = "hier")

#MCMC samples/output
s1.hier <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))

s2.hier <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))

s3.hier <- readRDS(paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))

# Breakdown of offspring for each parent
rents.hier <- readRDS(paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

#Breakdown of samples drawn from simulation
sample.info.hier <- readRDS(paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))


#Examine results quickly
results.all <- results.norm %>% bind_rows(results.uniform, results.hier)

results.all %>% group_by(prior, total_samples, parameter, in_interval) %>% 
  summarize(n()) %>% 
  View()


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


#----------------------Quick analysis of results----------------------------
#today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
jags_params <- c("Nf", "Nm", "surv", "lam") #Specify parameters

head(results)

#Temporary fix since results saved funny ... 
#Lambda and survival were switched, so switching them back ... 
HPD2.5_surv.norm <- results.norm %>% filter(parameter == "lam") %>% 
  pull(HPD2.5)

HPD2.5_lam.norm <- results.norm %>% filter(parameter == "surv") %>% 
  pull(HPD2.5)  


HPD97.5_surv.norm <- results.norm %>% filter(parameter == "lam") %>% 
  pull(HPD97.5)

HPD97.5_lam.norm <- results.norm %>% filter(parameter == "surv") %>% 
  pull(HPD97.5)  

results.norm[results.norm$parameter == "lam",]$HPD2.5 <- HPD2.5_lam.norm
results.norm[results.norm$parameter == "surv",]$HPD2.5 <- HPD2.5_surv.norm

results.norm[results.norm$parameter == "lam",]$HPD97.5 <- HPD97.5_lam.norm
results.norm[results.norm$parameter == "surv",]$HPD97.5 <- HPD97.5_surv.norm

results2.norm <- results.norm %>% dplyr::select(-c(in_interval, relative_bias))

#Calculate relative bias
 results2.norm <- results2.norm %>% 
   mutate(relative_bias = round(((Q50 - truth)/truth)*100,1)) %>%
   mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N")) 
 #%>% 
  # mutate(percent_sampled = round((total_samples/pop_size_mean) * 100, 0)

 
 #-----------Median Relative bias by sample size-------------------------
#Median relative bias by sample size - normal prior
 results2.norm %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(median = median(relative_bias), n = n())
 
#Median relative bias by sample size - uniform prior
 results.uniform %>% group_by(total_samples, parameter) %>% 
   filter(N.prior_max == 1000) %>% 
   dplyr::summarize(median = median(relative_bias), n = n())
 
 #Median relative bias by sample size - hierarchical
 results.hier %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(median = median(relative_bias), n = n())

 #-----------Within HPDI?-------------------------
 #Within HPD interval?
results2.norm %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

 results.uniform %>% group_by(total_samples, parameter) %>%
   filter(N.prior_max == 1000) %>% 
   dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)
 
 results.hier %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)


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
boxPlots.file <- paste0(results_plots_location, "BoxPlots_comparison.pdf")
pdf(file = boxPlots.file, width = 14, height = 4)
sim.samples.all <- c(200, 300, 400)
prior.max <- c(1000, 2000, 3000)
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
  
  RB_boxplot.uniform <- results.uniform %>% filter(N.prior_max == 1000 & total_samples == samps) %>% 
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
HPD.200.norm <- NULL
HPD.300.norm <- NULL
HPD.400.norm <- NULL
HPD.200.uniform <- NULL
HPD.300.uniform <- NULL
HPD.400.uniform <- NULL
HPD.200.hier <- NULL
HPD.300.hier <- NULL
HPD.400.hier <- NULL

for(i in 1:length(s1.norm)){
  for(j in 1:length(intervals)){
    #Normal prior
    #200 samples
    post.HPD.200.norm <- combine.mcmc(s1.norm[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.norm <- post.HPD.200.norm[row.names(post.HPD.200.norm) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200.norm <- rbind(HPD.200.norm, post.HPD.200.norm)
    
    #Normal prior
    #300 samples
    post.HPD.300.norm <- combine.mcmc(s2.norm[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.300.norm <- post.HPD.300.norm[row.names(post.HPD.300.norm) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.300.norm <- rbind(HPD.300.norm, post.HPD.300.norm)

    #Normal prior
    #400 samples
    post.HPD.400.norm <- combine.mcmc(s3.norm[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.400.norm <- post.HPD.400.norm[row.names(post.HPD.400.norm) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.400.norm <- rbind(HPD.400.norm, post.HPD.400.norm)
    

    #Hierarchical prior
    #200 samples
    post.HPD.200.hier <- combine.mcmc(s1.hier[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.200.hier <- post.HPD.200.hier[row.names(post.HPD.200.hier) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.200.hier <- rbind(HPD.200.hier, post.HPD.200.hier)
    
    #Hierarchical prior
    #300 samples
    post.HPD.300.hier <- combine.mcmc(s2.hier[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.300.hier <- post.HPD.300.hier[row.names(post.HPD.300.hier) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.300.hier <- rbind(HPD.300.hier, post.HPD.300.hier)
    
    #Hierarchical prior
    #400 samples
    post.HPD.400.hier <- combine.mcmc(s3.hier[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD.400.hier <- post.HPD.400.hier[row.names(post.HPD.400.hier) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")
    
    HPD.400.hier <- rbind(HPD.400.hier, post.HPD.400.hier)
    
  }
}


for(i in 1:length(s1.uniform)){
  for(j in 1:length(intervals)){
#Uniform
#Uniform prior
#200 samples
post.HPD.200.uniform <- combine.mcmc(s1.uniform[[i]]) %>% 
  HPDinterval(prob = intervals[j]) %>%
  data.frame() %>% 
  mutate(interval = intervals[j], iteration = i)

post.HPD.200.uniform <- post.HPD.200.uniform[row.names(post.HPD.200.uniform) %in% jags_params,] %>% #Remove deviance
  rownames_to_column(var = "parameter")

HPD.200.uniform <- rbind(HPD.200.uniform, post.HPD.200.uniform)

#Uniform prior
#300 samples
post.HPD.300.uniform <- combine.mcmc(s2.uniform[[i]]) %>% 
  HPDinterval(prob = intervals[j]) %>%
  data.frame() %>% 
  mutate(interval = intervals[j], iteration = i)

post.HPD.300.uniform <- post.HPD.300.uniform[row.names(post.HPD.300.uniform) %in% jags_params,] %>% #Remove deviance
  rownames_to_column(var = "parameter")

HPD.300.uniform <- rbind(HPD.300.uniform, post.HPD.300.uniform)

#Uniform prior
#400 samples
post.HPD.400.uniform <- combine.mcmc(s3.uniform[[i]]) %>% 
  HPDinterval(prob = intervals[j]) %>%
  data.frame() %>% 
  mutate(interval = intervals[j], iteration = i)

post.HPD.400.uniform <- post.HPD.400.uniform[row.names(post.HPD.400.uniform) %in% jags_params,] %>% #Remove deviance
  rownames_to_column(var = "parameter")

HPD.400.uniform <- rbind(HPD.400.uniform, post.HPD.400.uniform)
  }
}

# tail(HPD.200)
# levels(factor(HPD.200$interval))
# 
# #300 samples
# 
# for(i in 1:length(sim.samples.2_MCMC)){
#   for(j in 1:length(intervals)){
#     post.HPD <- combine.mcmc(sim.samples.2_MCMC[[i]]) %>% 
#       HPDinterval(prob = intervals[j]) %>%
#       data.frame() %>% 
#       mutate(interval = intervals[j], iteration = i)
#     
#     post.HPD <- post.HPD[row.names(post.HPD) %in% jags_params,] %>%  #Remove deviance
#       rownames_to_column(var = "parameter") 
#     
#     HPD.300 <- rbind(HPD.300, post.HPD)
#   }
# }
# 
# tail(HPD.300)
# levels(factor(HPD.300$interval))
# 
# #400 samples
# 
# for(i in 1:length(sim.samples.3_MCMC)){
#   for(j in 1:length(intervals)){
#     post.HPD <- combine.mcmc(sim.samples.3_MCMC[[i]]) %>% 
#       HPDinterval(prob = intervals[j]) %>%
#       data.frame() %>% 
#       mutate(interval = intervals[j], iteration = i)
#     
#     post.HPD <- post.HPD[row.names(post.HPD) %in% jags_params,] %>% #Remove deviance
#       rownames_to_column(var = "parameter") 
#       
#     HPD.400 <- rbind(HPD.400, post.HPD)
#   }
# }
# 
# tail(HPD.400)
# levels(factor(HPD.400$interval))


#--------------Create dataframes for viz and analysis----------------------

#--------------------------Normal prior-----------------------------------#
HPD.200.norm.4viz <- results.norm %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.norm, by = c("parameter", "iteration"))

HPD.300.norm.4viz <- results.norm %>% filter(total_samples == 300) %>% 
  right_join(HPD.300.norm, by = c("parameter", "iteration"))

HPD.400.norm.4viz <- results.norm %>% filter(total_samples == 400) %>% 
  right_join(HPD.400.norm, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval

HPD.200.norm.summary <- HPD.200.norm.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.300.norm.summary <- HPD.300.norm.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.300.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD400.norm.summary <- HPD.400.norm.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.400.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.norm.summary <- HPD.200.norm.summary %>% inner_join(HPD.300.norm.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD400.norm.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.norm.summary.tidy <- HPD.norm.summary %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
    ) %>% 
  mutate(prior = "normal")

date.of.simulation <- "01Jan2022"
purpose <- "estSurvLam"

write_csv(HPD.norm.summary, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose, "normalPrior.csv"))


#--------------------------Uniform prior-----------------------------------#
HPD.200.uniform.4viz <- results.uniform %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.uniform, by = c("parameter", "iteration"))

HPD.300.uniform.4viz <- results.uniform %>% filter(total_samples == 300) %>% 
  right_join(HPD.300.uniform, by = c("parameter", "iteration"))

HPD.400.uniform.4viz <- results.uniform %>% filter(total_samples == 400) %>% 
  right_join(HPD.400.uniform, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.uniform.summary <- HPD.200.uniform.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.300.uniform.summary <- HPD.300.uniform.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.300.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD400.uniform.summary <- HPD.400.uniform.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.400.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.uniform.summary <- HPD.200.uniform.summary %>% inner_join(HPD.300.uniform.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD400.uniform.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.uniform.summary.tidy <- HPD.uniform.summary %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(prior = "uniform")

date.of.simulation <- "03Jan2022"
purpose <- "testUniformPriors"
write_csv(HPD.uniform.summary, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose, "uniformPrior.csv"))


#--------------------------Hierarchical-----------------------------------#
HPD.200.hier.4viz <- results.hier %>% filter(total_samples == 200) %>% 
  right_join(HPD.200.hier, by = c("parameter", "iteration"))

HPD.300.hier.4viz <- results.hier %>% filter(total_samples == 300) %>% 
  right_join(HPD.300.hier, by = c("parameter", "iteration"))

HPD.400.hier.4viz <- results.hier %>% filter(total_samples == 400) %>% 
  right_join(HPD.400.hier, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.hier.summary <- HPD.200.hier.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.300.hier.summary <- HPD.300.hier.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.300.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD400.hier.summary <- HPD.400.hier.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.400.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.hier.summary <- HPD.200.hier.summary %>% inner_join(HPD.300.hier.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD400.hier.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.hier.summary.tidy <- HPD.hier.summary %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
  ) %>% 
  mutate(prior = "hier")


date.of.simulation <- "09Jan2022"
purpose <- "testHierarchical2"
write_csv(HPD.hier.summary, file = paste0(results_location, "HPD.summary_", date.of.simulation, "_", purpose, "hierarchical.csv"))

#create dataframes for plotting by sample size
norm.200.4viz <- HPD.norm.summary.tidy %>% filter(samples == "200.samples")
norm.300.4viz <- HPD.norm.summary.tidy %>% filter(samples == "300.samples")
norm.400.4viz <- HPD.norm.summary.tidy %>% filter(samples == "400.samples")

hier.200.4viz <- HPD.hier.summary.tidy %>% filter(samples == "200.samples")
hier.300.4viz <- HPD.hier.summary.tidy %>% filter(samples == "300.samples")
hier.400.4viz <- HPD.hier.summary.tidy %>% filter(samples == "400.samples")

uniform.200.4viz <- HPD.uniform.summary.tidy %>% filter(samples == "200.samples")
uniform.300.4viz <- HPD.uniform.summary.tidy %>% filter(samples == "300.samples")
uniform.400.4viz <- HPD.uniform.summary.tidy %>% filter(samples == "400.samples")


all.200.4viz <- norm.200.4viz %>% bind_rows(hier.200.4viz, uniform.200.4viz)
all.300.4viz <- norm.300.4viz %>% bind_rows(hier.300.4viz, uniform.300.4viz)
all.400.4viz <- norm.400.4viz %>% bind_rows(hier.400.4viz, uniform.400.4viz)

#-----------------------Viz----------------------------------------------------
#---------------------line plot -----------------------------------------------
#---------------------200 samples----------------------------------------------#
Nf.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "Nf"), 
                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "Nm"), 
       aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "surv"), 
                         aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.200samps.HPDI.plot <- ggplot(data = all.200.4viz %>% filter(parameter == "lam"), 
                        aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.200 <- ggarrange(Nf.200samps.HPDI.plot, Nm.200samps.HPDI.plot, surv.200samps.HPDI.plot, lam.200samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.200, fig.lab = "200 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")


#---------------------300 samples----------------------------------------------#
Nf.300samps.HPDI.plot <- ggplot(data = all.300.4viz %>% filter(parameter == "Nf"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.300samps.HPDI.plot <- ggplot(data = all.300.4viz %>% filter(parameter == "Nm"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.300samps.HPDI.plot <- ggplot(data = all.300.4viz %>% filter(parameter == "surv"), 
                                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.300samps.HPDI.plot <- ggplot(data = all.300.4viz %>% filter(parameter == "lam"), 
                                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.300 <- ggarrange(Nf.300samps.HPDI.plot, Nm.300samps.HPDI.plot, surv.300samps.HPDI.plot, lam.300samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.300, fig.lab = "300 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

#---------------------400 samples----------------------------------------------#
Nf.400samps.HPDI.plot <- ggplot(data = all.400.4viz %>% filter(parameter == "Nf"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("A) Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank())


Nm.400samps.HPDI.plot <- ggplot(data = all.400.4viz %>% filter(parameter == "Nm"), 
                                aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("B) Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

surv.400samps.HPDI.plot <- ggplot(data = all.400.4viz %>% filter(parameter == "surv"), 
                                  aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("C) Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

lam.400samps.HPDI.plot <- ggplot(data = all.400.4viz %>% filter(parameter == "lam"), 
                                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = prior, fill = prior, shape = prior), size = 2.5, alpha = .7) + 
  geom_smooth(aes(y = percent_in_interval, color = prior, linetype = ), se = FALSE, size = 1.5, show.legend = FALSE) + 
  geom_line(aes(y = interval.scale), size = 2) + 
  ggtitle("D) Lambda") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank(), 
        axis.title.y = element_blank())

plots.400 <- ggarrange(Nf.400samps.HPDI.plot, Nm.400samps.HPDI.plot, surv.400samps.HPDI.plot, lam.400samps.HPDI.plot, common.legend = TRUE)

annotate_figure(plots.400, fig.lab = "400 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")



head(results)

#Save to pdf
HPDI.file <- paste0(results_plots_location, "HPDI_prior_comparison.pdf")
pdf(file = HPDI.file, width = 11, height = 8)

#200 samples
annotate_figure(plots.200, fig.lab = "200 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")
#300 samples
annotate_figure(plots.300, fig.lab = "300 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")
#400 samples
annotate_figure(plots.400, fig.lab = "400 Samples", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")

dev.off()



















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
