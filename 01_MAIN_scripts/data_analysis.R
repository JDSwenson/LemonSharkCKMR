library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)

#----------------Quick analysis ------------------------------
#Check results from model diagnostics
# Results
results <- read_csv("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.results/CKMR_results_14Dec2021newseeds2.csv")

# Posterior samples
jags.model.200 <- readRDS("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/CKMR_modelout_14Dec2021_200_samples_thin15_draw15000_newSeeds2")

jags.model.300 <- readRDS("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/CKMR_modelout_14Dec2021_300_samples_thin15_draw15000_newSeeds2")

jags.model.400 <- readRDS("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/CKMR_modelout_14Dec2021_400_samples_thin15_draw15000_newSeeds2")

# Breakdown of offspring for each parent
rents <- readRDS("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/CKMR_parents.breakdown_27Dec2021_thin15_draw15000")

today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
jags_params <- c("Nf", "Nm", "surv") #Specify parameters


head(results)

#Calculate relative bias
results2 <- results %>% 
  mutate(relative_bias = round(((median - truth)/truth)*100,1)) %>%
  mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N")) %>% 
  mutate(percent_sampled = round((total_samples/pop_size_mean) * 100, 0))

#Median relative bias by sample size
results.temp2 %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(median = median(relative_bias), n = n())

#Check convergence
results %>% summarize(converged = sum(Rhat >= 1.1))

#Within HPD interval?
results.temp2 %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

funny.results <- results %>% filter(truth > `97.5` | truth < `2.5`)
View(funny.results)


#Violin plot
RB_violin <- ggplot(data=results, aes(x=factor(parameter))) +
  geom_violin(aes(y=relative_bias)) +
  #ylim(-50, 160) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Iteration", y="Relative bias", title="Relative bias by iteration") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")

RB_violin





#-------------------Calculate HPD intervals---------------------------
# How many estimates fall in different HPD intervals?
#Calculate HPD interval for each iteration
intervals <- c(seq(from = 0.95, to = 0.05, by = -0.05))
iterations = 100

#200 samples
HPD.200 <- NULL
for(i in 1:length(jags.model.200)){
  for(j in 1:length(intervals)){
    post.HPD <- combine.mcmc(jags.model.200[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD <- post.HPD[row.names(post.HPD) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter")

    HPD.200 <- rbind(HPD.200, post.HPD)
  }
}

tail(HPD.200)
levels(factor(HPD.200$interval))

#300 samples
HPD.300 <- NULL
for(i in 1:length(jags.model.300)){
  for(j in 1:length(intervals)){
    post.HPD <- combine.mcmc(jags.model.300[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD <- post.HPD[row.names(post.HPD) %in% jags_params,] %>%  #Remove deviance
      rownames_to_column(var = "parameter") 
    
    HPD.300 <- rbind(HPD.300, post.HPD)
  }
}

tail(HPD.300)
levels(factor(HPD.300$interval))

#400 samples
HPD.400 <- NULL
for(i in 1:length(jags.model.400)){
  for(j in 1:length(intervals)){
    post.HPD <- combine.mcmc(jags.model.400[[i]]) %>% 
      HPDinterval(prob = intervals[j]) %>%
      data.frame() %>% 
      mutate(interval = intervals[j], iteration = i)
    
    post.HPD <- post.HPD[row.names(post.HPD) %in% jags_params,] %>% #Remove deviance
      rownames_to_column(var = "parameter") 
      
    HPD.400 <- rbind(HPD.400, post.HPD)
  }
}

tail(HPD.400)
levels(factor(HPD.400$interval))

#Create dataframes for viz and analysis
HPD.200.4viz <- results2 %>% filter(total_samples == 200) %>% 
  right_join(HPD.200, by = c("parameter", "iteration"))

HPD.300.4viz <- results2 %>% filter(total_samples == 300) %>% 
  right_join(HPD.300, by = c("parameter", "iteration"))

HPD.400.4viz <- results2 %>% filter(total_samples == 400) %>% 
  right_join(HPD.400, by = c("parameter", "iteration"))


#How many iterations fall within the expected HPD interval
HPD.200.summary <- HPD.200.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.200.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.300.summary <- HPD.300.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.300.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD400.summary <- HPD.400.4viz %>% group_by(parameter, interval) %>% 
  mutate(in_HPD_interval = ifelse(lower < truth & truth < upper, "Y", "N")) %>% 
  dplyr::summarize(percent_in_interval.400.samples = sum(in_HPD_interval == "Y")/n() * 100)

HPD.summary <- HPD.200.summary %>% inner_join(HPD.300.summary, by = c("parameter", "interval")) %>% 
  inner_join(HPD400.summary, by = c("parameter", "interval")) %>% 
  mutate(interval.scale = interval*100)

HPD.summary.tidy <- HPD.summary %>% 
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "samples",
    names_prefix = "percent_in_interval.",
    values_to = "percent_in_interval"
    )

write_csv(HPD.summary, file = paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.results/HPD.summary", today, ".csv"))

#-----------------------Viz----------------------------------------------------
#---------------------line plot -----------------------------------------------
Nf.HPDI.plot <- ggplot(data = HPD.summary.tidy %>% filter(parameter == "Nf"), 
                 aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = samples)) + 
  geom_smooth(aes(y = percent_in_interval, color = samples), se = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("Female abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())


Nm.HPDI.plot <- ggplot(data = HPD.summary.tidy %>% filter(parameter == "Nm"), 
                       aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = samples)) + 
  geom_smooth(aes(y = percent_in_interval, color = samples), se = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("Male abundance") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())

surv.HPDI.plot <- ggplot(data = HPD.summary.tidy %>% filter(parameter == "surv"), 
                       aes(x = interval.scale)) +
  geom_point(aes(y = percent_in_interval, color = samples)) + 
  geom_smooth(aes(y = percent_in_interval, color = samples), se = FALSE) + 
  geom_line(aes(y = interval.scale), size = 1.5) + 
  ggtitle("Survival") + 
  labs(x = "HPDI",
       y = "Percent in HPDI") +
  theme(legend.title = element_blank())


ggarrange(Nf.HPDI.plot, Nm.HPDI.plot, surv.HPDI.plot)


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
