library(tidyverse)
library(RColorBrewer)

#Check results from model diagnostics
# Results
results <- read_csv("G://My Drive/Personal_Drive/R/CKMR/Model.results/Model.validation/neutralGrowth_estSurv_06Dec2021.csv")

head(results)

#Median relative bias by sample size
results %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(median = median(relative_bias), n = n())

#Check convergence
results %>% summarize(converged = sum(Rhat >= 1.1))

#Within HPD interval?
results %>% group_by(total_samples, parameter) %>% 
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