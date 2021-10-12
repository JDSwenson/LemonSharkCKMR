library(tidyverse)
library(AIC)


#Naive model
SB.NM.df <- read_csv("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/02_IBS/Dovi_IBS_model_validation/Lemon_sharks/results/skipped_breeding/Dovi_neutral_lambda_SB_NM_07.07.2021.csv") %>% 
  mutate(Sex = factor(Sex, levels = c("M", "F", "All")))

#Tailored model
SB.TM.df <- read_csv("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/02_IBS/Dovi_IBS_model_validation/Lemon_sharks/results/skipped_breeding/Dovi_neutral_lambda_SB_TM_07.07.2021.csv") %>% 
  mutate(Sex = factor(Sex, levels = c("M", "F", "All")))

#Function to calculate AIC
aic <- function(L, k){
  AIC <- 2*k - 2*log(L)
  return(AIC)
}

#Run AIC function on naive and tailored model likelihoods
SB.NM.df2 <- SB.NM.df %>% filter(Sex == "F") %>% 
  mutate(aic = aic(Likelihood, 1))

SB.TM.df2 <- SB.TM.df %>% filter(Sex == "F") %>% 
  mutate(aic = aic(Likelihood, 1))

SB.NM.aic <- SB.NM.df2$aic
SB.TM.aic <- SB.TM.df2$aic

#Store AIC values in separate dataframe and see which model is best
aic_scores <- data.frame(cbind(SB.NM.aic, SB.TM.aic))

aic_scores %>% mutate(best_fit = ifelse(SB.NM.aic > SB.TM.aic, "TM", "NM")) %>% 
  group_by(best_fit) %>% 
  summarize(n())
