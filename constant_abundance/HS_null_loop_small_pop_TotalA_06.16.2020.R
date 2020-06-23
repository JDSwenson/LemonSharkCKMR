#THIS IS THE MOST RECENT SCRIPT 6/16/2019

######## This code assumes abundance is constant through time and returns a single CKMR abundance estimate for the sampled time period ########

# Formats the data for a half-sibling model
# The code simulates kinship assignment with the same probabilities as the model i.e. it is code to test the null model
rm(list=ls())
library(tidyverse)
library(optimx)

#################Set up adult life history info ################
# Simulate and assign ages to adults based on life history
# Max age: two options: 
# 37 from https://www.saveourseasmagazine.com/lemon-sharks-old/
# 25 from White et. al. 2014
# Juvenile survival values from Dibattista 2007: 0.55 (mort = 0.45)
# Adult survival value from White et. al. 2014: 0.85 (mort=0.15)
# Litter size from Feldheim et. al. 2002 (avg=6, max=18)
# Age-at-maturity from Brown and Gruber 1988 (11.6 for M, 12.7 for F)
# About 77 juvenile lemon sharks inhabit the nursery at a given time (White et al 2014). So, 10sqrt(N) = 88 samples total

#Set population parameters for sampling
sample_ages <- 7 #specify how many cohorts are represented in the samples, since abundance estimates are given for the cohorts
Surv_age <- rep(0.85, sample_ages) #Set survival for adults
Surv_age[1] <- 0.55 #First year survival at 0.55; rest at 0.85
#Surv_age[1:2] <- c(rep(0.55,2)) #First two years, survival at 0.55
Prop_age <- rep(1, sample_ages) #Initialize variable Prop_age
for(iage in 2:length(Prop_age)) Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1] #Set proportion of animals from each cohort surviving to the next age class based on Surv_age
Prop_age=Prop_age/sum(Prop_age)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

n_yrs=sample_ages #Number of years from earliest individual to end of study

####Set up simulation parameters####
#Set Estimated_truth to set the sample size
Estimated_truth <- 200
#n_samples <- round(10*sqrt(Estimated_truth), 0) 
#Total number of samples taken over length of study
#n_samples_per_yr <- round(n_samples/length(samp_study_yrs),0)
iterations <- 200 #Set number of iterations to run in the loop

setwd("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/")
source("constant_abundance/models/get_P_lemon_HS_TotalA.R")
source("likelihood_functions/lemon_neg_log_like_HS_TotalA.R")

Total_abundance <- 200
Pars=c(log(Total_abundance))

P=get_P_lemon(Pars=Pars,P_Parent=P_Parent,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation started at ", Sys.time()))

sim_results <- data.frame(matrix(0, nrow = iterations, ncol = 5))
for(samps in 1:5){
  n_samples <- c(30, 60, 90, 120, 150)[samps]
  sampled_ages <- data.frame(matrix(0, nrow = n_samples, ncol = iterations))
  for(iter in 1:iterations) {
    #############################Sample population##############################
    Samples <- c()
    #Samples = matrix(0,n_samples,2) #matrix filled with 0s; each row is a sample(piece of tissue)
    
    #Samples=sample(c(7:1),n_samples,replace=TRUE) #birth year - probability of sampling is random
    Samples=sample(c(7:1),n_samples,replace=TRUE,prob=Prop_age) #birth year - probability of sampling based on Prop_age
    #head(Samples)
    
    Data <- t(combn(Samples, m=2))
    Data <- plyr::count(t(apply(Data, 1, sort)))
    colnames(Data)[1:2] <- c("Old_sib_birth", "Young_sib_birth")
    
    ##Add sex to third column, if helpful
    #np <- nrow(samp2)
    #samp2 <- data.frame(samp2)
    #samp2[,3] <- rbinom(n=np, size=1, prob=prob_m) #sex; 0=female
    #colnames(samp2)[1:3] <- c("Old_sib", "Young_sib", "Sex")
    
    #Remove same-cohort comparisons
    Data2 <- Data[which(Data$Old_sib_birth != Data$Young_sib_birth),]
    
    #Remove unused rows (where young sib birth is <4)
    Data2 <- Data2[which(Data2$Young_sib_birth>=4),]
    
    #Add probabilities of parentage to columns 4 & 5
    for(i in 1:nrow(Data2)){
      Data2$Parent_prob[i] <- P$P_Parent[Data2[i,1],Data2[i,2]]
    }
    #Data2
    
    #Randomly assign matches
    for(i in 1:nrow(Data2)){
      Data2$Matches[i] <- sum(rpois(n=Data2$freq[i], lambda=Data2$Parent_prob[i]))
    }
    
    Data2$Non_matches <- Data2$freq - Data2$Matches #Remove matches from non-matches
    ######################### Format pairwise comparison for CKMR ###########################
    ##Ultimately, want four dataframes:
    #1) Rows contain positive comparisons for mothers
    #2) Rows contain positive comparisons for fathers
    #3) Rows contain negative comparisons for mothers
    #4) Rows contain negtive comparisons for fathers
    
    ###################Create dataframes for CKMR model####################
    
    #head(Data2)
    Data_matches <- Data2[, c(1:2,5)]
    Data_non_matches <- Data2[, c(1:2,6)]
    
    ####LOOP####
    #if(nrow(Data_dad_yes)>=1 & nrow(Data_mom_yes)>=1){
    
    #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Parent=Data_non_matches, Pairs_Parent=Data_matches, P_Parent=P_Parent, n_yrs=n_yrs, t_start=t_start, t_end=t_end)
    
    #summary(CK_fit)
    
    #Store estimates in vector -- need to change based on number of parameters
    estimates <- exp(CK_fit$p1[1])
    
    #compute variance covariance matrix
    D=diag(length(Pars))*estimates #derivatives of transformations
    
    VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))
    
    #cat(paste("Estimated adult abundance: ",c(exp(CK_fit[1])),", SE = ",c(SE[1]),"\n"))
    #cat(paste("Estimated mature male abundance: ",c(exp(CK_fit[1,5:8])),", SE = ",c(SE[5:8]),"\n"))
    
    sim_results[iter, 1] <- round(exp(CK_fit$p1[1]),0) #N
    sim_results[iter, 2] <- round(SE[1],0) #N SE
    sim_results[iter, 3] <- sum(Data_matches[,3]) #Number of parents
    sim_results[iter, 4] <- exp(Pars[1]) #Truth
    sim_results[iter, 5] <- n_samples
    sim_results[iter, 6] <- "combined sex"
    
    save(CK_fit, file=paste0("../Results/model_objects/Lemon_CKModel_HS_TotalA_null_iteration_", iter))
    print(paste0("finished iteration", iter, " at: ", Sys.time()))
    }
  colnames(sim_results) <- c("N_est", "N_SE", "Parents_detected", "Truth", "Total_samples",  "Sex")
  colnames(sampled_ages) <- paste0("iter_", c(1:iterations))
  #f_results <- data.frame(f_results)
  #m_results <- data.frame(m_results)
  
  #f_results$Sex <- c(rep("F", times=nrow(f_results)))
  #m_results$Sex <- c(rep("M", times=nrow(m_results)))
  
  all_ests <- sim_results %>% 
    mutate(Relative_bias = round(((N_est - Truth)/Truth)*100,1))
  
  write.table(all_ests, file = paste0("Halfsib_sim_small_pop_TotalA_", n_samples, "_samps_06.16.2020.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}

sim_end_time <- Sys.time()
sim_end_time-sim_start_time
#Need to figure out which units for below
#print(paste0("Run time of simulation: ", round(sim_end_time - sim_start_time, 2), ""))
  print(paste0("Simulation finished at ", Sys.time()))
#head(sim_results, 30)