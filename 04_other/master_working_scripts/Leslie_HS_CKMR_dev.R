#rm(list=ls())
#Sex-specific estimates of N;
#Population growth fixed to observed population growth of the whole population.
#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(foreach)
library(parallel)
library(doParallel)
library(optimx)
#Load individual packages because tidyverse won't load on cluster
library(plyr)
library(dplyr)
library(tidyr)
library(popbio)
library(mpmtools)

set.seed(47)
#----------Simulation parameters-----------
t_start = 120 #First year we take samples
t_end = 125 #Last year we take samples
#min_est_cohort <- 42 #First year we are estimating abundance
#max_est_cohort <- t_end #Last year we are estimating abundance
#est_yrs <- min_est_cohort:max_est_cohort #All years for which we are estimating abundance
n_yrs <- t_end #Number of years of study

samp_yrs = t_start:t_end #Number of years being sampled

#-----------------Leslie Matrix--------------------
#Simulate a population with a Leslie Matrix
batchSize <- 4 #Average size of brood: 8 divided by 2 (for skipped-breeding)
fb <- batchSize/2 #Batch size for Leslie Matrix - divided by two for females only
maxAge <- 30 #max age of lemon sharks
surv <- 0.85 #adult survival

#Make Leslie Matrix
Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = c(0.6, 0.7, 0.81, rep(surv, times = 27), 0), #survival
  mx = c(rep(0, times = 12), rep(fb, times = 19)) #age-specific birth rates (female proportion of the population only)
)


A1_pre <- make_Leslie_matrix(Leslie_input) #pre-birth census ie no year-0 individuals
#View(A1)
A1_post <- pre_to_post(Amat = A1_pre, S0 = .6) #convert to post-birth census that has year 0 individuals
#View(A1_post)

#Calculate dominant eigenvalue (i.e. population growth rate) - should be the same for pre and post
(lam <- lambda1(A1_pre)) #Save lambda for use in CKMR model
lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)

stable_age <- mpmtools::stable_stage(A1_pre) #fishSim seems to use a pre-breeding census for the age distribution for the Founder population

#--------------Population simulation-------------------
#Project population forward in time
Abundance_year0 <- c(20000, rep(0, times = (maxAge-1)))

Year1 <- A1_pre %*% Abundance_year0

nYears <- n_yrs        # set the number of years to project
TMat <- A1_pre     # define the projection matrix
InitAbund <- Abundance_year0    # define the initial abundance

allYears <- matrix(0,nrow=nrow(TMat),ncol=nYears+1) # build a storage array for all abundances!
allYears[,1] <- InitAbund  # set the year 0 abundance                                    
for(t in 2:(nYears+1)){   # loop through all years and fill in Leslie Matrix
  allYears[,t] <-  TMat %*% allYears[,t-1]
}

#in allYears, columns are years, rows are ages
allYears <- data.frame(allYears)
colnames(allYears) <- paste("Yr_", c(1:(nYears+1)))
allYears <- allYears %>% mutate_all(round, digits=0)

#Adult population size by year
adult_truth <- allYears %>% slice(12:30) %>% 
  summarise_all(sum)

#Total population size by year
all_truth <- allYears %>% summarize_all()

#Create maturity vector
mat_m = mat_f = mat_a <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
#Knife-edge maturite for males, females, and all
mat_m[12:maxAge] =  mat_f[12:maxAge] = mat_a[12:maxAge] <- 1  

#Set ages at which males and females are mature for use in the CKMR model.
m_adult_age <- c(12:maxAge)
f_adult_age <- c(12:maxAge)

iterations <- 500


#------------------------Start loop-----------------
#Loop over sample size
for(samps in 1:4){
  results <- NULL #initialize results (main dataframe to store results from each iteration)
  n_samples <- c(300, 500, 750, 1000)[samps] #Set number of samples
  
  set.seed(47) #Set the seed so for each loop sample size is the only variable that changes
  for(iter in 1:iterations) {
    #index <- (iter*2) + 1 #Index for filling in the matrix at the end of each loop over n_samples -- not needed anymore

    Samples <- c() #initialize sample vector
    
    Samples <- sample(c(1:maxAge), n_samples, replace=TRUE, prob=stable_age) #assign age to samples based on the stable age distribution. 
    
    Samples2 <- t_end-Samples #Convert age to birth year
    min_est_cohort <- min(Samples2) #Reference year for abundance estimate is the birth year of the oldest individual
    
    #Starting parameters
    N_f = N_m <- adult_truth[,min_est_cohort]/2
    N_a <- adult_truth[,min_est_cohort]
    Pars1 <- c(log(N_f),log(N_m))
    Pars2 <- c(log(N_a))
    
    #Create pairwise comparison matrix of all sampled individuals (ie birth years)
    Data <- t(combn(Samples2, m=2))
    Data <- plyr::count(t(apply(Data, 1, sort)))
    colnames(Data)[1:2] <- c("Old_sib_birth", "Young_sib_birth")
    
    ##Add sex to third column, if helpful
    #np <- nrow(samp2)
    #samp2 <- data.frame(samp2)
    #samp2[,3] <- rbinom(n=np, size=1, prob=prob_m) #sex; 0=female
    #colnames(samp2)[1:3] <- c("Old_sib", "Young_sib", "Sex")
    
    #Remove same-cohort comparisons
    Data2 <- Data[which(Data$Old_sib_birth != Data$Young_sib_birth),]

#--------------------Source models and likelihood functions----------
    #For PC
    setwd(".")
    source("./functions/get_P_lemon_HS_sex-specific.R")
    source("./functions/lemon_neg_log_lik_HS_sex_specific.R")
    source("./functions/get_P_lemon_HS_sex-aggregated.R")
    source("./functions/lemon_neg_log_lik_HS_sex-aggregated.R")    
    
    # #For cluster
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS.R")
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_Sex_specific.R")
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS_TotalA.R")
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_TotalA.R")
    
    
    #Add probabilities of kinship from model to pairwise comparison dataframe 
    for(i in 1:nrow(Data2)){
      Data2$Mom_prob[i] <- P$P_Mother[Data2[i,1],Data2[i,2]]
      Data2$Dad_prob[i] <- P$P_Father[Data2[i,1],Data2[i,2]]
      Data2$Parent_prob[i] <- P_TotalA$P_Parent[Data2[i,1],Data2[i,2]]
    }
    #Data2
    
    #Assign kinship based on probabilities from model
    for(i in 1:nrow(Data2)){
      Data2$Mom_Matches[i] <- sum(rbinom(n=Data2$freq[i], size=1, p=Data2$Mom_prob[i]))
      Data2$Dad_Matches[i] <- sum(rbinom(n=Data2$freq[i], size=1, p=Data2$Dad_prob[i]))
      Data2$Total_Matches[i] <- sum(rbinom(n=Data2$freq[i], size=1, p=Data2$Parent_prob[i]))
    }
    
    
#Separate pairwise comparison dataframe into positive and negative comparisons
    ## Ultimately, want four dataframes for the sex-specific model:
    #1) Positive comparisons for mothers
    #2) Positive comparisons for fathers
    #3) Negative comparisons for mothers
    #4) Negtive comparisons for fathers
    ## And two for the sex-aggregated model
    #5) Positive comparisons for all adult
    #6) Negative comparisons for all adults
    mom_positives <- Data2[which(Data2$Mom_Matches >0), c(1:2,7)]
    dad_positives <- Data2[which(Data2$Dad_Matches >0), c(1:2,8)]
    dad_negatives = mom_negatives = parent_negatives <- Data2
    dad_negatives[,3] <- dad_negatives[,3] - dad_negatives[,8]
    dad_negatives <- dad_negatives[,c(1:3)]
    mom_negatives[,3] <- mom_negatives[,3] - mom_negatives[,7]
    mom_negatives <- mom_negatives[,c(1:3)]
    parent_positives <- Data2[which(Data2$Total_Matches >0), c(1:2,9)]
    parent_negatives[,3] <- parent_negatives[,3] - parent_negatives[,9]
    parent_negatives <- parent_negatives[,c(1:3)]
    
#--------------Fit models-----------------
    #Fit model = nlminb version
    # CK_fit <- nlminb(start = Pars, objective = lemon_neg_log_lik, hessian = TRUE, Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, t_start=t_start, t_end=t_end)
    # 
    # #exp(CK_fit$par)
    # #compute variance covariance matrix - nlminb
    # D=diag(length(Pars))*c(exp(CK_fit$p[1]),exp(CK_fit$p[2])) #derivatives of transformations
    # VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    # VC = (t(D)%*%VC_trans%*%D) #delta method
    # SE=round(sqrt(diag(VC)),0)    
    
    #Fit model - sex-specific
    CK_fit1 <- optimx(par=Pars1,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, t_start=t_start, t_end=t_end)
    
    #Fit model - sex-aggregated
    CK_fit2 <- optimx(par=Pars2,fn=lemon_neg_log_lik_TotalA,hessian=TRUE, method="BFGS", Negatives_Parent=parent_negatives, Pairs_Parent=parent_positives, P_Parent=P_Parent, t_start=t_start, t_end=t_end)
    
    #summary(CK_fit1)
    #exp(CK_fit1[1:2])
    
    #compute variance covariance matrix - sex-specific
    D1=diag(length(Pars1))*c(exp(CK_fit1$p1[1]),exp(CK_fit1$p2[1])) #derivatives of transformations
    VC_trans1 = solve(attr(CK_fit1, "details")["BFGS" ,"nhatend"][[1]])
    VC1 = (t(D1)%*%VC_trans1%*%D1) #delta method
    SE1=round(sqrt(diag(VC1)),0)
    
    
    #compute variance covariance matrix - sex-aggregated
    D2=diag(length(Pars2))*exp(CK_fit2$p1[1]) #derivatives of transformations
    VC_trans2 = solve(attr(CK_fit2, "details")["BFGS" ,"nhatend"][[1]])
    VC2 = (t(D2)%*%VC_trans2%*%D2) #delta method
    SE2 = round(sqrt(diag(VC2)),0)
    
    
    #Combine above to make dataframe with truth and estimates side-by-side
    #store years from youngest sibling in comparisons to end of study
    yrs <- c(min(mom_positives$Young_sib_birth, dad_positives$Young_sib_birth, parent_positives$Young_sib_birth):t_end)
    
    estimates <- data.frame(cbind(round(exp(c(CK_fit1$p1[1], CK_fit1$p2[1], CK_fit2$p1[1])),0)), c(SE1, SE2), c("F", "M", "All")) #CKMR estimates, SE, and sex
    estimates <- cbind(estimates, c(rep(as.numeric(adult_truth[min_est_cohort]/2), times=2), as.numeric(adult_truth[min_est_cohort]))) #True adult abundance in first year of data
    colnames(estimates) <- c("CKMR_estimate", "SE", "Sex", "Truth")
    estimates
    
    #Also store number of parents detected, lambda, and number of samples
    metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3]), sum(parent_positives[,3])), c(rep(lam, times = 3)), c(rep(n_samples, times=3)))
    colnames(metrics) <- c("Parents_detected", "Pop_growth", "Samples")
    
    
    #-----------------Loop end-----------------------------
    #Add results from this loop to results from previous loops
    results <- rbind(results, cbind(estimates, metrics))
    
    
    print(paste0("finished iteration", iter, " at: ", Sys.time()))
  }
  
  #Calculate relative bias
  results <- results %>% 
    mutate(Relative_bias = round(((CKMR_estimate - Truth)/Truth)*100,1))

  #Median relative bias by sex  
  #results %>% group_by(Sex) %>% 
  #  summarize(median = median(Relative_bias), n = n())
  
  write.table(results, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/Leslie.null_", n_samples, ".samples_02.03.2021_ages.correct.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}

#Quickly visualize results
library(ggpubr)
ggplot(data=results, aes(x=factor(Samples))) +
  geom_boxplot(aes(y=Relative_bias, fill=Sex)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="A) Relative Bias (Ages correct)") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")
