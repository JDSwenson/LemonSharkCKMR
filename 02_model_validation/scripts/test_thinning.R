#### LEMON SHARKS: DOVI'S BS MODEL

#Load packages
#devtools::install_github("BruceKendall/mpmtools")

#library(optimx)
library(tidyverse) # safe to ignore conflicts with filter() and lag()
#library(popbio) #not needed
#library(mpmtools) #not needed
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(postpack)

rm(list=ls())

#HP
source("~/R/working_directory/CKMR/LemonSharkCKMR_GitHub/01_MAIN_scripts/functions/Dovi_IBS.R")
source("~/R/working_directory/CKMR/LemonSharkCKMR_GitHub/01_MAIN_scripts/functions/pairwise_comparisons.R")

#Dell
source("~/R/working_directory/LemonSharkCKMR/01_MAIN_scripts/functions/Dovi_IBS.R")
source("~/R/working_directory/LemonSharkCKMR/01_MAIN_scripts/functions/pairwise_comparisons.R")
#######################################################################
### DATA-GENERATING MODEL


# Note on sequencing: Births happen at beginning of each year, followed by deaths 
# (i.e. a female who dies in year 10 can still give birth and have surviving pups in year 10)

####---------Simulation parameters---------####

init.adult.pop.size <- 1000 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M - has to add up to 1

YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival <- 0.825 # CHANGED FROM 0.9; Adult survival

repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 # CHANGED FROM 30; set the maximum age allowed in the simulation

mating.periodicity <- 1 # CHANGED FROM 2; number of years between mating; assigned to an individual 
                        # and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # CHANGED FROM c(1:3); vector of potential number of mates per mating
#avg.num.offspring <- 3 # NOT USED? CHANGED FROM 3; set the average number of offspring per mating (from a poisson distribution) #JDS Q

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation #JDS Q
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff

# stable age distribution - JDS Q
props <- rep(NA, max.age+1)
props[1] <- f
props[2] <- f * YOY.survival
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival
#props[repro.age+1] <- props[repro.age] * juvenile.survival + Adult.survival
for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * Adult.survival

prop.Adult <- sum(props[(repro.age+1):(max.age+1)])/sum(props)

Nages <- round(props[-1] * init.adult.pop.size) 
init.pop.size <- sum(Nages) # all ages except YOYs

burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs = t_end <- burn.in + Num.years

iterations <- 5 # CHANGED FROM 100; Number of iterations to loop over
#rseeds <- sample(1:1000000,iterations)
load("rseeds.rda")

#set.seed(885767)

####---------Sampling parameters---------####

#sample.years <- c(n_yrs - c(2:0)) #For three years of sampling
sample.years <- n_yrs #One year of sampling
#sample.size <- 300 #sample size per year
sample.vec <- c(400, 600, 800) #vector to sample over per year

#--------------Start simulation loop--------------

### MOVED SAMPLING BELOW, SO EXTRACT VARIOUS SAMPLE SIZES FROM THE SAME POPULATION

results <- NULL #initialize results array
post.samps.thin5 <- NULL #initialize array for samples from posterior

for(iter in 1:iterations) {
  set.seed(rseeds[iter])

  #Run individual based simulation.
  #Outputs a list, where the first element is the list of population parameters for each year, and the second is the population size
  out <- simulate.pop(init.pop.size = init.pop.size, 
               init.prop.female = init.prop.female,
               Nages = Nages,
               mating.periodicity = mating.periodicity,
               repro.age = repro.age,
               YOY.survival = YOY.survival,
               juvenile.survival = juvenile.survival,
               Adult.survival = Adult.survival,
               max.age = max.age,
               num.mates = num.mates,
               ff = ff,
               burn.in = burn.in,
               Num.years = Num.years)
  
  loopy.list <- out[[1]] #List of dataframes for each year of simulation
  pop.size <- out[[2]] #population parameters for each year of simulation
  
  ##===============================================================================================
  ####---------quick analysis of population growth (lambda)-----####
  
  plot(pop.size$population_size, pch=19, ylim=c(0.9*min(pop.size$population_size), 1.1*max(pop.size$population_size)))
  
  #} # end loop over iters
  
  ####---------Checking population parameters-------####
  
  #nrow(YOY.df)/length(mothers) # Average fecundity for last year; remember whether they've skipped breeding
  
  #True adult abundance
  pop.size <- pop.size %>% mutate(Total.adult.pop = Male.adult.pop + Female.adult.pop)
  
  #Calculate population growth for whole population
  total.lambda <- NULL
  for(l in 2:nrow(pop.size)){ 
    total.lambda.1 <- pop.size$population_size[l]/pop.size$population_size[l-1]
    total.lambda <- c(total.lambda, total.lambda.1)
  }
  
  #Calculate population growth for adults
  adult.lambda <- NULL
  for(l in 2:nrow(pop.size)){ 
    adult.lambda.1 <- pop.size$Total.adult.pop[l]/pop.size$Total.adult.pop[l-1]
    adult.lambda <- c(adult.lambda, adult.lambda.1)
  }
  
  total.lambda <- c(NA, total.lambda)
  adult.lambda <- c(NA, adult.lambda)
  pop.size$total.lambda <- total.lambda
  pop.size$adult.lambda <- adult.lambda
  
  #plot(total.lambda[(burn.in+1):n_yrs], pch=19)
  #abline(h=1, lty=3)
  
  mean.total.lam <- mean(pop.size$total.lambda[(burn.in+1):n_yrs], na.rm=T) # mean Lambda
  sd(pop.size$total.lambda[(burn.in+1):n_yrs], na.rm=T) # sd Lambda
  
  #Survival
sVec <- NULL
for(yr in 2:length(loopy.list)){
    sVec[yr] <- length(which(loopy.list[[yr]]$Survival=='S' & loopy.list[[yr]]$age.x>=repro.age))/length(which(loopy.list[[yr]]$age.x>=repro.age)) #survival of adults
}
  sVec <- c(NA, sVec) #Adult survival over all the simulation years
  
  length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age))/length(which(loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age)) #survival of juveniles
  length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x==0))/length(which(loopy.list[[n_yrs]]$age.x==0)) #survival of YOY
  

    #####################################################################################
  ### DATA SAMPLING
  
  for(samps in 1:length(sample.vec)){
  #try({  
     
    sample.size <- sample.vec[samps] #To loop over different sample sizes, draw a different number of samples each time
    
    ####------------------------Collect samples---------------------####
    #Initialize dataframes
    sample.df_all.info <- NULL
    sample.df_temp <- NULL
    
    #Sample population each year in sample.years and make dataframe of samples with all metadata
    for(i in sample.years){
      sample.df_temp <- loopy.list[[i]] %>% dplyr::slice_sample(n = sample.size)
      sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp)
    }
    
    #Keep just one instance of each individual (to avoid self-recapture) and sort by birth year so when we make the pairwise comparison matrix, Ind_1 is always older than Ind_2 (bc it comes first)
    sample.df_all.info <- sample.df_all.info %>% distinct(indv.name, .keep_all = TRUE) %>% 
      dplyr::arrange(birth.year, desc()) 
    
    
    ####-------------Construct pairwise comparison matrix--------------####
    pairwise.out <- build.pairwise(sample.df_all.info)
    
    #Separate list elements into different dataframes
    pairwise.df_all.info <- pairwise.out[[1]]
    positives <- pairwise.out[[2]]
    
    #Check that there are no repeat comparisons -- compare number of distinct comparisons to total number of comparisons.
    #Should return TRUE
    pairwise.df_all.info %>% distinct(Ind_1, Ind_2) %>% nrow() == nrow(pairwise.df_all.info)

    
    nrow(positives) #How many positive comparisons total?
    #head(mom_comps)
    #head(dad_comps)
    
    #Remove full sibs (and other unwanted individuals) from sample dataframe
    #Should incorporate full sibs into the model at some point ... 
    filtered.samples <- filter_indvs(sample.df_all.info, positives)
    
    #Reconstruct pairwise comparison dataframes from filtered sample list
    pairwise.out.filt <- build.pairwise(filtered.samples)
    
    #Separate list elements into different dataframes
    pairwise.df_all.info.filt <- pairwise.out.filt[[1]] 
    positives.filt <- pairwise.out.filt[[2]]
    mom_comps <- pairwise.out.filt[[3]] 
    dad_comps <- pairwise.out.filt[[4]]
        

        #-------------Kinship probabilities - Half-sib-------------------
    min_cohort <- n_yrs-40 # CHANGED THIS FROM MAX.AGE; set first year for calculating mean (arbitrary)
    mean.adult.lambda <- mean(pop.size$adult.lambda[min_cohort:n_yrs], na.rm=T) # mean Lambda over years of estimation for adults ### HOW DO WE KNOW THIS?
        
    #---------------Do I need anything in this section?--------------------#
    #m_adult_age <- f_adult_age <- c(repro.age:max.age) # Set ages at which males and females are mature. Called by kinship probability function.
    #pop_growth_all_mean <- mean(pop.size$Lambda[1:nrow(pop.size)], na.rm=T)
    
    
    #### Fit Bayesian model ####
    ####PICK UP HERE AFTER 10/27/2021####
    
    ################ STEP 1: PREPARE DATA #################
    #Need priors for:
    #Number of adults (uninformative)
    #maybe eventually survival (beta -- conjugate prior for binomial)
    
    #Define data
    jags_data = list(
      #Moms
      MHSP = mom_comps$yes,
      mom_n_comps = mom_comps$all,
      mom_ys_birth = mom_comps$Ind_2_birth,
      mom_os_birth = mom_comps$Ind_1_birth,
      mom_yrs = nrow(mom_comps),
      
      #Dads
      FHSP = dad_comps$yes,
      dad_n_comps = dad_comps$all,
      dad_ys_birth = dad_comps$Ind_2_birth,
      dad_os_birth = dad_comps$Ind_1_birth,
      dad_yrs = nrow(dad_comps),
      
      #Fix other potential parameters
      #surv = surv,
      lam = mean.adult.lambda,
      min_cohort = min_cohort
    )
    
    #Define initial values for priors
    
    ################### STEP 2: SPECIFY JAGS MODEL CODE ##################
    ##CKMR code
    HS_model = function(){
      #PRIORS
      Nf ~ dnorm(0, 1.0E-6) #Uninformative prior for female abundance
      Nm ~ dnorm(0, 1.0E-6) ##Uninformative prior for male abundance
      surv ~ dbeta(1 ,1) #Assumes constant survival of males and females
      
      #Likelihood
      for(i in 1:mom_yrs){
        MHSP[i] ~ dbin((surv^(mom_ys_birth[i] - mom_os_birth[i]))/(Nf*lam^(mom_ys_birth[i]-min_cohort)), mom_n_comps[i])
      }
      for(j in 1:dad_yrs){
        FHSP[j] ~ dbin((surv^(dad_ys_birth[j] - dad_os_birth[j]))/(Nm*lam^(dad_ys_birth[j]-min_cohort)), dad_n_comps[j])
      }
    }
    
    jags_file = paste0("~/R/working_directory/LemonSharkCKMR/models/HS_neutralGrowth_estNSurv_iteration_", iter, ".txt")
    write_model(HS_model, jags_file)
    
    
    ########### STEP 3: SPECIFY INITIAL VALUES ##################
    jags_inits = function(nc) {
      inits = list()
      for(c in 1:nc){
        inits[[c]] = list(
          surv = 0.8,
          Nf = rnorm(1, mean = 500, sd = 100),
          Nm = rnorm(1, mean = 500, sd = 100)
        )
      }
      return(inits)
    }
    
    ########## STEP 4: SET NODES TO MONITOR ################
    jags_params = c("Nf", "Nm", "surv")
    n_params = length(jags_params) #used to autofill dataframe later
    
    
    ########### STEP 5: SET MCMC DIMENSIONS ################
    jags_dims = c(
      ni = 5000,  # number of post-burn-in samples per chain
      nb = 5000,  # number of burn-in samples
      nt = 5,     # thinning rate
      nc = 2      # number of chains
    )
    
    
    ##### STEP 6: RUN THE MODEL WITH JAGS #####
    
    post = jagsUI::jags.basic(data = jags_data, #If using postpack from AFS workshop
                              
                              #post = rjags::jags(data = jags_data, #If wanting to use other diagnostics
                              model.file = jags_file,
                              inits = jags_inits(jags_dims["nc"]),
                              parameters.to.save = jags_params,
                              n.adapt = 1000,
                              n.iter = sum(jags_dims[c("ni", "nb")]),
                              n.thin = jags_dims["nt"],
                              n.burnin = jags_dims["nb"],
                              n.chains = jags_dims["nc"],
                              parallel = F
    )
    
    ##### STEP 7: CONVERGENCE DIAGNOSTICS #####
    # view convergence diagnostic summaries for all monitored nodes
    model_summary <- data.frame(t(post_summ(post, jags_params, Rhat = T, neff = T)))
    
    #Save samples from posterior
    post.samps.1 <- as_tibble(post[[1]]) %>% mutate(iter = iter, chain = 1)
    post.samps.2 <- as_tibble(post[[2]]) %>% mutate(iter = iter, chain = 2)

##### Compile and report results ####
    #Combine above to make dataframe with truth and estimates side-by-side
    #store years from youngest sibling in comparisons to end of study
    yrs <- c(min_cohort:t_end)
    
    #Extract true values from year of estimation (ie min_cohort)
    Mom_truth <- round(pop.size$Female.adult.pop[min_cohort],0)
    Dad_truth <- round(pop.size$Male.adult.pop[min_cohort], 0)
    surv_truth <- round(mean(sVec[min_cohort:n_yrs]), 4)
    #Adult_truth <- round(pop.size$Total.adult.pop[min_cohort], 0)
    Mom_min <- min(pop.size$Female.adult.pop[min_cohort:n_yrs])
    Mom_max <- max(pop.size$Female.adult.pop[min_cohort:n_yrs])
    Dad_min <- min(pop.size$Male.adult.pop[min_cohort:n_yrs])
    Dad_max <- max(pop.size$Male.adult.pop[min_cohort:n_yrs])
    #Adult_min <- min(pop.size$Total.adult.pop[min_cohort:n_yrs])
    #Adult_max <- max(pop.size$Total.adult.pop[min_cohort:n_yrs])
    surv_min <- min(sVec[min_cohort:n_yrs])
    surv_max <- max(sVec[min_cohort:n_yrs])
    
    #Create dataframe of estimates and truth
    estimates <- model_summary %>% rownames_to_column(var = "parameter") %>% 
      mutate(min = c(Mom_min, Dad_min, surv_min), max = c(Mom_max, Dad_max, surv_max), truth = c(Mom_truth, Dad_truth, surv_truth)) %>%
      as_tibble() %>% 
      rename(`50` = X50., `2.5` = X2.5., `97.5` = X97.5.)
    
    #Extract more metrics that can help with troubleshooting
    total_samples <- sample.size * length(sample.years)
    pop_size_mean <- round(mean(pop.size$population_size[min_cohort:n_yrs]),0)
    
    metrics <- cbind(c(sum(mom_comps[,3]), sum(dad_comps[,3]), sum(mom_comps[,3]) + sum(dad_comps[3])), 
                     c(rep(mean.adult.lambda, times = n_params)), 
                     c(rep(total_samples, times = n_params)),
                     c(rep(pop_size_mean, times = n_params)),
                     c(rep(iter, times = n_params)))
    colnames(metrics) <- c("parents_detected", "mean_adult_lambda", "total_samples", "pop_size_mean", "iteration")
    
    #-----------------Loop end-----------------------------    
    #Bind results from previous iterations with current iteration
    results <- rbind(results, cbind(estimates, metrics))
    post.samps.thin5 <- rbind(post.samps.thin5, post.samps.1, post.samps.2)
  
#  }) # end try clause    
  } # end loop over sample sizes
  
#  write.table(results, file = paste0("~/R/working_directory/temp_results/neutralGrowth_estSurv_iteration", iter, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
  
  print(paste0("finished iteration", iter, " at: ", Sys.time()))
} # end loop over iterations

##################################################################################
### SAVE AND CHECK RESULTS


#Calculate relative bias for all estimates
results2 <- results %>% 
  mutate(relative_bias = round(((mean - truth)/truth)*100,1)) %>% # CHANGED TABLE NAME SO CAN BUILD & CHECK RESULTS ITERATIVELY
  mutate(in_interval = ifelse(`2.5` < truth & truth < `97.5`, "Y", "N")) %>% 
  mutate(percent_sampled = round((total_samples/pop_size_mean) * 100, 0))

#Within credible interval?
results2 %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(median = median(relative_bias), n = n())

#Home computer: Dell Precision
write.table(results2, file = paste0("~/R/working_directory/LemonSharkCKMR/02_model_validation/results/Lemon_Shark_life_history/model_validation/test.thinning_thin1.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#MGHPCC
#write.table(results2, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/Dovi_lambdaModel_06_22.2021_neutralPopGrowth.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#write.table(age_dist, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_age.distributions_", total_samples, ".samples_02.10.2021_ages.correct_age.dist.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#write.table(survival_at_age, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_survival.at.age_", total_samples, ".samples_02.10.2021_ages.correct_surv.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
#}

####-------Quick viz of results--------####

#Box plot of relative bias
ggplot(data=results2, aes(x=factor(total_samples))) +
  geom_boxplot(aes(y=relative_bias, fill=parameter)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")



##################################################################################