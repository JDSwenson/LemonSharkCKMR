#### LEMON SHARKS: DOVI'S IBS MODEL

#Load packages
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

source("./01_MAIN_scripts/functions/Dovi_IBS.R")
source("./01_MAIN_scripts//functions/pairwise_comparisons.R")

########## DATA-GENERATING MODEL ##########
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
mating.periodicity <- 1 # CHANGED FROM 2; number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # CHANGED FROM c(1:3); vector of potential number of mates per mating
#avg.num.offspring <- 3 # NOT USED? CHANGED FROM 3; set the average number of offspring per mating (from a poisson distribution)

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff

# stable age distribution
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
n_yrs <- burn.in + Num.years #Total number of simulation years
min_cohort <- n_yrs - 5 # Set year of estimation

iterations <- 100 # Number of iterations to loop over
rseeds <- sample(1:1000000,iterations)
#load("rseeds.rda")

#set.seed(885767)


####-------------- Sampling parameters -------------####
#sample.years <- c(n_yrs - c(2:0)) #For three years of sampling
sample.years <- n_yrs #One year of sampling
#sample.size <- 300 #sample size per year
sample.vec <- c(400, 600, 800) #vector to sample over per year


####-------------- Start simulation loop -------------------####
# Moved sampling below so extract different sample sizes from same population
results <- NULL #initialize results array

for(iter in 1:iterations) {
  set.seed(rseeds[iter])

  #Run individual based simulation
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

#Outputs a list, where the first element is a list of various population parameters for each year, and the second is the population size for each year
  loopy.list <- out[[1]] #List of dataframes for each year of simulation
  pop.size <- out[[2]] #population parameters for each year of simulation
  

####------------- Examine simulation results ---------------####
  plot(pop.size$population_size, pch=19, ylim=c(0.9*min(pop.size$population_size), 1.1*max(pop.size$population_size))) #Plot population size through time
  
  #nrow(YOY.df)/length(mothers) # Average fecundity for last year; remember whether they've skipped breeding
  
  #Make dataframe of abundance for each year
  pop.size <- pop.size %>% mutate(Total.adult.pop = Male.adult.pop + Female.adult.pop)
  
  #Calculate population growth for whole population
  total.lambda <- NULL
  for(l in 2:nrow(pop.size)){ 
    total.lambda.1 <- pop.size$population_size[l]/pop.size$population_size[l-1]
    total.lambda <- c(total.lambda, total.lambda.1)
  }
  
  #Calculate population growth for adults only
  adult.lambda <- NULL
  for(l in 2:nrow(pop.size)){ 
    adult.lambda.1 <- pop.size$Total.adult.pop[l]/pop.size$Total.adult.pop[l-1]
    adult.lambda <- c(adult.lambda, adult.lambda.1)
  }
  
  #Add NA to first element of population growth vectors
  total.lambda <- c(NA, total.lambda)
  adult.lambda <- c(NA, adult.lambda)
  mean.adult.lambda <- mean(pop.size$adult.lambda[min_cohort:n_yrs], na.rm=T) # mean Lambda over years of estimation for adults ### HOW DO WE KNOW THIS?
  
  #Add population growth per year to pop.size dataframe
  pop.size$total.lambda <- total.lambda
  pop.size$adult.lambda <- adult.lambda
  
  #plot(total.lambda[(burn.in+1):n_yrs], pch=19)
  #abline(h=1, lty=3)
  
  (mean.total.lam <- mean(pop.size$total.lambda[(burn.in+1):n_yrs], na.rm=T)) # mean population growth for whole population
  sd(pop.size$total.lambda[(burn.in+1):n_yrs], na.rm=T) # sd Lambda
  
  #Calculate survival for each year
sVec <- NULL #Make empty vector to save yearly survival rates

#Store annual survival of adults
for(yr in 2:length(loopy.list)){
    sVec[yr] <- length(which(loopy.list[[yr]]$Survival=='S' & loopy.list[[yr]]$age.x>=repro.age))/length(which(loopy.list[[yr]]$age.x>=repro.age))
}
  sVec <- c(NA, sVec) #Add NA to survival vector for first year of simulation

  length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age))/length(which(loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age)) #survival of juveniles
  
  length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x==0))/length(which(loopy.list[[n_yrs]]$age.x==0)) #survival of YOY
  

  ####--------------------Collect samples---------------------####
  #Loop over sample sizes stored in sample.vec  
  for(samps in 1:length(sample.vec)){

        sample.size <- sample.vec[samps] #Specify sample size
    
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
    
    #Check that there are no repeat comparisons -- compare number of distinct comparisons to total number of comparisons. Should return TRUE
    pairwise.df_all.info %>% distinct(Ind_1, Ind_2) %>% nrow() == nrow(pairwise.df_all.info)

    
    nrow(positives) #How many positive comparisons total?
    
    #Remove full sibs (and other unwanted individuals) from sample dataframe, if present
    #Should incorporate full sibs into the model at some point ... 
    filtered.samples <- filter_indvs(sample.df_all.info, positives)
    
    #Reconstruct pairwise comparison dataframes from filtered sample list
    pairwise.out.filt <- build.pairwise(filtered.samples)
    
    #Separate list elements into different dataframes
    pairwise.df_all.info.filt <- pairwise.out.filt[[1]] 
    positives.filt <- pairwise.out.filt[[2]]
    mom_comps <- pairwise.out.filt[[3]] 
    dad_comps <- pairwise.out.filt[[4]]


    ########## Fit CKMR model ##########
    #-------------- STEP 1: PREPARE DATA ----------------#
    #Need priors for:
    #Number of adults (uninformative)
    #Survival (beta -- conjugate prior for binomial; uninformative)
    
    #Define data
    jags_data = list(
      #Moms
      MHSP = mom_comps$yes, # Positive maternal half-sibs; Y
      mom_n_comps = mom_comps$all, # Number of total maternal comparisons; R
      mom_ys_birth = mom_comps$Ind_2_birth, # birth year of younger sib; b
      mom_os_birth = mom_comps$Ind_1_birth, # birth year of older sib; a
      mom_yrs = nrow(mom_comps), # number of cohort comparisons to loop over 
      
      #Dads
      FHSP = dad_comps$yes, # Positive paternal half-sibs; Y
      dad_n_comps = dad_comps$all, # Number of total paternal comparisons; R
      dad_ys_birth = dad_comps$Ind_2_birth, # birth year of younger sib; b
      dad_os_birth = dad_comps$Ind_1_birth, # birth year of older sib; a
      dad_yrs = nrow(dad_comps), # number of cohort comparisons to loop over
      
      #Fix other potential parameters
      #surv = surv,
      lam = mean.adult.lambda, # fix lambda to mean lambda of adults over estimation period
      min_cohort = min_cohort # estimation year i.e. year the estimate will be focused on
    )
    

    #----------------- STEP 2: SPECIFY JAGS MODEL CODE ---------------#
    HS_model = function(){

      #PRIORS
      Nf ~ dnorm(0, 1.0E-6) # Uninformative prior for female abundance
      Nm ~ dnorm(0, 1.0E-6) # Uninformative prior for male abundance
      surv ~ dbeta(1 ,1) # Uninformative prior for adult survival
      
      #Likelihood
      for(i in 1:mom_yrs){ # Loop over maternal cohort comparisons
        MHSP[i] ~ dbin((surv^(mom_ys_birth[i] - mom_os_birth[i]))/(Nf*lam^(mom_ys_birth[i]-min_cohort)), mom_n_comps[i]) # Sex-specific CKMR model equation
      }
      for(j in 1:dad_yrs){ # Loop over paternal cohort comparisons
        FHSP[j] ~ dbin((surv^(dad_ys_birth[j] - dad_os_birth[j]))/(Nm*lam^(dad_ys_birth[j]-min_cohort)), dad_n_comps[j]) # Sex-specific CKMR model equation
      }
    }

    # Write model    
    jags_file = paste0("./models/HS_neutralGrowth_estSurv_iteration_", iter, ".txt")
    write_model(HS_model, jags_file)
    
    
    #------------ STEP 3: SPECIFY INITIAL VALUES ---------------#
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
    
    #------------ STEP 4: SET NODES TO MONITOR ---------------#
    jags_params = c("Nf", "Nm", "surv")
    n_params = length(jags_params) #used to autofill dataframe later
    
    
    #------------- STEP 5: SET MCMC DIMENSIONS ---------------#
    jags_dims = c(
      ni = 5000,  # number of post-burn-in samples per chain
      nb = 5000,  # number of burn-in samples
      nt = 1,     # thinning rate
      nc = 2      # number of chains
    )
    
    
    #---------------- STEP 6: RUN JAGS ---------------#
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
    
    #---------------- STEP 7: CONVERGENCE DIAGNOSTICS -----------------#
    # view convergence diagnostic summaries for all monitored nodes
    model_summary <- data.frame(t(post_summ(post, jags_params, Rhat = T, neff = T)))
  }
}
