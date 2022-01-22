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
library(runjags)
library(postpack)
library(coda)
library(ggmcmc)

rm(list=ls())

source("./01_MAIN_scripts/functions/Dovi_IBS.R")
source("./01_MAIN_scripts//functions/pairwise_comparisons.R")


#----------------Set output file locations ------------------------------
#Check results from model diagnostics
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

#rseeds <- sample(1:1000000,iterations)
#save(rseeds, file = "rseeds_12.27.rda")
load("rseeds_12.27.rda")

seeds <- "Seeds12.27"

purpose <- "prior_posterior_density_SDMax"

temp_location <- "~/R/working_directory/temp_results/"
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/"
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/models/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.results/"

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "HS_"
parents_prefix <- "CKMR_parents.breakdown"
sample.prefix <- "CKMR_sample.info"


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
est.year <- n_yrs - 5 # Set year of estimation

#rseeds <- sample(1:1000000,iterations)
load("rseeds_12.27.rda")

####------------- MCMC parameters ----------------####
ni <- 30000 # number of post-burn-in samples per chain
nb <- 40000 # number of burn-in samples
nt <- 15     # thinning rate
nc <- 2      # number of chains


####-------------- Sampling parameters -------------####
#sample.years <- c(n_yrs - c(1:0)) #For two years of sampling
sample.years <- n_yrs #One year of sampling
#sample.size <- 300 #sample size per year
sample.vec <- c(200, 300, 400) #vector to sample over per year

####-------------- Start simulation loop -------------------####
# Moved sampling below so extract different sample sizes from same population

# Initialize arrays for saving results
 results <- NULL
 sample.info <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL

sim.samples.1 <- paste0(sample.vec[1]*length(sample.years), ".samples")
sim.samples.2 <- paste0(sample.vec[2]*length(sample.years), ".samples")
sim.samples.3 <- paste0(sample.vec[3]*length(sample.years), ".samples")

####Initialize array from previous checkpoint
#Results
# results <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, "_iter_", iter, ".csv"))
# 
# #Model output for diagnostics
# sims.list.1 <- readRDS(file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))
# 
# sims.list.2 <- readRDS(file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))
# 
# sims.list.3 <- readRDS(file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))

# Detailed info on samples and parents to examine in more detail
# sample.info <- readRDS(file = paste0(temp_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))


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
  parents.tibble <- out[[3]] #%>% #Tibble for each parent for each year to check the distribution later
    #mutate(iteration = iter) 

  
  #Save results from the simulation
  source("./01_MAIN_scripts/functions/query_results.R")


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
    
    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    
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
    mom_comps <- pairwise.out.filt[[3]] %>% as_tibble()
    dad_comps <- pairwise.out.filt[[4]] %>% as_tibble()


    ########## Fit CKMR model ##########
    #-------------- STEP 1: PREPARE DATA ----------------#
    #Need priors for:
    #Number of adults (uninformative)
    #Survival (beta -- conjugate prior for binomial; uninformative)
    lam.sd <- 0.02277 #Value derived from Leslie matrix
    lam.tau <- 1/(lam.sd^2) 
    sd.vec <- c(100, 500, 1000, 2000, 5000, 10000)
  
      for(i in 1:length(sd.vec)){
    sd.max <- sd.vec[i]
    tau.max <- 1/(sd.max^2)
    
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
      est.year = est.year, # estimation year i.e. year the estimate will be focused on
      lam.tau = lam.tau,
      tau.max = tau.max
    )
    

    #----------------- STEP 2: SPECIFY JAGS MODEL CODE ---------------#
    #Convert tau to SD (for interpretation)
    #tau <- 1E-6
    #(sd <- sqrt(1/tau))
    #tau <- 
    
    HS_model = function(){

      #PRIORS
      mu ~ dunif(1, 10000)
      tau ~ dunif(.000000000000001, tau.max)
      Nf ~ dnorm(mu, tau) # Uninformative prior for female abundance
      Nm ~ dnorm(mu, tau) # Uninformative prior for male abundance
      surv ~ dbeta(1 ,1) # Uninformative prior for adult survival
      lam ~ dnorm(1, lam.tau)
      
      #Likelihood
      for(i in 1:mom_yrs){ # Loop over maternal cohort comparisons
        MHSP[i] ~ dbin((surv^(mom_ys_birth[i] - mom_os_birth[i]))/(Nf*lam^(mom_ys_birth[i]-est.year)), mom_n_comps[i]) # Sex-specific CKMR model equation
      }
      for(j in 1:dad_yrs){ # Loop over paternal cohort comparisons
        FHSP[j] ~ dbin((surv^(dad_ys_birth[j] - dad_os_birth[j]))/(Nm*lam^(dad_ys_birth[j]-est.year)), dad_n_comps[j]) # Sex-specific CKMR model equation
      }
    }

    # Write model    
    jags_file = paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/models/HS_neutralGrowth_testPriors.txt")
    write_model(HS_model, jags_file)
    
    
    #------------ STEP 3: SPECIFY INITIAL VALUES ---------------#
    jags_inits = function(nc) {
      inits = list()
      for(c in 1:nc){
        inits[[c]] = list(
          surv = 0.8,
          Nf = rnorm(1, mean = 500, sd = 100),
          Nm = rnorm(1, mean = 500, sd = 100),
          lam = 1
        )
      }
      return(inits)
    }
    
    #------------ STEP 4: SET NODES TO MONITOR ---------------#
    jags_params = c("Nf", "Nm", "surv", "lam")
    n_params = length(jags_params) #used to autofill dataframe later
    
    
    #------------- STEP 5: SET MCMC DIMENSIONS ---------------#
    jags_dims = c(
      ni = ni,  # number of post-burn-in samples per chain
      nb = nb,  # number of burn-in samples
      nt = nt,     # thinning rate
      nc = nc      # number of chains
    )
    
    MCMC.settings <- paste0("thin", jags_dims[names(jags_dims) == "nt"], "_draw", jags_dims[names(jags_dims) == "ni"], "_burn", jags_dims[names(jags_dims) == "nb"])
    
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
                              parallel = T
    )
    
    #---------------- STEP 7: CONVERGENCE DIAGNOSTICS -----------------#
    # view convergence diagnostic summaries for all monitored nodes
    # 2.5, 50, and 97.5 are quantiles in model.summary
    model.summary <- data.frame(t(post_summ(post, jags_params, Rhat = T, neff = T))) %>% 
      rownames_to_column(var = "parameter")
    
    #Calculate HPD intervals - 95%
    post.95 <- combine.mcmc(post) %>% 
      HPDinterval() %>% 
      data.frame() %>% 
      rownames_to_column(var = "parameter")
    post.95 <- post.95 %>% filter(parameter %in% jags_params) #Remove deviance
      
    
    #Combine into data.frame
    model.summary2 <- model.summary %>% left_join(post.95, by = "parameter") %>% 
      rename(HPD2.5 = lower, HPD97.5 = upper) %>% 
      dplyr::select(parameter, Q2.5 = X2.5., Q97.5 = X97.5., Q50 = X50., mean = mean, sd = sd, HPD2.5, HPD97.5, Rhat, neff)

    
    #---------------Plot prior vs posterior-----------------
    N.prior <- rnorm(n = 10000, mean = runif(1, min = 1, max = 10000), sd = runif(1, 1, sd.max)) %>% 
      as_tibble() %>% 
      mutate(Chain = "prior")
    
    surv.prior <- rbeta(n = 10000, shape1 = 1, shape2 = 2) %>% 
      as_tibble() %>% 
      mutate(Chain = "prior")
    
    lam.prior <- rnorm(n = 10000, mean = 1, sd = lam.sd) %>% 
      as_tibble() %>% 
      mutate(Chain = "prior")
    
    #Prepare posteriors for plotting
    Nf.post <- ggs(post) %>% 
      filter(Parameter == "Nf") %>% 
      mutate(Chain = ifelse(Chain == 1, "Chain 1", "Chain 2"))
    
    Nm.post <- ggs(post) %>% 
      filter(Parameter == "Nm") %>% 
      mutate(Chain = ifelse(Chain == 1, "Chain 1", "Chain 2"))
    
    surv.post <- ggs(post) %>% 
      filter(Parameter == "surv") %>% 
      mutate(Chain = ifelse(Chain == 1, "Chain 1", "Chain 2"))
    
    lam.post <- ggs(post) %>% 
      filter(Parameter == "lam") %>% 
      mutate(Chain = ifelse(Chain == 1, "Chain 1", "Chain 2"))
    
    Nf.density <- ggs_density(Nf.post) + 
      geom_density(data = N.prior, aes(x = value), alpha = 0.6) + 
      theme(legend.title = element_blank())
    
    Nm.density <- ggs_density(Nm.post) + 
      geom_density(data = N.prior, aes(x = value), alpha = 0.6) +
      theme(legend.title = element_blank())
    
    surv.density <- ggs_density(surv.post) + 
      geom_density(data = surv.prior, aes(x = value), alpha = 0.6) + 
      theme(legend.title = element_blank())
    
    lam.density <- ggs_density(lam.post) + 
      geom_density(data = lam.prior, aes(x = value), alpha = 0.6) + 
      theme(legend.title = element_blank())
    
    priorPostPlots <- ggarrange(Nf.density, Nm.density, surv.density, lam.density, common.legend = TRUE)
    
    ppp <- annotate_figure(priorPostPlots, top = text_grob(paste0("Density Plot: SD max = ", sd.max), face = "bold", size = 14))

ggsave(filename = paste0(purpose, "_", sd.max,"_", date.of.simulation, ".pdf"))
    }
    