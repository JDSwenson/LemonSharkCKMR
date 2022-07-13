#-------------- STEP 1: PREPARE DATA ----------------
yrs <- c(estimation.year:n_yrs)
ref.year.mom <- min(mom_comps.all$ref.year)
ref.year.dad <- min(dad_comps.all$ref.year)

#Create vectors of data for JAGS
#Mom
#HS
mom.mort.yrs <- mom_comps.all$mort.yrs
mom.popGrowth.yrs <- mom_comps.all$pop.growth.yrs
mom.n.comps <- mom_comps.all$all
mom.positives <- mom_comps.all$yes
mom.yrs <- nrow(mom_comps.all)

#Mom
#PO
# mom_comps.PO <- mom_comps.all %>% dplyr::filter(type == "PO")
# 
# mom.mort.yrs_PO <- mom_comps.PO$mort.yrs
# mom.popGrowth.yrs_PO <- mom_comps.PO$pop.growth.yrs
# mom.n.comps_PO <- mom_comps.PO$all
# mom.positives_PO <- mom_comps.PO$yes
# mom.yrs_PO <- nrow(mom_comps.PO)

#Dad
dad.mort.yrs <- dad_comps.all$mort.yrs
dad.popGrowth.yrs <- dad_comps.all$pop.growth.yrs
dad.n.comps <- dad_comps.all$all
dad.positives <- dad_comps.all$yes
dad.yrs <- nrow(dad_comps.all)
#dad.R0 <- dad_comps.all$R0

#Set mean and sd (precision) for lambda
#lam.tau <- 1/(lambda.prior.sd^2) #Value derived from Leslie matrix

#Calculate parameters for beta distribution from mean and variance for survival
 # surv.betaParams <- estBetaParams(survival.prior.mean, survival.prior.sd^2)
 # surv.alpha <- surv.betaParams[[1]]
 # surv.beta <- surv.betaParams[[2]]


  #Define data
  jags_data = list(
    #Mom
    #HS
    mom.mort.yrs = mom.mort.yrs,
    mom.popGrowth.yrs = mom.popGrowth.yrs,
    mom.n.comps = mom.n.comps,
    mom.positives = mom.positives,
    mom.yrs = mom.yrs,
    #mom.R0 = mom.R0,
    
    
    #Dad
    dad.mort.yrs = dad.mort.yrs,
    dad.popGrowth.yrs = dad.popGrowth.yrs,
    dad.n.comps = dad.n.comps,
    dad.positives = dad.positives,
    dad.yrs = dad.yrs
    #dad.R0 = dad.R0,
    

    #survival
     # surv.alpha = surv.alpha,
     # surv.beta = surv.beta
    
      )
  
  #------------ STEP 2: SPECIFY INITIAL VALUES ---------------#
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        #If estimating all parameters
        survival = runif(1, min=0.5, max=0.95),
        Nf = rnorm(1, mean = 500, sd = 100),
        Nm = rnorm(1, mean = 500, sd = 100),
        lambda = 1
      )
    }
    return(inits)
  }
  
  #----------------- STEP 3: SPECIFY JAGS MODEL CODE ---------------#
  #Convert tau to SD (for interpretation)
  #tau <- 1E-6
  #(sd <- sqrt(1/tau))
  #tau <- 
  
  HS.only_model = function(){
    
    #PRIORS - uninformative
    mu ~ dunif(1, 10000)
    sd ~ dunif(1, 10000)
    Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
    Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
    survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival
    lambda ~ dunif(0.95, 1.05) #Diffuse prior for lambda
    
    #PRIORS - informative
#    survival ~ dbeta(surv.alpha, surv.beta);T(0.001,0.999) #Informative prior
#    surv ~ dnorm(adult.survival, 1/(survival.prior.sd^2));T(0.5, 0.99) #Informative prior
    
    
    #Likelihood
    #Moms
    #HS + PO
    for(i in 1:mom.yrs){ # Loop over maternal cohort comparisons
      mom.positives[i] ~ dbin((survival^mom.mort.yrs[i])/(Nf*(lambda^mom.popGrowth.yrs[i])), mom.n.comps[i]) # Sex-specific CKMR model equation
    }
    
    #Dads
    #HS + PO
    for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
      dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
    }
  }
  
  # Write model
  jags_file = paste0(jags.model_location, purpose, "_iteration_", iter, ".txt")
  write_model(HS.only_model, jags_file)
  
  
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


if(s == 1){
  sims.list.1[[iter]] <- post
 } else if(s == 2){
  sims.list.2[[iter]] <- post
} else if(s == 3){
  sims.list.3[[iter]] <- post
} else if(s == 4){
  sims.list.4[[iter]] <- post
}

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
  dplyr::select(parameter, Q2.5 = X2.5., Q97.5 = X97.5., Q50 = X50., mean = mean, sd = sd, HPD2.5, HPD97.5, Rhat, neff) %>% 
  mutate(iteration = iter, seed = rseed)
