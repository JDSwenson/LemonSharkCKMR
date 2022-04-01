#-------------- STEP 1: PREPARE DATA ----------------#
#Create vectors of data for JAGS
#Mom
mom.mort.yrs <- mom_comps.all$mort.yrs
mom.popGrowth.yrs <- mom_comps.all$pop.growth.yrs
mom.n.comps <- mom_comps.all$all
mom.positives <- mom_comps.all$yes
mom.yrs <- nrow(mom_comps.all)
#mom.R0 <- mom_comps.all$R0

#Dad
dad.mort.yrs <- dad_comps.all$mort.yrs
dad.popGrowth.yrs <- dad_comps.all$pop.growth.yrs
dad.n.comps <- dad_comps.all$all
dad.positives <- dad_comps.all$yes
dad.yrs <- nrow(dad_comps.all)
#dad.R0 <- dad_comps.all$R0

#Set mean and sd (precision) for lambda
lam.tau <- 1/(0.02277^2) #Value derived from Leslie matrix
#N.tau <- 1/(sd^2)

#Define data
jags_data = list(
  #Mom
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
  dad.yrs = dad.yrs,
  #dad.R0 = dad.R0,
  
  #estimation.year = estimation.year, # estimation year i.e. year the estimate will be focused on
  #N.tau = 1E-6,
  lam.tau = lam.tau
)


#----------------- STEP 2: SPECIFY JAGS MODEL CODE ---------------#
#Convert tau to SD (for interpretation)
#tau <- 1E-6
#(sd <- sqrt(1/tau))
#tau <- 

HS.PO_model = function(){
  
  #PRIORS
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  surv ~ dnorm(0.825, 1/(0.005^2)) # Informative prior for adult survival
  #surv ~ dnorm(Adult.survival, 1/(.02)^2) #Informative prior
  lam ~ dnorm(1, lam.tau)
  
  #Likelihood
  #Moms
  for(i in 1:mom.yrs){ # Loop over maternal cohort comparisons
    mom.positives[i] ~ dbin((surv^mom.mort.yrs[i])/(Nf*(lam^mom.popGrowth.yrs[i])), mom.n.comps[i]) # Sex-specific CKMR model equation
  }

  #Dads
    for(j in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[j] ~ dbin((surv^dad.mort.yrs[j])/(Nm*(lam^dad.popGrowth.yrs[j])), dad.n.comps[j]) # Sex-specific CKMR model equation
  }
}

# Write model
jags_file = paste0(jags.model_location, purpose, "_iteration_", iter, ".txt")
write_model(HS.PO_model, jags_file)


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


if(samps == 1){
  sims.list.1[[iter]] <- post
} else if(samps == 2){
  sims.list.2[[iter]] <- post
} else if(samps == 3){
  sims.list.3[[iter]] <- post
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
  dplyr::select(parameter, Q2.5 = X2.5., Q97.5 = X97.5., Q50 = X50., mean = mean, sd = sd, HPD2.5, HPD97.5, Rhat, neff)