#-------------- STEP 1: PREPARE DATA ----------------
# yrs <- c(estimation.year:n_yrs)
# ref.year <- min(mom_comps.all$ref.year, dad_comps.all$ref.year)

#Create vectors of data for JAGS
#Mom
#HS
#mom_comps.HS <- mom_comps.all %>% dplyr::filter(type == "HS")
mom.oncycle <- mom_comps.HS$mom.oncycle
mom.mort.yrs_HS <- mom_comps.HS$year_gap
mom.popGrowth.yrs_HS <- mom_comps.HS$pop.growth.yrs
mom.n.comps_HS <- mom_comps.HS$all
mom.positives_HS <- mom_comps.HS$yes
mom.yrs_HS <- nrow(mom_comps.HS)


#Dad
dad.mort.yrs_HS <- dad_comps.HS$year_gap
dad.popGrowth.yrs_HS <- dad_comps.HS$pop.growth.yrs
dad.n.comps_HS <- dad_comps.HS$all
dad.positives_HS <- dad_comps.HS$yes
dad.yrs_HS <- nrow(dad_comps.HS)
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
    mom.mort.yrs_HS = mom.mort.yrs_HS,
    mom.oncycle = mom.oncycle,
    mom.popGrowth.yrs_HS = mom.popGrowth.yrs_HS,
    mom.n.comps_HS = mom.n.comps_HS,
    mom.positives_HS = mom.positives_HS,
    mom.yrs_HS = mom.yrs_HS,

    
    #Dad
    dad.mort.yrs_HS = dad.mort.yrs_HS,
    dad.popGrowth.yrs_HS = dad.popGrowth.yrs_HS,
    dad.n.comps_HS = dad.n.comps_HS,
    dad.positives_HS = dad.positives_HS,
    dad.yrs_HS = dad.yrs_HS,
    
    #Lambda
    # lambda.prior.mean = lambda.prior.mean,
    # lam.tau = lam.tau,

    #survival
     # surv.alpha = surv.alpha,
     # surv.beta = surv.beta,
    
    #In case I want to use the truncated normal prior
#    adult.survival = adult.survival, 
#    survival.prior.sd = survival.prior.sd,

    # #Breeding interval
    #psi = psi.truth,
    a = mating.periodicity
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
        lambda = 1,
        psi = runif(1, min=0, max=1)
        
      )
    }
    return(inits)
  }
  
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
post = jagsUI::jags(data = jags_data, #If using postpack from AFS workshop
                          
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
model.summary <- data.frame(t(post_summ(post$samples, jags_params, Rhat = T, neff = T))) %>% 
  rownames_to_column(var = "parameter")

#Calculate HPD intervals - 95%
post.95 <- combine.mcmc(post$samples) %>% 
  HPDinterval() %>% 
  data.frame() %>% 
  rownames_to_column(var = "parameter")
post.95 <- post.95 %>% filter(parameter %in% jags_params) #Remove deviance


#Combine into data.frame
model.summary2 <- model.summary %>% left_join(post.95, by = "parameter") %>% 
  rename(HPD2.5 = lower, HPD97.5 = upper) %>% 
  dplyr::select(parameter, Q2.5 = X2.5., Q97.5 = X97.5., Q50 = X50., mean = mean, sd = sd, HPD2.5, HPD97.5, Rhat, neff)
