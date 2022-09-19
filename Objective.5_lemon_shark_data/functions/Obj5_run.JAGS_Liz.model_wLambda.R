#-------------- STEP 1: PREPARE DATA ----------------
#yrs <- c(estimation.year:n_yrs)
#ref.year <- min(mom_comps.all$ref.year, dad_comps.all$ref.year)

#Create vectors of data for JAGS
#Mom
#HS - even years
mom_comps.HS_even <- mom_comps.HS %>% dplyr::filter(BI == "even")

mom.mort.yrs_HS.even <- mom_comps.HS_even$year_gap
mom.popGrowth.yrs_HS.even <- mom_comps.HS_even$pop.growth.yrs
mom.n.comps_HS.even <- mom_comps.HS_even$all
mom.positives_HS.even <- mom_comps.HS_even$yes
mom.yrs_HS.even <- nrow(mom_comps.HS_even)
#mom.R0 <- mom_comps.all$R0

#Mom
#HS - odd years
mom_comps.HS_odd <- mom_comps.HS %>% dplyr::filter(BI == "odd")

mom.mort.yrs_HS.odd <- mom_comps.HS_odd$year_gap
mom.popGrowth.yrs_HS.odd <- mom_comps.HS_odd$pop.growth.yrs
mom.n.comps_HS.odd <- mom_comps.HS_odd$all
mom.positives_HS.odd <- mom_comps.HS_odd$yes
mom.yrs_HS.odd <- nrow(mom_comps.HS_odd)

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
dad.mort.yrs <- dad_comps.HS$year_gap
dad.popGrowth.yrs <- dad_comps.HS$pop.growth.yrs
dad.n.comps <- dad_comps.HS$all
dad.positives <- dad_comps.HS$yes
dad.yrs <- nrow(dad_comps.HS)
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
    #HS: even years
    mom.mort.yrs_HS.even = mom.mort.yrs_HS.even,
    mom.popGrowth.yrs_HS.even = mom.popGrowth.yrs_HS.even,
    mom.n.comps_HS.even = mom.n.comps_HS.even,
    mom.positives_HS.even = mom.positives_HS.even,
    mom.yrs_HS.even = mom.yrs_HS.even,
    #mom.R0 = mom.R0,
    
    #Mom
    #HS: odd years
    mom.mort.yrs_HS.odd = mom.mort.yrs_HS.odd,
    mom.popGrowth.yrs_HS.odd = mom.popGrowth.yrs_HS.odd,
    mom.n.comps_HS.odd = mom.n.comps_HS.odd,
    mom.positives_HS.odd = mom.positives_HS.odd,
    mom.yrs_HS.odd = mom.yrs_HS.odd,
    #mom.R0 = mom.R0,
    
    
    #Dad
    dad.mort.yrs = dad.mort.yrs,
    dad.popGrowth.yrs = dad.popGrowth.yrs,
    dad.n.comps = dad.n.comps,
    dad.positives = dad.positives,
    dad.yrs = dad.yrs,
    #dad.R0 = dad.R0,
    

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
        psi = runif(1, min=0.75, max=1.0),
        lambda = 1
        
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
