#Define data
jags_data = list(
  #Mom
  #HS: even years
  mom.mort.yrs_HS.on = mom.mort.yrs_HS.on,
  mom.popGrowth.yrs_HS.on = mom.popGrowth.yrs_HS.on,
  mom.n.comps_HS.on = mom.n.comps_HS.on,
  mom.positives_HS.on = mom.positives_HS.on,
  mom.yrs_HS.on = mom.yrs_HS.on,
  
  #Mom
  #HS: odd years
  mom.mort.yrs_HS.off = mom.mort.yrs_HS.off,
  mom.popGrowth.yrs_HS.off = mom.popGrowth.yrs_HS.off,
  mom.n.comps_HS.off = mom.n.comps_HS.off,
  mom.positives_HS.off = mom.positives_HS.off,
  mom.yrs_HS.off = mom.yrs_HS.off,
  
  
  #Dad
  dad.mort.yrs = dad.mort.yrs,
  dad.popGrowth.yrs = dad.popGrowth.yrs,
  dad.n.comps = dad.n.comps,
  dad.positives = dad.positives,
  dad.yrs = dad.yrs,
  
  
  # #Breeding interval
  a = mating.periodicity
)

cat("JAGS data defined for multiennial HS only model\n")

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

    sims.list.1[[iter]] <- post
    

#---------------- STEP 5: CONVERGENCE DIAGNOSTICS -----------------#
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
  dplyr::select(parameter, Q2.5 = X2.5., Q97.5 = X97.5., Q50 = X50., mean = mean, sd = sd, HPD2.5, HPD97.5, Rhat, neff) %>% 
  mutate(iteration = iter, seed = rseed)
