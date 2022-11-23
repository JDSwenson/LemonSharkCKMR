######################### STEP 1: PREPARE DATA ########################
#------------------------ Annual model ------------------------
if(model == "annual"){
yrs <- c(estimation.year:n_yrs)
ref.year <- min(mom_comps.all$ref.year, dad_comps.all$ref.year)

#Create vectors of data for JAGS
#Mom
#HS - even years
mom_comps.HS_even <- mom_comps.all %>% dplyr::filter(type == "HS", BI == "even")
mom.mort.yrs_HS.even <- mom_comps.HS_even$mort.yrs
mom.popGrowth.yrs_HS.even <- mom_comps.HS_even$pop.growth.yrs
mom.n.comps_HS.even <- mom_comps.HS_even$all
mom.positives_HS.even <- mom_comps.HS_even$yes
mom.yrs_HS.even <- nrow(mom_comps.HS_even)
#mom.R0 <- mom_comps.all$R0

#Mom
#HS - odd years
mom_comps.HS_odd <- mom_comps.all %>% dplyr::filter(type == "HS", BI == "odd")
mom.mort.yrs_HS.odd <- mom_comps.HS_odd$mort.yrs
mom.popGrowth.yrs_HS.odd <- mom_comps.HS_odd$pop.growth.yrs
mom.n.comps_HS.odd <- mom_comps.HS_odd$all
mom.positives_HS.odd <- mom_comps.HS_odd$yes
mom.yrs_HS.odd <- nrow(mom_comps.HS_odd)


#Dad
dad.mort.yrs <- dad_comps.all$mort.yrs
dad.popGrowth.yrs <- dad_comps.all$pop.growth.yrs
dad.n.comps <- dad_comps.all$all
dad.positives <- dad_comps.all$yes
dad.yrs <- nrow(dad_comps.all)

#Calculate parameters for beta distribution from mean and variance for survival
 surv.betaParams <- estBetaParams(survival.prior.mean, survival.prior.sd^2)
 surv.alpha <- surv.betaParams[[1]]
 surv.beta <- surv.betaParams[[2]]
}


#------------------------ Multiennial model ------------------------
if(model == "multiennial"){
    #Create vectors of data for JAGS
    #Mom
    #HS - even years
    mom_comps.HS_even <- mom_comps.all %>% dplyr::filter(type == "HS", BI == "even")
    mom.mort.yrs_HS.even <- mom_comps.HS_even$mort.yrs
    mom.popGrowth.yrs_HS.even <- mom_comps.HS_even$pop.growth.yrs
    mom.n.comps_HS.even <- mom_comps.HS_even$all
    mom.positives_HS.even <- mom_comps.HS_even$yes
    mom.yrs_HS.even <- nrow(mom_comps.HS_even)
    
    #Mom
    #HS - odd years
    mom_comps.HS_odd <- mom_comps.all %>% dplyr::filter(type == "HS", BI == "odd")
    mom.mort.yrs_HS.odd <- mom_comps.HS_odd$mort.yrs
    mom.popGrowth.yrs_HS.odd <- mom_comps.HS_odd$pop.growth.yrs
    mom.n.comps_HS.odd <- mom_comps.HS_odd$all
    mom.positives_HS.odd <- mom_comps.HS_odd$yes
    mom.yrs_HS.odd <- nrow(mom_comps.HS_odd)
    
    #Dad
    dad.mort.yrs <- dad_comps.all$mort.yrs
    dad.popGrowth.yrs <- dad_comps.all$pop.growth.yrs
    dad.n.comps <- dad_comps.all$all
    dad.positives <- dad_comps.all$yes
    dad.yrs <- nrow(dad_comps.all)
    
    #If using a multiennial model that includes parent-offspring pairs, need to define those comparisons as well.
    if(HS.only != "yes"){ 
    #Mom
    #PO
    mom_comps.PO <- mom_comps.all %>% dplyr::filter(type == "PO")
    mom.mort.yrs_PO <- mom_comps.PO$mort.yrs
    mom.popGrowth.yrs_PO <- mom_comps.PO$pop.growth.yrs
    mom.n.comps_PO <- mom_comps.PO$all
    mom.positives_PO <- mom_comps.PO$yes
    mom.yrs_PO <- nrow(mom_comps.PO)
  }
}
 


######################### STEP 2: DEFINE DATA FOR JAGS ########################
#------------------------- Annual model; no lambda -------------------------
 #========================= Model validation =========================#
#Informed prior on survival#
if(model == "annual"){
  if(jags_file == paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation")){
    #Define data
    jags_data = list(
      #Mom
      mom.mort.yrs = mom.mort.yrs,
      mom.popGrowth.yrs = mom.popGrowth.yrs,
      mom.n.comps = mom.n.comps,
      mom.positives = mom.positives,
      mom.yrs = mom.yrs,
      
      #Dad
      dad.mort.yrs = dad.mort.yrs,
      dad.popGrowth.yrs = dad.popGrowth.yrs,
      dad.n.comps = dad.n.comps,
      dad.positives = dad.positives,
      dad.yrs = dad.yrs,
      
      #survival
      surv.alpha = surv.alpha,
      surv.beta = surv.beta
    
  )
}

 #------------------------- Annual model; Diffuse prior on survival -------------------------
   if(jags_file == paste0(jags.model_location, "HS.PO_noLambda_annual_model.txt") |
      jags_file == paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt") |
      jags_file == paste0(jags.model_location, "HS.PO_wideLambda_annual_model.txt")){
     
     #Define data
     jags_data = list(
       #Mom
       mom.mort.yrs = mom.mort.yrs,
       mom.popGrowth.yrs = mom.popGrowth.yrs,
       mom.n.comps = mom.n.comps,
       mom.positives = mom.positives,
       mom.yrs = mom.yrs,
       
       #Dad
       dad.mort.yrs = dad.mort.yrs,
       dad.popGrowth.yrs = dad.popGrowth.yrs,
       dad.n.comps = dad.n.comps,
       dad.positives = dad.positives,
       dad.yrs = dad.yrs,
   )
 }
}

#------------------------- Multiennial model -------------------------

if(model == "multiennial" & HS.only == "yes"){
  #Define data
  jags_data = list(
    #Mom
    #HS: even years
    mom.mort.yrs_HS.even = mom.mort.yrs_HS.even,
    mom.popGrowth.yrs_HS.even = mom.popGrowth.yrs_HS.even,
    mom.n.comps_HS.even = mom.n.comps_HS.even,
    mom.positives_HS.even = mom.positives_HS.even,
    mom.yrs_HS.even = mom.yrs_HS.even,
    
    #Mom
    #HS: odd years
    mom.mort.yrs_HS.odd = mom.mort.yrs_HS.odd,
    mom.popGrowth.yrs_HS.odd = mom.popGrowth.yrs_HS.odd,
    mom.n.comps_HS.odd = mom.n.comps_HS.odd,
    mom.positives_HS.odd = mom.positives_HS.odd,
    mom.yrs_HS.odd = mom.yrs_HS.odd,
    
    
    #Dad
    dad.mort.yrs = dad.mort.yrs,
    dad.popGrowth.yrs = dad.popGrowth.yrs,
    dad.n.comps = dad.n.comps,
    dad.positives = dad.positives,
    dad.yrs = dad.yrs,
    
    
    # #Breeding interval
    a = mating.periodicity
  )
  
}



if(model == "multiennial" & HS.only != "yes"){
  #Define data
  jags_data = list(
    #Mom
    #HS: even years
    mom.mort.yrs_HS.even = mom.mort.yrs_HS.even,
    mom.popGrowth.yrs_HS.even = mom.popGrowth.yrs_HS.even,
    mom.n.comps_HS.even = mom.n.comps_HS.even,
    mom.positives_HS.even = mom.positives_HS.even,
    mom.yrs_HS.even = mom.yrs_HS.even,
    
    #Mom
    #HS: odd years
    mom.mort.yrs_HS.odd = mom.mort.yrs_HS.odd,
    mom.popGrowth.yrs_HS.odd = mom.popGrowth.yrs_HS.odd,
    mom.n.comps_HS.odd = mom.n.comps_HS.odd,
    mom.positives_HS.odd = mom.positives_HS.odd,
    mom.yrs_HS.odd = mom.yrs_HS.odd,
    
    #Mom
    #PO
    mom.mort.yrs_PO = mom.mort.yrs_PO,
    mom.popGrowth.yrs_PO = mom.popGrowth.yrs_PO,
    mom.n.comps_PO = mom.n.comps_PO,
    mom.positives_PO = mom.positives_PO,
    mom.yrs_PO = mom.yrs_PO,
    
    
    #Dad
    dad.mort.yrs = dad.mort.yrs,
    dad.popGrowth.yrs = dad.popGrowth.yrs,
    dad.n.comps = dad.n.comps,
    dad.positives = dad.positives,
    dad.yrs = dad.yrs,
    
    # #Breeding interval
    a = mating.periodicity
  )   
}


######################### STEP 3: SPECIFY INITIAL VALUES #########################
#------------------------- Annual model -------------------------
if(model == "annual"){
  #If no lambda, then do not specify a starting value for this parameter (which isn't in the model)
  if(jags_file == paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation") | jags_file == paste0(jags.model_location, "HS.PO_noLambda_annual_model.txt")){
    jags_inits = function(nc) {
      inits = list()
      for(c in 1:nc){
       inits[[c]] = list(
         survival = runif(1, min=0.5, max=0.95),
         Nf = rnorm(1, mean = 500, sd = 100),
         Nm = rnorm(1, mean = 500, sd = 100)
     )
   }
   return(inits)
   }
  }
  #If lambda is in the model, then specify a starting value for this parameter
  if(jags_file == paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt") | jags_file == paste0(jags.model_location, "HS.PO_wideLambda_annual_model.txt")){
    jags_inits = function(nc) {
      inits = list()
      for(c in 1:nc){
        inits[[c]] = list(
          survival = runif(1, min=0.5, max=0.95),
          Nf = rnorm(1, mean = 500, sd = 100),
          Nm = rnorm(1, mean = 500, sd = 100),
          lambda = 1
        )
      }
      return(inits)
    }
  }
}
 
 
#------------------------- Multiennial model -------------------------
if(model == "multiennial"){
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        survival = runif(1, min=0.5, max=0.95),
        Nf = rnorm(1, mean = 500, sd = 100),
        Nm = rnorm(1, mean = 500, sd = 100),
        psi = runif(1, min=0, max=1),
        lambda = 1
      )
    }
    return(inits)
  }
}
 
 
#------------- STEP 4: SET MCMC DIMENSIONS ---------------#
n_params = length(jags_params) #used to autofill dataframe later

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


if(s == 1){
  sims.list.1[[iter]] <- post
 } else if(s == 2){
  sims.list.2[[iter]] <- post
} else if(s == 3){
  sims.list.3[[iter]] <- post
} else if(s == 4){
  sims.list.4[[iter]] <- post
}

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