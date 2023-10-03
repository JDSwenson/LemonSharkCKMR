######################### STEP 1: PREPARE DATA ########################
#------------------------ Annual model ------------------------
if(model == "annual.model"){
yrs <- c(estimation.year:n_yrs)
ref.year.mom <- min(mom_comps.all$ref.year)
ref.year.dad <- min(dad_comps.all$ref.year)

#Create vectors of data for JAGS
#Mom
mom.mort.yrs <- mom_comps.all$mort.yrs
mom.popGrowth.yrs <- mom_comps.all$pop.growth.yrs
mom.n.comps <- mom_comps.all$all
mom.positives <- mom_comps.all$yes
mom.yrs <- nrow(mom_comps.all)


#Dad
dad.mort.yrs <- dad_comps.all$mort.yrs
dad.popGrowth.yrs <- dad_comps.all$pop.growth.yrs
dad.n.comps <- dad_comps.all$all
dad.positives <- dad_comps.all$yes
dad.yrs <- nrow(dad_comps.all)

}


#------------------------ Multiennial model ------------------------
if(model == "multiennial.model"){
  yrs <- c(estimation.year:n_yrs)
  ref.year.mom <- min(mom_comps.all$ref.year)
  ref.year.dad <- min(dad_comps.all$ref.year)
  
    #Create vectors of data for JAGS
    #Mom
    #HS - even years
    mom_comps.HS_on <- mom_comps.all %>% dplyr::filter(type == "HS", BI == "on")
    mom.mort.yrs_HS.on <- mom_comps.HS_on$mort.yrs
    mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
    mom.n.comps_HS.on <- mom_comps.HS_on$all
    mom.positives_HS.on <- mom_comps.HS_on$yes
    mom.yrs_HS.on <- nrow(mom_comps.HS_on)
    
    #Mom
    #HS - off years
    mom_comps.HS_off <- mom_comps.all %>% dplyr::filter(type == "HS", BI == "off")
    mom.mort.yrs_HS.off <- mom_comps.HS_off$mort.yrs
    mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
    mom.n.comps_HS.off <- mom_comps.HS_off$all
    mom.positives_HS.off <- mom_comps.HS_off$yes
    mom.yrs_HS.off <- nrow(mom_comps.HS_off)
    
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
 #========================= Model validation =========================#
if(model == "annual.model"){
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
       
       #general
       estimation.year = estimation.year,
       est.year.calibrate = est.year.calibrate
   )
 }

#------------------------- Multiennial model -------------------------
if(model == "multiennial.model" & HS.only == "yes"){
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
    
    
    #general
    estimation.year = estimation.year,
    est.year.calibrate = est.year.calibrate,
    a = mating.periodicity
  )
  
}



if(model == "multiennial.model" & HS.only != "yes"){
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
    
    #general
    estimation.year = estimation.year,
    est.year.calibrate = est.year.calibrate,
    a = mating.periodicity
  )   
}


######################### STEP 3: SPECIFY INITIAL VALUES #########################
#------------------------- Annual model -------------------------
if(model == "annual.model"){
  #If lambda is in the model, then specify a starting value for this parameter
    jags_inits = function(nc) {
      inits = list()
      for(c in 1:nc){
        inits[[c]] = list(
          survival = runif(1, min=0.5, max=0.95),
          Nf0 = rnorm(1, mean = 500, sd = 100),
          Nm0 = rnorm(1, mean = 500, sd = 100),
          lambda = 1
        )
      }
      return(inits)
    }
  }
 
#------------------------- Multiennial model -------------------------
if(model == "multiennial.model"){
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        survival = runif(1, min=0.5, max=0.95),
        Nf0 = rnorm(1, mean = 500, sd = 100),
        Nm0 = rnorm(1, mean = 500, sd = 100),
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

# if(exists("s") == TRUE){
# if(s == 1){
#   sims.list.1[[iter]] <- post
#  } else if(s == 2){
#   sims.list.2[[iter]] <- post
# } else if(s == 3){
#   sims.list.3[[iter]] <- post
# } else if(s == 4){
#   sims.list.4[[iter]] <- post
# }
# }

#Save results from different estimation years.
if(est == 1){
  sims.list.1[[iter]] <- post
} else if(est == 2){
  sims.list.2[[iter]] <- post
} else if(est == 3){
  sims.list.3[[iter]] <- post
} else if(est == 4){
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