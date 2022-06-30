library(rjags)
library(jagsUI)
library(R2jags)
library(runjags)
library(postpack)

jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_Lemon.Shark.data/models/"
jags_params = c("Nfb", "psi", "Nm", "surv", "lam")

#----------------------- MCMC parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

#-------------- STEP 1: PREPARE DATA ----------------
#yrs <- c(estimation.year:n_yrs)

#Create vectors of data for JAGS
#Mom
#HS - even years
mom_comps.HS_even <- mom_comps.HS %>% dplyr::filter(BI == "even")

mom.year_gap_HS.even <- mom_comps.HS_even$year_gap
mom.popGrowth.yrs_HS.even <- mom_comps.HS_even$pop.growth.yrs
mom.n.comps_HS.even <- mom_comps.HS_even$all
mom.positives_HS.even <- mom_comps.HS_even$yes
mom.yrs_HS.even <- nrow(mom_comps.HS_even)
#mom.R0 <- mom_comps.all$R0

#Mom
#HS - odd years
mom_comps.HS_odd <- mom_comps.HS %>% dplyr::filter(BI == "odd")

mom.year_gap_HS.odd <- mom_comps.HS_odd$year_gap
mom.popGrowth.yrs_HS.odd <- mom_comps.HS_odd$pop.growth.yrs
mom.n.comps_HS.odd <- mom_comps.HS_odd$all
mom.positives_HS.odd <- mom_comps.HS_odd$yes
mom.yrs_HS.odd <- nrow(mom_comps.HS_odd)

#Mom
#PO
# mom_comps.PO <- mom_comps.all %>% dplyr::filter(type == "PO")
# 
# mom.year_gap_PO <- mom_comps.PO$year_gap
# mom.popGrowth.yrs_PO <- mom_comps.PO$pop.growth.yrs
# mom.n.comps_PO <- mom_comps.PO$all
# mom.positives_PO <- mom_comps.PO$yes
# mom.yrs_PO <- nrow(mom_comps.PO)

#Dad
dad.year_gap <- dad_comps.HS$year_gap
dad.popGrowth.yrs <- dad_comps.HS$pop.growth.yrs
dad.n.comps <- dad_comps.HS$all
dad.positives <- dad_comps.HS$yes
dad.yrs <- nrow(dad_comps.HS)
#dad.R0 <- dad_comps.all$R0

#Set mean and sd (precision) for lambda
lam.tau <- 1/(0.02277^2) #Value derived from Leslie matrix
#N.tau <- 1/(sd^2)



  #Define data
  jags_data = list(
    #Mom
    #HS: even years
    mom.year_gap_HS.even = mom.year_gap_HS.even,
    mom.popGrowth.yrs_HS.even = mom.popGrowth.yrs_HS.even,
    mom.n.comps_HS.even = mom.n.comps_HS.even,
    mom.positives_HS.even = mom.positives_HS.even,
    mom.yrs_HS.even = mom.yrs_HS.even,
    #mom.R0 = mom.R0,
    
    #Mom
    #HS: odd years
    mom.year_gap_HS.odd = mom.year_gap_HS.odd,
    mom.popGrowth.yrs_HS.odd = mom.popGrowth.yrs_HS.odd,
    mom.n.comps_HS.odd = mom.n.comps_HS.odd,
    mom.positives_HS.odd = mom.positives_HS.odd,
    mom.yrs_HS.odd = mom.yrs_HS.odd,
    #mom.R0 = mom.R0,
    
    
    #Dad
    dad.year_gap = dad.year_gap,
    dad.popGrowth.yrs = dad.popGrowth.yrs,
    dad.n.comps = dad.n.comps,
    dad.positives = dad.positives,
    dad.yrs = dad.yrs,
    #dad.R0 = dad.R0,
    

    #Lambda
    lam.tau = lam.tau,

    #survival
    Adult.survival = Adult.survival
      )
  
  #------------ STEP 2: SPECIFY INITIAL VALUES ---------------#
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        #If estimating all parameters
        surv = runif(1, min=.5, max=.95),
        #surv = .8,
        Nfb = rnorm(1, mean = 500, sd = 100),
        Nm = rnorm(1, mean = 500, sd = 100),
        lam = 1
        #psi = runif(1, min=0, max=1)
        
      )
    }
    return(inits)
  }
  
  #----------------- STEP 3: SPECIFY JAGS MODEL CODE ---------------#
  #Convert tau to SD (for interpretation)
  #tau <- 1E-6
  #(sd <- sqrt(1/tau))
  #tau <- 
  
  HS.only_LS.model = function(){
    
    #PRIORS - uninformative
    mu ~ dunif(1, 10000)
    sd ~ dunif(1, 10000)
    Nfb ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
    Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
    #surv ~ dbeta(1 ,1) # Uninformative prior for adult survival
    lam ~ dnorm(1, lam.tau)
    psi ~ dunif(0, 1) #Percent of animals breeding bi-ennially
    
    
    #PRIORS - informative
    surv ~ dnorm(Adult.survival, 1/(.000001)^2) #Informative prior
    
    #Likelihood
    #Moms
    #HS - even years
    for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
      mom.positives_HS.even[i] ~ dbin((surv^mom.year_gap_HS.even[i])/(Nfb*(lam^mom.popGrowth.yrs_HS.even[i])), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
    }
    
    #Moms
    #HS - odd years
    for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
      mom.positives_HS.odd[j] ~ dbin(((surv^mom.year_gap_HS.odd[j])*(1-psi))/(Nfb*(lam^mom.popGrowth.yrs_HS.odd[j])), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
    }
    
    #Dads
    #HS + PO
    for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
      dad.positives[f] ~ dbin((surv^dad.year_gap[f])/(Nm*(lam^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
    }
  }
  
  # Write model
  jags_file = paste0(jags.model_location, purpose, ".txt")
  write_model(HS.only_LS.model, jags_file)
  
  
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
