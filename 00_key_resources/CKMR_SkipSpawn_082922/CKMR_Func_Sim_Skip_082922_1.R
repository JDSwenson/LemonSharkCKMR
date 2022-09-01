CKMR.sim.pop.skipspawn = function(Ntarget, Prop_f, phi, mat, RunYears, StartYear = 0, a, psi)
{
  
  Prop_m = (1 - Prop_f) #Get prop of males
  maxage = which(phi^seq(1000) <= 0.01)[1] #set the maximum age for calculating init age props to where <= 1% of inds in cohort are left
  init.age.props = (phi^seq(0,maxage))/sum((phi^seq(0,maxage))) #Calc init age props 
  
  #Set up the array and first year
  ind.mat = matrix(NA, nrow = RunYears, ncol = Ntarget*2)
  #Set IDs that hold all the necessary info
  #This string manipulation is a terribly inefficient way of doing things but it's easy to think about
  #Will translate to a 3D array at some point..
  ind.mat[1,1:Ntarget] = paste(seq(Ntarget),
                               StartYear - sample(x = seq(0,maxage),size = Ntarget, prob = init.age.props, replace = T),
                               sample(x = c("F","M"),size = Ntarget, prob = c(Prop_f, Prop_m), replace = T),
                               sample(x = c("A", "S"),size = Ntarget, prob = c((1-psi), psi), replace = T), #Do they skip or are they annual?
                               "NA", "NA",
                               sep = "_")
  
  for(i in seq(2,RunYears))
  {
    #Report
    cat(paste0("Simulating year ",i,"\r"))
    #Get the pop for the previous year
    NprevYear = sum(!is.na(ind.mat[(i-1),]))
    #Kill some
    lived = rbinom(n = NprevYear,size = 1, prob = (phi))
    #Put the survivors into the start of this year
    ind.mat[i,(1:sum(lived))] = ind.mat[(i-1),][which(lived == 1)]
    #Figure out roughly how many new YOY we need for stability, then get the actual number
    births.target = Ntarget - sum(lived)
    births.actual = rpois(n = 1, lambda = births.target)
    #Get pools of potential parents, add skipped spawning for females
    potmoms = ind.mat[i,][which(sapply(strsplit(ind.mat[i,], "_"), `[`, 3) == "F" & 
                                  (StartYear + i) - as.numeric(sapply(strsplit(ind.mat[i,], "_"), `[`, 2)) >= mat)]
    potmoms.ann = potmoms[which((sapply(strsplit(potmoms, "_"), `[`, 4) == "A"))] #Annual spawners
    potmoms.skip = potmoms[which((sapply(strsplit(potmoms, "_"), `[`, 4) == "S"))] #Skip spawners
    potmoms.skip = potmoms.skip[which(((StartYear + i) - as.numeric(sapply(strsplit(potmoms.skip, "_"), `[`, 2)) / a) == 
                                        round((StartYear + i) - as.numeric(sapply(strsplit(potmoms.skip, "_"), `[`, 2)) / a))]
    potmoms = c(potmoms.ann, potmoms.skip)
    
    potdads = ind.mat[i,][which(sapply(strsplit(ind.mat[i,], "_"), `[`, 3) == "M" & 
                                  (StartYear + i) - as.numeric(sapply(strsplit(ind.mat[i,], "_"), `[`, 2)) >= mat)]
    
    #Assign parents to each actual birth
    moms = sample(potmoms, size = births.actual, replace = T)
    moms = sapply(strsplit(moms, "_"), `[`, 1)
    dads = sample(potdads, size = births.actual, replace = T)
    dads = sapply(strsplit(dads, "_"), `[`, 1)
    #Get new IDs
    IDstart = max(as.numeric(sapply(strsplit(ind.mat[(i-1),], "_"), `[`, 1)), na.rm = T)+1
    newIDs = paste(seq(IDstart, IDstart + births.actual),
                   StartYear + i,
                   sample(x = c("F","M"),size = births.actual, prob = c(Prop_f, Prop_m), replace = T),
                   sample(x = c("A", "S"),size = births.actual, prob = c((1-psi), psi), replace = T),
                   moms, dads,
                   sep = "_")
    #Add back to the matrix
    ind.mat[i,((sum(lived)+1):(sum(lived)+births.actual+1))] = newIDs
  }
  
  return(ind.mat)
}

CKMR.sim.samp.skipspawn = function(ind.mat, Nsample, SampYear, nYears = 1)
{
  #Sample randomly from pop
  samp.vec = NULL
  for(i in seq(nYears))
  {
    samp.vec = c(samp.vec,
               ind.mat[dim(ind.mat)[1],sample(seq(sum(!is.na(ind.mat[(SampYear - (i - 1)),]))), size = Nsample, replace = F)])
  }
  #Extract info for each from ID
  samp.df = data.frame(ID = sapply(strsplit(samp.vec, "_"), `[`, 1),
                       Cohort = sapply(strsplit(samp.vec, "_"), `[`, 2),
                       Sex = sapply(strsplit(samp.vec, "_"), `[`, 3),
                       Mom = sapply(strsplit(samp.vec, "_"), `[`, 5),
                       Dad = sapply(strsplit(samp.vec, "_"), `[`, 6))
  #Create set of pairwise comparisons
  samp.pairwise = expand.grid(Ind1 = samp.df$ID, Ind2 = samp.df$ID)
  samp.pairwise = samp.pairwise[which(samp.pairwise$Ind1 != samp.pairwise$Ind2),] #Eliminate self-self comparisons
  #Add cohort data
  samp.pairwise$Cohort1 = as.numeric(samp.df$Cohort[match(samp.pairwise$Ind1, samp.df$ID)]) 
  samp.pairwise$Cohort2 = as.numeric(samp.df$Cohort[match(samp.pairwise$Ind2, samp.df$ID)])
  samp.pairwise = samp.pairwise[which(samp.pairwise$Cohort1 > samp.pairwise$Cohort2 ),] #Remove same-age and reversed combos
  #Add parental data
  samp.pairwise$Mom1 = samp.df$Mom[match(samp.pairwise$Ind1, samp.df$ID)] 
  samp.pairwise$Mom2 = samp.df$Mom[match(samp.pairwise$Ind2, samp.df$ID)] 
  samp.pairwise$Dad1 = samp.df$Dad[match(samp.pairwise$Ind1, samp.df$ID)]
  samp.pairwise$Dad2 = samp.df$Dad[match(samp.pairwise$Ind2, samp.df$ID)]
  #Identify HSPs
  samp.pairwise$isMHSP = samp.pairwise$Mom1 == samp.pairwise$Mom2 
  samp.pairwise$isPHSP = samp.pairwise$Dad1 == samp.pairwise$Dad2
  
  samp.pairwise$AgeDiff = samp.pairwise$Cohort1 - samp.pairwise$Cohort2
  
  return(samp.pairwise)
}

CKMR.est.skipspawn = function(samp.pairwise, psi.fix, a.fix, jags.model_location)
{
  library(tidyverse)
  library(rjags)
  library(jagsUI)
  library(runjags)
  library(postpack)
  
samp.pairwise.sum.dad = suppressMessages(samp.pairwise %>%
  group_by(AgeDiff, Cohort2) %>%
  summarize(nComps = n(),
            nPHSP = sum(isPHSP)))

  samp.pairwise$OnCycle = FALSE
  samp.pairwise[which(samp.pairwise$AgeDiff/a == 
                        round(samp.pairwise$AgeDiff/a)),]$OnCycle = TRUE
  
  samp.pairwise.sum.mom = suppressMessages(samp.pairwise %>%
    group_by(AgeDiff,  OnCycle) %>%
      #filter(!OnCycle) %>% #To isolate one part of the MHSP prob for testing 
    summarize(nComps = n(),
              nMHSP = sum(isMHSP)))
  
  jags_data = list(
    dad.mort.yrs_HS = samp.pairwise.sum.dad$AgeDiff,
    dad.n.comps_HS = samp.pairwise.sum.dad$nComps,
    dad.positives_HS = samp.pairwise.sum.dad$nPHSP,
    dad.datalen = dim(samp.pairwise.sum.dad)[1],
    
    mom.mort.yrs_HS = samp.pairwise.sum.mom$AgeDiff,
    mom.n.comps_HS = samp.pairwise.sum.mom$nComps,
    mom.positives_HS = samp.pairwise.sum.mom$nMHSP,
    mom.datalen = dim(samp.pairwise.sum.mom)[1],
    mom.on.spawn.cycle = as.numeric(samp.pairwise.sum.mom$OnCycle),
    #psi = psi.fix, #Turn on to fix psi
    a = a.fix
  )
  
  #------------ STEP 2: SPECIFY INITIAL VALUES ---------------#
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        #If estimating all parameters
        survival = runif(1, min=0.05, max=0.95),
        Nf = rnorm(1, mean = 500, sd = 100),
        psi = runif(1, min=0.05, max=0.95), #Turn off if fixing psi
        Nm = rnorm(1, mean = 500, sd = 100)
        #lambda = 1,
      )
    }
    return(inits)
  }
  
  #----------------- STEP 3: SPECIFY JAGS MODEL CODE ---------------#
  #Convert tau to SD (for interpretation)
  #tau <- 1E-6
  #(sd <- sqrt(1/tau))
  #tau <- 
  
  HS.only_Skip_model = function(){
    
    #PRIORS - uninformative
    mu ~ dunif(1, 10000)
    sd ~ dunif(1, 10000)
    Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
    Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
    survival ~ dunif(0, 1) # Uninformative prior for adult survival
    psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi


    #Ben's skip-spawn pre-calcs
    #agedist <- survival^(1:1000) / sum(survival^(1:1000))#Get stable age dist

    pNf.ann <- (1 - psi) #Number of annual spawners is easy
    pNf.skip <- psi * ((1 / (1-(survival^a))) / (1 / (1-survival)))
    
    #Likelihood
    for(i in 1:mom.datalen){ # Loop over maternal cohort comparisons
      
      #Skip-spawn naive likelihood
       # mom.positives_HS[i] ~ dbin((survival^mom.mort.yrs_HS[i]) / Nf,
       #                            mom.n.comps_HS[i])
      
      #Patterson et al. preprint likelihood (only for a=2)
       # mom.positives_HS[i] ~ dbin((((survival^mom.mort.yrs_HS[i]) / ((Nf * (2 - psi)) / 2))* mom.on.spawn.cycle[i]) + #For on-cycle spawning (annual + psi/a)
       #                              ((((1-psi)*survival^mom.mort.yrs_HS[i]) / ((Nf * (2 - psi)) / 2))* (1-mom.on.spawn.cycle[i])), #For off-cycle spawning (annual only)
       #                            mom.n.comps_HS[i])

      #Ben's skip-spawn likelihood
      mom.positives_HS[i] ~ dbin((((survival^mom.mort.yrs_HS[i]) / (Nf * pNf.ann + Nf * pNf.skip)) * mom.on.spawn.cycle[i]) + #For on-cycle spawning (annual + psi/a)
                                   ((((survival^mom.mort.yrs_HS[i]) * ((pNf.ann / (pNf.ann + pNf.skip)))) / (Nf * pNf.ann)) * (1 - mom.on.spawn.cycle[i])), #For off-cycle spawning (annual only)
                                 mom.n.comps_HS[i])
    }
    
    
    for(j in 1:dad.datalen){ # Loop over paternal cohort comparisons
      dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/Nm), dad.n.comps_HS[j])
    }
    
  }
  
  # Write model
  jags_file = paste0(jags.model_location,"Skiptest", ".txt")
  write_model(HS.only_Skip_model, jags_file)
  
  jags_params = c("Nf", "Nm", "survival", "psi")
  #jags_params = c("Nf", "Nm", "survival")
  n_params = length(jags_params) #used to autofill dataframe later
  
  
  #------------- STEP 5: SET MCMC DIMENSIONS ---------------#
  jags_dims = c(
    ni = 20000,  # number of post-burn-in samples per chain
    nb = 10000,  # number of burn-in samples
    nt = 10,     # thinning rate
    nc = 3      # number of chains
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
  
  
  #Combine into data.frame
  model.summary2 <- model.summary %>% left_join(post.95, by = "parameter") %>% 
    rename(HPD2.5 = lower, HPD97.5 = upper) %>% 
    dplyr::select(parameter, Q2.5 = X2.5., Q97.5 = X97.5., Q50 = X50., mean = mean, sd = sd, HPD2.5, HPD97.5, Rhat, neff)
  
  return(model.summary2)
}



