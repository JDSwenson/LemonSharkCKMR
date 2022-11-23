HS.PO_noLambda_annual_model_validation = function(){
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  
  #PRIORS - informative
  survival ~ dbeta(surv.alpha, surv.beta);T(0.001,0.999) 
  
  #Likelihood
  #Moms
  #HS + PO
  for(i in 1:mom.yrs){ # Loop over maternal cohort comparisons
    mom.positives[i] ~ dbin((survival^mom.mort.yrs[i])/Nf, mom.n.comps[i]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/Nm, dad.n.comps[f]) # Sex-specific CKMR model equation
  }
}

####
HS.PO_noLambda_annual_model = function(){
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival

  #Likelihood
  #Moms
  #HS + PO
  for(i in 1:mom.yrs){ # Loop over maternal cohort comparisons
    mom.positives[i] ~ dbin((survival^mom.mort.yrs[i])/Nf, mom.n.comps[i]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/Nm, dad.n.comps[f]) # Sex-specific CKMR model equation
  }
}


####
HS.PO_narrowLambda_annual_model = function(){
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


####
HS.PO_wideLambda_annual_model = function(){
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival
  lambda ~ dunif(0.80, 1.20) #Diffuse prior for lambda
  

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



####
HS.only_noLambda_Skip_model = function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1.0) #Percent of animals breeding bi-ennially; CHANGED from dunif(0,1)
  

  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
    mom.positives_HS.even[i] ~ dbin((a*(survival^mom.mort.yrs_HS.even[i]))/((a + psi - (a*psi))*(Nf)), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
    mom.positives_HS.odd[j] ~ dbin(((survival^mom.mort.yrs_HS.odd[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf)), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
}


####

HSPOP_noLambda_Skip_model= function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1.0) #Percent of animals breeding bi-ennially; CHANGED from dunif(0,1)
  

  
  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
    mom.positives_HS.even[i] ~ dbin((a*(survival^mom.mort.yrs_HS.even[i]))/((a + psi - (a*psi))*(Nf)), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
    mom.positives_HS.odd[j] ~ dbin(((survival^mom.mort.yrs_HS.odd[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf)), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
  }

  #Moms
  #PO
  for(k in 1:mom.yrs_PO){
    mom.positives_PO[k] ~ dbin((survival^mom.mort.yrs_PO[k])/Nf, mom.n.comps_PO[k]) # Sex-specific CKMR model equation
  }
    
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/Nm, dad.n.comps[f]) # Sex-specific CKMR model equation
  }
}

####

HS.only_wideLambda_Skip_model = function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  lambda ~ dunif(0.8, 1.2)
  
  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
    mom.positives_HS.even[i] ~ dbin((a*(survival^mom.mort.yrs_HS.even[i]))/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.even[i]))), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
    mom.positives_HS.odd[j] ~ dbin(((survival^mom.mort.yrs_HS.odd[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.odd[j]))), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
  
}

####

HSPOP_wideLambda_Skip_model= function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  lambda ~ dunif(0.8, 1.2)
  
  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
    mom.positives_HS.even[i] ~ dbin((a*(survival^mom.mort.yrs_HS.even[i]))/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.even[i]))), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
    mom.positives_HS.odd[j] ~ dbin(((survival^mom.mort.yrs_HS.odd[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.odd[j]))), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #PO
  for(k in 1:mom.yrs_PO){
    mom.positives_PO[k] ~ dbin((survival^mom.mort.yrs_PO[k])/(Nf*(lambda^mom.popGrowth.yrs_PO[k])), mom.n.comps_PO[k]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
}
####

HS.only_narrowLambda_Skip_model = function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  lambda ~ dunif(0.95, 1.05)
  
  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
    mom.positives_HS.even[i] ~ dbin((a*(survival^mom.mort.yrs_HS.even[i]))/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.even[i]))), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
    mom.positives_HS.odd[j] ~ dbin(((survival^mom.mort.yrs_HS.odd[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.odd[j]))), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
  
}

####

HSPOP_narrowLambda_Skip_model= function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  lambda ~ dunif(0.95, 1.05)
  
  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.even){ # Loop over maternal cohort comparisons
    mom.positives_HS.even[i] ~ dbin((a*(survival^mom.mort.yrs_HS.even[i]))/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.even[i]))), mom.n.comps_HS.even[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.odd){ # Loop over maternal cohort comparisons
    mom.positives_HS.odd[j] ~ dbin(((survival^mom.mort.yrs_HS.odd[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.odd[j]))), mom.n.comps_HS.odd[j]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #PO
  for(k in 1:mom.yrs_PO){
    mom.positives_PO[k] ~ dbin((survival^mom.mort.yrs_PO[k])/(Nf*(lambda^mom.popGrowth.yrs_PO[k])), mom.n.comps_PO[k]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
}

# Write models
jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation")
write_model(HS.PO_noLambda_annual_model_validation, jags_file)

jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model.txt")
write_model(HS.PO_noLambda_annual_model, jags_file)

jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
write_model(HS.PO_narrowLambda_annual_model, jags_file)

jags_file = paste0(jags.model_location, "HS.PO_wideLambda_annual_model.txt")
write_model(HS.PO_wideLambda_annual_model, jags_file)

jags_file = paste0(jags.model_location, "HS.only_noLambda_Skip_model.txt")
write_model(HS.only_noLambda_Skip_model, jags_file)

jags_file = paste0(jags.model_location, "HSPOP_noLambda_Skip_model.txt")
write_model(HSPOP_noLambda_Skip_model, jags_file)

jags_file = paste0(jags.model_location, "HS.only_wideLambda_Skip_model.txt")
write_model(HS.only_wideLambda_Skip_model, jags_file)

jags_file = paste0(jags.model_location, "HSPOP_wideLambda_Skip_model.txt")
write_model(HSPOP_wideLambda_Skip_model, jags_file)

jags_file = paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt")
write_model(HS.only_narrowLambda_Skip_model, jags_file)

jags_file = paste0(jags.model_location, "HSPOP_narrowLambda_Skip_model.txt")
write_model(HSPOP_narrowLambda_Skip_model, jags_file)
