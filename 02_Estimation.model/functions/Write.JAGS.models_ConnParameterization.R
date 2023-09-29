#Wrote annual model but not biennial model
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models
  
####
HS.PO_narrowLambda_annual_model_Conn = function(){
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf0 ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm0 ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival
  lambda ~ dunif(0.95, 1.05) #Narrow prior for lambda
  
  #Likelihood
  #Moms
  #HS + PO
  for(i in 1:mom.yrs){ # Loop over maternal cohort comparisons
    mom.positives[i] ~ dbin((survival^mom.mort.yrs[i])/(Nf0*(lambda^mom.popGrowth.yrs[i])), mom.n.comps[i]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm0*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
  
  #derived quantities
  Nft = Nf0 * lambda^(estimation.year - est.year.calibrate)
  Nmt = Nm0 * lambda^(estimation.year - est.year.calibrate)
}

####
HS.PO_wideLambda_annual_model_Conn = function(){
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf0 ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm0 ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5, 0.95) # Uninformative prior for adult survival
  lambda ~ dunif(0.8, 1.2) #Narrow prior for lambda
  
  #Likelihood
  #Moms
  #HS + PO
  for(i in 1:mom.yrs){ # Loop over maternal cohort comparisons
    mom.positives[i] ~ dbin((survival^mom.mort.yrs[i])/(Nf0*(lambda^mom.popGrowth.yrs[i])), mom.n.comps[i]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm0*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
  
  #derived quantities
  Nft = Nf0 * lambda^(estimation.year - est.year.calibrate)
  Nmt = Nm0 * lambda^(estimation.year - est.year.calibrate)
}


####
HS.only_narrowLambda_Skip_model_Conn = function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf0 ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm0 ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.99) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  lambda ~ dunif(0.95, 1.05)
  
  #Likelihood
  #Moms
  #HS - even years
  for(i in 1:mom.yrs_HS.on){ # Loop over maternal cohort comparisons
    mom.positives_HS.on[i] ~ dbin((a*(survival^mom.mort.yrs_HS.on[i]))/((a + psi - (a*psi))*(Nf0*(lambda^mom.popGrowth.yrs_HS.on[i]))), mom.n.comps_HS.on[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.off){ # Loop over maternal cohort comparisons
    mom.positives_HS.off[j] ~ dbin(((survival^mom.mort.yrs_HS.off[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf0*(lambda^mom.popGrowth.yrs_HS.off[j]))), mom.n.comps_HS.off[j]) # Sex-specific CKMR model equation
  }
  
  #Dads
  #HS + PO
  for(f in 1:dad.yrs){ # Loop over paternal cohort comparisons
    dad.positives[f] ~ dbin((survival^dad.mort.yrs[f])/(Nm0*(lambda^dad.popGrowth.yrs[f])), dad.n.comps[f]) # Sex-specific CKMR model equation
  }
  
  #Derived quantities
  #derived quantities
  Nft = Nf0 * lambda^(estimation.year - est.year.calibrate)
  Nmt = Nm0 * lambda^(estimation.year - est.year.calibrate)
  Nfb0 = Nf0/a
  Nfbt = Nft/a
}

####

HS.PO_narrowLambda_Skip_model_Conn = function(){
  
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
  for(i in 1:mom.yrs_HS.on){ # Loop over maternal cohort comparisons
    mom.positives_HS.on[i] ~ dbin((a*(survival^mom.mort.yrs_HS.on[i]))/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.on[i]))), mom.n.comps_HS.on[i]) # Sex-specific CKMR model equation
  }
  
  #Moms
  #HS - odd years
  for(j in 1:mom.yrs_HS.off){ # Loop over maternal cohort comparisons
    mom.positives_HS.off[j] ~ dbin(((survival^mom.mort.yrs_HS.off[j])*(1-psi)*a)/((a + psi - (a*psi))*(Nf*(lambda^mom.popGrowth.yrs_HS.off[j]))), mom.n.comps_HS.off[j]) # Sex-specific CKMR model equation
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
  
  #Derived quantities
  Nfb1 <- ((a + psi - a*psi)/a)*Nf
  Nfb2 <- Nf/a
  
}



jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model_Conn.txt")
write_model(HS.PO_narrowLambda_annual_model_Conn, jags_file)

jags_file = paste0(jags.model_location, "HS.PO_wideLambda_annual_model_Conn.txt")
write_model(HS.PO_wideLambda_annual_model_Conn, jags_file)

jags_file = paste0(jags.model_location, "HS.only_narrowLambda_Skip_model_Conn.txt")
write_model(HS.only_narrowLambda_Skip_model_Conn, jags_file)

jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_Skip_model_Conn.txt")
write_model(HS.PO_narrowLambda_Skip_model_Conn, jags_file)