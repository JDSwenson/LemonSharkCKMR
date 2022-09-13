HS.only_noLambda_Skip_model = function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  
  #Ben's skip-spawn pre-calcs
  
  pNf.ann <- (1 - psi) #Number of annual spawners is easy
  
  pNf.skip <- ((psi * (1 - survival) / (1 - survival^(a-1))))
  pEff <- (pNf.ann + pNf.skip)
  
  #Likelihood
  for(i in 1:mom.yrs_HS){ # Loop over maternal cohort comparisons
    
    #Ben's skip-spawn likelihood
    mom.positives_HS[i] ~ dbin((((survival^(mom.mort.yrs_HS[i])) / (Nf * pEff)) * mom.oncycle[i]) + #For on-cycle spawning (annual + psi/a)
                                 ((((survival^mom.mort.yrs_HS[i]) * (pNf.ann / pEff)) / (Nf * pEff)) * (1 - mom.oncycle[i])), #For off-cycle spawning (annual only)
                               mom.n.comps_HS[i])
  }
  
  for(j in 1:dad.yrs_HS){ # Loop over paternal cohort comparisons
    dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/Nm), dad.n.comps_HS[j])
  }

}

####

HSPOP_noLambda_Skip_model= function(){
  
  #PRIORS - uninformative
  mu ~ dunif(1, 10000)
  sd ~ dunif(1, 10000)
  Nf ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for female abundance
  Nm ~ dnorm(mu, 1/(sd^2)) # Uninformative prior for male abundance
  survival ~ dunif(0.5,0.95) # Uninformative prior for adult survival
  psi ~ dunif(0, 1) #Percent of animals breeding non-annually #Turn off if fixing psi
  
  #Ben's skip-spawn pre-calcs
  
  pNf.ann <- (1 - psi) #Number of annual spawners is easy
  
  pNf.skip <- ((psi * (1 - survival) / (1 - survival^(a-1))))
  pEff <- (pNf.ann + pNf.skip)
  
  #Likelihood
  for(i in 1:mom.yrs_HS){ # Loop over maternal cohort comparisons
    
    #Ben's skip-spawn likelihood
    mom.positives_HS[i] ~ dbin((((survival^(mom.mort.yrs_HS[i])) / (Nf * pEff)) * mom.oncycle[i]) + #For on-cycle spawning (annual + psi/a)
                                 ((((survival^mom.mort.yrs_HS[i]) * (pNf.ann / pEff)) / (Nf * pEff)) * (1 - mom.oncycle[i])), #For off-cycle spawning (annual only)
                               mom.n.comps_HS[i])
  }
  
  for(j in 1:dad.yrs_HS){ # Loop over paternal cohort comparisons
    dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/Nm), dad.n.comps_HS[j])
  }

  for(k in 1:POP.yrs)
  {
    positives_POPs[k] ~ dbin(((survival^POP.mort.yrs[k]) / Nm * POP.isPaternal[k]) +
                               ((survival^POP.mort.yrs[k]) / Nf * (1 - POP.isPaternal[k])),n.comps_POPs[k])
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
  
  #Ben's skip-spawn pre-calcs
  
  pNf.ann <- (1 - psi) #Number of annual spawners is easy
  
  pNf.skip <- ((psi * (1 - survival) / (1 - survival^(a-1))))
  pEff <- (pNf.ann + pNf.skip)
  
  #Likelihood
  for(i in 1:mom.yrs_HS){ # Loop over maternal cohort comparisons
    
    #Ben's skip-spawn likelihood
    mom.positives_HS[i] ~ dbin((((survival^(mom.mort.yrs_HS[i])) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * mom.oncycle[i]) + #For on-cycle spawning (annual + psi/a)
                                 ((((survival^mom.mort.yrs_HS[i]) * (pNf.ann / pEff)) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * (1 - mom.oncycle[i])), #For off-cycle spawning (annual only)
                               mom.n.comps_HS[i])
  }
  
  for(j in 1:dad.yrs_HS){ # Loop over paternal cohort comparisons
    dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/(Nm * (lambda^dad.popGrowth.yrs_HS[j]))), dad.n.comps_HS[j])
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
  
  #Ben's skip-spawn pre-calcs
  
  pNf.ann <- (1 - psi) #Number of annual spawners is easy
  
  pNf.skip <- ((psi * (1 - survival) / (1 - survival^(a-1))))
  pEff <- (pNf.ann + pNf.skip)
  
  #Likelihood
  for(i in 1:mom.yrs_HS){ # Loop over maternal cohort comparisons
    
    #Ben's skip-spawn likelihood
    mom.positives_HS[i] ~ dbin((((survival^(mom.mort.yrs_HS[i])) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * mom.oncycle[i]) + #For on-cycle spawning (annual + psi/a)
                                 ((((survival^mom.mort.yrs_HS[i]) * (pNf.ann / pEff)) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * (1 - mom.oncycle[i])), #For off-cycle spawning (annual only)
                               mom.n.comps_HS[i])
  }
  
  for(j in 1:dad.yrs_HS){ # Loop over paternal cohort comparisons
    dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/(Nm * (lambda^dad.popGrowth.yrs_HS[j]))), dad.n.comps_HS[j])
  }
  
  for(k in 1:POP.yrs)
  {
    positives_POPs[k] ~ dbin((( survival^POP.mort.yrs[k]) / (Nm * (lambda^POP.popGrowth.yrs[k]))  * POP.isPaternal[k]) +
                               (( survival^POP.mort.yrs[k])/ (Nf * (lambda^POP.popGrowth.yrs[k])) * (1 - POP.isPaternal[k])),n.comps_POPs[k])
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
  
  #Ben's skip-spawn pre-calcs
  
  pNf.ann <- (1 - psi) #Number of annual spawners is easy
  
  pNf.skip <- ((psi * (1 - survival) / (1 - survival^(a-1))))
  pEff <- (pNf.ann + pNf.skip)
  
  #Likelihood
  for(i in 1:mom.yrs_HS){ # Loop over maternal cohort comparisons
    
    #Ben's skip-spawn likelihood
    mom.positives_HS[i] ~ dbin((((survival^(mom.mort.yrs_HS[i])) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * mom.oncycle[i]) + #For on-cycle spawning (annual + psi/a)
                                 ((((survival^mom.mort.yrs_HS[i]) * (pNf.ann / pEff)) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * (1 - mom.oncycle[i])), #For off-cycle spawning (annual only)
                               mom.n.comps_HS[i])
  }
  
  for(j in 1:dad.yrs_HS){ # Loop over paternal cohort comparisons
    dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/(Nm * (lambda^dad.popGrowth.yrs_HS[j]))), dad.n.comps_HS[j])
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
  
  #Ben's skip-spawn pre-calcs
  
  pNf.ann <- (1 - psi) #Number of annual spawners is easy
  
  pNf.skip <- ((psi * (1 - survival) / (1 - survival^(a-1))))
  pEff <- (pNf.ann + pNf.skip)
  
  #Likelihood
  for(i in 1:mom.yrs_HS){ # Loop over maternal cohort comparisons
    
    #Ben's skip-spawn likelihood
    mom.positives_HS[i] ~ dbin((((survival^(mom.mort.yrs_HS[i])) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * mom.oncycle[i]) + #For on-cycle spawning (annual + psi/a)
                                 ((((survival^mom.mort.yrs_HS[i]) * (pNf.ann / pEff)) / (Nf * pEff * (lambda^mom.popGrowth.yrs_HS[i]))) * (1 - mom.oncycle[i])), #For off-cycle spawning (annual only)
                               mom.n.comps_HS[i])
  }
  
  for(j in 1:dad.yrs_HS){ # Loop over paternal cohort comparisons
    dad.positives_HS[j] ~ dbin(((survival^dad.mort.yrs_HS[j])/(Nm * (lambda^dad.popGrowth.yrs_HS[j]))), dad.n.comps_HS[j])
  }
  
  for(k in 1:POP.yrs)
  {
    positives_POPs[k] ~ dbin(((survival^POP.mort.yrs[k]) / (Nm * (lambda^POP.popGrowth.yrs[k]))  * POP.isPaternal[k]) +
                               ((survival^POP.mort.yrs[k]) / (Nf * (lambda^POP.popGrowth.yrs[k])) * (1 - POP.isPaternal[k])),n.comps_POPs[k])
  }
}

# Write models

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


