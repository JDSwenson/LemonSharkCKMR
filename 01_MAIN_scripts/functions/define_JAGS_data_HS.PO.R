#Define data
jags_data = list(
  #Year 1
  #Moms - half-sib
  MHSP1 = HS.mom_comps[[1]]$yes, # Positive maternal half-sibs; Y
  mom_n_comps.HS1 = HS.mom_comps[[1]]$all, # Number of total maternal comparisons; R
  mom_ys_birth1 = HS.mom_comps[[1]]$Ind_2_birth, # birth year of younger sib; b
  mom_os_birth1 = HS.mom_comps[[1]]$Ind_1_birth, # birth year of older sib; a
  mom_yrs.HS1 = nrow(HS.mom_comps[[1]]), # number of cohort comparisons to loop over 
  
  #Dads - half-sib
  FHSP1 = HS.dad_comps[[1]]$yes, # Positive paternal half-sibs; Y
  dad_n_comps.HS1 = HS.dad_comps[[1]]$all, # Number of total paternal comparisons; R
  dad_ys_birth1 = HS.dad_comps[[1]]$Ind_2_birth, # birth year of younger sib; b
  dad_os_birth1 = HS.dad_comps[[1]]$Ind_1_birth, # birth year of older sib; a
  dad_yrs.HS1 = nrow(HS.dad_comps[[1]]), # number of cohort comparisons to loop over
  
  
  #Moms - PO
  MPOP1 = PO.mom_comps[[1]]$yes,
  mom_n_comps.PO1 = PO.mom_comps[[1]]$all,
  mom_off.birth1 = PO.mom_comps[[1]]$offspring.birth.year,
  mom.birth1 = PO.mom_comps[[1]]$parent.birth.year,
  mom.capture1 = PO.mom_comps[[1]]$parent.capture.year,
  mom_yrs.PO1 = nrow(PO.mom_comps[[1]]),
  
  #Dads - PO
  FPOP1 = PO.dad_comps[[1]]$yes,
  dad_n_comps.PO1 = PO.dad_comps[[1]]$all,
  dad_off.birth1 = PO.dad_comps[[1]]$offspring.birth.year,
  dad.birth1 = PO.dad_comps[[1]]$parent.birth.year,
  dad.capture1 = PO.dad_comps[[1]]$parent.capture.year,
  dad_yrs.PO1 = nrow(PO.dad_comps[[1]]),
  
  
  #Year 2
  #Moms - half-sib
  MHSP2 = HS.mom_comps[[2]]$yes, # Positive maternal half-sibs; Y
  mom_n_comps.HS2 = HS.mom_comps[[2]]$all, # Number of total maternal comparisons; R
  mom_ys_birth2 = HS.mom_comps[[2]]$Ind_2_birth, # birth year of younger sib; b
  mom_os_birth2 = HS.mom_comps[[2]]$Ind_2_birth, # birth year of older sib; a
  mom_yrs.HS2 = nrow(HS.mom_comps[[2]]), # number of cohort comparisons to loop over 
  
  #Dads - half-sib
  FHSP2 = HS.dad_comps[[2]]$yes, # Positive paternal half-sibs; Y
  dad_n_comps.HS2 = HS.dad_comps[[2]]$all, # Number of total paternal comparisons; R
  dad_ys_birth2 = HS.dad_comps[[2]]$Ind_2_birth, # birth year of younger sib; b
  dad_os_birth2 = HS.dad_comps[[2]]$Ind_2_birth, # birth year of older sib; a
  dad_yrs.HS2 = nrow(HS.dad_comps[[2]]), # number of cohort comparisons to loop over
  
  #Moms - PO
  MPOP2 = PO.mom_comps[[2]]$yes,
  mom_n_comps.PO2 = PO.mom_comps[[2]]$all,
  mom_off.birth2 = PO.mom_comps[[2]]$offspring.birth.year,
  mom.birth2 = PO.mom_comps[[2]]$parent.birth.year,
  mom.capture2 = PO.mom_comps[[2]]$parent.capture.year,
  mom_yrs.PO2 = nrow(PO.mom_comps[[2]]),
  
  #Dads - PO
  FPOP2 = PO.dad_comps[[2]]$yes,
  dad_n_comps.PO2 = PO.dad_comps[[2]]$all,
  dad_off.birth2 = PO.dad_comps[[2]]$offspring.birth.year,
  dad.birth2 = PO.dad_comps[[2]]$parent.birth.year,
  dad.capture2 = PO.dad_comps[[2]]$parent.capture.year,
  dad_yrs.PO2 = nrow(PO.dad_comps[[2]]),
  
  #Fix other potential parameters
  #surv = surv,
  yr1 = min_cohort, # estimation year i.e. year the estimate will be focused on
  yr2 = n_yrs,
  tau = 1E-6
)


#----------------- STEP 2: SPECIFY JAGS MODEL CODE ---------------#
#Convert tau to SD (for interpretation)
#tau <- 1E-6
#(sd <- sqrt(1/tau))

HS.PO_model = function(){
  
  #PRIORS
  Nf1 ~ dnorm(0, tau) # Uninformative prior for female abundance year 1
  Nm1 ~ dnorm(0, tau) # Uninformative prior for male abundance year 1
  Nf2 ~ dnorm(0, tau) # Uninformative prior for female abundance year 2
  Nm2 ~ dnorm(0, tau) # Uninformative prior for male abundance year 2
  surv ~ dbeta(1 ,1) # Uninformative prior for adult survival
  
  #Year 1
  #HS Likelihood - year 1
  for(i in 1:mom_yrs.HS1){ # Loop over maternal cohort comparisons
    MHSP1[i] ~ dbin((surv^(mom_ys_birth1[i] - mom_os_birth1[i]))/Nf1, mom_n_comps.HS1[i]) # Sex-specific CKMR model equation
  }
  for(j in 1:dad_yrs.HS1){ # Loop over paternal cohort comparisons
    FHSP1[j] ~ dbin((surv^(dad_ys_birth1[j] - dad_os_birth1[j]))/Nm1, dad_n_comps.HS1[j]) # Sex-specific CKMR model equation
  }
  
  #PO Likelihood - year 1
  for(k in 1:mom_yrs.PO1){ # Loop over maternal cohort comparisons
    MPOP1[k] ~ dbin(1/Nf1, mom_n_comps.PO1[k]) # Sex-specific CKMR model equation
  }
  for(l in 1:dad_yrs.PO1){ # Loop over paternal cohort comparisons
    FPOP1[l] ~ dbin(1/Nm1, dad_n_comps.PO1[l]) # Sex-specific CKMR model equation
  }

    #HS Likelihood - year 2
  for(m in 1:mom_yrs.HS2){ # Loop over maternal cohort comparisons
    MHSP2[m] ~ dbin((surv^(mom_ys_birth2[m] - mom_os_birth2[m]))/Nf2, mom_n_comps.HS2[m]) # Sex-specific CKMR model equation
  }
  for(n in 1:dad_yrs.HS2){ # Loop over paternal cohort comparisons
    FHSP2[n] ~ dbin((surv^(dad_ys_birth2[n] - dad_os_birth2[n]))/Nm2, dad_n_comps.HS2[n]) # Sex-specific CKMR model equation
  }
  
  #PO Likelihood - year 2: known alive
  for(q in 1:mom_yrs.PO2){ # Loop over maternal cohort comparisons
    MPOP2[q] ~ dbin(1/Nf2, mom_n_comps.PO2[q]) # Sex-specific CKMR model equation
  }
  for(r in 1:dad_yrs.PO2){ # Loop over paternal cohort comparisons
    FPOP2[r] ~ dbin(1/Nm2, dad_n_comps.PO2[r]) # Sex-specific CKMR model equation
  }
}