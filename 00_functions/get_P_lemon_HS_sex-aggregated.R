#Sex-aggregated model for calculating prior probability of kinship for half-siblings based on birth years

adult_mat_age <- c(7:max.age) #Adult age of females -males assumed to be the same
adult_age_at_mat <- min(adult_mat_age) #minimum age an adult could be mature

#Set up empty array that will be filled with function below
P_Parent = array(0, dim=c(n_yrs,n_yrs)) #Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

#CKMR model: populate array with kinship probabilities based on birth years
get_P_lemon_TotalA <- function(Pars2,P_Parent,t_start,t_end){
  N_A=exp(Pars2[1]) #number of mature adults

  for(os_birth in min_est_cohort:(n_yrs-1)){  #Loop over possible ages of older sibs
    for(ys_birth in max(os_birth+1):n_yrs){ #Loop over possible ages of younger sibs
      if((ys_birth - os_birth) <= ((maxAge+1) - adult_age_at_mat)){ #If the adult could have been mature in the birth year of both individual, then fill the array with the appropriate probability of kinship
        
        #Probability of kinship based on birth year
        #See Hillary et al (2018) equation (3)
        P_Parent[os_birth, ys_birth] <- (4/(N_A*lam^(ys_birth - min_est_cohort)))*(surv^(ys_birth - os_birth))
      } else P_Parent[os_birth, ys_birth] <- 0 #If it's not possible, set kinship probability to 0
    }
  }
  return(list(P_Parent=P_Parent)) #return makes sure this is moved out of the loop into the environment
}

#Store function as P_TotalA
P_TotalA=get_P_lemon_TotalA(Pars2=Pars2,P_Parent,t_start=t_start,t_end=t_end)

#View(P_TotalA$P_Parent)