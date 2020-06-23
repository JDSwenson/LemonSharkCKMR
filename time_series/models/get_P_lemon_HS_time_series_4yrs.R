#Model for calculating prior probability of kinship for half-siblings

P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1:est_ages]) #number of mature females
  index <- min_est_cohort-1
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
    for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
      P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F[ys_birth-index] #N_F is the number of females alive when the younger sibling was born
      } else P_Mother[os_birth, ys_birth] <- 0
        }
    }
  N_M=exp(Pars[(est_ages+1):(est_ages*2)])
  #N_M=exp(Pars[(n_yrs+1):(n_yrs*2)]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
      P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M[ys_birth-index]
      } else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}
#P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P$P_Mother
#P$P_Father
