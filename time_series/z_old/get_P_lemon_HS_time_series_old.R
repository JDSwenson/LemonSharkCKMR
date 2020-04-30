#Model for calculating prior probability of kinship for half-siblings

#Not assessing survival - keeping constant at 0.85
Surv <- 0.85

max_age = 25 #max age of lemon sharks

m_adult_age <- c(12:25) #Set ages at which males and females are mature.
f_adult_age <- c(13:25)
m_mat <- c(rep(0,11), rep(1,14)) #Set proportion of mature males and females at each age -- assumes knife-edge maturity.
f_mat <- c(rep(0,12), rep(1,13))

#If estimating survival, activate below code
#Pars=c(log(a_priori_abund/2),log(a_priori_abund/2), a_priori_surv)

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1:est_yrs]) #number of mature females
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
    for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
      P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F[ys_birth] #N_F is the number of females alive when the younger sibling was born
      } else P_Mother[os_birth, ys_birth] <- 0
        }
    }
  N_M=exp(Pars[(est_yrs+1):(est_yrs*2)])
  #N_M=exp(Pars[(n_yrs+1):(n_yrs*2)]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
      P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M[ys_birth]
      } else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}