#9/16/19 Just changed the maturity age in the likelihood functions to match the maturity age used to generate the data (different for males and females)

#Not assessing survival - keeping constant at 0.85
#Pars=c(log(a_priori_abund/2),log(a_priori_abund/2))
Surv <- 0.85
max_age = 25 #max age of lemon sharks
m_adult_age <- c(12:25) #Set ages at which males and females are mature.
f_adult_age <- c(13:25)
m_mat <- c(rep(0,11), rep(1,14)) #Set proportion of mature males and females at each age -- assumes knife-edge maturity.
f_mat <- c(rep(0,12), rep(1,13))

P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)
### Maybe come back and use a truncated distribution instead (package truncnorm)


get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant) - (total reproductive output from females)
    for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
      P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F
      } else P_Mother[os_birth, ys_birth] <- 0
        }
      }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
      P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M
      } else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}