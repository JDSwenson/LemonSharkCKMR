#Not assessing survival - keeping constant at 0.85
Surv <- 0.85
max_age = 30 #max age of lemon sharks
m_adult_age <- c(12:max_age) #Set ages at which males and females are mature.
f_adult_age <- c(13:max_age)

P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)
### Maybe come back and use a truncated distribution instead (package truncnorm)


get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant)
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

#P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

#View(P$P_Mother)