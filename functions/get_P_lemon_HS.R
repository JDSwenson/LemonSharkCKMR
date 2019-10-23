#9/16/19 Just changed the maturity age in the likelihood functions to match the maturity age used to generate the data (different for males and females)

##Assumes all offspring were born in the study years -- will need to change if expanding to offsring born before study years
P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)
### Maybe come back and use a truncated distribution instead (package truncnorm)

#Not assessing survival - keeping constant at 0.85
Pars=c(log(a_priori_abund/2),log(a_priori_abund/2))
Surv <- 0.85
#Estimating survival
#Pars=c(log(a_priori_abund/2),log(a_priori_abund/2), a_priori_surv)
#Surv <- Pars[3]
#P$P_Mother[20,29]

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant) - (total reproductive output from females)
    for(os_birth in 1:n_yrs){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F
        }
      }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:n_yrs){  #> = after, < = before
    for(ys_birth in 1:n_yrs){
      P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}
