#Model for calculating prior probability of kinship for half-siblings
#Not assessing survival - keeping constant at 0.85
Surv <- 0.85

max_age = 30 #max age of lemon sharks

adult_age <- c(12:30) #Adult age of females -males assumed to be the same
adult_mat <- c(rep(0,11), rep(1,19)) #knife-edge maturity

P_Parent = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

get_P_lemon_TotalA <- function(Pars2,P_Parent,t_start,t_end){
  N_A=exp(Pars2[1]) #number of mature females
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
    for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in max(os_birth+1):n_yrs){
        P_Parent[os_birth, ys_birth] <- (4/(N_A*lam^(ys_birth - min_est_cohort)))*(Surv^(ys_birth - os_birth)) #N_F is the number of females alive when the younger sibling was born
        }
    }
  return(list(P_Parent=P_Parent)) #return makes sure this is moved out of the loop into the environment
}

P_TotalA=get_P_lemon_TotalA(Pars2=Pars2,P_Parent,t_start=t_start,t_end=t_end)
#P_TotalA$P_Parent