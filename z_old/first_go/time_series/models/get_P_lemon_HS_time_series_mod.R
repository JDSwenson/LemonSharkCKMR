#Model for calculating prior probability of kinship for half-siblings
#Estimate years 10:17 of data 5/18/20
#Removed max age requirement from function because there are immortal lemon sharks in the dataset

#Not assessing survival - keeping constant at 0.85
Surv <- 0.85
max_age = 25 #max age of lemon sharks
m_adult_age <- c(12:25) #Set ages at which males and females are mature.
f_adult_age <- c(13:25)
m_mat <- c(rep(0,11), rep(1,14)) #Set proportion of mature males and females at each age -- assumes knife-edge maturity.
f_mat <- c(rep(0,12), rep(1,13))

P_Mother = P_Father = array(NA,dim=c(max_est_cohort,max_est_cohort)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1:n_yrs]) #number of mature females
  index <- min_est_cohort-1
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
    for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      #if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
      P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F[ys_birth-index] #N_F is the number of females alive when the younger sibling was born
      #} else P_Mother[os_birth, ys_birth] <- 0
        }
    }
  N_M=exp(Pars[(n_yrs+1):(n_yrs*2)])
  #N_M=exp(Pars[(n_yrs+1):(n_yrs*2)]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      #if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
      P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M[ys_birth-index]
      #} else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}
#P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P$P_Mother
#P$P_Father
