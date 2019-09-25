#9/16/19 Just changed the maturity age in the likelihood functions to match the maturity age used to generate the data (different for males and females)

#n_ages=25 #Max age of lemon sharks
#n_yrs=29
#t_start=26
#t_end=29
#mat_age=rep(0, n_ages)
#mat_age[12:n_ages]=1
P_Mother = P_Father = array(0,dim=c(n_yrs,n_yrs,n_yrs))  #creates two empty arrays, one for mother and one for father.  Dimensions are parent birth year, parent capture year, offspring birth year (all of which are specified by n_yrs)
### Maybe come back and use a truncated distribution instead (package truncnorm)
#a_priori_abund <- round(rnorm(n=1, mean=80, sd=80),0)
#a_priori_abund <- ifelse(a_priori_abund < 10, 10, a_priori_abund)

#Not assessing survival - keeping constant at 0.85
Pars=c(log(a_priori_abund/2),log(a_priori_abund/2))
Surv <- 0.85
#Estimating survival
#Pars=c(log(a_priori_abund/2),log(a_priori_abund/2), a_priori_surv)
#Surv <- Pars[3]

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant) - (total reproductive output from females)
    for(abirth in 1:n_yrs){  #> = after, < = before
    for(acapt in max(abirth, t_start):t_end){ 
      for(jbirth in 1:n_yrs){
          if(jbirth <= acapt & jbirth > abirth & ((jbirth - abirth) <= n_ages)){
          P_Mother[abirth,acapt,jbirth] <- f_mat[jbirth -abirth]/N_F
          } else if(jbirth > acapt & jbirth > abirth & ((jbirth - abirth) <= n_ages)) {
            P_Mother[abirth, acapt, jbirth] <- (Surv*(jbirth - acapt))/N_F
          }
        }
      }
    }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(abirth in 1:n_yrs){  #dad
    for(acapt in max(abirth, t_start):t_end){ 
      for(jbirth in 1:n_yrs){
        if(jbirth <= acapt & jbirth > abirth & ((jbirth - abirth) <= n_ages)) {
          P_Father[abirth, acapt, jbirth] <- m_mat[jbirth - abirth]/N_M
        } else if(jbirth > acapt & jbirth > abirth & ((jbirth - abirth) <= n_ages)) {
          P_Father[abirth, acapt, jbirth] <- (Surv*(jbirth - acapt))/N_M
        }
      }
    } 
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}