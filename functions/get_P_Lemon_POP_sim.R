#n_ages=25 #Max age of lemon sharks
#n_yrs=29
#t_start=26
#t_end=29
#mat_age=rep(0, n_ages)
#mat_age[12:n_ages]=1
n_yrs=t_end
n_ages <- maxAge
P_Mother = P_Father = array(0,dim=c(n_yrs,n_yrs,n_yrs))  #creates two empty arrays, one for mother and one for father.  Dimensions are parent birth year, parent capture year, offspring birth year (all of which are specified by n_yrs)
### Maybe come back and use a truncated distribution instead (package truncnorm)
#a_priori_abund <- round(rnorm(n=1, mean=80, sd=80),0)
#a_priori_abund <- ifelse(a_priori_abund < 10, 10, a_priori_abund)
f_adult_age <- min(which(femaleCurve>0))
m_adult_age <- min(which(maleCurve>0))

#Not assessing survival - keeping constant at 0.85
#Pars=c(log(a_priori_abund/2),log(a_priori_abund/2))
Surv <- 0.85
#Estimating survival
#Pars=c(log(a_priori_abund/2),log(a_priori_abund/2), a_priori_surv)
#Surv <- Pars[3]

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant) - (total reproductive output from females)
  
  #Dimensions are parent birth year, parent capture year, offspring birth year
   
   for(abirth in 1:n_yrs){  #a < b = a was born after b (is older), > = opposite
     for(acapt in max(abirth, t_start):t_end){
       for(jbirth in 1:n_yrs){
         if(jbirth <= acapt & jbirth > abirth & ((acapt - abirth) <= n_ages) & ((jbirth - abirth) <= n_ages)){
          P_Mother[abirth,acapt,jbirth] <- mat_f[jbirth -abirth]/N_F #If juvenile was born before the adult was captured AND after the adult was born AND the adult wasn't too old to have been captured OR to have given birth to the juvenile, then probability equals f_mat at the age the adult was when the juvenile was born over N_f
          } else if(jbirth > acapt & jbirth - abirth >= f_adult_age & ((acapt - abirth) <= n_ages) & ((jbirth - abirth) <= n_ages)) {
            P_Mother[abirth, acapt, jbirth] <- (Surv^(jbirth - acapt))/N_F #If juvenile was born after the adult was captured AND the adult was old enough to have given birth AND the adult wasn't too old to have been captured OR to have given birth to the juvenile, then the probability equals prob of survival from adult capture year to juvenile birth
            #Age distribution was set with survival = 0.85, so this should work
          }
        }
      }
    }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(abirth in 1:n_yrs){  #dad
    for(acapt in max(abirth, t_start):t_end){ 
      for(jbirth in 1:n_yrs){
        if(jbirth <= acapt & jbirth > abirth & ((acapt - abirth) <= n_ages) & ((jbirth - abirth) <= n_ages)){
          P_Father[abirth,acapt,jbirth] <- mat_m[jbirth -abirth]/N_M #If juvenile was born before the adult was captured AND after the adult was born AND the adult wasn't too old, then probability equals f_mat at the age the adult was when the juvenile was born over N_f
        } else if(jbirth > acapt & jbirth - abirth >= m_adult_age & ((acapt - abirth) <= n_ages) & ((jbirth - abirth) <= n_ages)) {
          P_Father[abirth, acapt, jbirth] <- Surv^((jbirth - acapt))/N_M #If juvenile was born after the adult was captured AND the adult was old enough to have given birth AND the adult wasn't too old, then the probability equals prob of survival from adult capture year to juvenile birth
          #Age distribution was set with survival = 0.85, so this should work
        }
      }
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}