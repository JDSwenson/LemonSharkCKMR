P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)
### Maybe come back and use a truncated distribution instead (package truncnorm)

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars1[1]) #number of mature females (assume time constant)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge+1) - firstBreed)){
        P_Mother[os_birth, ys_birth] <- (surv^(ys_birth - os_birth))/N_F
      } else P_Mother[os_birth, ys_birth] <- 0
    }
  }
  N_M=exp(Pars1[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge+1) - firstBreed)){
        P_Father[os_birth,ys_birth] <- (surv^(ys_birth - os_birth))/N_M
      } else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)



#Sex-aggregated model for calculating prior probability of kinship for half-siblings based on birth years

#Set up empty array that will be filled with function below
P_Parent = array(0, dim=c(n_yrs,n_yrs)) #Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

#CKMR model: populate array with kinship probabilities based on birth years
get_P_lemon_TotalA <- function(Pars2,P_Parent,t_start,t_end){
  N_A=exp(Pars2[1]) #number of mature adults
  
  for(os_birth in 1:(n_yrs-1)){  #Loop over possible ages of older sibs
    for(ys_birth in max(os_birth+1):n_yrs){ #Loop over possible ages of younger sibs
      if((ys_birth - os_birth) <= ((maxAge+1) - firstBreed)){ #If the adult could have been mature in the birth year of both individual, then fill the array with the appropriate probability of kinship
        
        #Probability of kinship based on birth year
        #See Hillary et al (2018) equation (3)
        P_Parent[os_birth, ys_birth] <- (4/(N_A^(ys_birth - os_birth)))*(surv^(ys_birth - os_birth))
      } else P_Parent[os_birth, ys_birth] <- 0 #If it's not possible, set kinship probability to 0
    }
  }
  return(list(P_Parent=P_Parent)) #return makes sure this is moved out of the loop into the environment
}

#Store function as P_TotalA
P_TotalA=get_P_lemon_TotalA(Pars2=Pars2,P_Parent,t_start=t_start,t_end=t_end)


#------------------Likelihood functions--------------------------
#For sex-specific model
lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,n_yrs,t_start,t_end){
  
  P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){
    loglik=loglik+Negatives_Mother[irow,3]*log(1-P$P_Mother[Negatives_Mother[irow,1],Negatives_Mother[irow,2]])
  } 
  for(irow in 1:nrow(Negatives_Father)){
    loglik=loglik+Negatives_Father[irow,3]*log(1-P$P_Father[Negatives_Father[irow,1],Negatives_Father[irow,2]])
  }  
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    loglik=loglik+Pairs_Mother[irow,3]*log(P$P_Mother[Pairs_Mother[irow,1],Pairs_Mother[irow,2]])
  }
  for(irow in 1:nrow(Pairs_Father)){
    loglik=loglik+Pairs_Father[irow,3]*log(P$P_Father[Pairs_Father[irow,1],Pairs_Father[irow,2]])
  }  
  -loglik
}



#Likelihood for sex-aggregated model
lemon_neg_log_lik_TotalA <- function(Pars2, Negatives_Parent, Pairs_Parent, P_Parent, t_start, t_end){
  
  P_TotalA=get_P_lemon_TotalA(Pars2=Pars2,P_Parent,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Parent)){
    loglik=loglik+Negatives_Parent[irow,3]*log(1-P_TotalA$P_Parent[Negatives_Parent[irow,1],Negatives_Parent[irow,2]])
  } 
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Parent)){
    loglik=loglik+Pairs_Parent[irow,3]*log(P_TotalA$P_Parent[Pairs_Parent[irow,1],Pairs_Parent[irow,2]])
  }
  -loglik
}
