P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs))

surv <- 0.85

get_P_lemon <- function(Pars1,P_Mother,P_Father,t_start,t_end){
  N_F=exp(Pars1[1]) #number of mature females (assume time constant)
  
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      
      #Basic way of doing it - one pop_growth
      P_Mother[os_birth, ys_birth] <- (surv^(ys_birth - os_birth))/(N_F*lam^(ys_birth-min_est_cohort))
      
      #More complicated way of doing it using observed pop_growth
      #   if(ys_birth == min_est_cohort_F){
      #   P_Mother[os_birth, ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/(N_F*pop_growth_F[min_est_cohort_F])
      # }
      # else  { #Multiply by the observed population growth rate
      #   NtempF <- N_F
      #   for(z in 1:length(min_est_cohort_F:(ys_birth-1))){
      #   popGF <- pop_growth_F[min_est_cohort_F:(ys_birth-1)]
      #   NtempF <- NtempF*popGF[z]
      #   }
      #   P_Mother[os_birth, ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/NtempF
      #   }
    }
  }
  N_M=exp(Pars1[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      
      #Basic way of doing it: one pop growth
      P_Father[os_birth,ys_birth] <- (surv^(ys_birth - os_birth))/(N_M*lam^(ys_birth-min_est_cohort))
      
      #More complicated way using observed pop growth
      # if(ys_birth == min_est_cohort_M){
      #   P_Father[os_birth,ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/(N_M*pop_growth_M[min_est_cohort_M])
      # 
      # }
      # else{ #Multiply by the observed population growth rate
      #   NtempM <- N_M
      #   for(z in 1:length(min_est_cohort_M:(ys_birth-1))){
      #     popGM <- pop_growth_M[min_est_cohort_M:(ys_birth-1)]
      #     NtempM <- NtempM*popGM[z]
      #   }
      #   P_Father[os_birth,ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/NtempM
      #}
    }
  }
  
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}

P=get_P_lemon(Pars1=Pars1,P_Mother=P_Mother,P_Father=P_Father,t_start=t_start,t_end=t_end)