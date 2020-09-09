#Model for calculating prior probability of kinship for half-siblings
#Estimate years 3:6 of data 5/26/20
#Removed max age requirement from function because there are immortal lemon sharks in the dataset
#Confirmed that below function returns expected values!! 5/20/2020


#Not assessing survival - keeping constant at 0.85
Surv <- 0.85
max_age = 25 #max age of lemon sharks
m_adult_age <- c(12:25) #Set ages at which males and females are mature.
f_adult_age <- c(13:25)
m_mat <- c(rep(0,11), rep(1,14)) #Set proportion of mature males and females at each age -- assumes knife-edge maturity.
f_mat <- c(rep(0,12), rep(1,13))

unique_moms <- unique(as.character(mom_positives$Mother))
unique_dads <- unique(as.character(dad_positives$Father))

#Pars
#est_yrs

#Dimensions are older sib birth year, younger sib birth year, and number of unique moms
P_Mother_pos <- array(NA,dim=c(max_est_cohort,max_est_cohort, length(unique_moms)))
P_Father_pos <- array(NA,dim=c(max_est_cohort,max_est_cohort, length(unique_dads))) 

#Dimensions are older sib birth year and younger sib birth year
P_Mother_neg = P_Father_neg <- array(NA, dim=c(max_est_cohort, max_est_cohort))

#Function for positive comparisons
get_P_lemon_pos <- function(Pars,P_Mother_pos,P_Father_pos,n_yrs,t_start,t_end){
  N_F=exp(Pars[1:n_yrs]) #number of mature females - intialize
  index <- min_est_cohort-1 #Connects year of estimation with appropriate element in Pars when used below
    for(os_birth in 1:(max_est_cohort-1)){  #Loop through older sib birth years
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){ #Loop through younger sib birth years
      for(mommy in 1:length(unique_moms)) { #Loop through unique moms
        mom_i <- unique_moms[mommy] #Mom being evaluated in present iteration
        
        #Did mom have pups in birth years of older and younger sibs being compared? 1 if yes.
        os_indicator <- length(sum_table_mom$pups[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == os_birth)]) 
        ys_indicator <- length(sum_table_mom$pups[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == ys_birth)])
        
        #If the mother had offspring in both the older sib and younger sib birth years, then ...
      if(os_indicator == 1 & ys_indicator == 1) {
        mom_pups_os <- sum_table_mom$pups[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == os_birth)] #Observed number of pups mom had during birth year of older sibling
        mean_pups_os <- sum_table_mom$mean_pups_yr[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == os_birth)] #Mean number of pups during birth year of older sibling
        mom_pups_ys <- sum_table_mom$pups[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == ys_birth)] #Observed number of pups from mom during birth year of younger sibling
        mean_pups_ys <- sum_table_mom$mean_pups_yr[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == ys_birth)] #Mean number of pups during birth year of younger sibling
          
      P_Mother_pos[os_birth, ys_birth, mommy] <- ((mom_pups_os / mean_pups_os) * 
        (mom_pups_ys / (mean_pups_ys * N_F[ys_birth-index]))) * (Surv^(ys_birth - os_birth)) #Store probability of kinship in P_Mother_pos
      
      } else P_Mother_pos[os_birth, ys_birth, mommy] <- 0 #If mother didn't have sampled offspring during both of these years, give a 0 to the associated cell (it won't be used, since negative comparisons use a different dataframe)
      }
    }
    }
  N_M=exp(Pars[(n_yrs+1):(n_yrs*2)])
  for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      for(daddy in 1:length(unique_dads)) {
        dad_i <- unique_dads[daddy] #dad being evaluated in present iteration
        
        #Did dad have pups in birth years of older and younger sibs being compared? 1 if yes.
        os_indicator <- length(sum_table_dad$pups[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == os_birth)]) 
        ys_indicator <- length(sum_table_dad$pups[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == ys_birth)]) 
        
        #If the Father had offspring in both the older sib and younger sib birth years, then store the 
        if(os_indicator == 1 & ys_indicator == 1) {
          dad_pups_os <- sum_table_dad$pups[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == os_birth)] #Observed number of pups dad had during birth year of older sibling
          mean_pups_os <- sum_table_dad$mean_pups_yr[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == os_birth)] #Mean number of pups during birth year of older sibling
          dad_pups_ys <- sum_table_dad$pups[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == ys_birth)] #Observed number of pups from dad during birth year of younger sibling
          mean_pups_ys <- sum_table_dad$mean_pups_yr[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == ys_birth)] #Mean number of pups during birth year of younger sibling
          
          P_Father_pos[os_birth, ys_birth, daddy] <- ((dad_pups_os / mean_pups_os) * 
                                                        (dad_pups_ys / (mean_pups_ys * N_M[ys_birth-index]))) * (Surv^(ys_birth - os_birth)) #Store probability of kinship in P_Father_pos
          
        } else P_Father_pos[os_birth,ys_birth, daddy] <- 0
    }
    }
  }
  return(list(P_Mother_pos=P_Mother_pos, P_Father_pos=P_Father_pos)) #return makes sure this is moved out of the loop into the environment
}
#P_pos=get_P_lemon_pos(Pars=Pars,P_Mother_pos=P_Mother_pos,P_Father_pos=P_Father_pos,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P_pos$P_Mother_pos[23:45,,21]
#unique_moms
#P$P_Father

#In case there are issues, check out below -- R was weird once and was registering the same values as different
#unique_moms[26]
#head(mom_positives)
#mom_positives[mom_positives$Mother=="FEMALE#051",]
#mom_positives
#test <- mom_positives[mom_positives$Mother=="FEMALE#051",]
#test
#test[1,] == test[2,]

get_P_lemon_neg <- function(Pars,P_Mother_neg,P_Father_neg,n_yrs,t_start,t_end){
  N_F=exp(Pars[1:n_yrs]) #number of mature females
  index <- min_est_cohort-1
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
  for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      
      #mean_pups_os <- sum_table_mom$mean_pups_yr[which(sum_table_mom$BirthY == os_birth)][1] #Mean number of pups during birth year of older sibling
      
      #mean_pups_ys <- sum_table_mom$mean_pups_yr[which(ssum_table_mom$BirthY == ys_birth)][1] #Mean number of pups during birth year of younger sibling
      
      P_Mother_neg[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F[ys_birth-index] #If expected values cancel out
      
      #P_Mother_neg[os_birth, ys_birth] <- ( / mean_pups_os) * (mom_pups_ys / (mean_pups_ys * N_F[ys_birth-index]))) * (Surv^(ys_birth - os_birth)) #Store probability of kinship in P_Father_pos
      
      #} else P_Mother[os_birth, ys_birth] <- 0
    }
  }
  N_M=exp(Pars[(n_yrs+1):(n_yrs*2)])
  #N_M=exp(Pars[(n_yrs+1):(n_yrs*2)]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      #if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
      
      #P_Father_neg[os_birth,ys_birth] <-((dad_pups_os / mean_pups_os) * (dad_pups_ys / (mean_pups_ys * N_F[ys_birth-index]))) * (Surv^(ys_birth - os_birth))
      P_Father_neg[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M[ys_birth-index] #If expected values cancel out
      
      #} else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother_neg=P_Mother_neg, P_Father_neg=P_Father_neg)) #return makes sure this is moved out of the loop into the environment
}
#P_neg=get_P_lemon_neg(Pars=Pars,P_Mother_neg=P_Mother_neg,P_Father_neg=P_Father_neg,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P_neg$P_Mother_neg
