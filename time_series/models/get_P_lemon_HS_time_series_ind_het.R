#Model for calculating prior probability of kinship for half-siblings
#Estimate years 3:6 of data 5/19/20
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

P_Mother_pos <- array(NA,dim=c(max_est_cohort,max_est_cohort, length(unique_moms))) #Dimensions are older sib birth year, younger sib birth year, and number of unique moms
P_Father_pos <- array(NA,dim=c(max_est_cohort,max_est_cohort, length(unique_dads))) #Dimensions are older sib birth year, younger sib birth year, and number of unique moms
P_Mother_neg = P_Father_neg <- array(NA, dim=c(max_est_cohort, max_est_cohort))

get_P_lemon_pos <- function(Pars,P_Mother_pos,P_Father_pos,n_yrs,t_start,t_end){
  N_F=exp(Pars[1:n_yrs]) #number of mature females
  index <- min_est_cohort-1
    for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      for(mommy in 1: length(unique_moms)) {
      if(length(sum_table_mom$pups[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$Juv_study_DOB == os_birth)]) == 1 & length(sum_table_mom$pups[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$Juv_study_DOB == ys_birth)]) == 1){
      P_Mother_pos[os_birth, ys_birth, mommy] <- (sum_table_mom$pups[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$Juv_study_DOB == os_birth)] / 
                                                    sum_table_mom$mean_pups_yr[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$Juv_study_DOB == os_birth)]) * 
        (sum_table_mom$pups[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$Juv_study_DOB == ys_birth)] /
           (sum_table_mom$mean_pups_yr[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$Juv_study_DOB == ys_birth)] * N_F[ys_birth-index])) *
           (Surv^(ys_birth - os_birth)) #N_F is the number of females alive when the younger sibling was born
      } else P_Mother_pos[os_birth, ys_birth, mommy] <- 0
      }
    }
    }
  N_M=exp(Pars[(n_yrs+1):(n_yrs*2)])
  for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      for(daddy in 1: length(unique_dads)) {
        if(length(sum_table_dad$pups[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$Juv_study_DOB == os_birth)]) == 1 & length(sum_table_dad$pups[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$Juv_study_DOB == ys_birth)]) == 1){
        #if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
        P_Father_pos[os_birth, ys_birth, daddy] <- (sum_table_dad$pups[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$Juv_study_DOB == os_birth)] / 
                                                      sum_table_dad$mean_pups_yr[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$Juv_study_DOB == os_birth)]) * 
          (sum_table_dad$pups[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$Juv_study_DOB == ys_birth)] /
             (sum_table_dad$mean_pups_yr[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$Juv_study_DOB == ys_birth)] * N_M[ys_birth-index])) *
          (Surv^(ys_birth - os_birth))
        } else P_Father_pos[os_birth,ys_birth, daddy] <- NA
    }
    }
  }
  return(list(P_Mother_pos=P_Mother_pos, P_Father_pos=P_Father_pos)) #return makes sure this is moved out of the loop into the environment
}
#P_pos=get_P_lemon_pos(Pars=Pars,P_Mother_pos=P_Mother_pos,P_Father_pos=P_Father_pos,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P_pos$P_Mother_pos[2,4,28]
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
      #if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
      P_Mother_neg[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F[ys_birth-index] #N_F is the number of females alive when the younger sibling was born
      #} else P_Mother[os_birth, ys_birth] <- 0
    }
  }
  N_M=exp(Pars[(n_yrs+1):(n_yrs*2)])
  #N_M=exp(Pars[(n_yrs+1):(n_yrs*2)]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(max_est_cohort-1)){  #> = after, < = before
    for(ys_birth in max(min_est_cohort,(os_birth+1)):max_est_cohort){
      #if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
      P_Father_neg[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M[ys_birth-index]
      #} else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother_neg=P_Mother_neg, P_Father_neg=P_Father_neg)) #return makes sure this is moved out of the loop into the environment
}
#P_neg=get_P_lemon_neg(Pars=Pars,P_Mother_neg=P_Mother_neg,P_Father_neg=P_Father_neg,n_yrs=n_yrs,t_start=t_start,t_end=t_end)