#TRoubleshoot model

################### MOTHERS #################################
#Assign values to variables in loop for troubleshooting
os_birth = min_est_cohort
ys_birth = min_est_cohort+1
mommy=1
mom_i = unique_moms[mommy]

#Look at vector of unique moms
unique_moms

#Look at table of pups per mom per year
sum_table_mom

#Look at the row corresponding to the values assigned above
sum_table_mom[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$BirthY == os_birth),]
sum_table_mom[which(sum_table_mom$Mother == unique_moms[mommy] & sum_table_mom$BirthY == ys_birth),]

#Check specific values being assigned to variables
(mom_pups_os <- sum_table_mom$pups[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == os_birth)]) #Observed number of pups mom had during birth year of older sibling
(mean_pups_os <- sum_table_mom$mean_pups_yr[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == os_birth)]) #Mean number of pups during birth year of older sibling
(mom_pups_ys <- sum_table_mom$pups[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == ys_birth)]) #Observed number of pups from mom during birth year of younger sibling
(mean_pups_ys <- sum_table_mom$mean_pups_yr[which(sum_table_mom$Mother == mom_i & sum_table_mom$BirthY == ys_birth)]) #Mean number of pups during birth year of younger sibling

#Check that value being assigned based on equation/model is appropriate
(P_Mother_pos[os_birth, ys_birth, mommy] <- ((mom_pups_os / mean_pups_os) * 
                                              (mom_pups_ys / (mean_pups_ys * N_F[ys_birth-index]))) * (Surv^(ys_birth - os_birth)))


##############################  FATHERS ###########################################
#Assign values to variables in loop for troubleshooting
os_birth = min_est_cohort
ys_birth = min_est_cohort+1
daddy=1
dad_i = unique_dads[daddy]
index=min_est_cohort-1

#Look at vector of unique dads
unique_dads

#Look at table of pups per dad per year
sum_table_dad

#Look at the row corresponding to the values assigned above
sum_table_dad[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$BirthY == os_birth),]
sum_table_dad[which(sum_table_dad$Father == unique_dads[daddy] & sum_table_dad$BirthY == ys_birth),]

#Check specific values being assigned to variables
(dad_pups_os <- sum_table_dad$pups[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == os_birth)]) #Observed number of pups dad had during birth year of older sibling
(mean_pups_os <- sum_table_dad$mean_pups_yr[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == os_birth)]) #Mean number of pups during birth year of older sibling
(dad_pups_ys <- sum_table_dad$pups[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == ys_birth)]) #Observed number of pups from dad during birth year of younger sibling
(mean_pups_ys <- sum_table_dad$mean_pups_yr[which(sum_table_dad$Father == dad_i & sum_table_dad$BirthY == ys_birth)]) #Mean number of pups during birth year of younger sibling

#Check that value being assigned based on equation/model is appropriate
(P_Father_pos[os_birth, ys_birth, daddy] <- ((dad_pups_os / mean_pups_os) * 
                                               (dad_pups_ys / (mean_pups_ys * N_M[ys_birth-index]))) * (Surv^(ys_birth - os_birth)))

N_M=exp(Pars[(n_yrs+1):(n_yrs*2)])
N_F=exp(Pars[1:n_yrs])


#Model function
P_pos=get_P_lemon_pos(Pars=Pars,P_Mother_pos=P_Mother_pos,P_Father_pos=P_Father_pos,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

P_pos$P_Mother_pos[24:45,,1] #show years 24:45 (will stop at 23 by default)
P_pos$P_Father_pos[24:45,,1]
sum_table_dad[which(sum_table_dad$Father == unique_dads[daddy]),]



#TRoubleshoot likelihood function
#Dimensions are parent birth year, parent capture year, offspring birth year
P$P_Mother[10, 22, 23]


lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother_pos,P_Father_pos, P_Mother_neg, P_Father_neg, n_yrs,t_start,t_end){
  P_pos=get_P_lemon_pos(Pars=Pars,P_Mother_pos=P_Mother_pos,P_Father_pos=P_Father_pos,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  P_neg=get_P_lemon_neg(Pars=Pars,P_Mother_neg=P_Mother_neg,P_Father_neg=P_Father_neg,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  Negatives_Mother <- mom_negatives
  Negatives_Father <- dad_negatives
  Pairs_Mother <- mom_positives
  Pairs_Father <- dad_positives
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){
    loglik=loglik+Negatives_Mother[irow,3]*log(1-P_neg$P_Mother_neg[Negatives_Mother[irow,1],Negatives_Mother[irow,2]])
  } 
  loglik
  
  for(irow in 1:nrow(Negatives_Father)){
    loglik=loglik+Negatives_Father[irow,3]*log(1-P_neg$P_Father_neg[Negatives_Father[irow,1],Negatives_Father[irow,2]])
  }
  loglik
  
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    (loglik=loglik + Pairs_Mother[irow,3]*log(P_pos$P_Mother_pos[Pairs_Mother[irow,1],Pairs_Mother[irow,2], which(unique_moms==Pairs_Mother[irow, 4])]))
  }
  loglik
  for(irow in 1:nrow(Pairs_Father)){
    (loglik=loglik+Pairs_Father[irow,3]*log(P_pos$P_Father_pos[Pairs_Father[irow,1],Pairs_Father[irow,2], which(unique_dads==Pairs_Father[irow, 4])]))
  }  
  -loglik
}
irow=1
loglik=0

Pairs_Mother
P_pos$P_Mother_pos[32:45,,1]

unique_dads
P_pos$P_Father_pos[,,4]
