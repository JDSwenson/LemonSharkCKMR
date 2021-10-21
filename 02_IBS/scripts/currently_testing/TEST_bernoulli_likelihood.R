mom_comps <- mom_positives %>% 
  rename(yes = freq) %>% 
  full_join(mom_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(all = yes + no)

dad_comps <- dad_positives %>% 
  rename(yes = freq) %>% 
  full_join(dad_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(all = yes + no)


#Need to figure out how to connect the probability below with the likelihood function
get_P_lemon_bern <- function(Pars, M_comps, P_comps){
  M_prob <- (surv^(ys_birth - os_birth))/(N_F*lam^(ys_birth-min_cohort))
}




lemon_neg_log_lik_bern <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,t_start,t_end){
  
  P=get_P_lemon(Pars=Pars, M_comps = )
  
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