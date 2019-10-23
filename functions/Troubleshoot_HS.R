lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,n_yrs,t_start,t_end){
  
P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

(Negatives_Father <- Data_dad_no)
(Negatives_Mother <- Data_mom_no)
(Pairs_Mother <- Data_mom_yes)
(Pairs_Father <- Data_dad_yes)
P$P_Mother[25,29]

head(Juv_lems)

  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){
    loglik=loglik+Negatives_Mother[irow,3]*log(1-P$P_Mother[Negatives_Mother[irow,1],Negatives_Mother[irow,2]])
  }
  loglik
  for(irow in 1:nrow(Negatives_Father)){
    loglik=loglik+Negatives_Father[irow,3]*log(1-P$P_Father[Negatives_Father[irow,1],Negatives_Father[irow,2]])
  }
  loglik
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    loglik=loglik+Pairs_Mother[irow,3]*log(P$P_Mother[Pairs_Mother[irow,1],Pairs_Mother[irow,2]])
  }
  loglik
  for(irow in 1:nrow(Pairs_Father)){
    loglik=loglik+Pairs_Father[irow,3]*log(P$P_Father[Pairs_Father[irow,1],Pairs_Father[irow,2]])
  }  
  -loglik
}

lemon_neg_log_lik(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,n_yrs,t_start,t_end)

test_mat <- array(0, dim=c(10,10))
for(i in 1:10){
  for(j in max(i, 6):10){
    test_mat[i,j] <- j*i
  }
}
test_mat
