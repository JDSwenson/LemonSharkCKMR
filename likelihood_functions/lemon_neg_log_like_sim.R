lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,n_yrs,t_start,t_end){
  
P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){
    loglik=loglik+Negatives_Mother[irow,4]*log(1-P$P_Mother[Negatives_Mother[irow,1],Negatives_Mother[irow,2],Negatives_Mother[irow,3]])
  }
  for(irow in 1:nrow(Negatives_Father)){
    loglik=loglik+Negatives_Father[irow,4]*log(1-P$P_Father[Negatives_Father[irow,1],Negatives_Father[irow,2],Negatives_Father[irow,3]])
  }  
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    loglik=loglik+Pairs_Mother[irow,4]*log(P$P_Mother[Pairs_Mother[irow,1],Pairs_Mother[irow,2],Pairs_Mother[irow,3]])
  }
  for(irow in 1:nrow(Pairs_Father)){
    loglik=loglik+Pairs_Father[irow,4]*log(P$P_Father[Pairs_Father[irow,1],Pairs_Father[irow,2],Pairs_Father[irow,3]])
  }  
  -loglik
}
