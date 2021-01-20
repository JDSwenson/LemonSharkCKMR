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