lemon_neg_log_lik <- function(Pars, Negatives_Parent, Pairs_Parent, P_Parent, n_yrs, t_start, t_end){
  
  P=get_P_lemon(Pars=Pars,P_Parent,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Parent)){
    loglik=loglik+Negatives_Parent[irow,3]*log(1-P$P_Parent[Negatives_Parent[irow,1],Negatives_Parent[irow,2]])
  } 
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Parent)){
    loglik=loglik+Pairs_Parent[irow,3]*log(P$P_Parent[Pairs_Parent[irow,1],Pairs_Parent[irow,2]])
  }
  -loglik
}
