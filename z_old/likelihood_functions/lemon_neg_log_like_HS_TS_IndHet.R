#P_pos$P_Mother_pos[mom_positives[110,1],mom_positives[110,2], which(unique_moms==mom_positives[110, 4])]

#P_pos$P_Mother_pos
#head(mom_positives)
#unique_moms
#dim(P_Mother_pos)
#which(unique_moms==mom_positives[110, 4])
#unique_moms

lemon_neg_log_lik <- function(Pars, Negatives_Mother, Negatives_Father, Pairs_Mother, Pairs_Father, P_Mother_pos, P_Father_pos, P_Mother_neg, P_Father_neg, n_yrs, t_start, t_end){
  
  P_pos=get_P_lemon_pos(Pars=Pars,P_Mother_pos=P_Mother_pos,P_Father_pos=P_Father_pos,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  P_neg=get_P_lemon_neg(Pars=Pars,P_Mother_neg=P_Mother_neg,P_Father_neg=P_Father_neg,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){
    loglik=loglik+Negatives_Mother[irow,3]*log(1-P_neg$P_Mother_neg[Negatives_Mother[irow,1],Negatives_Mother[irow,2]])
  } 
  for(irow in 1:nrow(Negatives_Father)){
    loglik=loglik+Negatives_Father[irow,3]*log(1-P_neg$P_Father_neg[Negatives_Father[irow,1],Negatives_Father[irow,2]])
  }  
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    loglik=loglik + Pairs_Mother[irow,3]*log(P_pos$P_Mother_pos[Pairs_Mother[irow,1],Pairs_Mother[irow,2], which(unique_moms==Pairs_Mother[irow, 4])])
  }
  for(irow in 1:nrow(Pairs_Father)){
    loglik=loglik+Pairs_Father[irow,3]*log(P_pos$P_Father_pos[Pairs_Father[irow,1],Pairs_Father[irow,2], which(unique_dads==Pairs_Father[irow, 4])])
  }  
  -loglik
}
