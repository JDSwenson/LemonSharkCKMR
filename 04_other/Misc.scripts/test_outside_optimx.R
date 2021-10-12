Negatives_Mother <- mom_negatives   
Negatives_Father <- dad_negatives
Pairs_Mother <- mom_positives
Pairs_Father <- dad_positives

values <- NULL
lik <- NULL

#Test parameter space for abundance and lambda
## Remove inner loop if only testing abundance
#Gives error when parameter value is 0
for(i in seq(50, 10000, by = 10)){
Pars[1] = Pars[2] <- log(i)

for(j in seq(0.95, 1.05, by = 0.005)){
  Pars[3] <- log(j)
  P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,t_start=t_start,t_end=t_end)

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
  
  print(i)
  print(j)
  print(-loglik)
  
  values[i] <- i
  lik[i] <- -loglik
}
}

min(lik)

head(values)
lik
sum(is.nan(lik))
df <- data.frame(cbind(values, lik))



#Test parameter space for lambda
for(i in seq(0.95, 1.05, by = 0.005)){
  Pars[3] <- log(i)
  
  P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,t_start=t_start,t_end=t_end)
  
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
  
  print(i)
  print(-loglik)
  
  values[i] <- i
  lik[i] <- -loglik
}
