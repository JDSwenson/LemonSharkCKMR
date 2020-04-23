#Datasets have the following columns: birth of parent, death of parent, birth of offspring, frequency of comparison
beluga_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,Maturity_age,n_yrs,t_start,t_end){
  P=get_P_beluga(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,Maturity_age=Maturity_age,n_yrs=n_yrs,t_start=t_start,t_end=t_end) #Define the function that generates the priors for the neg_log_lik function
  
  #Pars = VECTOR of initial parameter values 
    #c(mature females alive at t1, mature males alive at t1, male maturity age, slope for population trend)
  #Negatives_Mother = 4 column MATRIX of pairwise comparisons between mothers and young     that did NOT yield a POP
    #[,1] = mother birth year
    #[,2] = mother death year
    #[,3] = young birth year
    #[,4] = frequency with which that combination appears
  #Negatives_Father = same as above, but with Father
  #Pairs_Mother = 4 column MATRIX of pairwise comparisons between mothers and young     that DID yield a POP; columns same as above
  #Pairs_Father = same as above, but with Father
  #P_Mother = 3 dimensional ARRAY with prior probabilities of kinship, dimensions are: mother birth year, mother death year, and offspring birth year.
  #P_Father = same as above, but with father
  #Maturity_age = VECTOR with length n_ages
  #n_yrs = 60 (total time period examined)
  #t_start = 41 (first year of study)
  #t_end = 60 (last year of study)
  
  loglik=0 #set initial value of loglik to 0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){ #Loop through each row of Negative_Mother - each row is a unique combination of parent birth year, parent death year, and offspring birth year
    loglik=loglik+Negatives_Mother[irow,4]*log(1-P$P_Mother[Negatives_Mother[irow,1],Negatives_Mother[irow,2],Negatives_Mother[irow,3]]) #add to loglik the frequency for each combination times the log of 1 minus the prior probability of each combination (specified by finding the prior value in the associated cell of P_Mother). 1 minus the prior prob is the prior prob of NOT being kin. 
    #So, frequency of each occurence times the probability of the occurrence. Pretty straightforward.
  }
  for(irow in 1:nrow(Negatives_Father)){ #Same as above, but with Fathers
    loglik=loglik+Negatives_Father[irow,4]*log(1-P$P_Father[Negatives_Father[irow,1],Negatives_Father[irow,2],Negatives_Father[irow,3]])
  }  

    #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    loglik=loglik+Pairs_Mother[irow,4]*log(P$P_Mother[Pairs_Mother[irow,1],Pairs_Mother[irow,2],Pairs_Mother[irow,3]]) #Same as above, except use positive priors instead of negative
  }
  for(irow in 1:nrow(Pairs_Father)){#Same as above, but with Fathers
    loglik=loglik+Pairs_Father[irow,4]*log(P$P_Father[Pairs_Father[irow,1],Pairs_Father[irow,2],Pairs_Father[irow,3]])
  }  
  -loglik #Take the negative value of the log likelihood
}