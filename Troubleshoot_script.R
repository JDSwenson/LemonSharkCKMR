P$P_Father[36, 40]
ys_birth <- 90
os_birth <- 80


P_Mother = P_Father = array(0,dim=c(n_yrs,n_yrs))

surv <- 0.9
m_age_at_mat <- min(m_adult_age)
f_age_at_mat <- min(f_adult_age)

get_P_cownose <- function(Pars1,P_Mother,P_Father,t_start,t_end){
  N_F=exp(Pars1[1]) #number of mature females (assume time constant)
  
  for(os_birth in min_est_cohort:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge+1) - f_age_at_mat)){
      #Basic way of doing it - one pop_growth
      P_Mother[os_birth, ys_birth] <- (surv^(ys_birth - os_birth))/(N_F*lam^(ys_birth-min_est_cohort))
      } else P_Mother[os_birth, ys_birth] <- 0
      #More complicated way of doing it using observed pop_growth
      #   if(ys_birth == min_est_cohort_F){
      #   P_Mother[os_birth, ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/(N_F*pop_growth_F[min_est_cohort_F])
      # }
      # else  { #Multiply by the observed population growth rate
      #   NtempF <- N_F
      #   for(z in 1:length(min_est_cohort_F:(ys_birth-1))){
      #   popGF <- pop_growth_F[min_est_cohort_F:(ys_birth-1)]
      #   NtempF <- NtempF*popGF[z]
      #   }
      #   P_Mother[os_birth, ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/NtempF
      #   }
    }
  }
  N_M=exp(Pars1[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in min_est_cohort:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge+1) - m_age_at_mat)){
      #Basic way of doing it: one pop growth
      P_Father[os_birth,ys_birth] <- (surv^(ys_birth - os_birth))/(N_M*lam^(ys_birth-min_est_cohort))
      } else P_Father[os_birth,ys_birth] <- 0
      #More complicated way using observed pop growth
      # if(ys_birth == min_est_cohort_M){
      #   P_Father[os_birth,ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/(N_M*pop_growth_M[min_est_cohort_M])
      # 
      # }
      # else{ #Multiply by the observed population growth rate
      #   NtempM <- N_M
      #   for(z in 1:length(min_est_cohort_M:(ys_birth-1))){
      #     popGM <- pop_growth_M[min_est_cohort_M:(ys_birth-1)]
      #     NtempM <- NtempM*popGM[z]
      #   }
      #   P_Father[os_birth,ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/NtempM
      #}
    }
  }
  
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}

P=get_P_cownose(Pars1=Pars1,P_Mother=P_Mother,P_Father=P_Father,t_start=t_start,t_end=t_end)

#Model for calculating prior probability of kinship for half-siblings
#Not assessing survival - keeping constant at 0.85
Surv <- 0.85

max_age = 30 #max age of lemon sharks

adult_mat_age <- c(12:30) #Adult age of females -males assumed to be the same
adult_mat <- c(rep(0,11), rep(1,19)) #knife-edge maturity
adult_age_at_mat <- min(adult_mat_age)

P_Parent = array(0, dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

get_P_lemon_TotalA <- function(Pars2,P_Parent,t_start,t_end){
  N_A=exp(Pars2[1]) #number of mature females
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
  for(os_birth in min_est_cohort:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in max(os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge+1) - adult_age_at_mat)){
      P_Parent[os_birth, ys_birth] <- (4/(N_A*lam^(ys_birth - min_est_cohort)))*(Surv^(ys_birth - os_birth)) #N_F is the number of females alive when the younger sibling was born
    } else P_Parent[os_birth, ys_birth] <- 0
    }
  }
  return(list(P_Parent=P_Parent)) #return makes sure this is moved out of the loop into the environment
}

P_TotalA=get_P_lemon_TotalA(Pars2=Pars2,P_Parent,t_start=t_start,t_end=t_end)
#P_TotalA$P_Parent


View(P$P_Mother)



Negatives_Mother <- mom_negatives   
Negatives_Father <- dad_negatives
Pairs_Mother <- mom_positives
Pairs_Father <- dad_positives
    
lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,n_yrs,t_start,t_end){
      
      P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
      
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    pairs <- findRelativesPar(indiv = indiv, sampled = TRUE)
    POPs <- pairs[pairs$OneTwo == 1,1:2] ##Parent-Offspring pairs
    HSPs <- pairs[pairs$TwoTwo == 1,1:2] ##Half-Sibling pairs - verified from fishSim vignette
    non_POPs <- pairs[pairs$OneTwo == 0,1:2] ##pairs that are not POPs
    non_HSPs <- pairs[pairs$TwoTwo == 0,1:2] ##pairs that are not half-sibs
    
    # think of this dataframe as storing only information about the younger fish. Stores ID and birth year only (at present) for every individual in indiv, but renames the ID column to younger so it can be joined with HSPs_tbl below)
    youngerbirthyears <- indiv %>%
      select(Me, BirthY, Mum, Dad) %>% 
      rename("younger" = Me, "Young_sib_birth" = BirthY, "Young_sib_mom" = Mum, "Young_sib_dad" = Dad)
    
    #Create dataframe with IDs of sampled half-sibs and columns named appropriately for joins
    HSPs_tbl <- HSPs %>% 
      rename(
        Me = Var1,
        younger = Var2) %>%
      mutate_all(as.character)
    
    #Extract the values of indiv for the older individual in each comparison.
    #Join all the information from indiv with the IDs of the sampled individuals in HSPs_2_tbl. 
    HSPs_2_tbl <- inner_join(indiv, HSPs_tbl, by = "Me") %>%
      left_join(youngerbirthyears, by = "younger")  %>%
      rename("Old_sib_birth" = BirthY, "Old_sib_mom" = Mum, "Old_sib_dad" = Dad) %>%   select(c(Old_sib_birth, Young_sib_birth, Old_sib_mom, Young_sib_mom, Old_sib_dad, Young_sib_dad))
    
    #Split HSPs dataframe into MHS and PHS pairs
    #Filter out intra-cohort comparisons
    mom_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,] %>% 
      select(Old_sib_birth, Young_sib_birth) %>% 
      filter(Young_sib_birth != Old_sib_birth) %>% 
      plyr::count() %>% 
      select(Old_sib_birth, Young_sib_birth, freq)
    
    
    dad_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad,] %>% 
      select(Old_sib_birth, Young_sib_birth) %>% 
      filter(Young_sib_birth != Old_sib_birth) %>% 
      plyr::count() %>% 
      select(Old_sib_birth, Young_sib_birth, freq)
    
    min_est_cohort_F <- min(mom_positives$Young_sib_birth)
    min_est_cohort_M <- min(dad_positives$Young_sib_birth)
    min_est_cohort_all <- min(min_est_cohort_F, min_est_cohort_M)
    
    # Rename the columns so that the join functions know what to work with. The column `Me` is still the actual id column. However, we also need the column `younger` so that we can look up the 
    #birth years stored in `youngerbirthyears`
    non_HSPs_tbl <- non_HSPs %>% 
      rename(
        Me = Var1,
        younger = Var2) %>%
      mutate_all(as.character) # convert columns to characters since not all levels of indiv$Me are present in non_POPs_tbl
    #head(non_HSPs_tbl)
    
    ##Extract the values of indiv for the older individual in each comparison.
    #Join all the information from indiv with the IDs of the sampled individuals in non_HSPs_2_tbl. 
    non_HSPs_2_tbl <- inner_join(indiv, non_HSPs_tbl, by = "Me") %>%
      left_join(youngerbirthyears, by = "younger")  %>%
      rename(Old_sib_birth = BirthY) %>% 
      select(c(Old_sib_birth, Young_sib_birth))
    
    #Create table of all negative comparisons, count occurences, and filter intracohort comparisons
    mom_negatives = dad_negatives <- non_HSPs_2_tbl %>% 
      plyr::count() %>% 
      filter(Young_sib_birth != Old_sib_birth)
    
    #-------------Kinship probabilities - Half-sib-------------------
    #Not assessing survival - keeping constant at 0.85
    Surv <- 0.85
    max_age = 30 #max age of lemon sharks
    m_adult_age <- c(12:max_age) #Set ages at which males and females are mature.
    f_adult_age <- c(12:max_age)
    pop_growth_F_mean <- mean(pop_growth_F[min_est_cohort_F:length(pop_growth_F)])
    pop_growth_M_mean <- mean(pop_growth_M[min_est_cohort_M:length(pop_growth_M)])
    pop_growth_all_mean <- mean(pop_growth_all[min_est_cohort_all:length(pop_growth_all)])
    
    #Calculate mean survival-at-age
    mean_alive_at_age <- alive_at_age[1:30,] %>% rowMeans() %>% 
      round(digits = 0)
    #store mean adult survival
    mean_adult_surv <- mean(mean_alive_at_age[13:30]/mean_alive_at_age[12:29])
    
    
    P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)
    ### Maybe come back and use a truncated distribution instead (package truncnorm)
    
    get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
      N_F=exp(Pars[1]) #number of mature females (assume time constant)
      for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
        for(ys_birth in (os_birth+1):n_yrs){
          if(ys_birth == min_est_cohort_F){
            P_Mother[os_birth, ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/(N_F*pop_growth_F[min_est_cohort_F])
          }
          else  { #Multiply by the observed population growth rate
            NtempF <- N_F 
            for(z in 1:length(min_est_cohort_F:(ys_birth-1))){
              popGF <- pop_growth_F[min_est_cohort_F:(ys_birth-1)]
              NtempF <- NtempF*popGF[z]
            }
            P_Mother[os_birth, ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/NtempF
          }
        }
      }
      N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
      for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
        for(ys_birth in (os_birth+1):n_yrs){
          if(ys_birth == min_est_cohort_M){
            P_Father[os_birth,ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/(N_M*pop_growth_M[min_est_cohort_M])
            
          } 
          else{ #Multiply by the observed population growth rate
            NtempM <- N_M 
            for(z in 1:length(min_est_cohort_M:(ys_birth-1))){
              popGM <- pop_growth_M[min_est_cohort_M:(ys_birth-1)]
              NtempM <- NtempM*popGM[z]
            }
            P_Father[os_birth,ys_birth] <- (mean_adult_surv^(ys_birth - os_birth))/NtempM
          }
        }
      }
      return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
    }
  }
  
  P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  #------------------Likelihood function--------------------------
  lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,n_yrs,t_start,t_end){
    
    P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
    
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
  
  #Fit model
  
  #Fit model
  CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, n_yrs=n_yrs, t_start=t_start, t_end=t_end)
  
  #summary(CK_fit)
  #exp(CK_fit[1:2])
  
  #compute variance covariance matrix
  D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1])) #derivatives of transformations
  VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
  VC = (t(D)%*%VC_trans%*%D) #delta method
  SE=round(sqrt(diag(VC)),0)
  
  #Rename columns
  colnames(CK_fit)[1:2] <- c("N_est_F", "N_est_M")
  
  #Combine above to make dataframe with truth and estimates side-by-side
  #store years from youngest sibling in comparisons to end of study
  yrs <- c(min(mom_positives$Young_sib_birth, dad_positives$Young_sib_birth):t_end)
  
  estimates <- data.frame(cbind(t(round(exp(CK_fit[1:2]),0)), SE, c("F", "M")))
  colnames(estimates) <- c("CKMR_estimate", "SE", "Sex")
  (estimates <- cbind(estimates, truth = c(round(Mom_truth[min_est_cohort_F],0), round(Dad_truth[min_est_cohort_M],0))))
  
  #-----------------Loop end-----------------------------
  results[index,] <- c(round(exp(CK_fit[1]),0), round(SE[1],0), round(Mom_truth[min_est_cohort_F],0), "F", sum(mom_positives[,3]), pop_growth_all_mean, n_samples)
  results[index+1,] <- c(round(exp(CK_fit[2]),0), round(SE[2],0), round(Dad_truth[min_est_cohort_M],0), "M", sum(dad_positives[,3]), pop_growth_all_mean, n_samples)
  colnames(results) <- c("N_est", "SE", "Truth", "Sex", "Parents_detected", "Pop_growth", "Samples")
  
  ##Transpose proportion-at-age so ages are different columns
  age_pT <- data.frame(t(age_p[,-1]))
  colnames(age_pT) <- paste0("Age_", age_p[,1])
  age_pT$iter <- iter
  
  #Store proportion-at-age for each year of the fishSim simulation for each iteration of the loop.
  age_dist <- rbind(age_dist, age_pT)
  
  
  ##calculate survival for each age
  alive_at_age2 <- data.frame(alive_at_age[-31,-1])
  
  obs_survival <- data.frame()
  for(i in 1:(nrow(alive_at_age2)-1)){
    for(j in 1:(ncol(alive_at_age2)-1)){
      obs_survival[i, j] <- round(alive_at_age2[i+1, j+1]/alive_at_age2[i,j], digits=2)
    }
  }
  
  obs_survivalT <- data.frame(t(obs_survival))
  colnames(obs_survivalT) <- paste0("Age_", c(1:29))
  obs_survivalT$iter <- iter
  survival_at_age <- rbind(survival_at_age, obs_survivalT)
  rownames(survival_at_age) <- paste0("sim_yr_", c(1:n_yrs))
  
  print(paste0("finished iteration", iter, " at: ", Sys.time()))
}
results <- results %>% 
  mutate(Relative_bias = round(((N_est - Truth)/Truth)*100,1))

results %>% group_by(Sex) %>% 
  summarize(median = median(Relative_bias), n = n())

write.table(results, file = paste0("results/sex.combined.popgrowth_sex.specific.N", n_samples, "_samples.01.12.2021.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

write.table(age_dist, file = paste0("results/age_distributions_", n_samples, "_samples.01.12.2021.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

write.table(survival_at_age, file = paste0("results/survival_at_age_", n_samples, "_samples.01.14.2021.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}








#TRoubleshoot model
os_birth <- 52
ys_birth <- 56

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(f_adult_age)+1)){
        P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/(N_F*pop_growth^(ys_birth-min_est_cohort_F))
      } else P_Mother[os_birth, ys_birth] <- 0
    }
  }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= (max_age - min(m_adult_age)+1)){
        P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/(N_M*pop_growth^(ys_birth-min_est_cohort_M))
      } else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
P$P_Mother[41,59]

#------------------Likelihood function--------------------------
lemon_neg_log_lik <- function(Pars,mom_negatives,dad_negatives,mom_positives,dad_positives,P_Mother,P_Father,n_yrs,t_start,t_end){
  
  P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(mom_negatives)){
    loglik=loglik+mom_negatives[irow,3]*log(1-P$P_Mother[mom_negatives[irow,1],mom_negatives[irow,2]])
  }
  loglik
  
  for(irow in 1:nrow(dad_negatives)){
    loglik=loglik+dad_negatives[irow,3]*log(1-P$P_Father[dad_negatives[irow,1],dad_negatives[irow,2]])
  loglik
    }  
  loglik
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(mom_positives)){
    loglik=loglik+mom_positives[irow,3]*log(P$P_Mother[mom_positives[irow,1],mom_positives[irow,2]])
  }
  loglik
  for(irow in 1:nrow(dad_positives)){
    loglik=loglik+dad_positives[irow,3]*log(P$P_Father[dad_positives[irow,1],dad_positives[irow,2]])
  }  
  -loglik
}




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

###################Total Abundance ##########################
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
