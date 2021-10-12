#Work on skipped-breeding model
ys_birth <- 85
os_birth <- 73

mom_positives %>% mutate(x = Ind_1_birth %% 2)
  
if(ys_birth %% 2 == 0){
  print("Hi")
}

#Try running model outside of optimization algorithm
lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,t_start,t_end){
  
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
  -loglik
}

lemon_neg_log_lik(Pars=Pars,Negatives_Mother = mom_negatives,Negatives_Father = dad_negatives,Pairs_Mother = mom_positives_sub,Pairs_Father = dad_positives_sub,P_Mother = P_Mother,P_Father = P_Father,t_start = t_start,t_end = t_end)

mom_positives_sub <- mom_positives %>% sample_n(size = 5)
dad_positives_sub <- dad_positives %>% sample_n(size = 5)

#Check oldest and youngest parents in fishSim
mom_indiv <- indiv %>% select(Mum, BirthY)
dad_indiv <- indiv %>% select(Dad, BirthY)

mom_indiv %>% rename(Me = Mum, offspring_birth = BirthY) %>% 
  left_join(indiv, by = "Me") %>% 
  mutate(age_m = offspring_birth - BirthY) %>% 
  arrange(age_m) %>% 
  head(10)

dad_indiv %>% rename(Me = Dad, offspring_birth = BirthY) %>% 
  left_join(indiv, by = "Me") %>% 
  mutate(age_m = offspring_birth - BirthY) %>% 
  arrange(age_m) %>% 
  head(10)






#Check oldest and youngest parents
head(mom_positives)
head(indiv)
mom_positives %>% rename(Me = Mother) %>% 
  left_join(indiv, by = "Me") %>% 
  select(Me, Old_sib_birth, Young_sib_birth, BirthY) %>% 
  mutate(age_young = Young_sib_birth - BirthY, age_old = Old_sib_birth - BirthY) %>% 
  summarize(min_age = min(age_old), max_age = max(age_young))

dad_positives %>% rename(Me = Father) %>% 
  left_join(indiv, by = "Me") %>% 
  select(Me, Old_sib_birth, Young_sib_birth, BirthY) %>% 
  mutate(age_young = Young_sib_birth - BirthY, age_old = Old_sib_birth - BirthY) %>% 
  summarize(min_age = min(age_old), max_age = max(age_young))


#Tweak the function to see if it changes the estimate ... 
get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      #if((ys_birth - os_birth) <= (max_age - min(f_adult_age))){
        P_Mother[os_birth, ys_birth] <- (Surv^(ys_birth - os_birth))/N_F
      #} else P_Mother[os_birth, ys_birth] <- 0
    }
  }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      #if((ys_birth - os_birth) <= (max_age - min(m_adult_age))){
        P_Father[os_birth,ys_birth] <- (Surv^(ys_birth - os_birth))/N_M
      #} else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)


(mom_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,] %>% 
    select(Old_sib_birth, Young_sib_birth) %>% 
    filter(Young_sib_birth != Old_sib_birth) %>% 
    plyr::count() %>% 
    select(Old_sib_birth, Young_sib_birth, freq))



#Check to make sure truth calcualtions are correctly placed
Dad_truth_before <- c()
Mom_truth_before <- c()

Dad_truth_after <- c()
Mom_truth_after <- c()
for (y in 1:length(sim_yrs)) {
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "truncPoisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve)
  
  Dad_truth_before[y] <- indiv %>% 
    filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    summarize(n())
  
  Mom_truth_before[y] <- indiv %>% 
    filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    summarize(n())
  
  indiv <- mort(indiv, type = death_type, year=y, ageMort = ageMort, maxAge = maxAge)
  indiv <- birthdays(indiv)
  #Store true values for each year
  Dad_truth_after[y] <- indiv %>% 
    filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    summarize(n())
  
  Mom_truth_after[y] <- indiv %>% 
    filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    summarize(n())
}

Mom_truth_before <- unlist(Mom_truth_before) 
Dad_truth_before <- unlist(Dad_truth_before)
Mom_truth_after <- unlist(Mom_truth_after)
Dad_truth_after <- unlist(Dad_truth_after)

Mom_truth_before
Mom_truth_after

Dad_truth_before
Dad_truth_after
