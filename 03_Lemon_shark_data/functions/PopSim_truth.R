########## Compile and report results #########

#Calculate geometric mean of survival for each reference year
#See: https://www.r-bloggers.com/2021/08/calculate-geometric-mean-in-r/
#ref.vec <- ref.tibble %>% distinct(ref.yr) %>% pull(ref.yr)


  (surv_mean <- exp(mean(log(sVec[est.year.calibrate:estimation.year]))))
  
  #Arithemetic mean of survival
  (surv_mean_AM <- round(mean(sVec[ref.year:(n_yrs-1)]), 4))
  


##Equivalent way of calculating geometric mean
# sVec2 <- NULL
# 
# #Loop over the survival values for each year; note that survival here represents survival FROM the year it's recorded e..g sVec[2] represents survival from year 2 to year 3
# for(s in ref.value:(n_yrs - 1)){
#   sTemp <- sVec[s]
#   sVec2 <- c(sVec2, sTemp)
# }
# 
# delta <- n_yrs - ref.value
# 
# surv_mean <- (prod(sVec2))^(1/delta)

#surv_mean <- (sVec[n_yrs]/sVec[ref.value])^(1/(n_yrs - ref.value))
surv_min <- min(sVec[est.year.calibrate:(n_yrs-1)]) #Minimum survival over sample period
surv_max <- max(sVec[est.year.calibrate:(n_yrs-1)]) #Maximum survival over sample period


#psi
psi_truth.vec <- NULL
#Check proportion of conformists and non-conformists that are breeding each year
for(i in ref.year:n_yrs){
  Newmoms.vec <- loopy.list[[i]] %>% dplyr::filter(age.x == 0) %>% pull(mother.x)
  
  moms.all.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec)
  
  psi_truth.vec[i] <- moms.all.df %>% group_by(repro.strategy) %>% 
    summarize(n = n()) %>% 
    mutate(perc.conf = round(n/sum(n), 2)) %>% #percent conformists
    dplyr::filter(repro.strategy == "conformist") %>% 
    pull(perc.conf)
}

psi_truth <- round(mean(psi_truth.vec[ref.year:n_yrs]), 2)

#Add the calculated values to ref.tibble2
ref.tibble2 <- tibble(survival = surv_mean,
                      surv_AM = surv_mean_AM,
                      surv_min = surv_min,
                      surv_max = surv_max,
                      psi = psi_truth)

true.values <- ref.tibble2 %>% pivot_longer(cols = c(survival, psi), names_to = "parameter", values_to = "all.truth") %>% 
  mutate(iteration = iter,
         seed = rseed,
         estimation.year = estimation.year,
         T0 = est.year.calibrate)