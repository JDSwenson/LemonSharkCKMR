#Fills the arrays with probabilities of POP kinship based on ages. ideath1 starts at 41 because that is when the study starts. This means every time the loop completes, the first 40 columns are not given any values. Seems to produce the same values with line 8 is adjusted to for(ideath1 in t_start:t_end)
n_ages=40
get_P_beluga <- function(Pars,P_Mother,P_Father,Maturity_age,n_yrs, t_start, t_end){
  n_ages=length(Maturity_age)
  Log_trend = (0.01/(1+exp(-Pars[4]))-0.005)*(0:(t_end-1))  #trend on log scale limited to +/- 0.01/yr -- assuming population isn't growing
  N_F_yr=exp(Pars[1]+Log_trend)  #number of mature females (assume time constant) - (total reproductive output from females)
  for(ibirth1 in 1:n_yrs){  #n_yrs number of iterations, ibirth1 = mother's birth?
    for(ideath1 in max(ibirth1, t_start):t_end){ #ideath1 is mother's death year
      for(ibirth2 in 1:n_yrs){ #ibirth2 is birth year of offspring. Loop over all study years.
        if(ibirth2<=ideath1 & ibirth2>ibirth1 & ((ibirth2-ibirth1)<=n_ages)) #If offspring was born before mother died & after mother was born & if the mother was not older than the max age when the offspring was born, then ...
          P_Mother[ibirth1,ideath1,ibirth2]=Maturity_age[ibirth2-ibirth1]/N_F_yr[ibirth2] #Probability of kinship between a mother that was born in ibirth1 and died in ideath1 being the parent of offspring that was born in ibirth2 = the probability of fecundity for the age of the mother at the offspring's birth divided by the total number of females alive at the offspring's birth.
      }
    } 
  }
  N_M_yr=exp(Pars[2]+Log_trend) #number of mature males (time constant) - (total reproductive output from males)
  Paternity_age=rep(0,n_ages)
  age_50 = n_ages/(1+exp(-Pars[3]))  #logistic location parameter has range (0,n_ages)
  Paternity_age=1/(1+exp(-10*(c(1:n_ages)-age_50)))  #"almost" knife edged maturity for males (continuous logistic function to be kind to optimizer)
  #cat(Paternity_age)
  for(ibirth1 in 1:n_yrs){  #dad
    for(ideath1 in max(ibirth1, t_start):t_end){ 
      for(ibirth2 in 1:n_yrs){
        if(ibirth2<=ideath1 & ibirth2>ibirth1 & ((ibirth2-ibirth1)<=n_ages))
          P_Father[ibirth1,ideath1,ibirth2]=Paternity_age[ibirth2-ibirth1]/N_M_yr[ibirth2]
      }
    } 
  }
  return(list(P_Mother=P_Mother,P_Father=P_Father))
}