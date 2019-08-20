n_ages=25 #Max age of lemon sharks
n_yrs=28
t_start=25
t_end=28
mat_age=rep(0, n_ages)
mat_age[12:n_ages]=1
P_Mother = P_Father = array(0,dim=c(n_yrs,n_yrs,n_yrs))  #creates two empty arrays, one for mother and one for father.  Dimensions are parent birth year, parent death year, offspring birth year (all of which are specified by n_yrs)

Pars=c(log(40),log(40)) #Based on the avg # juveniles in Bimini, not adults, and assumes equal sex ratio

get_P_lemon <- function(Pars,P_Mother,P_Father,n_yrs,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant) - (total reproductive output from females)
  for(ibirth1 in 1:n_yrs){  #
    for(ideath1 in max(ibirth1, t_start):t_end){ 
      for(ibirth2 in 1:n_yrs){
        if(ibirth2<=ideath1 & ibirth2>ibirth1 & ((ibirth2-ibirth1)<=n_ages))
          P_Mother[ibirth1,ideath1,ibirth2]=mat_age[ibirth2-ibirth1]/N_F
      }
    } 
  }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(ibirth1 in 1:n_yrs){  #dad
    for(ideath1 in max(ibirth1, t_start):t_end){ 
      for(ibirth2 in 1:n_yrs){
        if(ibirth2<=ideath1 & ibirth2>ibirth1 & ((ibirth2-ibirth1)<=n_ages))
          P_Father[ibirth1,ideath1,ibirth2]=mat_age[ibirth2-ibirth1]/N_M
      }
    } 
  }
  return(list(P_Mother=P_Mother,P_Father=P_Father))
}
