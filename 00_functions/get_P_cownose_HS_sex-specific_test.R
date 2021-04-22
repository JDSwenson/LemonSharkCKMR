#Sex-specific model for calculating prior probability of kinship for half-siblings based on birth years

#Set up empty array that will be filled with function below
P_Mother = P_Father = array(0,dim=c(n_yrs,n_yrs)) #Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

#CKMR model: populate array with kinship probabilities based on birth years
get_P_cownose <- function(Pars1,P_Mother,P_Father,t_start,t_end){

  N_F=exp(Pars1[1]) #number of mature females
  for(os_birth in 1:(n_yrs-1)){  #This loops over the possible ages of older sibs
    for(ys_birth in (os_birth+1):n_yrs){ #This loops over possible ages of younger sibs
      if((ys_birth - os_birth) <= (max.age - repro.age)){ #Make sure the difference in birth years allows for an animal to be mature in both years. 
        P_Mother[os_birth, ys_birth] <- (surv^(ys_birth - os_birth))/(N_F*lam^(ys_birth-year_of_est))
      } else P_Mother[os_birth, ys_birth] <- 0 #If the animal could not have been mature in the birth year of each individual, give a prob of 0

    }
  }
##Repeat the above with males

  N_M=exp(Pars1[2]) #number of mature males
  for(os_birth in 1:(n_yrs-1)){  #loop over possbile ages of older sib
    for(ys_birth in (os_birth+1):n_yrs){ #loop over possible ages of younger sib
      if((ys_birth - os_birth) <= ((max.age+1) - repro.age)){ #Make sure adult could have been mature in both birth years
        
        #Fill array with kinship probability from half-sib CKMR equation
        #See Bravington 2016 equation 3.10
        P_Father[os_birth,ys_birth] <- (surv^(ys_birth - os_birth))/(N_M*lam^(ys_birth-year_of_est))
      } else P_Father[os_birth,ys_birth] <- 0 #If the animal could not have been mature in the birth year of each individual, give a prob of 0
    }
  }
  
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}

#Assign function to P
P=get_P_cownose(Pars1=Pars1,P_Mother=P_Mother,P_Father=P_Father,t_start=t_start,t_end=t_end)

#Look at array
#View(P$P_Mother)
#View(P$P_Father)