P_Mother_HS = P_Father_HS = array(NA,dim=c(n_yrs,n_yrs))
P_Mother_PO = P_Father_PO = array(0,dim=c(n_yrs,n_yrs,n_yrs))  #creates two empty arrays, one for mother and one for father.  Dimensions are parent birth year, parent capture year, offspring birth year (all of which are specified by n_yrs)

f_mat_age <- which(mat_f == 1)
m_mat_age <- which(mat_m == 1)

surv <- 0.85

get_P_lemon <- function(Pars1,P_Mother_HS,P_Father_HS, P_Mother_PO, P_Father_PO){
  #Females
  N_F=exp(Pars1[1]) #number of mature females
  
  for(older_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(younger_birth in (older_birth+1):n_yrs){

      #Half-sib model      
      #Basic way of doing it - one pop_growth
      P_Mother_HS[older_birth, younger_birth] <- (surv^(younger_birth - older_birth))/(N_F*lam^(younger_birth-min_est_cohort))
      
      #Parent-offspring model
      for(older_capt in max(older_birth, t_start):t_end){
        if(younger_birth <= older_capt & younger_birth > older_birth & ((older_capt - older_birth) <= max_age) & ((younger_birth - older_birth) <= max_age)){
          P_Mother_PO[older_birth, older_capt, younger_birth] <- mat_f[younger_birth - older_birth]/N_F #If juvenile was born before the adult was captured AND after the adult was born AND the adult wasn't too old to have been captured OR to have given birth to the juvenile, then probability equals f_mat at the age the adult was when the juvenile was born over N_f
        } else if(younger_birth > older_capt & younger_birth - older_birth >= min(f_mat_age) & ((older_capt - older_birth) <= max_age) & ((younger_birth - older_birth) <= max_age)) {
          P_Mother_PO[older_birth, older_capt, younger_birth] <- (surv^(younger_birth - older_capt))/N_F
        }
        }
    }
  }
  
  #Males
  N_M=exp(Pars1[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(older_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(younger_birth in (older_birth+1):n_yrs){
      
      #Half-sib model
      #Basic way of doing it: one pop growth
      P_Father_HS[older_birth,younger_birth] <- (surv^(younger_birth - older_birth))/(N_M*lam^(younger_birth-min_est_cohort))
      
      #Parent-offspring model
      for(older_capt in max(older_birth, t_start):t_end){
        if(younger_birth <= older_capt & younger_birth > older_birth & ((older_capt - older_birth) <= max_age) & ((younger_birth - older_birth) <= max_age)){
          P_Father_PO[older_birth, older_capt, younger_birth] <- mat_m[younger_birth - older_birth]/N_F #If juvenile was born before the adult was captured AND after the adult was born AND the adult wasn't too old to have been captured OR to have given birth to the juvenile, then probability equals f_mat at the age the adult was when the juvenile was born over N_f
        } else if(younger_birth > older_capt & younger_birth - older_birth >= min(m_mat_age) & ((older_capt - older_birth) <= max_age) & ((younger_birth - older_birth) <= max_age)) {
          P_Father_PO[older_birth, older_capt, younger_birth] <- (surv^(younger_birth - older_capt))/N_F
      }
    }
    }
  }
  
  return(list(P_Mother_HS=P_Mother_HS, P_Father_HS=P_Father_HS, P_Mother_PO = P_Mother_PO, P_Father_PO = P_Father_PO)) #return makes sure this is moved out of the loop into the environment
}

P=get_P_lemon(Pars1=Pars1, P_Mother_HS=P_Mother_HS, P_Father_HS=P_Father_HS, P_Mother_PO=P_Mother_PO, P_Father_PO=P_Father_PO)