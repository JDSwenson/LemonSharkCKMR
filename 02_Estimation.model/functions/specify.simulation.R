##########Specify which model to use ###################################

specify.model <- function(){
noLambda_model.validation.scenarios <- c("scenario_1_model.validation") #Includes a tight prior for survival; no lambda
noLambda_annual.model.scenarios <- c("scenario_1.2.1", "scenario_1.2.2", "scenario_1.2.3", "scenario_2.1.1", "scenario_2.1.2", "scenario_2.1.3", "scenario_2.1.4") #Includes a uniform prior for survival; no lambda
narrowLambda_annual.model.scenarios <- c("scenario_2.2.1", "scenario_2.2.2", "scenario_2.2.3", "scenario_2.2.4", "scenario_3.1.1", "scenario_3.2.1", "scenario_3.3.1", "scenario_3.4.1", "scenario_3.5.1", "scenario_3.6.1", "scenario_3.7.1", "scenario_4.4", "scenario_4.5", "scenario_4.6")
wideLambda_annual.model.scenarios <- c("scenario_2.3.1", "scenario_2.3.2", "scenario_2.3.3", "scenario_2.3.4")
narrowLambda_skip.model.scenarios <- c("scenario_3.1.2", "scenario_3.2.2", "scenario_3.3.2", "scenario_3.4.2", "scenario_3.5.2", "scenario_3.6.2", "scenario_3.7.2", "scenario_4.1", "scenario_4.2", "scenario_4.3")

if(scenario %in% noLambda_model.validation.scenarios){
  jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation.txt") #Annual model w/o lambda but with informed prior on survival
} else if(scenario %in% noLambda_annual.model.scenarios){
  jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model.txt") #Annual JAGS model w/o lambda and with diffuse prior on survival
} else if(scenario %in% narrowLambda_annual.model.scenarios){
  jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
} else if(scenario %in% wideLambda_annual.model.scenarios){
  jags_file = paste0(jags.model_location, "HS.PO_wideLambda_annual_model.txt")
} else if(scenario %in% narrowLambda_skip.model.scenarios){
  jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") 
}

return(jags_file)
}

if(objective == 1){
  ########################## Objective 1 #########################
  #------------------------- Set input file parameters from population simulation -------------------------#
  PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
  PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
  mating.periodicity <- 1 #number of years between mating for females
  
  #------------------------- Objective 1 common settings -------------------------#
  jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
  source("./02_Estimation.model/functions/Obj123.functions.R")
  model <- "annual.model" #For naming output files
  

  if(scenario %in% c("scenario_1_model.validation")){
    #========================= Scenario 2.1 =========================
    
  #Set up informed prior for survival
  survival.prior.mean <- adult.survival
  survival.prior.cv <- 0.05
  survival.prior.sd <- survival.prior.mean * survival.prior.cv
  scenario <- "scenario_1_model.validation" #For naming output files and calculating truth.
  
  } else if(scenario %in% c("scenario_1.2.1", "scenario_1.2.2", "scenario_1.2.3")){
    
    cat(paste0("Testing base-case CKMR model with diffuse survival prior"))
    
    if(scenario == "scenario_1.2.2"){
      
      cat(paste0(" and include aunt|uncle niece|nephew relationships as imposter HSPs"))
      test.decoys <- "yes" #Do we include aunt|uncle niece|nephew pairs as HSPs? Only include as "yes" if explicitly testing; otherwise blank or no.
      filter.decoys <- "no" #Do we filter aunt|uncle niece|nephew pairs based on a year-gap criterion? Should only be "yes" when test.decoys is also "yes". Will also require a threshold to be set below.
      
    
  } else if(scenario == "scenario_1.2.3"){
    
    test.decoys <- "yes" #Do we include aunt|uncle niece|nephew pairs as HSPs? Only include as "yes" if explicitly testing; otherwise blank or no.
    filter.decoys <- "yes" #Do we filter aunt|uncle niece|nephew pairs based on a year-gap criterion? Should only be "yes" when test.decoys is also "yes". Will also require a threshold to be set below.
    year.gap.threshold <- repro.age #Set a maximum number of years HSPs can be separated by to be included as HSPs. Goal is to have an empirical threshold for removing potential aunt|uncle niece|nephew pairs. 
    
    cat(paste0(" and include aunt|uncle niece|nephew relationships as imposter HSPs but filter by age gap: cannot have more than ", year.gap.threshold, " between birth years to be included as HSPs."))
    
    
  }
  }
  
} else if(objective ==2){
  ########################## Objective 2 #########################
  #------------------------- Set input file locations -------------------------#
  PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
  mating.periodicity <- 1 #number of years between mating for females

  #------------------------- Objective 2 common settings -------------------------#
  cat(paste0("Testing population growth (Objective 2)"))
  model <- "annual.model" #For naming output files
  
  source("./02_Estimation.model/functions/Obj123.functions.R")
  
  if(scenario %in% c("scenario_2.1.1", "scenario_2.1.2", "scenario_2.1.3", "scenario_2.1.4")){
    #========================= Scenario 2.1 =========================
    
    cat(paste0(" with a naive model"))
    jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
    
    if(scenario == "scenario_2.1.1"){
      #------------------------- Scenario 2.1.1: Small population decline; no lambda in model
      cat(paste0(" and a slightly decreasing population."))
      
      PopSim.lambda <- "lambda.slight.decrease"

      
    } else if(scenario == "scenario_2.1.2"){
      
      #------------------------- Scenario 2.1.2: Small population growth; no lambda in model
      cat(paste0(" and a slightly increasing population."))
      
      PopSim.lambda <- "lambda.slight.increase"
      
      
    } else if(scenario == "scenario_2.1.3"){
      
      #------------------------- Scenario 2.1.3: Substantial population decline; no lambda in model
      cat(paste0(" and a severely decreasing population."))
      
      PopSim.lambda <- "lambda.extreme" 

    } else if(scenario == "scenario_2.1.4"){
      
      #------------------------- Scenario 2.1.3: Substantial population decline; no lambda in model
      cat(paste0(" and a stable population."))

      PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
      
    }
  } else if(scenario %in% c("scenario_2.2.1", "scenario_2.2.2", "scenario_2.2.3", "scenario_2.2.4")){
    
    #========================= Scenario 2.2 =========================
    cat(paste0(" with an adapted model and narrow prior"))
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
    if(scenario == "scenario_2.2.1"){
      #------------------------- Scenario 2.2.1: Small population decline; lambda in model w/ tight prior
      cat(paste0(" and a slightly decreasing population."))
      
      PopSim.lambda <- "lambda.slight.decrease"
      
      
    } else if(scenario == "scenario_2.2.2"){
      #------------------------- Scenario 2.2.2: Small population growth; lambda in model w/ tight prior
      cat(paste0(" and a slightly increasing population."))
      
      PopSim.lambda <- "lambda.slight.increase"
      
      
    } else if(scenario == "scenario_2.2.3"){
      
      #------------------------- Scenario 2.2.3: Substantial population decline; lambda in model w/ tight prior
      cat(paste0(" and a severely decreasing population."))
      
      PopSim.lambda <- "lambda.extreme" 

      
    } else if(scenario == "scenario_2.2.4"){
      
      #------------------------- Scenario 2.2.3: Substantial population decline; lambda in model w/ tight prior
      cat(paste0(" and a stable population."))
      
      PopSim.lambda <- "lambda.1" 

            
    }
  } else if(scenario %in% c("scenario_2.3.1", "scenario_2.3.2", "scenario_2.3.3", "scenario_2.3.4")){
    
    
    #========================= Scenario 2.3 =========================
    cat(paste0(" with an adapted model and diffuse prior"))
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
    if(scenario == "scenario_2.3.1"){
      #------------------------- Scenario 2.3.1: Small population decline; lambda in model w/ wide prior
      cat(paste0(" and a slightly decreasing population."))
      
      PopSim.lambda <- "lambda.slight.decrease"
      
      
      
    } else if(scenario == "scenario_2.3.2"){
      
      #------------------------- Scenario 2.3.2: Small population growth; lambda in model w/ wide prior
      cat(paste0(" and a slightly increasing population."))
      
      PopSim.lambda <- "lambda.slight.increase"
      
      

      
    } else if(scenario == "scenario_2.3.3"){
      
      #------------------------- Scenario 2.3.3: Substantial population decline; lambda in model w/ wide prior
      cat(paste0(" and a severely decreasing population."))
      
      PopSim.lambda <- "lambda.extreme" 
      
      
    } else if(scenario == "scenario_2.3.4"){
      
      #------------------------- Scenario 2.3.3: Substantial population decline; lambda in model w/ wide prior
      cat(paste0(" and a stable population."))
      
      PopSim.lambda <- "lambda.1" 
      
      
    }
  }
  
} else if(objective == 3){
  ########################## Objective 3 #########################

  #------------------------- Objective 3 common settings -------------------------#
  PopSim.lambda <- "lambda.1"
  source("./02_Estimation.model/functions/Obj123.functions.R") #Changed name of script that includes pairwise comparison and other functions
  mating.periodicity <- 2 #Overwrite below if triennial breeding
  
  cat(paste0("Testing multiennial breeding (Objective 3)"))
  
  if(scenario == "scenario_3.1.1"){
    #========================= Scenario 3.1 =========================
    #------------------------- Scenario 3.1.1: Biennial breeding; psi = 1; annual model
    cat(paste0(" with 100% biennial breeders and a naive model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated

    
  } else if(scenario == "scenario_3.1.2"){
    #------------------------- Scenario 3.1.2: Biennial breeding; psi = 1; multiennial model
    cat(paste0(" with 100% biennial breeders and a multiennial model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated

    
  } else if(scenario == "scenario_3.2.1"){
    #========================= Scenario 3.2 =========================
    #------------------------- Scenario 3.2.1: Biennial breeding; psi = 0.9; annual model
    cat(paste0(" with 90% biennial breeders and a naive model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
  } else if(scenario == "scenario_3.2.2"){
    
    #------------------------- Scenario 3.2.2: Biennial breeding; psi = 0.9; multiennial model
    cat(paste0(" with 90% biennial breeders and a multiennial model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
  } else if(scenario == "scenario_3.3.1"){
    #========================= Scenario 3.3 =========================
    #------------------------- Scenario 3.3.1: Biennial breeding; psi = 0.75; annual model
    cat(paste0(" with 75% biennial breeders and a naive model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
  } else if(scenario == "scenario_3.3.2"){
    #------------------------- Scenario 3.3.2: Biennial breeding; psi = 0.75; multiennial model
    cat(paste0(" with 75% biennial breeders and a multiennial model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
  } else if(scenario == "scenario_3.4.1"){
    
    #========================= Scenario 3.4 =========================
    #------------------------- Scenario 3.4.1: Biennial breeding; psi = 0.50; annual model
    cat(paste0(" with 50% biennial breeders and a naive model."))

    PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
  } else if(scenario == "scenario_3.4.2"){
    
    #------------------------- Scenario 3.4.2: Biennial breeding; psi = 0.50; multiennial model
    cat(paste0(" with 50% biennial breeders and a multiennial model."))

    PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated

    
  } else if(scenario == "scenario_3.5.1"){
    #========================= Scenario 3.5: Triennial breeding =========================
    #------------------------- Scenario 3.5.1: Triennial breeding; psi = 1; annual model
    cat(paste0(" with 100% triennial breeders and a naive model."))
    
    PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    mating.periodicity <- 3

  } else if(scenario == "scenario_3.5.2"){
    
    #------------------------- Scenario 3.5.2: Triennial breeding; psi = 1; multiennial model
    cat(paste0(" with 100% triennial breeders and a multiennial model."))
    
    PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
    mating.periodicity <- 3
    
    
  } else if(scenario == "scenario_3.6.1"){
    
    #========================= Scenario 3.6: Biennial breeding w/ stochasticity
    #------------------------- Scenario 3.6.1: Biennial breeding; psi = 1; annual model; 10% on-cycle breeders fail to breed; 10% off-cycle breeders do breed.
    cat(paste0(" with stochastic biennial breeders and a naive model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi1_offcycle" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
    
  } else if(scenario == "scenario_3.6.2"){
    #------------------------- Scenario 3.6.2: Biennial breeding; psi = 1; multiennial model; 10% on-cycle breeders fail to breed; 10% off-cycle breeders do breed.
    cat(paste0(" with stochastic biennial breeders and a multiennial model."))
    
    PopSim.breeding.schedule <- "biennial.breeding_psi1_offcycle" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated

  } else if(scenario == "scenario_3.7.1"){
    #========================= Scenario 3.7: Annual breeding =========================
    #------------------------- Scenario 3.7.1: Annual breeding; psi = 1; annual model
    cat(paste0(" with 100% annual breeders and an annual model."))
    
    PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    mating.periodicity <- 1
    
  } else if(scenario == "scenario_3.7.2"){
    
    #------------------------- Scenario 3.7.2: Annual breeding; psi = 1; multiennial model
    cat(paste0(" with 100% annual breeders and a multiennial model."))
    
    PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
    mating.periodicity <- 1
    
  }
  
} else if(objective ==4){
  ########################## Objective 4 #########################

  #------------------------- Objective 4 common settings -------------------------
  cat(paste0("Testing aging error (Objective 4)"))
  
  PopSim.lambda <- "lambda.1"
  source("./02_Estimation.model/functions/Obj4.functions.R") #Changed name of script that includes pairwise comparison and other functions
  
  
  if(scenario == "scenario_4.1"){
    
    #------------------------- Scenario 4.1: Biennial breeding; psi = 1; aging error 5% CV
    cat(paste0(" with biennial breeders and 5% CV."))
    
    age.cv <- 0.05
    mating.periodicity <- 2
    PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated

    
  }else if(scenario == "scenario_4.2"){
    #------------------------- Scenario 4.2: Biennial breeding; psi = 1; aging error 10% CV
    cat(paste0(" with biennial breeders and 10% CV."))
    
    age.cv <- 0.10
    mating.periodicity <- 2
    PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated

    
  }else if(scenario == "scenario_4.3"){
    
    #------------------------- Scenario 4.3: Biennial breeding; psi = 1; aging error 20% CV
    cat(paste0(" with biennial breeders and 20% CV."))
    
    age.cv <- 0.20
    mating.periodicity <- 2
    PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
    model <- "multiennial.model"
    jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
    
  }else if(scenario == "scenario_4.4"){
    
    #------------------------- Scenario 4.4: Annual breeding; aging error 5% CV
    cat(paste0(" with annual breeders and 5% CV."))
    
    age.cv <- 0.05
    mating.periodicity <- 1
    PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
  }else if(scenario == "scenario_4.5"){
    
    #------------------------- Scenario 4.5: Biennial breeding; psi = 1; aging error 10% CV
    cat(paste0(" with annual breeders and 10% CV."))
    
    age.cv <- 0.10
    mating.periodicity <- 1
    PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated

  }else if(scenario == "scenario_4.6"){
    
    #------------------------- Scenario 4.6: Biennial breeding; psi = 1; aging error 20% CV
    cat(paste0(" with annual breeders and 20% CV."))
    
    age.cv <- 0.20
    mating.periodicity <- 1
    PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
    model <- "annual.model"
    jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
    
  }
}