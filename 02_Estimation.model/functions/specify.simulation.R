##########Specify which model to use ###################################

specify.model <- function(){
noLambda_model.validation.scenarios <- c("scenario_1_model.validation")
noLambda_annual.model.scenarios <- c("scenario_2.1.1", "scenario_2.1.2", "scenario_2.1.3")
narrowLambda_annual.model.scenarios <- c("scenario_2.2.1", "scenario_2.2.2", "scenario_2.2.3", "scenario_3.1.1", "scenario_3.2.1", "scenario_3.3.1", "scenario_3.4.1", "scenario_3.5.1", "scenario_3.6.1", "scenario_4.4", "scenario_4.5", "scenario_4.6")
wideLambda_annual.model.scenarios <- c("scenario_2.3.1", "scenario_2.3.2", "scenario_2.3.3")
narrowLambda_skip.model.scenarios <- c("scenario_3.1.2", "scenario_3.2.2", "scenario_3.3.2", "scenario_3.4.2", "scenario_3.5.2", "scenario_3.6.2", "scenario_4.1", "scenario_4.2", "scenario_4.3")

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
  #------------------------- Set input file locations -------------------------#
  PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
  PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
  mating.periodicity <- 1 #number of years between mating for females
  
  #------------------------- Set output file locations -------------------------#
  MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.output/"
  results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.results/"
  
  
  #------------------------- Objective 1 common settings -------------------------#
  jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
  
  #Set up informed prior for survival
  survival.prior.mean <- adult.survival
  survival.prior.cv <- 0.05
  survival.prior.sd <- survival.prior.mean * survival.prior.cv
  scenario <- "scenario_1_model.validation" #For naming output files and calculating truth.
  model <- "annual.model" #For naming output files
  
  source("./02_Estimation.model/functions/Obj123.functions.R")
  
  cat(paste0("Testing model validation (objective 1)"))
  
} else if(objective ==2){
  ########################## Objective 2 #########################
  #------------------------- Set input file locations -------------------------#
  PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
  mating.periodicity <- 1 #number of years between mating for females
  
  #------------------------- Set output file locations -------------------------# 
  MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.output/"
  results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.results/"
  
  #------------------------- Objective 2 common settings -------------------------#
  cat(paste0("Testing population growth (Objective 2)"))
  model <- "annual.model" #For naming output files
  
  source("./02_Estimation.model/functions/Obj123.functions.R")
  
  if(scenario %in% c("scenario_2.1.1", "scenario_2.1.2", "scenario_2.1.3")){
    #========================= Scenario 2.1 =========================
    
    cat(paste0(" with a naive model"))
    
    if(scenario == "scenario_2.1.1"){
      #------------------------- Scenario 2.1.1: Small population decline; no lambda in model
      cat(paste0(" and a slightly decreasing population."))
      
      PopSim.lambda <- "lambda.slight.decrease"
      jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
      
      
    } else if(scenario == "scenario_2.1.2"){
      
      #------------------------- Scenario 2.1.2: Small population growth; no lambda in model
      cat(paste0(" and a slightly increasing population."))
      
      PopSim.lambda <- "lambda.slight.increase"
      jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated

      
    } else if(scenario == "scenario_2.1.3"){
      
      #------------------------- Scenario 2.1.3: Substantial population decline; no lambda in model
      cat(paste0(" and a severely decreasing population."))
      
      PopSim.lambda <- "lambda.extreme" 
      jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated

    }
  } else if(scenario %in% c("scenario_2.2.1", "scenario_2.2.2", "scenario_2.2.3")){
    
    #========================= Scenario 2.2 =========================
    cat(paste0(" with an adapted model and narrow prior"))
    
    if(scenario == "scenario_2.2.1"){
      #------------------------- Scenario 2.2.1: Small population decline; lambda in model w/ tight prior
      cat(paste0(" and a slightly decreasing population."))
      
      PopSim.lambda <- "lambda.slight.decrease"
      jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
      
    } else if(scenario == "scenario_2.2.2"){
      #------------------------- Scenario 2.2.2: Small population growth; lambda in model w/ tight prior
      cat(paste0(" and a slightly increasing population."))
      
      PopSim.lambda <- "lambda.slight.increase"
      jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
      
    } else if(scenario == "scenario_2.2.3"){
      
      #------------------------- Scenario 2.2.3: Substantial population decline; lambda in model w/ tight prior
      cat(paste0(" and a severely decreasing population."))
      
      PopSim.lambda <- "lambda.extreme" 
      jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
      
      
    }
  } else if(scenario %in% c("scenario_2.3.1", "scenario_2.3.2", "scenario_2.3.3")){
    
    #========================= Scenario 2.3 =========================
    cat(paste0(" with an adapted model and diffuse prior"))
    
    if(scenario == "scenario_2.3.1"){
      #------------------------- Scenario 2.3.1: Small population decline; lambda in model w/ wide prior
      cat(paste0(" and a slightly decreasing population."))
      
      PopSim.lambda <- "lambda.slight.decrease"
      jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
      
      
    } else if(scenario == "scenario_2.3.2"){
      
      #------------------------- Scenario 2.3.2: Small population growth; lambda in model w/ wide prior
      cat(paste0(" and a slightly increasing population."))
      
      PopSim.lambda <- "lambda.slight.increase"
      jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
      

      
    } else if(scenario == "scenario_2.3.3"){
      
      #------------------------- Scenario 2.3.3: Substantial population decline; lambda in model w/ wide prior
      cat(paste0(" and a severely decreasing population."))
      
      PopSim.lambda <- "lambda.extreme" 
      jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
      
    }
  }
  
} else if(objective == 3){
  ########################## Objective 3 #########################
  #------------------------- Set output file locations -------------------------#
  MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Nov2022/Model.output/"
  results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Nov2022/Model.results/"

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

  }
  
} else if(objective ==4){
  ########################## Objective 4 #########################
  #------------------------- Set output file locations ------------------------- 
  MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.output/"
  results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.results/"

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