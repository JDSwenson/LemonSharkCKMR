#################### Breeding schedule ######################
if(b.schedule == "annual"){
  #------------------------------ Annual ------------------------------
  breeding.schedule <- "annual.breeding"
  mating.periodicity <- 1 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
  non.conformists <- 0
  ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
  ff
  print("breeding schedule is annual")
  
} else if(b.schedule == "biennial"){
  
  #------------------------------ Biennial ------------------------------
  mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
  print("breeding schedule is biennial")
  
  if(input.psi == 1){
    #============================== psi 1 ==============================
    breeding.schedule <- "biennial.breeding_psi1"
    non.conformists <- 0
    
    print("psi equals 1")
    
    if(offcycle.breeding == "yes"){
      
      breeding.schedule <- "biennial.breeding_psi1_offcycle"
      
      percent.breed_off.cycle <- 0.1 #Percent of off-cycle mothers that will breed each year
      percent.skip_on.cycle <- 0.1 #Percent of on-cycle mothers that will skip breeding each year
      
      print(paste0("There are ", percent.breed_off.cycle*100, "% of females that will breed off-cycle and ", percent.skip_on.cycle*100, "% of females that will fail when on-cycle."))
      
    }
  } else if(input.psi != 1){
    
    breeding.schedule <- paste0("biennial.breeding_psi", input.psi)
    non.conformists <- 1 - input.psi #proportion of off-year breeders to randomly include off their breeding cycle
    
    print(paste0("psi equals ", 1 - input.psi))
  }
  
} else if(b.schedule == "triennial"){
  
  mating.periodicity <- 3 #number of years between mating; assigned to an individual and sticks with them through their life.
  breeding.schedule <- "triennial.breeding_psi1"
  non.conformists <- 0
  
  print("breeding schedule is triennial")
}

if(b.schedule == "biennial" | b.schedule == "triennial"){
  # Adjust fecundity ==============================================================
  ## for effective number of breeders each year, mating cycle, number of mates ====
  #(from liz) ====
  psi <- 1-non.conformists
  ff <- mating.periodicity/(mating.periodicity-psi*mating.periodicity+psi)*f/init.prop.female/mean(num.mates)
  ff
}

#################### Population growth ####################
if(input.popGrowth == "stable"){
  #------------------------------Stable------------------------------
  population.growth <- "lambda.1"
  
  print("population growth is stable")
  
} else if(input.popGrowth == "slight increase"){
  
  #------------------------------Slight increase------------------------------
  population.growth <- "lambda.slight.increase"
  # ff.shift <- ff+0.5 #Increase fecundity to slightly increase population growth - only works for annual breeding
  
  (Ma <- -log(Adult.survival)) #Mortality of adults
  (Mj <- -log(juvenile.survival)) #Mortality of juveniles
  (Sa <- exp(-Ma)) #Survival of adults
  (Sj <- exp(-Mj)) #Survival of juveniles
  
  Mx <- -0.01 #Extra mortality (make negative so it reduces mortality)
  (Sa <- exp(-Ma - Mx)) #Survival of adults
  (Sj <- exp(-Mj - Mx)) #Survival of juveniles
  
  Adult.survival <- Sa
  juvenile.survival <- Sj
  
  cat(paste0("Adult survival is ", round(Adult.survival, 0), ";\nJuvenile survival is ", round(juvenile.survival,0), ";\nPopulation will increase in size"))
  
} else if(input.popGrowth == "slight decline"){
  
  #------------------------------Slight decline------------------------------
  population.growth <- "lambda.slight.decrease"
  # ff.shift <- ff-0.5 #Decrease fecundity to slightly decrease population growth - only works for annual breeding
  (Ma <- -log(Adult.survival)) #Mortality of adults
  (Mj <- -log(juvenile.survival)) #Mortality of juveniles
  (Sa <- exp(-Ma)) #Survival of adults
  (Sj <- exp(-Mj)) #Survival of juveniles
  
  Mx <- 0.01 #Extra mortality
  (Sa <- exp(-Ma - Mx)) #Survival of adults
  (Sj <- exp(-Mj - Mx)) #Survival of juveniles
  
  Adult.survival <- Sa
  juvenile.survival <- Sj
  
  cat(paste0("Adult survival is ", round(Adult.survival, 0), ";\nJuvenile survival is ", round(juvenile.survival,0),";\nPopulation will decrease in size"))
  
} else if(input.popGrowth == "severe decline"){
  
  #------------------------------Substantial decrease------------------------------
  population.growth <- "lambda.extreme" #Mortality will change later inside the loop
  
  cat(paste0("Adult survival is ", round(Adult.survival,0), ";\nJuvenile survival is ", round(juvenile.survival,0),";\nPopulation will decline rapidly"))
}