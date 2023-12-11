#rm(list=ls())
library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(MASS)
library(popbio)
library(mpmtools)
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(runjags)
library(postpack)
library(coda)
library(FSA)

rm(list=ls())

#----------------Set output file locations ------------------------------
temp_location <- "~/R/working_directory/temp_results/"
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.output/"
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models
results_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims"

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"

source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5.functions.R")

########################## MCMC parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

######################### Simulation settings ######################
PopSim.breeding.schedule <- "biennial.breeding" #Can be annual.breeding or biennial.breeding
model <- "multiennial.model"

jags_params <- c("Nft", "Nfbt", "Nf0", "survival", "lambda", "psi") #List the parameters to be estimated
HS.only <- "yes"
downsample <- "no"
down.perc <- .3
filter.full.sibs <- "no"

jags_file <- paste0(jags.model_location, "MHS.only_narrowLambda_Skip_model_Conn.txt")

#---------------------------Read in and organize data------------------------------------
purpose <- "CKMR_LS_data"
rm.arb <- TRUE #Do we want to remove juvniles with arbitrary birth years (there are 48/2198 in the Bimini dataset)
#estimation.year <- 2000 #Will only be used if using the model with lambda, but the script will throw an error if this object isn't defined (could fix this with a conditional statement, but don't feel like it now :-p).

#Set path to data file
input.file <- "G://My Drive/Personal_Drive/R/CKMR/z_old_files/Objective.5_lemon_shark_data/data/Main_lemon_shark.csv"

#Save von bertalanffy growth function
vonBert <- vbFuns(param = "Typical")
max.age <- 50
repro.age = 12
mating.periodicity <- 2 #Will be translated to a

#Read in lemon shark data and specify column types with short string representation (c is character and d is double; D is date, but is difficult to work with)
lemon_data <- read_csv(input.file, col_types = "ccdccccddccc")
lemon_data <- lemon_data[which(lemon_data$Island=="BIM"),] #Subset for bimini

lemon_data %>% separate(`Tube Label`, sep=";", into = c("Juvenile_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4")) -> lemon_data2 #Separate Tube labels into capture history (there is one tube label per sample instance)

colnames(lemon_data2) <- c("PIT_tag","Capture_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4", "Capture_Year", "Capture_Date", "Island", "Site", "Sex", "PCL_cm", "TL_cm", "DOB", "Father", "Mother") #relabel columns
#head(lemon_data2)

lemon_ref <- lemon_data2[which(lemon_data2$DOB!=""),] #Subset for juveniles with known birth dates

#Double-check that all individuals have a DOB (they do)
lemon_ref %>% group_by(DOB) %>% summarize(n())

#head(lemon_ref)
nrow(lemon_ref)


juv_ref <- lemon_ref %>% dplyr::select(Capture_ID, 
                                       Capture_Year,
                                       Sex,
                                       DOB,
                                       Father,
                                       Mother,
                                       Site)

juv_ref %>% separate(DOB, sep="/", into=c("Yr1", "Yr2")) -> juv_ref #Separate uncertain birth years (e.g. 1993/1994) into separate columns
colnames(juv_ref) <- c("indv.name", "capture.year","sex","Yr1","Yr2","father.x","mother.x", "Site") #Rename columns to match expectations of pairwise comparison function


if(rm.arb == TRUE){
  #Remove samples with uncertain birth date
  juv_ref <- juv_ref %>% dplyr::filter(is.na(Yr2) == TRUE) %>% 
    dplyr::rename(birth.year = Yr1) %>% 
    dplyr::select(-Yr2)
  
} else {
#Randomly assign year if more than one specified
assign_yr <- function(yr1,yr2) {
  if(is.na(yr2)==TRUE) c <- yr1
  else sample(c(yr1, yr2), 1) -> c
  return(c)
  }

DOB = c()
for(i in 1:nrow(juv_ref)) {
DOB[i] <- assign_yr(juv_ref$Yr1[i], juv_ref$Yr2[i])
}

juv_ref$birth.year = DOB #Add randomly assigned birth year to juv_ref
}

#There are three instances where the indv.name (a.k.a. tube ID) is the same for multiple individuals. I manually checked these and they are different individuals in each case, so not recaptures.
juv_ref %>% group_by(indv.name) %>% summarize(num = n()) %>% dplyr::filter(num > 1)

#Which sites are represented in the dataset
juv_ref %>% group_by(Site) %>% summarize(num = n())

#Filter for North Island sites (where sampling was more structured). NOS and SHL are the north island sites
juv_ref.North <- juv_ref %>% dplyr::filter(Site %in% c("NOS", "SHL"))
juv_ref.North %>% group_by(Site) %>% summarize(num = n())

#Make reference table for adults and frequency of offspring by year
father_ref <- juv_ref.North %>% dplyr::select(birth.year, father.x) %>% 
  plyr::count() %>% 
  as_tibble()

mother_ref <- juv_ref.North %>% dplyr::select(birth.year, mother.x) %>% 
  plyr::count() %>% 
  as_tibble

#Look at breeding schedule of known mothers
mother_ref %>% drop_na(mother.x) %>% pivot_wider(names_from = birth.year, values_from = freq)

#Check out what percentage of females bred off-cycle
even.temp <- NULL
rogue.breeders <- NULL

unique.moms <- mother_ref %>% drop_na(mother.x) %>% 
  distinct(mother.x) %>% 
  pull(mother.x)

for(i in 1:length(unique.moms)){
  
  mom <- unique.moms[i]
  
  mother_ref2 <- mother_ref %>% drop_na(mother.x) %>% 
    dplyr::filter(mother.x == mom)
  
  birth.vec <- mother_ref2 %>%
    mutate(birth.year = as.numeric(birth.year)) %>% 
    pull(birth.year)
  
  if(length(birth.vec) > 1){
  
    rogue.temp <- 0
    
  for(j in 2:length(birth.vec)){
    if((birth.vec[j] - birth.vec[j-1]) %% 2 == 0){
      even.temp[j-1] <- 1
    } else if((birth.vec[j] - birth.vec[j-1]) %% 2 != 0){
      rogue.temp <- 1
  }
} 
  if(rogue.temp > 0){
    
  rogue.breeders[i] <- 1 
} else rogue.breeders[i] <- 0
} else next 
}

#I THINK this vector has a 1 for each mom that bred off-cycle at least once, and 0 for each mom that did not. NA if there was just one instance of offspring.
rogue.breeders

length(rogue.breeders[which(rogue.breeders == 1)])/length(rogue.breeders)

#Look at breeding schedule of known (we think) fathers
father_ref %>% drop_na(father.x) %>% pivot_wider(names_from = birth.year, values_from = freq)


#Calculate percent of sampled individuals with unknown mother
num.unknown.mom <- juv_ref.North %>% group_by(indv.name) %>% 
  summarize(unknown.mom = sum(is.na(mother.x) == TRUE)) %>% 
  dplyr::filter(unknown.mom > 0) %>% 
  nrow()

num.known.mom <- juv_ref.North %>% group_by(indv.name) %>% 
  summarize(unknown.mom = sum(is.na(mother.x) == TRUE)) %>% 
  dplyr::filter(unknown.mom == 0) %>% 
  nrow()

#Percent of sampled individuals with unknown mothers
round((num.unknown.mom/num.known.mom)*100,0)


#Replace NAs for parents with a random string of numbers (so assume for now that samples with NA for parents are unrelated)
# for(i in 1:nrow(juv_ref.North)){
#   if(is.na(juv_ref.North$father.x[i]) == TRUE){
#     juv_ref.North$father.x[i] <- paste(sample(letters, size = 20, replace = T), collapse="")
#   }
#   if(is.na(juv_ref.North$mother.x[i]) == TRUE){
#     juv_ref.North$mother.x[i] <- paste(sample(letters, size = 20, replace = T), collapse="")
# }
# }

head(juv_ref.North)

min(juv_ref.North$birth.year, na.rm = TRUE) #Earliest capture
max(juv_ref.North$birth.year, na.rm = TRUE) #Latest capture

#Remove unknown mothers, then for each group of full sibs keep just one
NoFullSibs.df <- juv_ref.North %>% drop_na(mother.x) %>%
  group_by(mother.x, father.x) %>%
  slice_sample(n = 1) %>%
  ungroup()

#Look at number of full sibs
total.sibs <- nrow(juv_ref.North) - nrow(NoFullSibs.df) #Calculate number of full sibs

#Number of individuals contributing to full sibs -- need to work on this
juv_ref.North %>% drop_na(mother.x) %>% 
  group_by(mother.x, father.x) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::summarize(sum(n))

#Different cohort sibs
juv_ref.North %>% group_by(mother.x, father.x) %>%
  filter(n() > 1) %>% 
  distinct(birth.year) %>%
  ungroup() %>% 
  group_by(mother.x, father.x) %>% 
  summarize(distinct.birth.years = n()) %>% 
  filter(distinct.birth.years > 1) %>% 
  nrow()

#To keep the same individual each time
  # distinct(mother.x, father.x, .keep_all = TRUE) %>%
  # as_tibble() #If there is more than one individual with the same mother AND father, then only keep one.

#save the years over which samples were taken
sample.years <- c(min(juv_ref.North$birth.year):max(juv_ref.North$birth.year))
#First year of estimation will be three years after the first year of sampling - that's the second year with solid data
first.est.year <- min(sample.years) + 4
n_yrs <- max(sample.years)

#Save samples
sample.info <- juv_ref.North %>% mutate(birth.year = as.integer(birth.year))
results <- NULL
mom.comps.tibble <- NULL
dad.comps.tibble <- NULL
post.samps_list.block3 <- list()
post.samps_list.block5 <- list()
post.samps_list.all <- list()

if(downsample == "yes"){
  iterations <- 100
} else iterations <- 1

#If need to pick up from a stopping point:
# results <- read_csv(paste0(results_location, "/LS_data_ConnModel_allSibsDown_2023.10.29_temp.csv"))
# 
# mom.comps.tibble <- readRDS(paste0(results_location, "/comps/mom.comps_LS_data_ConnModel_allSibsDown_2023.10.29_temp"))
# 
# dad.comps.tibble <- readRDS(paste0(results_location, "/comps/dad.comps_LS_data_ConnModel_allSibsDown_2023.10.29_temp"))

#------------------------Start loop-----------------------------------------
for(iter in 1:iterations){
  rep <- 1
for(block in first.est.year:n_yrs){
  
  rseed <- sample(c(1:10000000), size = 1)
  
  #Estimate Nt in the present
  estimation.year <- block
  
  #-------------Split samples--------------
  #Blocked  samples - five year
  samples.df.block5 <- sample.info %>% dplyr::filter(birth.year >= block - 4 & birth.year <= block)
  
  #Blocked  samples - three year
  samples.df.block3 <- sample.info %>% dplyr::filter(birth.year >= block - 2 & birth.year <= block)
  
  #all samples
  samples.df.all <- sample.info %>% dplyr::filter(birth.year <= block)
  
  
  #--------------Calculate reference year--------------
  #Five year block
  oldest.block5 <- min(samples.df.block5$birth.year)
  est.year.calibrate.block5 <- oldest.block5 + 1
  
  #Three year block
  oldest.block3 <- min(samples.df.block3$birth.year)
  est.year.calibrate.block3 <- oldest.block3 + 1
  
  #All samples
  oldest.all <- min(samples.df.all$birth.year)
  est.year.calibrate.all <- oldest.all + 1
  
  ####------------------------ Split block and all ----------------
  #----------------Start with blocked data ---------------------#
  #---------------Three year block------------------------
  #-------------Construct pairwise comparison matrix--------------
  cat(paste0("\nBuilding pairwise comparison matrix for three year sample blocks for iteration ", iter, ".\n"))

    filtered.samples.HS.df.block3 <- samples.df.block3 %>% dplyr::arrange(birth.year) %>% 
      ungroup()
  
    if(downsample != "yes"){
      
    if(filter.full.sibs == "yes"){
      filtered.samples.HS.df.block3 <- filtered.samples.HS.df.block3 %>%
        drop_na(mother.x) %>%
        distinct(mother.x, father.x, .keep_all = TRUE)
        # group_by(mother.x, father.x) %>%
        # slice_sample(n = 1) %>%
        # ungroup()
    }

  pairwise.df.HS.block3 <- data.frame(t(combn(filtered.samples.HS.df.block3$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
  colnames(pairwise.df.HS.block3) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(pairwise.df.HS.block3)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.HS.block3 <- filtered.samples.HS.df.block3 %>%
    dplyr::select(older.sib = indv.name, 
                  older.sib.birth = birth.year, 
                  older.sib.mom = mother.x,
                  older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.df_all.info.HS.block3 <- pairwise.df.HS.block3 %>% left_join(OlderSib_birthyears.HS.block3, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(filtered.samples.HS.df.block3, by = "indv.name") %>% 
    dplyr::rename("younger.sib" = indv.name, 
                  "younger.sib.birth" = birth.year, 
                  "younger.sib.mom" = mother.x, 
                  "younger.sib.dad" = father.x) %>%
    dplyr::select(older.sib, 
                  older.sib.birth, 
                  younger.sib, 
                  younger.sib.birth, 
                  older.sib.mom, 
                  younger.sib.mom, 
                  older.sib.dad, 
                  younger.sib.dad)
  
  pairwise.df_HS.filt.block3 <- pairwise.df_all.info.HS.block3 %>% 
    dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
    dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
  
  
  #Extract positive half-sib comparisons
  positives.all.block3 <- pairwise.df_HS.filt.block3 %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
    mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
  #nrow(positives.HS)
  
  ####----------------Construct pairwise comparison matrix for model----------####
  #Sex-specific half-sib
  mom_positives.all.block3 <- positives.all.block3 %>% filter(shared.parent == "mother") %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  dad_positives.all.block3 <- positives.all.block3 %>% filter(shared.parent == "father")  %>%
    dplyr::select(older.sib.birth, younger.sib.birth) %>%
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  
  #Make dataframes for negative comparisons
  hsp.negs.block3 <- pairwise.df_HS.filt.block3 %>%
    dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    as_tibble()
  
  mom_comps.all.block3 <- hsp.negs.block3 %>% 
    dplyr::rename(no = freq) %>% 
    left_join(mom_positives.all.block3, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - est.year.calibrate.block3) %>% 
    mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
  
  dad_comps.all.block3 <- hsp.negs.block3 %>% 
    dplyr::rename(no = freq) %>% 
    left_join(dad_positives.all.block3, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - est.year.calibrate.block3)
  
  num.samps_all.block3 <- nrow(filtered.samples.HS.df.block3) #Save the total number of samples
  mom.pos_all.block3 <- sum(mom_positives.all.block3$yes)
  
  #----------------------------Downsample----------------------------------------------#
  }else if(downsample == "yes"){
    #Downsample for HSPs
    #Calculate proportion of samples born in each birth year
    HS.samps.props.block3 <- filtered.samples.HS.df.block3 %>% group_by(birth.year) %>%
      summarize(num = n()) %>%
      ungroup() %>%
      mutate(num.samps = ceiling(num * down.perc))
    
    #Save vector of proportions corresponding to the correct index in the HS sample dataframe
    # HS.props.vec <- filtered.samples.HS.df %>% left_join(HS.samps.props, by = "birth.year") %>%
    #   pull(num.samps)
    
    HS.samps.df.down.block3 <- NULL #Initialize so I don't get an error
    HS.samps.temp.block3 <- NULL
    
    set.seed(rseed)
    
    for(n in 1:nrow(HS.samps.props.block3)){
      
      down.year.block3 <- HS.samps.props.block3$birth.year[n]
      
      HS.samps.temp.block3 <- filtered.samples.HS.df.block3 %>% 
        dplyr::filter(birth.year == down.year.block3) %>% 
        slice_sample(n = HS.samps.props.block3$num.samps[n]) %>% 
        arrange(birth.year) #arrange by birth year so they're in the correct order for combn
      
      HS.samps.df.down.block3 <- bind_rows(HS.samps.df.down.block3, HS.samps.temp.block3)
      
    }

    #----------------Filter full sibs------------------
    if(filter.full.sibs == "yes"){
      HS.samps.df.down.block3 <- HS.samps.df.down.block3 %>%
        drop_na(mother.x) %>%
        distinct(mother.x, father.x, .keep_all = TRUE)
      # group_by(mother.x, father.x) %>%
      # slice_sample(n = 1) %>%
      # ungroup()
    }
    
    num.samps_down.block3 <- nrow(HS.samps.df.down.block3) 
       
    ####----------------Construct pairwise comparison matrix after downsampling----------####
    
    pairwise.df.HS.block3 <- data.frame(t(combn(HS.samps.df.down.block3$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
    colnames(pairwise.df.HS.block3) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
    
    #Create dataframe of pairwise comparisons with just individual IDs
    head(pairwise.df.HS.block3)
    #head(pairwise.df)
    
    #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
    OlderSib_birthyears.HS.block3 <- filtered.samples.HS.df.block3 %>%
      dplyr::select(older.sib = indv.name, 
                    older.sib.birth = birth.year, 
                    older.sib.mom = mother.x,
                    older.sib.dad = father.x)
    
    #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
    #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
    pairwise.df_all.info.HS.block3 <- pairwise.df.HS.block3 %>% left_join(OlderSib_birthyears.HS.block3, by = "older.sib") %>% 
      as_tibble() %>%  
      left_join(filtered.samples.HS.df.block3, by = "indv.name") %>% 
      dplyr::rename("younger.sib" = indv.name, 
                    "younger.sib.birth" = birth.year, 
                    "younger.sib.mom" = mother.x, 
                    "younger.sib.dad" = father.x) %>%
      dplyr::select(older.sib, 
                    older.sib.birth, 
                    younger.sib, 
                    younger.sib.birth, 
                    older.sib.mom, 
                    younger.sib.mom, 
                    older.sib.dad, 
                    younger.sib.dad)
    
    pairwise.df_HS.filt.block3 <- pairwise.df_all.info.HS.block3 %>% 
      dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
      dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
    
    
    #Extract positive half-sib comparisons
    positives.HS.block3 <- pairwise.df_HS.filt.block3 %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
      mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
    #nrow(positives.HS)
    
    #Sex-specific half-sib
    mom_positives.all.block3 <- positives.HS.block3 %>% filter(shared.parent == "mother") %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    dad_positives.all.block3 <- positives.HS.block3 %>% filter(shared.parent == "father")  %>%
      dplyr::select(older.sib.birth, younger.sib.birth) %>%
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    
    #Make dataframes for negative comparisons
    hsp.negs.block3 <- pairwise.df_HS.filt.block3 %>%
      dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      as_tibble()
    
    mom_comps.all.block3 <- hsp.negs.block3 %>% 
      dplyr::rename(no = freq) %>% 
      left_join(mom_positives.all.block3, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.block3) %>% 
      mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
    
    dad_comps.all.block3 <- hsp.negs.block3 %>% 
      dplyr::rename(no = freq) %>% 
      left_join(dad_positives.all.block3, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.block3) 
    
    mom.pos_down.block3 <- sum(mom_positives.all.block3$yes)
    dad.pos_down.block3 <- sum(dad_positives.all.block3$yes)
  } #End downsample
  
  #Adult.survival <- 0.85 #If wanting to fix survival, or give it a tighter prior
  
  #---------------------------Prep data for JAGS----------------------------
  mom.HSPs.block3 <- mom_comps.all.block3 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  dad.HSPs.block3 <- dad_comps.all.block3 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  
  #yrs <- c(estimation.year:n_yrs)
  ref.year.mom.block3 <- min(mom_comps.all.block3$ref.year)
  ref.year.dad.block3 <- min(dad_comps.all.block3$ref.year)
  
  #Create vectors of data for JAGS
  #Mom
  #HS - even years
  mom_comps.HS_on <- mom_comps.all.block3 %>% dplyr::filter(BI == "even")
  
  mom.mort.yrs_HS.on <- mom_comps.HS_on$year_gap
  mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
  mom.n.comps_HS.on <- mom_comps.HS_on$all
  mom.positives_HS.on <- mom_comps.HS_on$yes
  mom.yrs_HS.on <- nrow(mom_comps.HS_on)
  
  #Mom
  #HS - off years
  mom_comps.HS_off <- mom_comps.all.block3 %>% dplyr::filter(BI == "odd")
  
  mom.mort.yrs_HS.off <- mom_comps.HS_off$year_gap
  mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
  mom.n.comps_HS.off <- mom_comps.HS_off$all
  mom.positives_HS.off <- mom_comps.HS_off$yes
  mom.yrs_HS.off <- nrow(mom_comps.HS_off)
  
  #Dad
  dad.mort.yrs <- dad_comps.all.block3$year_gap
  dad.popGrowth.yrs <- dad_comps.all.block3$pop.growth.yrs
  dad.n.comps <- dad_comps.all.block3$all
  dad.positives <- dad_comps.all.block3$yes
  dad.yrs <- nrow(dad_comps.all.block3)
  
  cat("Running model using three year blocks of samples for iteration ", iter, ".\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", block - 2, " to ", block, "\n"))
  
  est.year.calibrate <- est.year.calibrate.block3
  
  
  #----------------------------Fit JAGS model-------------------------------#
  if(sum(mom_comps.all.block3$yes) > 0){
    
    
    #seed <- sample(c(1:10000000), size = 1)
    #MHSP only model
    set.seed(rseed)

    tryCatch({
        source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
    
    if(downsample == "yes"){
      results.temp.block3 <- model.summary2 %>% mutate(T0 = est.year.calibrate, 
                                                estimation.year = block) %>% 
        mutate(years_sampled = (block - oldest.block3) + 1,
               mom_HSPs = sum(mom_comps.all.block3$yes),
               estimation_year = estimation.year,
               num.samps_down = num.samps_down.block3,
               mom_HSPs.down = mom.pos_down.block3,
               seed = rseed,
               time_window = "three year block") %>% 
        as_tibble()
      
    } else{
    
        results.temp.block3 <- model.summary2 %>% 
        mutate(estimation.year = block,
               T0 = est.year.calibrate.block3,
               num.samps_all = num.samps_all.block3) %>% 
        dplyr::select(parameter,
                      Q50,
                      HPD2.5,
                      HPD97.5,
                      mean,
                      sd,
                      Q2.5,
                      Q97.5,
                      Rhat,
                      estimation.year,
                      T0,
                      neff,
                      seed,
                      num.samps_all) %>% 
        mutate(HSPs_detected = mom.HSPs.block3) %>% 
        mutate(time_window = "three year block")
      
    }
    
    post.samps_list.block3[[rep]] <- post
    names(post.samps_list.block3)[[rep]] <- paste0("yrs_", oldest.block3, "_to_", block)
    
    print(paste0("Finished comparison ", oldest.block3, " to ", block))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  } else {next}

#### End Three year block ####
  

  #---------------Five year block------------------------
  sim <- "five year block"
  
  #-------------Construct pairwise comparison matrix--------------
  cat(paste0("\nBuilding pairwise comparison matrix for five year sample blocks for iteration ", iter, ".\n"))
  
  filtered.samples.HS.df.block5 <- samples.df.block5 %>% dplyr::arrange(birth.year) %>% 
    ungroup()
  
  if(downsample != "yes"){
  
  if(filter.full.sibs == "yes"){
    filtered.samples.HS.df.block5 <- filtered.samples.HS.df.block5 %>%
      drop_na(mother.x) %>%
      distinct(mother.x, father.x, .keep_all = TRUE)
    # group_by(mother.x, father.x) %>%
    # slice_sample(n = 1) %>%
    # ungroup()
  }
  
  pairwise.df.HS.block5 <- data.frame(t(combn(filtered.samples.HS.df.block5$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
  colnames(pairwise.df.HS.block5) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(pairwise.df.HS.block5)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.HS.block5 <- filtered.samples.HS.df.block5 %>%
    dplyr::select(older.sib = indv.name, 
                  older.sib.birth = birth.year, 
                  older.sib.mom = mother.x,
                  older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.df_all.info.HS.block5 <- pairwise.df.HS.block5 %>% left_join(OlderSib_birthyears.HS.block5, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(filtered.samples.HS.df.block5, by = "indv.name") %>% 
    dplyr::rename("younger.sib" = indv.name, 
                  "younger.sib.birth" = birth.year, 
                  "younger.sib.mom" = mother.x, 
                  "younger.sib.dad" = father.x) %>%
    dplyr::select(older.sib, 
                  older.sib.birth, 
                  younger.sib, 
                  younger.sib.birth, 
                  older.sib.mom, 
                  younger.sib.mom, 
                  older.sib.dad, 
                  younger.sib.dad)
  
  pairwise.df_HS.filt.block5 <- pairwise.df_all.info.HS.block5 %>% 
    dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
    dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
  
  
  #Extract positive half-sib comparisons
  positives.all.block5 <- pairwise.df_HS.filt.block5 %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
    mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
  #nrow(positives.HS)
  
  ####----------------Construct pairwise comparison matrix for model----------####
  #Sex-specific half-sib
  mom_positives.all.block5 <- positives.all.block5 %>% filter(shared.parent == "mother") %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  dad_positives.all.block5 <- positives.all.block5 %>% filter(shared.parent == "father")  %>%
    dplyr::select(older.sib.birth, younger.sib.birth) %>%
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  
  #Make dataframes for negative comparisons
  hsp.negs.block5 <- pairwise.df_HS.filt.block5 %>%
    dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    as_tibble()
  
  mom_comps.all.block5 <- hsp.negs.block5 %>% 
    dplyr::rename(no = freq) %>% 
    left_join(mom_positives.all.block5, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - est.year.calibrate.block5) %>% 
    mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
  
  dad_comps.all.block5 <- hsp.negs.block5 %>% 
    dplyr::rename(no = freq) %>% 
    left_join(dad_positives.all.block5, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - est.year.calibrate.block5)
  
  num.samps_all.block5 <- nrow(filtered.samples.HS.df.block5) #Save the total number of samples
  mom.pos_all.block5 <- sum(mom_positives.all.block5$yes)
  
  #----------------------------Downsample----------------------------------------------#
} else if(downsample == "yes"){
    #Downsample for HSPs
    #Calculate proportion of samples born in each birth year
    HS.samps.props.block5 <- filtered.samples.HS.df.block5 %>% group_by(birth.year) %>%
      summarize(num = n()) %>%
      ungroup() %>%
      mutate(num.samps = ceiling(num * down.perc))
    
    #Save vector of proportions corresponding to the correct index in the HS sample dataframe
    # HS.props.vec <- filtered.samples.HS.df %>% left_join(HS.samps.props, by = "birth.year") %>%
    #   pull(num.samps)
    
    HS.samps.df.down.block5 <- NULL #Initialize so I don't get an error
    HS.samps.temp.block5 <- NULL
    
    set.seed(rseed)
    
    for(n in 1:nrow(HS.samps.props.block5)){
      
      down.year.block5 <- HS.samps.props.block5$birth.year[n]
      
      HS.samps.temp.block5 <- filtered.samples.HS.df.block5 %>% 
        dplyr::filter(birth.year == down.year.block5) %>% 
        slice_sample(n = HS.samps.props.block5$num.samps[n]) %>% 
        arrange(birth.year) #arrange by birth year so they're in the correct order for combn
      
      HS.samps.df.down.block5 <- bind_rows(HS.samps.df.down.block5, HS.samps.temp.block5)
      
    }
  
    #----------------Filter full sibs------------------
    if(filter.full.sibs == "yes"){
      HS.samps.df.down.block5 <- HS.samps.df.down.block5 %>%
        drop_na(mother.x) %>%
        distinct(mother.x, father.x, .keep_all = TRUE)
      # group_by(mother.x, father.x) %>%
      # slice_sample(n = 1) %>%
      # ungroup()
    }
  
    num.samps_down.block5 <- nrow(HS.samps.df.down.block5)
    
    
    ####----------------Construct pairwise comparison matrix after downsampling----------####
    pairwise.df.HS.block5 <- data.frame(t(combn(HS.samps.df.down.block5$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
    colnames(pairwise.df.HS.block5) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
    
    #Create dataframe of pairwise comparisons with just individual IDs
    head(pairwise.df.HS.block5)
    #head(pairwise.df)
    
    #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
    OlderSib_birthyears.HS.block5 <- filtered.samples.HS.df.block5 %>%
      dplyr::select(older.sib = indv.name, 
                    older.sib.birth = birth.year, 
                    older.sib.mom = mother.x,
                    older.sib.dad = father.x)
    
    #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
    #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
    pairwise.df_all.info.HS.block5 <- pairwise.df.HS.block5 %>% left_join(OlderSib_birthyears.HS.block5, by = "older.sib") %>% 
      as_tibble() %>%  
      left_join(filtered.samples.HS.df.block5, by = "indv.name") %>% 
      dplyr::rename("younger.sib" = indv.name, 
                    "younger.sib.birth" = birth.year, 
                    "younger.sib.mom" = mother.x, 
                    "younger.sib.dad" = father.x) %>%
      dplyr::select(older.sib, 
                    older.sib.birth, 
                    younger.sib, 
                    younger.sib.birth, 
                    older.sib.mom, 
                    younger.sib.mom, 
                    older.sib.dad, 
                    younger.sib.dad)
    
    pairwise.df_HS.filt.block5 <- pairwise.df_all.info.HS.block5 %>% 
      dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
      dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
    
    
    #Extract positive half-sib comparisons
    positives.HS.block5 <- pairwise.df_HS.filt.block5 %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
      mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
    #nrow(positives.HS)
    
    #Sex-specific half-sib
    mom_positives.all.block5 <- positives.HS.block5 %>% filter(shared.parent == "mother") %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    dad_positives.all.block5 <- positives.HS.block5 %>% filter(shared.parent == "father")  %>%
      dplyr::select(older.sib.birth, younger.sib.birth) %>%
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    
    #Make dataframes for negative comparisons
    hsp.negs.block5 <- pairwise.df_HS.filt.block5 %>%
      dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      as_tibble()
    
    mom_comps.all.block5 <- hsp.negs.block5 %>% 
      dplyr::rename(no = freq) %>% 
      left_join(mom_positives.all.block5, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.block5) %>% 
      mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
    
    dad_comps.all.block5 <- hsp.negs.block5 %>% 
      dplyr::rename(no = freq) %>% 
      left_join(dad_positives.all.block5, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.block5) 
    
    mom.pos_down.block5 <- sum(mom_positives.all.block5$yes)
    dad.pos_down.block5 <- sum(dad_positives.all.block5$yes)
  } #End downsample
  
  #Adult.survival <- 0.85 #If wanting to fix survival, or give it a tighter prior
  
  #---------------------------Prep data for JAGS----------------------------
  mom.HSPs.block5 <- mom_comps.all.block5 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  dad.HSPs.block5 <- dad_comps.all.block5 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  
  #yrs <- c(estimation.year:n_yrs)
  ref.year.mom.block5 <- min(mom_comps.all.block5$ref.year)
  ref.year.dad.block5 <- min(dad_comps.all.block5$ref.year)
  
  #Create vectors of data for JAGS
  #Mom
  #HS - even years
  mom_comps.HS_on <- mom_comps.all.block5 %>% dplyr::filter(BI == "even")
  
  mom.mort.yrs_HS.on <- mom_comps.HS_on$year_gap
  mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
  mom.n.comps_HS.on <- mom_comps.HS_on$all
  mom.positives_HS.on <- mom_comps.HS_on$yes
  mom.yrs_HS.on <- nrow(mom_comps.HS_on)
  
  #Mom
  #HS - off years
  mom_comps.HS_off <- mom_comps.all.block5 %>% dplyr::filter(BI == "odd")
  
  mom.mort.yrs_HS.off <- mom_comps.HS_off$year_gap
  mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
  mom.n.comps_HS.off <- mom_comps.HS_off$all
  mom.positives_HS.off <- mom_comps.HS_off$yes
  mom.yrs_HS.off <- nrow(mom_comps.HS_off)
  
  #Dad
  dad.mort.yrs <- dad_comps.all.block5$year_gap
  dad.popGrowth.yrs <- dad_comps.all.block5$pop.growth.yrs
  dad.n.comps <- dad_comps.all.block5$all
  dad.positives <- dad_comps.all.block5$yes
  dad.yrs <- nrow(dad_comps.all.block5)
  
  cat("Running model using five year blocks of samples for iteration ", iter, ".\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", block - 4, " to ", block, "\n"))
  
  est.year.calibrate <- est.year.calibrate.block5
  
  
  #----------------------------Fit JAGS model-------------------------------#
  if(sum(mom_comps.all.block5$yes) > 0){
    
    
    #seed <- sample(c(1:10000000), size = 1)
    #MHSP only model
    set.seed(rseed)
    
    tryCatch({
    source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
    
    if(downsample == "yes"){
      results.temp.block5 <- model.summary2 %>% mutate(T0 = est.year.calibrate, 
                                                estimation.year = block) %>% 
        mutate(years_sampled = (block - oldest.block5) + 1,
               mom_HSPs = sum(mom_comps.all.block5$yes),
               estimation_year = estimation.year,
               num.samps_down = num.samps_down.block5,
               mom_HSPs.down = mom.pos_down.block5,
               seed = rseed,
               time_window = "five year block") %>% 
        as_tibble()
      
    } else{
      
      results.temp.block5 <- model.summary2 %>% 
        mutate(estimation.year = block,
               T0 = est.year.calibrate.block5,
               num.samps_all = num.samps_all.block5) %>% 
        dplyr::select(parameter,
                      Q50,
                      HPD2.5,
                      HPD97.5,
                      mean,
                      sd,
                      Q2.5,
                      Q97.5,
                      Rhat,
                      estimation.year,
                      T0,
                      neff,
                      seed,
                      num.samps_all) %>% 
        mutate(HSPs_detected = mom.HSPs.block5) %>% 
        mutate(time_window = "five year block")
      
    }
    
    post.samps_list.block5[[rep]] <- post
    names(post.samps_list.block5)[[rep]] <- paste0("yrs_", oldest.block5, "_to_", block)
    
    print(paste0("Finished comparison ", oldest.block5, " to ", block))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  } else {next}
  

  gc()
  
  #----------------Use all data ---------------------
  #-------------Construct pairwise comparison matrix--------------
  sim <- "all"
  #All samples - set up so we don't have to build the pairwise comparison matrix with every shift in estimation year
  cat(paste0("\nBuilding pairwise comparison matrix for all samples for iteration ", iter, ".\n"))

  filtered.samples.HS.df.all <- samples.df.all %>% dplyr::arrange(birth.year) %>% 
    ungroup()
  
  if(downsample != "yes"){
    
  if(filter.full.sibs == "yes"){
    filtered.samples.HS.df.all <- filtered.samples.HS.df.all %>%
      drop_na(mother.x) %>%
      distinct(mother.x, father.x, .keep_all = TRUE)
    # group_by(mother.x, father.x) %>%
    # slice_sample(n = 1) %>%
    # ungroup()
  }
  
  pairwise.df.HS.all <- data.frame(t(combn(filtered.samples.HS.df.all$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
  colnames(pairwise.df.HS.all) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(pairwise.df.HS.all)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.HS.all <- filtered.samples.HS.df.all %>%
    dplyr::select(older.sib = indv.name, 
                  older.sib.birth = birth.year, 
                  older.sib.mom = mother.x,
                  older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.df_all.info.HS.all <- pairwise.df.HS.all %>% left_join(OlderSib_birthyears.HS.all, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(filtered.samples.HS.df.all, by = "indv.name") %>% 
    dplyr::rename("younger.sib" = indv.name, 
                  "younger.sib.birth" = birth.year, 
                  "younger.sib.mom" = mother.x, 
                  "younger.sib.dad" = father.x) %>%
    dplyr::select(older.sib, 
                  older.sib.birth, 
                  younger.sib, 
                  younger.sib.birth, 
                  older.sib.mom, 
                  younger.sib.mom, 
                  older.sib.dad, 
                  younger.sib.dad)
  
  pairwise.df_HS.filt.all <- pairwise.df_all.info.HS.all %>% 
    dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
    dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
  
  
  #Extract positive half-sib comparisons
  positives.all.all <- pairwise.df_HS.filt.all %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
    mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
  #nrow(positives.HS)
  
  ####----------------Construct pairwise comparison matrix for model----------####
  #Sex-specific half-sib
  mom_positives.all.all <- positives.all.all %>% filter(shared.parent == "mother") %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  dad_positives.all.all <- positives.all.all %>% filter(shared.parent == "father")  %>%
    dplyr::select(older.sib.birth, younger.sib.birth) %>%
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  
  #Make dataframes for negative comparisons
  hsp.negs.all <- pairwise.df_HS.filt.all %>%
    dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    as_tibble()
  
  mom_comps.all.all <- hsp.negs.all %>% 
    dplyr::rename(no = freq) %>% 
    left_join(mom_positives.all.all, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - est.year.calibrate.all) %>% 
    mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
  
  dad_comps.all.all <- hsp.negs.all %>% 
    dplyr::rename(no = freq) %>% 
    left_join(dad_positives.all.all, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - est.year.calibrate.all)
  
  num.samps_all.all <- nrow(filtered.samples.HS.df.all) #Save the total number of samples
  mom.pos_all.all <- sum(mom_positives.all.all$yes)
  
  #----------------------------Downsample----------------------------------------------#
} else if(downsample == "yes"){
    #Downsample for HSPs
    #Calculate proportion of samples born in each birth year
    HS.samps.props.all <- filtered.samples.HS.df.all %>% group_by(birth.year) %>%
      summarize(num = n()) %>%
      ungroup() %>%
      mutate(num.samps = ceiling(num * down.perc))
    
    #Save vector of proportions corresponding to the correct index in the HS sample dataframe
    # HS.props.vec <- filtered.samples.HS.df %>% left_join(HS.samps.props, by = "birth.year") %>%
    #   pull(num.samps)
    
    HS.samps.df.down.all <- NULL #Initialize so I don't get an error
    HS.samps.temp.all <- NULL
    
    set.seed(rseed)
    
    for(n in 1:nrow(HS.samps.props.all)){
      
      down.year.all <- HS.samps.props.all$birth.year[n]
      
      HS.samps.temp.all <- filtered.samples.HS.df.all %>% 
        dplyr::filter(birth.year == down.year.all) %>% 
        slice_sample(n = HS.samps.props.all$num.samps[n]) %>% 
        arrange(birth.year) #arrange by birth year so they're in the correct order for combn
      
      HS.samps.df.down.all <- bind_rows(HS.samps.df.down.all, HS.samps.temp.all)
      
    }
    
    #----------------Filter full sibs------------------
    if(filter.full.sibs == "yes"){
      HS.samps.df.down.all <- HS.samps.df.down.all %>%
        drop_na(mother.x) %>%
        distinct(mother.x, father.x, .keep_all = TRUE)
#        group_by(mother.x, father.x) %>%
 #       slice_sample(n = 1) %>%
 #       ungroup()
    }
    
    num.samps_down.all <- nrow(HS.samps.df.down.all)
    
    
    ####----------------Construct pairwise comparison matrix after downsampling----------####
    pairwise.df.HS.all <- data.frame(t(combn(HS.samps.df.down.all$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
    colnames(pairwise.df.HS.all) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
    
    #Create dataframe of pairwise comparisons with just individual IDs
    head(pairwise.df.HS.all)
    #head(pairwise.df)
    
    #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
    OlderSib_birthyears.HS.all <- filtered.samples.HS.df.all %>%
      dplyr::select(older.sib = indv.name, 
                    older.sib.birth = birth.year, 
                    older.sib.mom = mother.x,
                    older.sib.dad = father.x)
    
    #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
    #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
    pairwise.df_all.info.HS.all <- pairwise.df.HS.all %>% left_join(OlderSib_birthyears.HS.all, by = "older.sib") %>% 
      as_tibble() %>%  
      left_join(filtered.samples.HS.df.all, by = "indv.name") %>% 
      dplyr::rename("younger.sib" = indv.name, 
                    "younger.sib.birth" = birth.year, 
                    "younger.sib.mom" = mother.x, 
                    "younger.sib.dad" = father.x) %>%
      dplyr::select(older.sib, 
                    older.sib.birth, 
                    younger.sib, 
                    younger.sib.birth, 
                    older.sib.mom, 
                    younger.sib.mom, 
                    older.sib.dad, 
                    younger.sib.dad)
    
    pairwise.df_HS.filt.all <- pairwise.df_all.info.HS.all %>% 
      dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
      dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
    
    
    #Extract positive half-sib comparisons
    positives.HS.all <- pairwise.df_HS.filt.all %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
      mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
    #nrow(positives.HS)
    
    #Sex-specific half-sib
    mom_positives.all.all <- positives.HS.all %>% filter(shared.parent == "mother") %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    dad_positives.all.all <- positives.HS.all %>% filter(shared.parent == "father")  %>%
      dplyr::select(older.sib.birth, younger.sib.birth) %>%
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    
    #Make dataframes for negative comparisons
    hsp.negs.all <- pairwise.df_HS.filt.all %>%
      dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      as_tibble()
    
    mom_comps.all.all <- hsp.negs.all %>% 
      dplyr::rename(no = freq) %>% 
      left_join(mom_positives.all.all, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.all) %>% 
      mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
    
    dad_comps.all.all <- hsp.negs.all %>% 
      dplyr::rename(no = freq) %>% 
      left_join(dad_positives.all.all, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.all) 
    
    mom.pos_down.all <- sum(mom_positives.all.all$yes)
    dad.pos_down.all <- sum(dad_positives.all.all$yes)
  } #End downsample
  
  #Adult.survival <- 0.85 #If wanting to fix survival, or give it a tighter prior
  
  #---------------------------Prep data for JAGS----------------------------
  mom.HSPs.all <- mom_comps.all.all %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  dad.HSPs.all <- dad_comps.all.all %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  
  #yrs <- c(estimation.year:n_yrs)
  ref.year.mom.all <- min(mom_comps.all.all$ref.year)
  ref.year.dad.all <- min(dad_comps.all.all$ref.year)
  
  #Create vectors of data for JAGS
  #Mom
  #HS - even years
  mom_comps.HS_on <- mom_comps.all.all %>% dplyr::filter(BI == "even")
  
  mom.mort.yrs_HS.on <- mom_comps.HS_on$year_gap
  mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
  mom.n.comps_HS.on <- mom_comps.HS_on$all
  mom.positives_HS.on <- mom_comps.HS_on$yes
  mom.yrs_HS.on <- nrow(mom_comps.HS_on)
  
  #Mom
  #HS - off years
  mom_comps.HS_off <- mom_comps.all.all %>% dplyr::filter(BI == "odd")
  
  mom.mort.yrs_HS.off <- mom_comps.HS_off$year_gap
  mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
  mom.n.comps_HS.off <- mom_comps.HS_off$all
  mom.positives_HS.off <- mom_comps.HS_off$yes
  mom.yrs_HS.off <- nrow(mom_comps.HS_off)
  
  #Dad
  dad.mort.yrs <- dad_comps.all.all$year_gap
  dad.popGrowth.yrs <- dad_comps.all.all$pop.growth.yrs
  dad.n.comps <- dad_comps.all.all$all
  dad.positives <- dad_comps.all.all$yes
  dad.yrs <- nrow(dad_comps.all.all)
  
  cat("Running model using all available samples for iteration ", iter, ".\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", oldest.all, " to ", block, "\n"))
  
  est.year.calibrate <- est.year.calibrate.all
  
  
  #----------------------------Fit JAGS model-------------------------------#
  if(sum(mom_comps.all.all$yes) > 0){
    
    
    #seed <- sample(c(1:10000000), size = 1)
    #MHSP only model
    set.seed(rseed)
    
    tryCatch({
    source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
    
    if(downsample == "yes"){
      results.temp.all <- model.summary2 %>% mutate(T0 = est.year.calibrate, 
                                                estimation.year = block) %>% 
        mutate(years_sampled = (block - oldest.all) + 1,
               mom_HSPs = sum(mom_comps.all.all$yes),
               estimation_year = estimation.year,
               num.samps_down = num.samps_down.all,
               mom_HSPs.down = mom.pos_down.all,
               seed = rseed,
               time_window = "all available samples") %>% 
        as_tibble()
      
    } else{
      
      results.temp.all <- model.summary2 %>% 
        mutate(estimation.year = block,
               T0 = est.year.calibrate.all,
               num.samps_all = num.samps_all.all) %>% 
        dplyr::select(parameter,
                      Q50,
                      HPD2.5,
                      HPD97.5,
                      mean,
                      sd,
                      Q2.5,
                      Q97.5,
                      Rhat,
                      estimation.year,
                      T0,
                      neff,
                      seed,
                      num.samps_all) %>% 
        mutate(HSPs_detected = mom.HSPs.all) %>% 
        mutate(time_window = "all available samples")
      
    }
    
    post.samps_list.all[[rep]] <- post
    names(post.samps_list.all)[[rep]] <- paste0("yrs_", oldest.all, "_to_", block)
    
    print(paste0("Finished comparison ", oldest.all, " to ", block))
    
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  } else {next}
  
  
  gc()
  
  
  
    ####------------------------ Combine block and all ----------------####
  #Save mom and dad pairwise comparison dataframes
  #Five year block
  mom_comps.all.block5 <- mom_comps.all.block5 %>% mutate(seed = rseed,
                                                          time_window = "five year block")
  
  dad_comps.all.block5 <- dad_comps.all.block5 %>% mutate(seed = rseed,
                                                          time_window = "five year block")
  
  
  #Three year block  
  mom_comps.all.block3 <- mom_comps.all.block3 %>% mutate(seed = rseed,
                                                          time_window = "three year block")
  
  dad_comps.all.block3 <- dad_comps.all.block3 %>% mutate(seed = rseed,
                                                          time_window = "three year block")
  
  
  #All samples
  mom_comps.all.all <- mom_comps.all.all %>% mutate(seed = rseed,
                                                    time_window = "all available samples")
  
  
  dad_comps.all.all <- dad_comps.all.all %>% mutate(seed = rseed,
                                                    time_window = "all available samples")
  
  mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all.block3, mom_comps.all.block5, mom_comps.all.all)
  
  dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all.block3, dad_comps.all.block5, dad_comps.all.all)
  
  
  #-----------------Store output files iteratively--------------------
  #Results
  if(downsample == "yes"){
    results.temp.block3 <- results.temp.block3 %>% mutate(iteration = iter)
    results.temp.block5 <- results.temp.block5 %>% mutate(iteration = iter)
    results.temp.all <- results.temp.all %>% mutate(iteration = iter)
  }
  results <- rbind(results, results.temp.block3, results.temp.block5, results.temp.all)
  
rep <- rep+1

}
  
  #Save as temporary files so they save and overwrite every iteration. Also will ask Shelley to run the code at the end.
  # write_csv(results, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/LS_data_ConnModel_allSibs_2023.12.08_temp.csv")
  # 
  # #Save final pairwise comparison matrices
  # saveRDS(mom.comps.tibble, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/mom.comps_LS_data_ConnModel_allSibs_2023.12.08_temp")
  # 
  # saveRDS(dad.comps.tibble, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/dad.comps_LS_data_ConnModel_allSibs_2023.12.08_temp")
  
}



#### End loop over blocks ####

results %>% dplyr::filter(estimation.year >= 1997) %>% 
  as_tibble() %>% 
  arrange(time_window, estimation.year)

results %>% dplyr::filter(parameter == "Nfbt", time_window == "all available samples", estimation.year == 2002) %>% distinct(.keep_all = TRUE)

#JAGS seemed to struggle, and kept throwing errors. Not sure why, since the model has been working fine: when filtering for one sibling it worked fine, and when downsampling it was fine, but not when downsampling AND filtering for one sibling. For some reason, some model runs got duplicated. But we can remove them.
results %>% distinct(.keep_all = TRUE) %>% dplyr::filter(parameter == "Nfbt", time_window == "all available samples") %>% dplyr::count(estimation.year)

results %>% distinct(.keep_all = TRUE) %>% 
  write_csv(file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/peer_review_rd2/LS_data/LS_data_ConnModel_filterSibs_2023.12.08.csv")

#Save final pairwise comparison matrices
saveRDS(mom.comps.tibble, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/peer_review_rd2/LS_data/mom.comps_LS_data_ConnModel_filterSibs_2023.12.08")

saveRDS(dad.comps.tibble, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/peer_review_rd2/LS_data/dad.comps_LS_data_ConnModel_filterSibs_2023.12.08")


#Save model output
saveRDS(post.samps_list.block3, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/peer_review_rd2/LS_data/Model.output/block3_model_results_filterSibs_2023.12.08")

saveRDS(post.samps_list.block5, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/peer_review_rd2/LS_data/Model.output/block5_model_results_filterSibs_2023.12.08")

saveRDS(post.samps_list.all, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lemon_shark_sims/peer_review_rd2/LS_data/Model.output/all_model_results_filterSibs_2023.12.08")
