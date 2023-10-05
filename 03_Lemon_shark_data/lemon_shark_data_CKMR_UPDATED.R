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

jags_params <- c("Nft", "Nfbt", "Nf0", "Nfb0", "survival", "lambda", "psi") #List the parameters to be estimated
HS.only <- "yes"
downsample <- "yes"
down.perc <- .3

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

#To keep the same individual each time
  # distinct(mother.x, father.x, .keep_all = TRUE) %>%
  # as_tibble() #If there is more than one individual with the same mother AND father, then only keep one.

#save the years over which samples were taken
sample.years <- c(min(NoFullSibs.df$birth.year):max(NoFullSibs.df$birth.year))
#First year of estimation will be three years after the first year of sampling - that's the second year with solid data
first.est.year <- sample.years + 3
n_yrs <- max(sample.years)

#Save samples
sample.info <- NoFullSibs.df %>% mutate(birth.year = as.integer(birth.year))
results <- NULL
mom.comps.tibble <- NULL
dad.comps.tibble <- NULL
post.samps_list.block3 <- list()
post.samps_list.block5 <- list()
post.samps_list.all <- list()
rseed <- sample(c(1:10000000), size = 1)
iterations <- 50


#------------------------Start loop-----------------------------------------
for(iter in 1:iterations){
for(block in first.est.year:n_yrs){
  
  rep <- 1
  
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
  cat(paste0("\nBuilding pairwise comparison matrix for three year sample blocks\n"))

    filtered.samples.HS.df.block3 <- samples.df.block3 %>% dplyr::arrange(birth.year) %>% 
      ungroup()

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
  if(downsample == "yes"){
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
    mom_positives.HS.block3 <- positives.HS.block3 %>% filter(shared.parent == "mother") %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    dad_positives.HS.block3 <- positives.HS.block3 %>% filter(shared.parent == "father")  %>%
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
      left_join(mom_positives.HS.block3, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.block3) %>% 
      mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
    
    dad_comps.all.block3 <- hsp.negs.block3 %>% 
      dplyr::rename(no = freq) %>% 
      left_join(dad_positives.HS.block3, by = c("older.sib.birth", "younger.sib.birth")) %>% 
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
  
  cat("Running model using three year blocks of samples\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", block - 2, " to ", block, "\n"))
  
  est.year.calibrate <- est.year.calibrate.block3
  
  
  #----------------------------Fit JAGS model-------------------------------#
  if(sum(mom_comps.all.block3$yes) > 0){
    
    
    #seed <- sample(c(1:10000000), size = 1)
    #MHSP only model
    set.seed(rseed)

        source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
    
    if(downsample == "yes"){
      results.temp.block3 <- model.summary2 %>% mutate(T0 = est.year.calibrate, 
                                                estimation.year = block) %>% 
        mutate(years_sampled = (block - oldest.block3) + 1,
               mom_HSPs = sum(mom_comps.all.block3$yes),
               estimation_year = estimation.year,
               num.samps_all = num.samps_all.block3,
               num.samps_down = num.samps_down.block3,
               mom_HSPS.all = mom.pos_all.block3,
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
                      seed) %>% 
        mutate(HSPs_detected = mom.HSPs.block3) %>% 
        mutate(time_window = "three year block")
      
    }
    
    post.samps_list.block3[[rep]] <- post
    names(post.samps_list.block3)[[rep]] <- paste0("yrs_", oldest.block3, "_to_", block)
    
    print(paste0("Finished comparison ", oldest.block3, " to ", block))
  } else {next}

#### End Three year block ####
  
####PICK UP HERE ON 10/5/2023: 
#Copy and paste above to five year window and all samples and edit ".block3" suffix
#Start testing
#Run on all samples (for main text)
#Run on downsampled samples (for supplementary)
#


  
  
    
  #---------------Five year block------------------------
  sim <- "five year block"
  
  #-------------Construct pairwise comparison matrix--------------
  cat(paste0("\nBuilding pairwise comparison matrix for five year sample blocks\n"))
  
  filtered.samples.HS.df.block5 <- samples.df.block5 %>% dplyr::arrange(birth.year) %>% 
    ungroup()
  
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
  if(downsample == "yes"){
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
    mom_positives.HS.block5 <- positives.HS.block5 %>% filter(shared.parent == "mother") %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    dad_positives.HS.block5 <- positives.HS.block5 %>% filter(shared.parent == "father")  %>%
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
      left_join(mom_positives.HS.block5, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.block5) %>% 
      mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
    
    dad_comps.all.block5 <- hsp.negs.block5 %>% 
      dplyr::rename(no = freq) %>% 
      left_join(dad_positives.HS.block5, by = c("older.sib.birth", "younger.sib.birth")) %>% 
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
  
  cat("Running model using five year blocks of samples\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", block - 4, " to ", block, "\n"))
  
  est.year.calibrate <- est.year.calibrate.block5
  
  
  #----------------------------Fit JAGS model-------------------------------#
  if(sum(mom_comps.all.block5$yes) > 0){
    
    
    #seed <- sample(c(1:10000000), size = 1)
    #MHSP only model
    set.seed(rseed)
    
    source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
    
    if(downsample == "yes"){
      results.temp.block5 <- model.summary2 %>% mutate(T0 = est.year.calibrate, 
                                                estimation.year = block) %>% 
        mutate(years_sampled = (block - oldest.block5) + 1,
               mom_HSPs = sum(mom_comps.all.block5$yes),
               estimation_year = estimation.year,
               num.samps_all = num.samps_all.block5,
               num.samps_down = num.samps_down.block5,
               mom_HSPS.all = mom.pos_all.block5,
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
                      seed) %>% 
        mutate(HSPs_detected = mom.HSPs.block5) %>% 
        mutate(time_window = "five year block")
      
    }
    
    post.samps_list.block5[[rep]] <- post
    names(post.samps_list.block5)[[rep]] <- paste0("yrs_", oldest.block5, "_to_", block)
    
    print(paste0("Finished comparison ", oldest.block5, " to ", block))
  } else {next}
  

  gc()
  
  #----------------Use all data ---------------------
  #-------------Construct pairwise comparison matrix--------------
  sim <- "all"
  #All samples - set up so we don't have to build the pairwise comparison matrix with every shift in estimation year
  cat(paste0("\nBuilding pairwise comparison matrix for all samples\n"))


  filtered.samples.HS.df.all <- samples.df.all %>% dplyr::arrange(birth.year) %>% 
    ungroup()
  
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
  if(downsample == "yes"){
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
    mom_positives.HS.all <- positives.HS.all %>% filter(shared.parent == "mother") %>% 
      dplyr::select(older.sib.birth, younger.sib.birth) %>% 
      plyr::count() %>% 
      dplyr::rename(yes = freq)
    
    dad_positives.HS.all <- positives.HS.all %>% filter(shared.parent == "father")  %>%
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
      left_join(mom_positives.HS.all, by = c("older.sib.birth", "younger.sib.birth")) %>% 
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
      mutate(all = yes + no) %>% 
      mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
      dplyr::rename(ref.year = younger.sib.birth) %>% 
      mutate(pop.growth.yrs = ref.year - est.year.calibrate.all) %>% 
      mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
    
    dad_comps.all.all <- hsp.negs.all %>% 
      dplyr::rename(no = freq) %>% 
      left_join(dad_positives.HS.all, by = c("older.sib.birth", "younger.sib.birth")) %>% 
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
  
  cat("Running model using all available samples\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", oldest.all, " to ", block, "\n"))
  
  est.year.calibrate <- est.year.calibrate.all
  
  
  #----------------------------Fit JAGS model-------------------------------#
  if(sum(mom_comps.all.all$yes) > 0){
    
    
    #seed <- sample(c(1:10000000), size = 1)
    #MHSP only model
    set.seed(rseed)
    
    source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
    
    if(downsample == "yes"){
      results.temp.all <- model.summary2 %>% mutate(T0 = est.year.calibrate, 
                                                estimation.year = block) %>% 
        mutate(years_sampled = (block - oldest.all) + 1,
               mom_HSPs = sum(mom_comps.all.all$yes),
               estimation_year = estimation.year,
               num.samps_all = num.samps_all.all,
               num.samps_down = num.samps_down.all,
               mom_HSPS.all = mom.pos_all.all,
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
                      seed) %>% 
        mutate(HSPs_detected = mom.HSPs.all) %>% 
        mutate(time_window = "all available samples")
      
    }
    
    post.samps_list.all[[rep]] <- post
    names(post.samps_list.all)[[rep]] <- paste0("yrs_", oldest.all, "_to_", block)
    
    print(paste0("Finished comparison ", oldest.all, " to ", block))
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
}

#### End loop over blocks ####
results %>% dplyr::filter(estimation.year >= 1997) %>% 
  as_tibble() %>% 
  arrange(time_window, estimation.year)

write_csv(results, file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/LS_data_ConnModel_2023.10.05.csv")

#Save final pairwise comparison matrices
saveRDS(mom.comps.tibble, file = paste0(results_location, "mom.comps_LS_data_ConnModel_2023.10.05"))

saveRDS(dad.comps.tibble, file = paste0(results_location, "dad.comps_LS_data_ConnModel_2023.10.05"))





#######################################################################################-
#--------------Run simulations over different subsets of data--------------------------
#######################################################################################-

estimation.year <- 2000

years <- c(min(NoFullSibs.df$birth.year):max(NoFullSibs.df$birth.year))

results <- NULL
post.samps_list <- list()
rep = 0
downsample <- "no"
down.perc <- 0.3 #Percentage of samples to retain from each year
#max.HSPs <- 150 #Not using anymore

rseeds <- sample(c(1:10000000), size = 200)

#Regular loop over different periods of estimation
#for(i in min(years)){
#  for(j in (i + 2): max(years)){ #To iteratively add to all years
for(i in min(years):max(years-3)){ #For four year intervals
  for(j in (i + 3)){ #To run in four year intervals

    estimation.year = j
    
    #When an error and need to start at a certain number
# for(i in 2012:max(years-1)){
#   for(j in (i + 2): max(years)){

        rep = rep+1
        set.seed(seed)
        
        filtered.samples.HS.df <- NoFullSibs.df %>% 
        #  group_by(mother.x, birth.year) %>% #Can comment out this line and the next to include all HSPs; uncomment to just use one indv per litter.
        #  slice_sample(n=1) %>% 
          dplyr::arrange(birth.year) %>% 
          mutate(birth.year = as.numeric(birth.year)) %>%  #Arrange so older sib always comes first in pairwise comparison matrix
          ungroup() %>% #Uncomment the pipe to include filtering
          dplyr::filter(birth.year >= i & birth.year <= j) #Comment out this line and the pipe in the previous line to use the full dataset to estimate abundance for each year

pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
  dplyr::select(older.sib = indv.name, 
         older.sib.birth = birth.year, 
         older.sib.mom = mother.x,
         older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(filtered.samples.HS.df, by = "indv.name") %>% 
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

pairwise.df_HS.filt <- pairwise.df_all.info.HS %>% 
  dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
  dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))


#Extract positive half-sib comparisons
positives.HS <- pairwise.df_HS.filt %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
  mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
#nrow(positives.HS)

####----------------Construct pairwise comparison matrix for model----------####
#Sex-specific half-sib
mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  dplyr::rename(yes = freq)

dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
  dplyr::select(older.sib.birth, younger.sib.birth) %>%
  plyr::count() %>% 
  dplyr::rename(yes = freq)


#Make dataframes for negative comparisons
hsp.negs <- pairwise.df_HS.filt %>%
  dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  as_tibble()

mom_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  mutate(mom.oncycle = ifelse(BI == "even", 1, 0))

dad_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) 

num.samps_all <- nrow(filtered.samples.HS.df) #Save the total number of samples
mom.pos_all <- sum(mom_positives.HS$yes)

#----------------------------Downsample----------------------------------------------#
if(downsample == "yes"){
#Downsample for HSPs
#Calculate proportion of samples born in each birth year
HS.samps.props <- filtered.samples.HS.df %>% group_by(birth.year) %>%
  summarize(num = n()) %>%
  ungroup() %>%
  mutate(num.samps = ceiling(num * down.perc))

#Save vector of proportions corresponding to the correct index in the HS sample dataframe
# HS.props.vec <- filtered.samples.HS.df %>% left_join(HS.samps.props, by = "birth.year") %>%
#   pull(num.samps)

HS.samps.df.down <- NULL #Initialize so I don't get an error
HS.samps.temp <- NULL

set.seed(seed)
for(n in 1:nrow(HS.samps.props)){

    down.year <- HS.samps.props$birth.year[n]
  
  HS.samps.temp <- filtered.samples.HS.df %>% dplyr::filter(birth.year == down.year) %>% 
    slice_sample(n = HS.samps.props$num.samps[n]) %>% 
    arrange(birth.year) #arrange by birth year so they're in the correct order for combn

HS.samps.df.down <- bind_rows(HS.samps.df.down, HS.samps.temp)

}

num.samps_down <- nrow(HS.samps.df.down)

####----------------Construct pairwise comparison matrix after downsampling----------####
pairwise.df.HS <- data.frame(t(combn(HS.samps.df.down$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
  dplyr::select(older.sib = indv.name, 
                older.sib.birth = birth.year, 
                older.sib.mom = mother.x,
                older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(filtered.samples.HS.df, by = "indv.name") %>% 
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

pairwise.df_HS.filt <- pairwise.df_all.info.HS %>% 
  dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
  dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))


#Extract positive half-sib comparisons
positives.HS <- pairwise.df_HS.filt %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
  mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
#nrow(positives.HS)

#Sex-specific half-sib
mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  dplyr::rename(yes = freq)

dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
  dplyr::select(older.sib.birth, younger.sib.birth) %>%
  plyr::count() %>% 
  dplyr::rename(yes = freq)


#Make dataframes for negative comparisons
hsp.negs <- pairwise.df_HS.filt %>%
  dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  as_tibble()

mom_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  mutate(mom.oncycle = ifelse(BI == "even", 1, 0))

dad_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) 

mom.pos_down <- sum(mom_positives.HS$yes)
dad.pos_down <- sum(dad_positives.HS$yes)
}


#----------------------------Fit JAGS model to Lemon Shark Data------------------------------
#--------------------Specify which JAGS model to use---------------------------------#
#jags_file <- paste0(jags.model_location, "HS.only_noLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "HSPOP_noLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "HS.only_wideLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "HSPOP_wideLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "HSPOP_narrowLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "LizModel_noLambda.txt")
#jags_file <- paste0(jags.model_location, "LizModel_Lambda.txt")
#jags_file <- paste0(jags.model_location, "MHSP.only_noLambda_Skip_model.txt")
jags_file <- paste0(jags.model_location, "MHSP.only_narrowLambda_Skip_model.txt")

#JAGS parameters
jags_params = c("Nf", "Nfb2", "survival", "psi", "lambda") #w lambda
#jags_params = c("Nf", "survival", "psi") #w/o lambda
mating.periodicity <- 2 #Will be translated to a

#Adult.survival <- 0.85 #If wanting to fix survival, or give it a tighter prior

#----------------------- MCMC & model parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

#----------------------------Fit JAGS model-------------------------------#
if(sum(mom_comps.HS$yes) > 0){


  set.seed(rseed.iter)
  #MHSP only model
  source("~/R/working_directory/LemonSharkCKMR/03_Lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")

  if(downsample == "yes"){
results.temp <- model.summary2 %>% mutate(min.year = i, max.year = j) %>% 
  mutate(years_sampled = j - i,
         mom_HSPs = sum(mom_comps.HS$yes),
         estimation_year = estimation.year,
         num.samps_all = num.samps_all,
         num.samps_down = num.samps_down,
         mom_HSPS.all = mom.pos_all,
         mom_HSPs.down = mom.pos_down,
         seed = rseed.iter)
  } else{
    results.temp <- model.summary2 %>% mutate(min.year = i, max.year = j) %>% 
      mutate(years_sampled = j - i,
             mom_HSPs = sum(mom_comps.HS$yes),
             estimation_year = estimation.year,
             num.samps_all = num.samps_all,
             seed = rseed.iter)
  }

results <- bind_rows(results, results.temp)
post.samps_list[[rep]] <- post
names(post.samps_list)[[rep]] <- paste0("yrs_", i, "_to_", j)

print(paste0("Finished comparison ", i, " to ", j))
} else {next}
}
}

results2 <- results %>% mutate(years_sampled = years_sampled+1)

write_csv(results2, file = "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/Model.results/CKMR_results_2023.06.22_FourYearIntervals_wLambd.csv")

#write_rds(post.samps_list, file = "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/Model.results/post_samples")







#######################################################################################-
#--------------Fit CKMR model to full dataset--------------------------
#######################################################################################-
filtered.samples.HS.df <- NoFullSibs.df %>% 
  group_by(mother.x, birth.year) %>% 
  slice_sample(n=1) %>% 
  dplyr::arrange(birth.year) %>% 
  mutate(birth.year = as.numeric(birth.year)) %>%  #Arrange so older sib always comes first in pairwise comparison matrix
  ungroup()
  
#To keep the same individual each time  
  # dplyr::distinct(mother.x, birth.year, .keep_all = TRUE) %>% #Keep just one individual per litter
  # dplyr::arrange(birth.year) %>% 
  # mutate(birth.year = as.numeric(birth.year)) #Arrange so older sib always comes first in pairwise comparison matrix

pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
  dplyr::select(older.sib = indv.name, 
                older.sib.birth = birth.year, 
                older.sib.mom = mother.x,
                older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(filtered.samples.HS.df, by = "indv.name") %>% 
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

pairwise.df_HS.filt <- pairwise.df_all.info.HS %>% 
  dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
  dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))


#Extract positive half-sib comparisons
positives.HS <- pairwise.df_HS.filt %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
  mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
#nrow(positives.HS)

#Sex-specific half-sib
mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  dplyr::rename(yes = freq)

dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
  dplyr::select(older.sib.birth, younger.sib.birth) %>%
  plyr::count() %>% 
  dplyr::rename(yes = freq)


#Make dataframes for negative comparisons
hsp.negs <- pairwise.df_HS.filt %>%
  dplyr::filter(younger.sib.mom != older.sib.mom) %>% #Only include maternal negatives
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  as_tibble()

mom_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  mutate(mom.oncycle = ifelse(BI == "even", 1, 0))

dad_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) 

(num.samps_all <- nrow(filtered.samples.HS.df)) #Save the total number of samples
(mom.pos_all <- sum(mom_positives.HS$yes))


jags_params = c("Nf", "survival", "psi") #w/o lambda
#jags_params = c("Nf", "survival", "psi", "lambda") #w/ lambda
mating.periodicity <- 2 #Will be translated to a

#----------------------- MCMC & model parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

jags_file <- paste0(jags.model_location, "MHSP.only_noLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "MHSP.only_narrowLambda_Skip_model.txt")

source("~/R/working_directory/LemonSharkCKMR/Objective.5_lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")

model.summary2








#######################################################################################-
#--------------Downsample loop--------------------------
#######################################################################################-
estimation.year <- 2000

years <- c(min(NoFullSibs.df$birth.year):max(NoFullSibs.df$birth.year))

results <- NULL
downsample <- "yes"
max.HSPs <- 150

rseeds <- sample(c(1:10000000), size = 200)
  
for(i in 1:200){
  rseed.iter <- rseeds[i]
  set.seed(rseed.iter)
  
  filtered.samples.HS.df <- NoFullSibs.df %>% 
    group_by(mother.x, birth.year) %>% 
    slice_sample(n=1) %>% 
    dplyr::arrange(birth.year) %>% 
    mutate(birth.year = as.numeric(birth.year)) %>%  #Arrange so older sib always comes first in pairwise comparison matrix
    ungroup()
  
  #To keep the same individual each time  
  # dplyr::distinct(mother.x, birth.year, .keep_all = TRUE) %>% #Keep just one individual per litter
  # dplyr::arrange(birth.year) %>% 
  # mutate(birth.year = as.numeric(birth.year)) #Arrange so older sib always comes first in pairwise comparison matrix
  
  pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
  colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(pairwise.df.HS)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
    dplyr::select(older.sib = indv.name, 
                  older.sib.birth = birth.year, 
                  older.sib.mom = mother.x,
                  older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(filtered.samples.HS.df, by = "indv.name") %>% 
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
  
  pairwise.df_HS.filt <- pairwise.df_all.info.HS %>% 
    dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
    dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
  
  
  #Extract positive half-sib comparisons
  positives.HS <- pairwise.df_HS.filt %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
    mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
  #nrow(positives.HS)
  
  #Sex-specific half-sib
  mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
    dplyr::select(older.sib.birth, younger.sib.birth) %>%
    plyr::count() %>% 
    dplyr::rename(yes = freq)
  
  
  #Make dataframes for negative comparisons
  hsp.negs <- pairwise.df_HS.filt %>%
    dplyr::filter(younger.sib.mom != older.sib.mom) %>% #Only include maternal negatives
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    as_tibble()
  
  mom_comps.HS <- hsp.negs %>% 
    dplyr::rename(no = freq) %>% 
    left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
    mutate(mom.oncycle = ifelse(BI == "even", 1, 0))
  
  dad_comps.HS <- hsp.negs %>% 
    dplyr::rename(no = freq) %>% 
    left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no) %>% 
    dplyr::rename(ref.year = younger.sib.birth) %>% 
    mutate(pop.growth.yrs = ref.year - estimation.year) 
  
  (num.samps_all <- nrow(filtered.samples.HS.df)) #Save the total number of samples
  (mom.pos_all <- sum(mom_positives.HS$yes))
  
  
  jags_params = c("Nf", "survival", "psi") #w/o lambda
  #jags_params = c("Nf", "survival", "psi", "lambda") #w/ lambda
  mating.periodicity <- 2 #Will be translated to a
  
  #----------------------- MCMC & model parameters ----------------------#
  ni <- 40000 # number of post-burn-in samples per chain
  nb <- 50000 # number of burn-in samples
  nt <- 20     # thinning rate
  nc <- 2      # number of chains
  
  jags_file <- paste0(jags.model_location, "MHSP.only_noLambda_Skip_model.txt")
  #jags_file <- paste0(jags.model_location, "MHSP.only_narrowLambda_Skip_model.txt")
  
#----------------Downsample the dataset--------------------------------------
  
#Downsample for HSPs
#Calculate proportion of samples born in each birth year
HS.samps.props <- filtered.samples.HS.df %>% group_by(birth.year) %>%
  summarize(num = n()) %>%
  ungroup() %>%
  mutate(prop = num/sum(num))

#Save vector of proportions corresponding to the correct index in the HS sample dataframe
HS.props.vec <- filtered.samples.HS.df %>% left_join(HS.samps.props, by = "birth.year") %>%
  pull(prop)

#Calculate proportion of positive HS comparisons for mom
HS_comps.mom.yes <- mom_comps.HS %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

HS_comps.mom.all <- mom_comps.HS %>%
  summarize(all = sum(all)) %>% 
  pull(all)

HS_prop.mom.yes <- HS_comps.mom.yes/HS_comps.mom.all #Calculate proportion of positive comparisons

HS_target.comps = round(max.HSPs/HS_prop.mom.yes, 0)
HS_target.samples <- round(sqrt(HS_target.comps), 0)

HS.samps.df.down <- NULL #Initialize so I don't get an error

if(HS_comps.mom.yes > max.HSPs){
  HS.samps.df.down <- filtered.samples.HS.df %>% slice_sample(n = HS_target.samples, weight_by = HS.props.vec) %>% arrange(birth.year) #arrange by birth year so they're in the correct order for combn
} else{
  HS.samps.df.down <- filtered.samples.HS.df
}

(num.samps_down <- nrow(HS.samps.df.down))

pairwise.df.HS <- data.frame(t(combn(HS.samps.df.down$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- HS.samps.df.down %>%
  dplyr::select(older.sib = indv.name, 
                older.sib.birth = birth.year, 
                older.sib.mom = mother.x,
                older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(filtered.samples.HS.df, by = "indv.name") %>% 
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

pairwise.df_HS.filt <- pairwise.df_all.info.HS %>% 
  dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
  dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))


#Extract positive half-sib comparisons
positives.HS <- pairwise.df_HS.filt %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
  mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
#nrow(positives.HS)

#Sex-specific half-sib
mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  dplyr::rename(yes = freq)

dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
  dplyr::select(older.sib.birth, younger.sib.birth) %>%
  plyr::count() %>% 
  dplyr::rename(yes = freq)


#Make dataframes for negative comparisons
hsp.negs <- pairwise.df_HS.filt %>%
  dplyr::filter(younger.sib.mom != older.sib.mom) %>% #Only include maternal negatives
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  as_tibble()

mom_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  mutate(BI = ifelse(year_gap %% 2 == 0, "even", "odd")) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  mutate(mom.oncycle = ifelse(BI == "even", 1, 0))

dad_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no) %>% 
  dplyr::rename(ref.year = younger.sib.birth) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) 

jags_params = c("Nf", "survival", "psi") #w/o lambda
#jags_params = c("Nf", "survival", "psi", "lambda") #w/ lambda
mating.periodicity <- 2 #Will be translated to a

#----------------------- MCMC & model parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

jags_file <- paste0(jags.model_location, "MHSP.only_noLambda_Skip_model.txt")
#jags_file <- paste0(jags.model_location, "MHSP.only_narrowLambda_Skip_model.txt")

#----------------------------Fit JAGS model-------------------------------#
if(sum(mom_comps.HS$yes) > 0){
  
  #MHSP only model
  source("~/R/working_directory/LemonSharkCKMR/Objective.5_lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")
  
  (mom.pos_down <- sum(mom_positives.HS$yes))
  
  results.temp <- model.summary2 %>% mutate(iteration = i) %>% 
    mutate(mom_HSPs = sum(mom_comps.HS$yes),
           num.samps_down = num.samps_down,
           mom_HSPs.down = mom.pos_down,
           seed = rseed.iter)
  
  results <- bind_rows(results, results.temp)

  print(paste0("Finished iteration ", i))
} else {next}
}


write_csv(results, file = "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/Model.results/CKMR_results_2022.12.15_downsampleFullDatasetLoop.csv")

results %>% as_tibble() %>% group_by(parameter) %>% 
  summarize(mean.median = mean(Q50),
            median.median = median(Q50),
            sd.media = sd(Q50),
            mean.HPD2.5 = mean(HPD2.5),
            meanHPD97.5 = mean(HPD97.5),
            mean.samps = mean(num.samps_down),
            sd.samps = sd(num.samps_down),
            mean.HSPs = mean(mom_HSPs.down),
            sd.HSPs = sd(mom_HSPs.down)
            )












#----------------Minimum number of moms and dads-------------------------------
#Calculate minimum number of mothers per year
mother_ref %>% 
  drop_na(mother.x) %>% 
  group_by(birth.year) %>% 
  summarize(min.moms = n()) %>% 
  dplyr::arrange(desc(min.moms)) %>% 
  dplyr::arrange(birth.year) %>% 
  View()
  
#Calculate average number of moms per year
mother_ref %>% 
  drop_na(mother.x) %>% 
  dplyr::filter(birth.year > 1994) %>% #Remove years before which we estimate abundance
  group_by(birth.year) %>% 
  summarize(min.moms = n()) %>% 
  dplyr::arrange(desc(min.moms)) %>% 
  dplyr::arrange(birth.year) %>% 
  summarize(mean.moms.per.yr = mean(min.moms))

mother_ref %>% dplyr::filter(is.na(mother.x) == TRUE)


father_ref %>% 
  group_by(birth.year) %>% 
  summarize(min.dads = n()) %>% 
  dplyr::arrange(desc(min.dads)) %>% 
  dplyr::filter(birth.year == estimation.year)

father_ref %>% dplyr::filter(is.na(father.x) == TRUE)

filtered.samples.HS.df %>% dplyr::distinct(father.x) %>% nrow()
filtered.samples.HS.df %>% dplyr::distinct(mother.x) %>% nrow()


#----------------Examine things more closely----------------
mom_comps.HS %>% group_by(BI) %>% 
  summarize(sum = sum(yes))








################# Simulate adult info ################
# Simulate and assign ages to adults based on life history
# Max age: two options: 
# 37 from https://www.saveourseasmagazine.com/lemon-sharks-old/
# 25 from White et. al. 2014
# Juvenile survival values from Dibattista 2007: 0.55 (mort = 0.45)
# Adult survival value from White et. al. 2014: 0.85 (mort=0.15)
# Litter size from Feldheim et. al. 2002 (avg=6, max=18)
# Age-at-maturity from Brown and Gruber 1988 (11.6 for M, 12.7 for F)
# About 77 juvenile lemon sharks inhabit the nursery at a given time (White et al 2014). So, 10sqrt(N) = 88 samples total

####Simulate ages for parents####
#Establish stable age distribution FOR ADULTS ONLY (We have ages for juveniles)
#VARIABLES
max.age <- 25
m.repro.age <- 12
mating.periodicity <- 2
m_adult_age <- c(m.repro.age:max.age) #Set ages at which males are mature -- used to simulate ages
m_mat <- c(rep(0,(m.repro.age-1)), rep(1,(max.age - m.repro.age + 1))) #Set proportion of mature males at each age -- assumes knife-edge maturity
Surv_age_m <- rep(0.85,length(m_adult_age)) #Set survival for adults
Prop_age_m <- rep(1, length(m_adult_age)) #Initialize variable Prop_age_m
for(iage in 2:length(Prop_age_m)) Prop_age_m[iage]=Prop_age_m[iage-1]*Surv_age_m[iage-1] #Set proportion of animals from each cohort surviving to the next age class based on Surv_age
Prop_age_m=Prop_age_m/sum(Prop_age_m)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

Uniq_Father <- unique(Father_ref$Father) #Save unique fathers

#Same as above but with mothers (which mature one year later)
f.repro.age <- 13
f_adult_age <- c(f.repro.age:max.age)
f_mat <- c(rep(0,(f.repro.age-1)), rep(1,(max.age - f.repro.age + 1)))
Surv_age_f <- rep(0.85,length(f_adult_age))
Prop_age_f <- rep(1, length(f_adult_age))
for(iage in 2:length(Prop_age_f)) Prop_age_f[iage]=Prop_age_f[iage-1]*Surv_age_f[iage-1] 
Prop_age_f=Prop_age_f/sum(Prop_age_f)

Uniq_Mother <- unique(Mother_ref$Mother)

#Simulate ages of mothers
#Commented out below 9/19
#Mother_age <- as.data.frame(cbind(Uniq_Mother,Mother_age), stringsAsFactors=FALSE)
#colnames(Mother_age) <- c("Mother", "Age")

# Want to base estimate on animals caught between 2004-2007 (the most heavily sampled years)

####Set up dataframe for population we will sample####
#Subset for animals captured in 2004-2007
t_start <- 2004
t_end <- 2007
lems <- subset(juv_ref, subset=Capture_Year %in% c(t_start:t_end), select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex))
#head(juv_ref)
#View(lems)
#Assign ages to parents in subset of data
#Mothers
#DELETE one of the for loops below -- I think the second will work fine b/c I replaced blanks with Unknown, but testing
#for(i in 1:length(lems[,1])) {
#  if(lems$Mother[i]!=""){
#    lems$Mother_age[i] <- Mother_age[which(Mother_age$Mother==lems$Mother[i]),2]
#  }
#}

#for(i in 1:length(lems[,1])) {
#    lems$Mother_age[i] <- Mother_age[which(Mother_age$Mother==lems$Mother[i]),2]
#}

#Change age into DOB (will be used in CKMR model)
#lems$Mother_DOB <- t_start-as.numeric(lems$Mother_age)
#min(lems$Mother_DOB)

#Fathers
#DELETE one of the below for loops (I think the second will work, but testing)
#for(i in 1:length(lems[,1])) {
 # if(lems$Father[i]!=""){
#    lems$Father_age[i] <- Father_age[which(Father_age$Father==lems$Father[i]),2]
#  }
#}

#for(i in 1:length(lems[,1])) {
#    lems$Father_age[i] <- Father_age[which(Father_age$Father==lems$Father[i]),2]
#}

#lems$Father_DOB <- t_start-as.numeric(lems$Father_age)
#head(lems)

#add the +1 so the first and last years of sampling are counted
yr_1 <- t_start-max_age
yrs <- c(seq(yr_1,t_end))
samp_yrs <- c(t_start:t_end)
study_yrs <- c(1:((t_end+1)-yr_1))
samp_study_yrs = c(26:29)
#study_yrs[which(yrs == lems$DOB[1])]


#Change birth year for juveniles to birth year relative to study
Juv_study_DOB <- c()
for(i in 1:length(lems[,1])){
  Juv_study_DOB[i] <- study_yrs[which(yrs == lems$DOB[i])]
}

Juv_lems <- as.data.frame(cbind(lems$Indiv_ID, Juv_study_DOB, lems$Capture_Year, lems$Sex), stringsAsFactors = FALSE)
colnames(Juv_lems) <- c("Indiv", "Juv_study_DOB", "Capture_Year", "Sex")
Juv_lems2 <- Juv_lems[,c(1,3:4)]

Fat_lems <- as.data.frame(cbind(lems$Father, lems$Capture_Year, rep("M", length(lems[,1]))), stringsAsFactors = FALSE)

Mot_lems <- as.data.frame(cbind(lems$Mother, lems$Capture_Year, rep("F", length(lems[,1]))), stringsAsFactors = FALSE)

colnames(Fat_lems) = colnames(Mot_lems) = colnames(Juv_lems2) = c("Indiv", "Capture_Year", "Sex")

#Add column for capture year that corresponds to study year (e.g. 26 instead of 2004)
for(i in 1:length(Juv_lems2[,1])){ 
  Juv_lems2$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Juv_lems2$Capture_Year[i])]
}
for(i in 1:length(Fat_lems[,1])){
  Fat_lems$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Fat_lems$Capture_Year[i])]
}
for(i in 1:length(Mot_lems[,1])){
  Mot_lems$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Mot_lems$Capture_Year[i])]
}

#colnames(Juv_lems) <- c("Indiv", "Capture_Year", "Sex")
#head(Fat_lems)
#head(Mot_lems)
#length(Mot_lems[,1])
#Remove duplicate entries (There are no duplicate juveniles - double-checked 9/13/19)
#Troubleshoot - Below will only keep the first instance of each, so more captures in first year (because assigned capture year to all individuals)
Mot_lems <- Mot_lems[!duplicated(Mot_lems$Indiv),]
#Mot_lems$uniqID <- NULL
Mot_lems <- Mot_lems[-c(which(Mot_lems[,1]=="Unknown")),] #Remove Unknown individuals
#length(Mot_lems[,1]) #Check number of unique mothers

#Fat_lems$uniqID <- paste0(Fat_lems$Indiv,Fat_lems$DOB)
#length(Fat_lems[,1])
Fat_lems <- Fat_lems[!duplicated(Fat_lems$Indiv),]
#Fat_lems$uniqID <- NULL
Fat_lems <- Fat_lems[-c(which(Fat_lems[,1]=="Unknown")),] #Remove unknown individuals
#length(Fat_lems[,1])

#Combine data frames
pop_df <- c()
pop_df <- rbind(Juv_lems2, Fat_lems, Mot_lems)
#head(pop_df)
#tail(pop_df)

####Sample from population####
#Want 88 samples total, spread over 4 years, so 22 samples per year
#White estimates 77 juveniles in the population
#VARIABLES to move to top
Estimated_truth <- 77
n_samples <- round(10*sqrt(Estimated_truth), 0)
samps=c()

#Make vector called samps that randomly samples the population
for(i in 1:4){
  cp_yr = c(2004:2007)[i]
  samps <- c(samps, sample(pop_df[which(pop_df$Capture_Year==cp_yr),1], size=n_samples/length(samp_study_yrs), replace=FALSE))
}

#length(samps)
#head(samps)
#head(pop_df)

#Create sample dataframe that has year of capture(death) year of birth, and sex
Samp_birth = c() #Initialize vectors
Samp_cap_yr = c()
Samp_sex = c()
  
#head(Juv_lems)
#head(Fat_lems)
#study_yrs
#samp_yrs

#samp_study_yrs[which(samp_yrs==Fat_lems$Capture_Year[1])]
#head(pop_df)
#Add capture year that corresponds to study year (e.g. 26 instead of 2004)
#head(Juv_lems)
#Simulate birth year (based on study year)
#as.numeric(Juv_lems$DOB[which(Juv_lems==samps[2])])
for(i in 1:length(samps)){
  if(samps[i] %in% Juv_lems$Indiv){
    Samp_birth[i] <- as.numeric(Juv_lems$Juv_study_DOB[which(Juv_lems==samps[i])])
  } else if(samps[i] %in% Fat_lems$Indiv){
    Samp_birth[i] <- Fat_lems$Capt_study_yr[which(Fat_lems$Indiv==samps[i])] - sample(m_adult_age, replace=TRUE, prob=Prop_age_m, size=1)
  } else if(samps[i] %in% Mot_lems$Indiv){
    Samp_birth[i] <- Mot_lems$Capt_study_yr[which(Mot_lems$Indiv==samps[i])] - sample(f_adult_age, replace=TRUE, prob=Prop_age_f, size=1)
  }
}
#Samp_birth   

#head(pop_df)
#Samp_birth[i] = as.numeric(pop_df[which(pop_df$Indiv==samps[i]),2]) #search pop_df for the Indiv ID and extract the associated birth year
for(i in 1:length(samps))Samp_cap_yr[i] = as.numeric(pop_df[which(pop_df$Indiv==samps[i]),4]) #search pop_df for the Indiv ID and extract the associated capture year
for(i in 1:length(samps))Samp_sex[i] = pop_df[which(pop_df$Indiv==samps[i]),3] #search pop_df for the Indiv ID and extract the associated sex
#Combine vectors into dataframe and rename columns
Samples <- as.data.frame(cbind(samps,Samp_birth,Samp_cap_yr, Samp_sex), stringsAsFactors = FALSE)
colnames(Samples) <- c("Indiv_ID", "Birth", "Death", "Sex")
#head(Samples)

Samples[,2] <- as.numeric(Samples[,2]) #Change columns to numeric; changing earlier doesn't work
Samples[,3] <- as.numeric(Samples[,3])
#Samples

#Assign father and mother for each sampled individual
#Can below be run with foreach? 
for(i in 1:length(Samples[,1])){
  for(j in 1:length(juv_ref[,1])){
    if(Samples[i,1]==juv_ref[j,1]){
      Samples$Father[i]=juv_ref[j,6]
      Samples$Mother[i]=juv_ref[j,7]
    } #else if(Samples[i,1]==juv_ref[j,6] | Samples[i,1]==juv_ref[j,7]){
     #Samples$Offspring[i]=juv_ref[j,1] 
    #}
  }
}

Samples$Dad_birth=Samples$Dad_death=Samples$Mom_birth=Samples$Mom_death=NA #Set birth and death year of parents to NA so missing data are easier to work with downstream
Samples$Mom=Samples$Dad=0 #Set to 0 so we can add a 1 to rows in which samples have IDed parents

#If the sampled animal also is a parent to another sampled animal, assign it a birth and death date
for(i in 1:length(Samples[,1])){
  if(Samples$Indiv_ID[i] %in% Samples$Father){
    Samples$Dad_birth[which(Samples$Father==Samples$Indiv_ID[i])] <- Samples$Birth[i]
    Samples$Dad_death[which(Samples$Father==Samples$Indiv_ID[i])] <- Samples$Death[i]
  }
  if(Samples$Indiv_ID[i] %in% Samples$Mother){
    Samples$Mom_birth[which(Samples$Mother==Samples$Indiv_ID[i])] <- Samples$Birth[i]
    Samples$Mom_death[which(Samples$Mother==Samples$Indiv_ID[i])] <- Samples$Death[i]
  }
}
#View(Samples)
#Set Mom and Data column to 1 if the sampled individual in the row has a sampled parent
for(i in 1:length(Samples[,1])){
  if(is.na(Samples$Dad_birth[i])==FALSE) Samples$Dad[i]=1
  if(is.na(Samples$Mom_birth[i])==FALSE) Samples$Mom[i]=1
}

#VARIABLE
####NOTE TO SELF - Streamline script so we skip the year and just give age according to study####
#Set first year from which we might see individuals in the data
#tail(Samples)
#for(i in 1:length(Samples[,1])) Samples$Birth[i] <- max(Samples$Birth[i]-yr_1,1) #Some samples have a birth year of 0 if not for max
#Samples$Death=Samples$Death-yr_1+1 #Add one to make sure dates of sampling are within specified range
#for(i in 1:length(Samples[,1])) Samples$Mom_birth[i] <- max(Samples$Mom_birth[i]-yr_1)
#Samples$Mom_death=Samples$Mom_death-yr_1+1
#for(i in 1:length(Samples[,1])) Samples$Dad_birth[i] <- max(Samples$Dad_birth[i]-yr_1)
#Samples$Dad_death=Samples$Dad_death-yr_1+1
#head(Samples)
#View(Samples)

####Set up pairwise comparison matrix####
#We don't want to compare individuals with themselves (hence -1), and we don't want to do each comparison twice (divide by 2)
Data <- data.frame(matrix(0,nrow=n_samples*(n_samples-1)/2,ncol=5)) # massive array for pairwise comparisons

#columns are: adult birth; adult death year; young birth; adult sex; probability of POP
counter=1
#View(Samples)
#brute force - could probably be improved w/ expand.grid or something. Fills first four columns of Data
for(iind1 in 1:(n_samples-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(iind2 in (iind1+1):n_samples){#iind2 starts at 2 and counts to n_samples
    if(Samples[iind1,2]>Samples[iind2,2]){ #compare birth years of each combination of samples. If birth year is greater for iind1 i.e. if iind2 was born before iind1. Lethal sampling, so iind2 starts at the number after iind1, and there is no need to go back before.
      #Fill Data array:
      #Data[,1] = adult birth,
      #Data[,2] = adult death year
      #Data[,3] = young birth
      #Data[,4] = adult sex
      Data[counter,1]=Samples[iind2,2] #put birth year of iind2 in adult birth year
      Data[counter,2]=Samples[iind2,3] #put year of capture of iind2 in adult death year (remember, lethal sampling)
      Data[counter,3]=Samples[iind1,2] #put birth year of iind1 in offspring birth year
      Data[counter,4]=Samples[iind2,4] #put sex of iind2 in adult sex
    }
    else{ #birth year greater for second value - reverse above
      Data[counter,1]=Samples[iind1,2]
      Data[counter,2]=Samples[iind1,3]
      Data[counter,3]=Samples[iind2,2]
      Data[counter,4]=Samples[iind1,4]
    }
    counter=counter+1
  }
}

colnames(Data) <- c("Adult_birth","Adult_death","Offspring_birth","Adult_sex", "Matches")
#length(Data[,1])
#head(Data)

############### Assign parentage ################
#Create separate dataframes for sampled moms and sampled dads
#DELETE below?
#SM2 <- aggregate(Sampled_Moms, by=list("Birth", "Mom_death", "Mom_birth"), FUN = "sum")
#head(Data)
dm <- Data[Data$Adult_sex == "F",] #Subset pairwise matrix for mother comparisons
dd <- Data[Data$Adult_sex == "M",] #Subset pairwise matrix for father comparisons

#Create column for unique combinations of adult birth, adult death, and offspring birth year
dm$MasterID <- paste(dm$Adult_birth, dm$Adult_death, dm$Offspring_birth, sep="-")
dd$MasterID <- paste(dd$Adult_birth, dd$Adult_death, dd$Offspring_birth, sep="-")
#head(dd)

#Count unique combos of adult birth, adult death, and offspring birth AND Master ID (should be the same)
dm2 <- plyr::count(dm[,c(1:3,6)])
dd2 <- plyr::count(dd[,c(1:3,6)])
colnames(dm2)[5] = colnames(dd2)[5] <- "No_matches" #Change frequency column to No_matches (assume all comparisons are no matches - will subtract matches later)
#head(dd2)

#Subset sampled moms and sampled dads
Sampled_Moms <- Samples[Samples$Mom == 1,]
####Troubleshoot####
#Correct for instances where adult is sampled and mother's age is less

Sampled_Dads <- Samples[Samples$Dad == 1,]
#View(Samples)

#Create dataframe of sampled moms (sm) that includes a column that counts the occurrances of specific combinations of offspring birth, adult birth, and adult death
sm <- Sampled_Moms %>% 
  #complete(f_year, month, age) %>%  #Expands each years data to include all months and ages
  group_by(Birth, Mom_birth, Mom_death) %>%  #ID which columns to group by
  summarise(Matches = sum(!is.na(Mom)))

#Create dataframe of sampled dads (sd) that includes a column that counts the occurrances of specific combinations of offspring birth, adult birth, and adult death
sd <- Sampled_Dads %>% 
  #complete(f_year, month, age) %>%  #Expands each years data to include all months and ages
  group_by(Birth, Dad_birth, Dad_death) %>%  #ID which columns to group by
  summarise(Matches = sum(!is.na(Dad)))

#Create MasterID column to merge dataframes with matches (sm & sd) with dataframes with all comparisons (dm2, dd2)
sd$MasterID <- paste(sd$Dad_birth, sd$Dad_death, sd$Birth, sep="-")
sm$MasterID <- paste(sm$Mom_birth, sm$Mom_death, sm$Birth, sep="-")

all_dads <- merge(dd2, sd[,4:5], by="MasterID", all=TRUE) #Merge father dataframes
all_dads$Matches[is.na(all_dads$Matches)==TRUE] <- 0 #Put 0 in all Matches cells with NA
all_dads$No_matches <- all_dads$No_matches-all_dads$Matches #Calculate number of no matches

#sum(all_dads$Matches)

all_moms <- merge(dm2, sm[,4:5], by="MasterID", all=TRUE) #Merge mother dataframes
all_moms$Matches[is.na(all_moms$Matches)==TRUE] <- 0 #Put 0 in all Matches cells with NA
all_moms$No_matches <- all_moms$No_matches-all_moms$Matches #Calculate number of no matches

#sum(all_moms$Matches, na.rm=TRUE)
#which(is.na(all_moms[,2])==TRUE) #Shouldn't be any
#which(is.na(all_dads[,2])==TRUE) #Shouldn't be any
#which(is.na(sm[,2])==TRUE)
#which(is.na(dm[,2])==TRUE)

#head(all_moms)
#Create dataframes for CKMR model!
Data_mom_yes <- all_moms[all_moms$Matches > 0, c(2:4,6)] 
Data_dad_yes <- all_dads[all_dads$Matches > 0, c(2:4,6)]
Data_dad_no <- all_dads[all_dads$Matches == 0, c(2:4,5)]
Data_mom_no <- all_moms[all_moms$Matches == 0, c(2:4,5)]

#In case offspring are born after parents die, band-aid fix to allow the model to run
#Change either adult birth year or death year so 1) adult does not die before offspring's birth, and 2) adult is not older than max age

for(i in 1:length(Data_mom_yes[,1])){
  if(Data_mom_yes$Offspring_birth[i] > Data_mom_yes$Adult_death[i] & Data_mom_yes$Offspring_birth[i]-Data_mom_yes$Adult_birth[i] > max_age){
  Data_mom_yes$Adult_birth[i] <- Data_mom_yes$Offspring_birth[i]-max_age
  Data_mom_yes$Adult_death[i] <- Data_mom_yes$Offspring_birth[i]
} else if(Data_mom_yes$Offspring_birth[i] > Data_mom_yes$Adult_death[i] & Data_mom_yes$Offspring_birth[i]-Data_mom_yes$Adult_birth[i] <= max_age) {
  Data_mom_yes$Adult_death[i] <- Data_mom_yes$Offspring_birth[i]
} else Data_mom_yes$Adult_death[i] <- Data_mom_yes$Adult_death[i]
}

for(i in 1:length(Data_dad_yes[,1])){
  if(Data_dad_yes$Offspring_birth[i] > Data_dad_yes$Adult_death[i] & Data_dad_yes$Offspring_birth[i]-Data_dad_yes$Adult_birth[i] > max_age){
    Data_dad_yes$Adult_birth[i] <- Data_dad_yes$Offspring_birth[i]-max_age
    Data_dad_yes$Adult_death[i] <- Data_dad_yes$Offspring_birth[i]
  } else if(Data_dad_yes$Offspring_birth[i] > Data_dad_yes$Adult_death[i] & Data_dad_yes$Offspring_birth[i]-Data_dad_yes$Adult_birth[i] <= max_age) {
    Data_dad_yes$Adult_death[i] <- Data_dad_yes$Offspring_birth[i]
  } else Data_dad_yes$Adult_death[i] <- Data_dad_yes$Adult_death[i]
}
#In case adult is older than max_age, band-aid fix to allow the model to run
#DELETE - think it's fixed?
#Data_dad_yes$Adult_birth <- ifelse(Data_dad_yes$Adult_death-Data_dad_yes$Adult_birth > max_age, Data_dad_yes$Adult_death-max_age, Data_dad_yes$Adult_birth)
#Data_mom_yes$Adult_birth <- ifelse(Data_mom_yes$Adult_death-Data_mom_yes$Adult_birth > max_age, Data_mom_yes$Adult_death-max_age, Data_mom_yes$Adult_birth)

#head(Data_mom)
#Data_mom2 <- Data_mom[!duplicated(Data_mom),]
#Data_dad2 <- Data_dad[!duplicated(Data_dad),]
# nrow(Data_dad2)
# 
# Data_mom3 <- merge(Data_mom2[,-5],SM2[,4:5], by = "MasterID", all = TRUE)
# Data_mom3$Matches[is.na(Data_mom3$Matches)==TRUE] <- 0
# table(Data_mom3$Matches)
# 
# Data_dad3 <- merge(Data_dad2[,-5],SD2[,4:5], by = "MasterID", all = TRUE)
# Data_dad3$Matches[is.na(Data_dad3$Matches)==TRUE] <- 0
# table(Data_dad3$Matches)



####Old code
# for(i in 1:length(Samples[,1])){
#   if(Samples$Mom[i]==1){
#     Data[which(Data[,1]==Samples$Mom_birth[i] & Data[,2]==Samples$Mom_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="F"),5] = Data[which(Data[,1]==Samples$Mom_birth[i] & Data[,2]==Samples$Mom_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="F"),5] + 1
#   } 
#   if(Samples$Dad[i]==1){
#     Data[which(Data[,1]==Samples$Dad_birth[i] & Data[,2]==Samples$Dad_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="M"),5]=Data[which(Data[,1]==Samples$Dad_birth[i] & Data[,2]==Samples$Dad_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="M"),5]+1
#   }
# }
# head(Data)
# Data_dad_no = Data[which(Data[,4]=="M" & Data[,5]==0),1:4]
# Data_dad_yes = Data[which(Data[,4]=="M" & Data[,5]>0),1:4]
# Data_mom_no = Data[which(Data[,4]=="F" & Data[,5]==0),1:4]
# Data_mom_yes = Data[which(Data[,4]=="F" & Data[,5]>0),1:4]
# Data_dad_no=count(Data_dad_no[,1:3])
# Data_mom_no=count(Data_mom_no[,1:3])
# Data_dad_yes=count(Data_dad_yes[,1:3])
# Data_mom_yes=count(Data_mom_yes[,1:3])
# sum(Data_mom_yes$freq)