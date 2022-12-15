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
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/Model.output/"
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/Model.results/"

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"

#---------------------------Read in and organize data------------------------------------
purpose <- "CKMR_LS_data"
rm.arb <- TRUE #Do we want to remove juvniles with arbitrary birth years (there are 48/2198 in the Bimini dataset)
estimation.year <- 2000 #Will only be used if using the model with lambda, but the script will throw an error if this object isn't defined (could fix this with a conditional statement, but don't feel like it now :-p).

#Set path to data file
input.file <- "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/data/Main_lemon_shark.csv"

#Save von bertalanffy growth function
vonBert <- vbFuns(param = "Typical")
max.age <- 50
repro.age = 12

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

source("~/R/working_directory/LemonSharkCKMR/Objective.5_lemon_shark_data/functions/Obj5.functions.R")

#Remove unknown mothers, then for each group of full sibs keep just one
NoFullSibs.df <- juv_ref.North %>% drop_na(mother.x) %>% 
  group_by(mother.x, father.x) %>% 
  slice_sample(n = 1)

#To keep the same individual each time
  # distinct(mother.x, father.x, .keep_all = TRUE) %>%
  # as_tibble() #If there is more than one individual with the same mother AND father, then only keep one.


#######################################################################################-
#--------------Run simulations over different subsets of data--------------------------
#######################################################################################-

estimation.year <- 2000

years <- c(min(NoFullSibs.df$birth.year):max(NoFullSibs.df$birth.year))

results <- NULL
post.samps_list <- list()
rep = 0
downsample <- "no"
max.HSPs <- 150

#Regular loop over different periods of estimation
  for(i in min(years):max(years-2)){ #Change to -1 if doing full loop
    for(j in (i + 2): max(years)){
 #   j= i + 4

    estimation.year = j
    
    #When an error and need to start at a certain number
# for(i in 2012:max(years-1)){
#   for(j in (i + 2): max(years)){

        rep = rep+1
    
        filtered.samples.HS.df <- NoFullSibs.df %>% 
          # group_by(mother.x, birth.year) %>% #Can comment out this line and the next to include all HSPs
          # slice_sample(n=1) %>% 
          dplyr::arrange(birth.year) %>% 
          mutate(birth.year = as.numeric(birth.year)) %>%  #Arrange so older sib always comes first in pairwise comparison matrix
          ungroup() %>%
          dplyr::filter(birth.year >= i & birth.year <= j)

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
jags_params = c("Nf", "survival", "psi", "lambda") #w lambda
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

  #MHSP only model
  source("~/R/working_directory/LemonSharkCKMR/Objective.5_lemon_shark_data/functions/Obj5_run.JAGS_MHSP.only.R")

results.temp <- model.summary2 %>% mutate(min.year = i, max.year = j) %>% 
  mutate(years_sampled = j - i,
         mom_HSPs = sum(mom_comps.HS$yes),
         estimation_year = estimation.year,
         num.samps_all = num.samps_all,
       #  num.samps_down = num.samps_down,
         mom_HSPS.all = mom.pos_all)
       #  mom_HSPs.down = mom.pos_down)
  

results <- bind_rows(results, results.temp)
post.samps_list[[rep]] <- post
names(post.samps_list)[[rep]] <- paste0("yrs_", i, "_to_", j)

print(paste0("Finished comparison ", i, " to ", j))
} else {next}
}
}

results2 <- results %>% mutate(years_sampled = years_sampled+1)

write_csv(results2, file = "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/Model.results/CKMR_results_2022.12.15_FullLitter_wLambda.csv")

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