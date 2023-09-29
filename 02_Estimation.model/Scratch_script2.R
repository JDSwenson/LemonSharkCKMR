rents.df <- readRDS(file = paste0(PopSim.location, "parents.breakdown_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

rents.df %>% dplyr::filter(parent.sex == "father") %>% 
  group_by(year) %>% 
  summarize(mean.offspring = mean(num.off),
            min.offspring = min(num.off),
            max.offspring = max(num.off)) %>% 
  summarize(min.off.overall = min(min.offspring),
            max.off.overall = max(max.offspring))


samps.test <- samples.df_all %>% dplyr::filter(iteration == 200, sampling.scheme == "sample.ALL.ages", sample.prop == 2)

imposters.test <- imposters.df %>% dplyr::filter(iteration == 200, sampling.scheme == "sample.ALL.ages", sample.prop == 2) %>% 
  distinct(.keep_all = TRUE)

which(duplicated(imposters.test) | duplicated(imposters.test, fromLast = TRUE))
imposters.test %>% arrange(aunt.unc)








imp.df <- imposters.test %>% mutate(lineage = ifelse(is.na(niece.nephew_maternal.grandmother) == TRUE & is.na(niece.nephew_maternal.grandfather) == TRUE, "maternal",
                                         ifelse(is.na(niece.nephew_paternal.grandmother) == TRUE & is.na(niece.nephew_paternal.grandfather) == TRUE, "paternal", NA))) %>% 
  dplyr::select(aunt.unc_relation, aunt.unc_sex, niece.nephew_relation, niece.nephew_sex, lineage, niece.nephew_paternal.grandmother, niece.nephew_paternal.grandfather, niece.nephew_maternal.grandmother, niece.nephew_maternal.grandfather, aunt.unc, niece.nephew, shared.relation, iteration, seed, sample.prop)

imp.df %>% View()

imp.df %>% dplyr::filter(is.na(lineage) == TRUE)

#Seems like there might be a correlation between aunt/uncle and being related through the maternal or paternal line ... 
imp.df %>% dplyr::count(aunt.unc_relation, lineage)



charlatan_HSPs %>% dplyr::count(sampling.scheme, iteration, sample.prop)
charlatan_HSPs %>% dplyr::arrange(shared.relation) %>% View()
charlatan_HSPs %>% distinct(.keep_all = TRUE)
charlatan_HSPs %>% duplicated() %>% sum()


build.pairwise <- function(filtered.samples.PO.list, filtered.samples.HS.df){
  
  #Different birth years for abundance estimate
  # OffBirth.years.PO.vec <- filtered.samples.PO.list %>% 
  #   bind_rows() %>%
  #   distinct(birth.year) %>% 
  #   pull(birth.year) 
  
  
  #initialize lists for output
  PO.mom.pairwise.list <- list()
  HS.mom.pairwise.list <- list()
  PO.dad.pairwise.list <- list()
  HS.dad.pairwise.list <- list()
  
  #Double-checked PO pairwise comparison code 08/06/2023
  
  #Loop over each dataframe corresponding to the offspring birth year
  for(y in 1:length(filtered.samples.PO.list)){
    #PO.year <- OffBirth.years.PO.vec[y]
    
    #Make vector of potential parents for the year
    PO.parents <- filtered.samples.PO.list[[y]] %>% dplyr::filter(relation == "parent") %>% 
      pull(indv.name)
    
    #Make vector of potential offspring for the year
    PO.offspring <- filtered.samples.PO.list[[y]] %>% dplyr::filter(relation == "offspring") %>% 
      pull(indv.name)
    
    #generate pairwise comparison matrix for all potential parents and offspring for the year
    pairwise.df.PO <- expand.grid(PO.parents, PO.offspring, stringsAsFactors = FALSE) %>% as_tibble()
    colnames(pairwise.df.PO) <- c("indv.name", "offspring.name") #Rename columns so they can easily be joined
    
    head(pairwise.df.PO)
    
    #Prep a second dataframe that reverses the comparison columns
    pwdf.check <- pairwise.df.PO %>% dplyr::rename(offspring.name = indv.name, indv.name = offspring.name) 
    pwdf.check2 <- pwdf.check %>% bind_rows(pairwise.df.PO)
    
    #Check that there are no duplicate comparisons - should be TRUE
    nrow(pairwise.df.PO) == pwdf.check2 %>% distinct(indv.name, offspring.name) %>% nrow()/2
    
    #Check that no individuals are being compared against themselves - should be TRUE
    nrow(pairwise.df.PO %>% dplyr::filter(indv.name == offspring.name)) == 0
    
    #Create dataframe that will be used to extract the birth years for the potential offspring from each comparison using joins.
    Offspring_birthyears.PO <- filtered.samples.PO.list[[y]] %>%
      dplyr::select(offspring.name = indv.name, 
                    offspring.birth = birth.year, 
                    offspring.mom = mother.x, 
                    offspring.dad = father.x)
    
    #Extract metadata for each individual in each comparison
    #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
    pairwise.df_all.info.PO <- pairwise.df.PO %>% left_join(Offspring_birthyears.PO, by = "offspring.name") %>% 
      as_tibble() %>% 
      left_join(filtered.samples.PO.list[[y]], by = "indv.name") %>% 
      dplyr::rename("parent.name" = indv.name, 
                    "parent.birth" = birth.year, 
                    #"parent.age" = age.x, 
                    "parent.sex"= sex, 
                    "parent.capture.year" = capture.year) %>%
      dplyr::select(parent.name, 
                    offspring.name, 
                    offspring.birth, 
                    offspring.mom, 
                    offspring.dad, 
                    parent.sex,
                    parent.capture.year, 
                    parent.age.in.offspring.birth.year = age.in.OffBirth.year)
    
    head(pairwise.df_all.info.PO)
    
    #Ok. We have a pairwise comparison matrix with all the info we need to assign kinship.
    
    #Don't think I need the commented code below, since I can just compare the assigned parent of the offspring with the potential parent in the comparison
    
    #Make dataframe of positive matches
    positives.PO <- pairwise.df_all.info.PO %>% 
      dplyr::filter(parent.name == offspring.mom | parent.name == offspring.dad) %>% 
      mutate(parent = ifelse(parent.name == offspring.mom, "mother", "father"))
    
    nrow(positives.PO)
    
    #Separate into sex-specific dataframes
    #Assumes we've previously filtered for all potential parents i.e. there should not be any potential parents here that were younger than reproductive age during the offspring birth year. Check here; should all be repro.age or older.
    pairwise.df_all.info.PO %>% arrange(parent.age.in.offspring.birth.year) %>% 
      pull(parent.age.in.offspring.birth.year) %>% 
      head()
    
    
    #Assuming we've already filtered for age, now we are concerned with counting the number of comparisons for each combination of offspring birth and parent capture year
    mom_positives.PO <- positives.PO %>% dplyr::filter(parent == "mother") %>% 
      dplyr::select(offspring.birth, 
                    parent.capture.year) %>% 
      plyr::count()
    
    dad_positives.PO <- positives.PO %>% dplyr::filter(parent == "father") %>% 
      dplyr::select(offspring.birth, 
                    parent.capture.year) %>% 
      plyr::count()
    
    #Confirmed that building the positive and negative dataframes separately gives the same result as the process below.
    
    #Create dataframe of all, negative, and positive comparisons for each combination of parent capture year and offspring birth year
    mom_comps.PO <- pairwise.df_all.info.PO %>%
      dplyr::filter(parent.sex == "F") %>% 
      dplyr::select(offspring.birth, parent.capture.year) %>% 
      plyr::count() %>% 
      rename(all = freq) %>% 
      left_join(mom_positives.PO, by = c("offspring.birth", "parent.capture.year")) %>% 
      rename(yes = freq) %>%
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(no = all - yes) %>% 
      mutate(no = replace_na(no, 0)) %>% 
      mutate(mort.yrs = ifelse(parent.capture.year < offspring.birth, offspring.birth - parent.capture.year, 0)) #Need a separate equation w/ survival for instances where the parent was sampled before the offspring. So flag here. Later, can sort based on whether there's a 0 here or not.
    
    dad_comps.PO <- pairwise.df_all.info.PO %>%
      dplyr::filter(parent.sex == "M") %>% 
      dplyr::select(offspring.birth, parent.capture.year) %>% 
      plyr::count() %>% 
      rename(all = freq) %>% 
      left_join(dad_positives.PO, by = c("offspring.birth", "parent.capture.year")) %>% 
      rename(yes = freq) %>%
      mutate(yes = replace_na(yes, 0)) %>% 
      mutate(no = all - yes) %>% 
      mutate(no = replace_na(no, 0)) %>% 
      mutate(mort.yrs = ifelse(parent.capture.year < offspring.birth, offspring.birth - parent.capture.year, 0)) #Need a separate equation w/ survival for instances where the parent was sampled before the offspring. So flag here. Later, can sort based on whether there's a 0 here or not.
    
    #Add the dataframe for this year to the list to be output
    PO.mom.pairwise.list[[y]] <- mom_comps.PO
    #names(PO.mom.pairwise.list)[[y]] <- paste0("PO.mom.pairwise.year_", PO.year)
    
    PO.dad.pairwise.list[[y]] <- dad_comps.PO
    #names(PO.dad.pairwise.list)[[y]] <- paste0("PO.dad.pairwise.year_", est.year)
  }
  
  
  
  #--------------Half sibling--------------------------
  filtered.samples.HS.df <- filtered.samples.HS.df %>% dplyr::arrange(birth.year) #Arrange so older sib always comes first in pairwise comparison matrix
  
  #Double-check that numbers are the same
  pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) %>%
    as_tibble() %>% 
    dplyr::filter(X1 != X2) %>% 
    distinct(X1, X2)
  
  colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  pwdf.test.HS <- pairwise.df.HS %>% dplyr::rename(indv.name = older.sib, older.sib = indv.name) %>% 
    bind_rows(pairwise.df.HS)
  
  pwdf.test.HS %>% distinct(indv.name, older.sib) %>% nrow()/2
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(pairwise.df.HS)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
    dplyr::select(older.sib = indv.name, 
                  older.sib.birth = birth.year, 
                  older.sib.age = age.x, 
                  older.sib.mom = mother.x,
                  older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(filtered.samples.HS.df, by = "indv.name") %>% 
    dplyr::rename("younger.sib" = indv.name, 
                  "younger.sib.birth" = birth.year, 
                  "younger.sib.age" = age.x, 
                  "younger.sib.mom" = mother.x, 
                  "younger.sib.dad" = father.x) %>%
    dplyr::select(older.sib, 
                  older.sib.birth, 
                  older.sib.age, 
                  younger.sib, 
                  younger.sib.birth, 
                  younger.sib.age, 
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
  #nrow(positives)
  
  ####----------------Split dataframes into final form for model----------####
  #Sex-specific half-sib
  mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    rename(yes = freq)
  
  dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
    dplyr::select(older.sib.birth, younger.sib.birth) %>%
    plyr::count() %>% 
    rename(yes = freq)
  
  
  #Make dataframes for negative comparisons
  #Confirmed the final total number of pairwise comparisons for moms and dads equals the number of rows in the pairwise comparison matrix
  mom.hsp.negs <- pairwise.df_HS.filt %>%
    dplyr::filter(younger.sib.mom != older.sib.mom) %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    as_tibble()
  
  mom_comps.HS <- mom.hsp.negs %>% 
    dplyr::rename(no = freq) %>% 
    left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no)
  
  dad.hsp.negs <- pairwise.df_HS.filt %>%
    dplyr::filter(younger.sib.dad != older.sib.dad) %>% 
    dplyr::select(older.sib.birth, younger.sib.birth) %>% 
    plyr::count() %>% 
    as_tibble()
  
  dad_comps.HS <- dad.hsp.negs %>% 
    dplyr::rename(no = freq) %>% 
    left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
    mutate(yes = replace_na(yes, 0)) %>% 
    mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(all = yes + no)
  
  
  # HS.mom.pairwise.list[[y]] <- mom_comps.HS
  # names(HS.mom.pairwise.list)[[y]] <- paste0("HS.mom.pairwise.year_", est.year)
  # 
  # HS.dad.pairwise.list[[y]] <- dad_comps.HS
  # names(HS.dad.pairwise.list)[[y]] <- paste0("HS.dad.pairwise.year_", est.year)
  
  
  #--------------Rename columns and combine PO and HS dataframes----------
  PO.mom.pairwise.df <- PO.mom.pairwise.list %>% 
    bind_rows() %>% 
    dplyr::rename(ref.year = offspring.birth) %>% 
    dplyr::select(ref.year, all, yes, mort.yrs) %>% 
    mutate(type = "PO", 
           parent = "mother")
  
  PO.dad.pairwise.df <- PO.dad.pairwise.list %>% 
    bind_rows() %>% 
    dplyr::rename(ref.year = offspring.birth) %>% 
    dplyr::select(ref.year, all, yes, mort.yrs) %>% 
    mutate(type = "PO",
           parent = "father")
  
  HS.mom.pairwise.df <- mom_comps.HS %>%
    dplyr::rename(ref.year = younger.sib.birth,
                  mort.yrs = year_gap) %>% 
    dplyr::select(ref.year, all, yes, mort.yrs) %>% 
    #dplyr::filter(mort.yrs < (2*repro.age)) %>% #Added 04/01/2022 to remove comparisons that could be confused for grand-parents
    mutate(type = "HS",
           parent = "mother")
  
  HS.dad.pairwise.df <- dad_comps.HS %>%
    dplyr::rename(ref.year = younger.sib.birth,
                  mort.yrs = year_gap) %>% 
    dplyr::select(ref.year, all, yes, mort.yrs) %>% 
    #dplyr::filter(mort.yrs < (2*repro.age)) %>% #Added 04/01/2022 to remove comparisons that could be confused for grand-parents
    mutate(type = "HS",
           parent = "father")
  
  #Combine PO and HS dataframes for each sex
  mom_comps.all <- bind_rows(HS.mom.pairwise.df, PO.mom.pairwise.df) %>% 
    mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
    arrange(desc(ref.year), mort.yrs) %>% 
    mutate(BI = ifelse(mort.yrs %% mating.periodicity == 0, "on", "off"))
  
  dad_comps.all <- bind_rows(HS.dad.pairwise.df, PO.dad.pairwise.df) %>% 
    mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
    arrange(desc(ref.year), mort.yrs)
  
  return(list(mom_comps.all, dad_comps.all, positives.HS))
}

samples %>% group_by(mother.x, father.x) %>%
  filter(n() > 1) %>% 
  arrange(mother.x)
  group_by(birth.year)
  summarize(n())
  distinct(birth.year) %>%
  ungroup() %>% 
  group_by(mother.x, father.x) %>% 
  summarize(distinct.birth.years = n()) %>% 
  filter(distinct.birth.years > 1) %>% 
  nrow()

#The vector of years for which we need to split samples into potential parents and offspring i.e. offspring birth years.
NoFullSibs.df %>% 
  dplyr::filter(age.x <= repro.age) %>% 
  distinct(birth.year) %>% 
  arrange(birth.year) %>%
  pull() 








sample.df_all.info %>% group_by(sample.prop, sampling.scheme, ref.yr) %>% 
  summarize(n())


samples.df.YOY <- samples.df
samples.df.juvs <- samples.df
samples.df.ALL <- samples.df
str(samples.df.YOY)

samples.df.YOY %>% group_by(iteration, sampling.scheme, sample.prop) %>% 
  summarize(samp.size = n()) %>% 
  group_by(sample.prop) %>% 
  summarize(mean_samp.size = mean(samp.size))

samples.df.juvs %>% group_by(iteration, sampling.scheme, sample.prop) %>% 
  summarize(samp.size = n()) %>% 
  group_by(sample.prop) %>% 
  summarize(mean_samp.size = mean(samp.size))

samples.df.ALL %>% group_by(iteration, sampling.scheme, sample.prop) %>% 
  summarize(samp.size = n()) %>% 
  group_by(sample.prop) %>% 
  summarize(mean_samp.size = mean(samp.size))





samples.df <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

samples.df %>% group_by(age.x) %>% summarize(n())

pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

(truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme)))



(a*(survival^2))/(((a + psi - (a*psi))*(Nf*(lambda^3))))



results.1 <- results.all %>% dplyr::filter(purpose == "HS.PO_downsample")

results.1
pop.size.1
rents.1.1 <- rents.1 %>% dplyr::distinct(iteration, year, parent, .keep_all = TRUE)
samples1

#Matches pop.size.1
rents.iter.yr <- rents.1.1 %>% group_by(iteration, year, parent.sex) %>% 
  summarize(n())

#Get unique number of parents sampled for each iteration
results.1$unique_parents_in_sample


#-----------Make pairwise comparison for each set of samples---------------#
init.adult.pop.size <- 1000 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M - has to add up to 1
YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival <- 0.825 # CHANGED FROM 0.9; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 #set the maximum age allowed in the simulation
mating.periodicity <- 1 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) #vector of potential number of mates per mating
#avg.num.offspring <- 3 # NOT USED? CHANGED FROM 3; set the average number of offspring per mating (from a poisson distribution)

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff

#Stable age distribution
props <- rep(NA, max.age+1)
props[1] <- f
props[2] <- f * YOY.survival
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival
#props[repro.age+1] <- props[repro.age] * juvenile.survival + Adult.survival

for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * Adult.survival
prop.Adult <- sum(props[(repro.age+1):(max.age+1)])/sum(props)
Nages <- round(props[-1] * init.adult.pop.size) 
init.pop.size <- sum(Nages) # all ages except YOYs

#Set length of simulation and estimation year
burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs <- burn.in + Num.years #Total number of simulation years
estimation.year <- n_yrs - 5 # Set year of estimation
source("./01_MAIN_scripts/functions/remove_dups.R")
source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_downsample.R")

mom_comps.samples1 <- dad_comps.samples1 <- NULL

options(dplyr.summarise.inform = FALSE)

for(iter in 1:100){
  print(paste0("working on iteration ", iter))
  
  #50 samples per year
  samples1.temp1 <- samples1 %>% filter(sample.size == 50 & iteration == iter)
  NoDups1.temp1 <- samples1.temp1 %>% split.dups()
  save.last1 <- NoDups1.temp1[[2]]
  filter1.temp1 <- save.last1 %>% filter.samples()
  PO.samples1.list1 <- filter1.temp1[[1]]
  HS.samples1.df1 <- filter1.temp1[[2]]
  pairwise.samples1 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list1, filtered.samples.HS.df = HS.samples1.df1)
  
  mom_comps.temp1 <- pairwise.samples1[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 50)
  dad_comps.temp1 <- pairwise.samples1[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 50)
  
  #150 samples per year
  samples1.temp2 <- samples1 %>% filter(sample.size == 150 & iteration == iter)
  NoDups1.temp2 <- samples1.temp2 %>% split.dups()
  save.last2 <- NoDups1.temp2[[2]]
  filter1.temp2 <- save.last2 %>% filter.samples()
  PO.samples1.list2 <- filter1.temp2[[1]]
  HS.samples1.df2 <- filter1.temp2[[2]]
  pairwise.samples2 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list2, filtered.samples.HS.df = HS.samples1.df2)
  
  mom_comps.temp2 <- pairwise.samples2[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 150)
  dad_comps.temp2 <- pairwise.samples2[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 150)
  
  #250 samples per year
  samples1.temp3 <- samples1 %>% filter(sample.size == 250 & iteration == iter)
  NoDups1.temp3 <- samples1.temp3 %>% split.dups()
  save.last3 <- NoDups1.temp3[[2]]
  filter1.temp3 <- save.last3 %>% filter.samples()
  PO.samples1.list3 <- filter1.temp3[[1]]
  HS.samples1.df3 <- filter1.temp3[[2]]
  pairwise.samples3 <- build.pairwise(filtered.samples.PO.list = PO.samples1.list3, filtered.samples.HS.df = HS.samples1.df3)
  
  mom_comps.temp3 <- pairwise.samples3[[1]] %>% 
    mutate(iteration = iter, samples.per.yr = 250)
  dad_comps.temp3 <- pairwise.samples3[[2]] %>% 
    mutate(iteration = iter, samples.per.yr = 250)
  
  
  mom_comps.samples1 <- rbind(mom_comps.samples1, mom_comps.temp1, mom_comps.temp2, mom_comps.temp3) %>% 
    as_tibble()
  dad_comps.samples1 <- rbind(dad_comps.samples1, dad_comps.temp1, dad_comps.temp2, dad_comps.temp3) %>% 
    as_tibble()
  
  
}

head(mom_comps.samples1)
head(dad_comps.samples1)

nrow(mom_comps.samples1)

#Split pairwise comparisons by number of samples so we can plot
mom_comps1.50samps <- mom_comps.samples1 %>% dplyr::filter(samples.per.yr == 50)
mom_comps1.150samps <- mom_comps.samples1 %>% dplyr::filter(samples.per.yr == 150)
mom_comps1.250samps <- mom_comps.samples1 %>% dplyr::filter(samples.per.yr == 250)

dad_comps1.50samps <- dad_comps.samples1 %>% dplyr::filter(samples.per.yr == 50)
dad_comps1.150samps <- dad_comps.samples1 %>% dplyr::filter(samples.per.yr == 150)
dad_comps1.250samps <- dad_comps.samples1 %>% dplyr::filter(samples.per.yr == 250)


rents.1
results.1 %>% arrange(desc(relative_bias)) %>% View()


#-----------Calculate C and expected HSPs-------------------
#Subset by iteration and model type
mom_comps1.50samps.iter1 <- mom_comps1.50samps %>% dplyr::filter(iteration == 1) %>% 
  dplyr::filter(type == "HS")

E.comps <- samples^2


#Summarize umber of distinct mortality years
age.gap.denom <- mom_comps1.50samps.iter1 %>% distinct(mort.yrs) %>% 
  summarize(s = sum(mort.yrs)) %>% 
  pull(s)

#What's the true abundance?
N.truth <- pop.size.1 %>% filter(year == 85 & iter == 1) %>% 
  distinct(Female.adult.pop) %>% 
  pull(Female.adult.pop)

#Number of comparisons for each year gap (mort.yrs)
comps.yr.df <- mom_comps1.50samps.iter1 %>% group_by(mort.yrs) %>% 
  summarize(comps.yr = sum(all))

#Calculate expected HSPs 
mom_comps1.50sampsC <- mom_comps1.50samps.iter1 %>% 
  group_by(mort.yrs) %>% 
  summarize(cross.cohort.instances = n()) %>% 
  mutate(cum.surv = Adult.survival ^ mort.yrs) %>% 
  mutate(Ct_y = cross.cohort.instances * cum.surv) %>% 
  mutate(C_yr = Ct_y/age.gap.denom)



weighted.mean.C <- sum(mom_comps1.50sampsC$C_yr)
total.comps <- sum(mom_comps1.50samps.iter1$all)



(Exp_HS <- (weighted.mean.C * total.comps)/N.truth)
  
Exp_HS = Y*(C/Nf)

#NEXT: calculate number of comparisons.
#sum across age gap
#Calculate weighted mean C
#Want to get value of expected number of HSPs