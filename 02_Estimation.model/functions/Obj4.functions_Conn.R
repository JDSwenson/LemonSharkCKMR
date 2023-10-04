####--------------------Split individual recaptures into first and last occasion------------------####
split.dups <- function(samples){
  
  #Identify duplicates
  dups.first.capture <- samples %>% group_by(indv.name) %>% 
    filter(n() > 1) %>% 
    distinct(indv.name, .keep_all = TRUE) %>% 
    ungroup() 
  
  dups.later.capture <- samples %>% group_by(indv.name) %>% 
    filter(n() > 1) %>% 
    anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>% 
    ungroup()
  
  #Save the first capture instance of recaptured individuals and sort so older individual comes first (for pairwise comparison matrix)
  samples.save.first <- samples %>% anti_join(dups.later.capture, by = c("indv.name", "capture.year")) %>% 
    dplyr::arrange(birth.year) %>% 
    as_tibble()
  
  #Save the last capture instance of recaptured individuals and sort so older individual comes first (for pairwise comparison matrix)
  samples.save.second <- samples %>% anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>%
    dplyr::arrange(birth.year) %>% 
    as_tibble()
  
  return(list(samples.save.first, samples.save.second))  
}

####---------------Filter samples -------------####
#Filters for full sibs, but also identifies appropriate sample sets for each offspring birth year (for PO) and each younger sib birth year (for half-sib)
filter.samples <- function(samples){

  #------------------Filter 1: full siblings----------------------
    NoFullSibs.df <- samples %>% distinct(mother.x, father.x, .keep_all = TRUE) %>%
    as_tibble() #If there is more than one individual with the same mother AND father, then only keep one.
    total.sibs <- nrow(samples) - nrow(NoFullSibs.df)
      
    dif.cohort.sibs <- samples %>% group_by(mother.x, father.x) %>%
      filter(n() > 1) %>% 
      distinct(birth.year.miss) %>%
      ungroup() %>% 
      group_by(mother.x, father.x) %>% 
      summarize(distinct.birth.years = n()) %>% 
      filter(distinct.birth.years > 1) %>% 
      nrow()

    #The vector of years for which we need to split samples into potential parents and offspring i.e. offspring birth years.
  OffBirth.years <- NoFullSibs.df %>% 
    dplyr::filter(age.miss == 0) %>% 
    distinct(birth.year.miss) %>% 
    arrange(birth.year.miss) %>%
    pull() 
  
  #Initialize lists to save output
  PO.samps.list <- list()
  HS.samps.list <- list()
  parents <- list()
  offspring <- list()
  

  #---------Filter 2: potential parents and offspring for each year------
  #For each offspring birth year, create a dataframe of potential parents i.e. sampled individuals that are old enough to have been a parent in that year.
  for(y in 1:length(OffBirth.years)){
    OffBirth.year <- OffBirth.years[y] #year of focus for this iteration
    
    #Dataframe of reproductively mature individuals during offspring birth year
    parents[[y]] <- NoFullSibs.df %>%
      mutate(age.in.OffBirth.year = age.miss + (OffBirth.year - capture.year)) %>% 
      filter(age.in.OffBirth.year >= repro.age) %>% 
      mutate(relation = "parent")
    
    #Dataframe of offspring born in this year
    offspring[[y]] <- NoFullSibs.df %>% 
      mutate(age.in.OffBirth.year = age.miss + (OffBirth.year - capture.year)) %>% #unnecessary line, I think
      filter(birth.year.miss == OffBirth.year) %>%
      mutate(relation = "offspring")
      
    #Add dataframe of all potential offspring and parents to a list where each dataframe corresponds to an offspring birth year
    PO.samps.list[[y]] <- rbind(parents[[y]], offspring[[y]]) %>% 
      dplyr::select(indv.name, 
                    birth.year.miss, 
                    capture.year,
                    age.in.OffBirth.year,
                    mother.x, 
                    father.x,
                    sex,
                    relation)
    
    #Not sure if this is needed, but add dataframe of potential half-sib samples for the year to a list where each dataframe corresponds to a specific birth year
    # HS.samps.list[[y]] <- NoFullSibs.df %>%
    #   filter(birth.year.miss <= OffBirth.year) #Filters for offspring that were born before or in Offbirth.year.miss
  }
  
  #Name list elements to correspond to sample year
  names(PO.samps.list) <- paste0("PO.samples.year_", c(OffBirth.years))
  #names(HS.samps.list) <- paste0("HS.samples.year_", c(OffBirth.years))
  
  #Not sure if this is needed, but create dataframes of samples for PO and HS comparisons. Think only need this for half-sibs
#  PO.samps.df <- bind_rows(PO.samps.list) #Confirmed that it's the correct number of rows
  HS.samps.df <- NoFullSibs.df %>% 
    dplyr::filter(age.miss < repro.age) #Only use juveniles for half-sib model
    
  print(paste0("There are ", total.sibs, " pairs of full siblings in the dataset, ", dif.cohort.sibs, " of these pairs were born in different years."))
  return(list(PO.samps.list, HS.samps.df, NoFullSibs.df, total.sibs, dif.cohort.sibs)) #Return list of possible parents and offspring for each year, and dataframe of potential half-sibs
  
}



####--------------Build pairwise comparison matrices------------####
#Input is the filtered list from above for POs, and the filtered dataframe for HS
build.pairwise <- function(filtered.samples.PO.list, filtered.samples.HS.df){

#Different birth years for abundance estimate
  # OffBirth.years.PO.vec <- filtered.samples.PO.list %>% 
  #   bind_rows() %>%
  #   distinct(birth.year.miss) %>% 
  #   pull(birth.year.miss) 

  
  #initialize lists for output
  PO.mom.pairwise.list <- list()
  HS.mom.pairwise.list <- list()
  PO.dad.pairwise.list <- list()
  HS.dad.pairwise.list <- list()
  
  
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
    pairwise.df.PO <- expand.grid(PO.parents, PO.offspring) 
    colnames(pairwise.df.PO) <- c("indv.name", "offspring.name") #Rename columns so they can easily be joined
  
head(pairwise.df.PO)

#Create dataframe that will be used to extract the birth years for the potential offspring from each comparison using joins.
Offspring_birthyears.PO <- filtered.samples.PO.list[[y]] %>%
  dplyr::select(offspring.name = indv.name, 
                offspring.birth = birth.year.miss, 
                offspring.mom = mother.x, 
                offspring.dad = father.x)

#Extract metadata for each individual in each comparison
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.PO <- pairwise.df.PO %>% left_join(Offspring_birthyears.PO, by = "offspring.name") %>% 
  as_tibble() %>% 
  left_join(filtered.samples.PO.list[[y]], by = "indv.name") %>% 
  dplyr::rename("parent.name" = indv.name, 
                "parent.birth" = birth.year.miss, 
                #"parent.age" = age.miss, 
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
filtered.samples.HS.df <- filtered.samples.HS.df %>% dplyr::arrange(birth.year.miss) #Arrange so older sib always comes first in pairwise comparison matrix
pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
  dplyr::select(older.sib = indv.name, 
         older.sib.birth = birth.year.miss, 
         older.sib.age = age.miss, 
         older.sib.mom = mother.x,
         older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(filtered.samples.HS.df, by = "indv.name") %>% 
  dplyr::rename("younger.sib" = indv.name, 
                "younger.sib.birth" = birth.year.miss, 
                "younger.sib.age" = age.miss, 
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
hsp.negs <- pairwise.df_HS.filt %>%
  dplyr::filter(younger.sib.mom != older.sib.mom & younger.sib.dad != older.sib.dad) %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  as_tibble()

mom_comps.HS <- hsp.negs %>% 
  dplyr::rename(no = freq) %>% 
  left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>% 
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(year_gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(all = yes + no)

dad_comps.HS <- hsp.negs %>% 
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
  mutate(pop.growth.yrs = ref.year - est.year.calibrate) %>% 
  arrange(desc(ref.year), mort.yrs) %>% 
  mutate(BI = ifelse(mort.yrs %% mating.periodicity == 0, "on", "off"))

dad_comps.all <- bind_rows(HS.dad.pairwise.df, PO.dad.pairwise.df) %>% 
  mutate(pop.growth.yrs = ref.year - est.year.calibrate) %>% 
  arrange(desc(ref.year), mort.yrs)

return(list(mom_comps.all, dad_comps.all, positives.HS))
}


#----------------Calculate number of HSPs for downsampling---------------
downsample <- function(mom_comps.all, dad_comps.all, HS.samps.df, PO.samps.list){

  #Downsample for HSPs
#Calculate proportion of samples born in each birth year
  HS.samps.props <- HS.samps.df %>% group_by(birth.year.miss) %>%
    summarize(num = n()) %>%
  ungroup() %>%
  mutate(prop = num/sum(num))

  #Save vector of proportions corresponding to the correct index in the HS sample dataframe
HS.props.vec <- HS.samps.df %>% left_join(HS.samps.props, by = "birth.year.miss") %>%
  pull(prop)

  #Calculate proportion of positive HS comparisons for mom
HS_comps.mom.yes <- mom_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

HS_comps.mom.all <- mom_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

HS_prop.mom.yes <- HS_comps.mom.yes/HS_comps.mom.all #Calculate proportion of positive comparisons

#Calculate proportion of positive HS comparisons for dad
HS_comps.dad.yes <- dad_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

HS_comps.dad.all <- dad_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

HS_prop.dad.yes <- HS_comps.dad.yes/HS_comps.dad.all

HS_prop.mean.yes <- mean(c(HS_prop.mom.yes, HS_prop.dad.yes)) 
HS_target.comps = round(max.HSPs/HS_prop.mean.yes, 0)
HS_target.samples <- round(sqrt(HS_target.comps), 0)

HS.samps.df.down <- NULL #Initialize so I don't get an error
if(HS_comps.mom.yes + HS_comps.dad.yes > max.HSPs){
  HS.samps.df.down <- HS.samps.df %>% slice_sample(n = HS_target.samples, weight_by = HS.props.vec)
} else{
  HS.samps.df.down <- HS.samps.df
}


#-----------------------Downsample for POPs-------------------------------------#
PO_comps.mom.yes <- mom_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

PO_comps.mom.all <- mom_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

PO_prop.mom.yes <- PO_comps.mom.yes/PO_comps.mom.all #Calculate proportion of positive comparisons

#Calculate proportion of positive PO comparisons for dad
PO_comps.dad.yes <- dad_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

PO_comps.dad.all <- dad_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

PO_prop.dad.yes <- PO_comps.dad.yes/PO_comps.dad.all

PO.all.yes <- PO_comps.dad.yes + PO_comps.mom.yes

#This doesn't work for POPs the same way it does for half sibs. So, just streamline it and cut number of samples down proportionally
# PO_prop.mean.yes <- mean(c(PO_prop.mom.yes, PO_prop.dad.yes)) 
# PO_target.comps = round(max.POPs/PO_prop.mean.yes, 0)
# PO_target.samples <- round(sqrt(PO_target.comps), 0)
# PO_target.samples.per.yr <- round(PO_target.samples/length(sample.years), 0)

#sample.size.rents.reduced <- (max.POPs/PO.all.yes)*sample.size.rents
#sample.size.juvs.reduced <- (max.POPs/PO.all.yes)*sample.size.juvs
#PO.down.prop <- max.POPs/PO.all.yes

PO.samps.list.down <- list()
if(PO_comps.mom.yes + PO_comps.dad.yes > max.POPs){

  #---------Filter 2: potential parents and offspring for each year------
#For each offspring birth year, create a dataframe of potential parents i.e. sampled individuals that are old enough to have been a parent in that year.
for(y in 1:length(PO.samps.list)){
  
  off.down <- PO.samps.list[[y]] %>% dplyr::filter(relation == "offspring") %>% 
    slice_sample(prop = PO.down.prop)
  
  rents.down <- PO.samps.list[[y]] %>% dplyr::filter(relation == "parent") %>% 
    slice_sample(prop = PO.down.prop)
  
  PO.samps.list.down[[y]] <- rbind(off.down, rents.down)
  
  } #End creation of PO sample list

names(PO.samps.list.down) <- paste0("PO.samples.year_", c(sample.years))

} else {PO.samps.list.down <- PO.samps.list} #End PO downsample

return(list(HS.samps.df.down,
            PO.samps.list.down))
}



calc.psi <- function(loopy.list){
  #Calculate breeding interval
  #Initialize sample dataframes
  BI.df <- NULL
  BI.df_temp <- NULL
  
  ref.years <- c(50:n_yrs)
  
  for(i in ref.years){ #Extract all YOY for each year sampled
    BI.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
      dplyr::filter(age.miss == 0)
    
    #Combine all
    BI.df <- bind_rows(BI.df, BI.df_temp)
  }
  
  #Extract unique mothers for each birth year
  BI.df <- BI.df %>% as_tibble() %>% 
    arrange(birth.year.miss) %>% 
    distinct(mother.x, birth.year.miss, .keep_all = TRUE)
  
  #-----------------Identify all relatives in population since estimation year and find breeding interval--------------------#
  pairwise.BI.df <- data.frame(t(combn(BI.df$indv.name, m=2)))  # generate all combinations of the elements of x, taken m at a time.
  colnames(pairwise.BI.df) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(BI.df)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.BI <- BI.df %>%
    dplyr::select(older.sib = indv.name, 
           older.sib.birth = birth.year.miss, 
           older.sib.age = age.miss, 
           older.sib.mom = mother.x,
           older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.BI.all.info <- pairwise.BI.df %>% left_join(OlderSib_birthyears.BI, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(BI.df, by = "indv.name") %>% 
    dplyr::rename("younger.sib" = indv.name, 
                  "younger.sib.birth" = birth.year.miss, 
                  "younger.sib.age" = age.miss, 
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
                  younger.sib.dad) %>% 
    dplyr::filter(older.sib.birth != younger.sib.birth)
  
  
  #Extract all positive maternal half-sib comparisons from sampled years
  positives.mom.BI <- pairwise.BI.all.info %>% filter(older.sib.mom == younger.sib.mom) %>% 
    mutate(yr.gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(BI = ifelse(yr.gap %% mating.periodicity == 0, "on", "off")) %>% #If the year gap between ha;f-sibs is divisible by 2, then call the breeding interval (BI) "even"
    dplyr::distinct(older.sib.birth, younger.sib.birth, younger.sib.mom, .keep_all = TRUE) #Duplicative of earlier; making sure cohort size isn't biasing things here
  
  #Everything matches up re: positives and 
  positives.mom.BI %>% group_by(BI) %>% summarize(n()) #How many positive comparisons for 
  
  skipped.moms <- positives.mom.BI %>% dplyr::filter(BI == "on") %>% 
    dplyr::distinct(younger.sib.mom) %>% 
    pull(younger.sib.mom) #Should be the same as the older sib mom
  
  annual.moms <- positives.mom.BI %>% dplyr::filter(BI == "off") %>% 
    dplyr::distinct(younger.sib.mom) %>% 
    pull(younger.sib.mom) #Should be the same as the older sib mom
  
  #Which moms ONLY reproduced in even years?
  skipped.only.moms <- skipped.moms[which(!skipped.moms %in% annual.moms)]
  all.moms <- unique(positives.mom.BI$younger.sib.mom)
  
  #Percent of individuals that only breed when years are evenly spaced
  (psi.truth <- length(skipped.only.moms)/length(all.moms))
  
}







#----------------Calculate expected matches-------------------
calc.Exp <- function(mom_comps.all, dad_comps.all){
  #Make dataframe of truth, for calculations of expectations
  #ADDED from main script on 10/04/2023
  truth.iter <- pop_size.df %>% dplyr::filter(iteration == iter,
                                              year == estimation.year) %>% 
    mutate(Nfb1 = Num.mothers,
           Nmb1 = Num.fathers) %>%
    dplyr::select(Nf = Female.adult.pop,
                  Nfb1,
                  Nfb2 = Num.mothers,
                  Nm = Male.adult.pop,
                  Nmb1,
                  Nmb2 = Num.fathers,
                  estimation.year = year,
                  iteration = iteration,
                  seed = seed) %>% 
    pivot_longer(cols = starts_with("N"),
                 names_to = "parameter",
                 values_to = "all.truth") %>% 
    #bind_rows(lambda.truth.df) %>% 
    mutate(sampling.scheme = s.scheme, 
           sample.prop = sample.proportion,
           T0 = est.year.calibrate)
  
  #Calculate expectations
  #Define Nf and Nm as objects for calculations
  Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
  
  Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
  
  #-----END ADDED
  
#-----------Calculate C and expected HSPs-------------------
#Calculate Exp(R) for MHSPs
mom_C.HS <- mom_comps.all %>% dplyr::filter(type == "HS") #Filter for HS relationships

#Summarize number of distinct mortality years
mom_age.gap.denom <- mom_C.HS %>% distinct(mort.yrs) %>% 
  summarize(s = sum(mort.yrs)) %>% 
  pull(s)

#What's the true abundance?
# mom_N.truth <- pop.size.tibble %>% filter(year == estimation.year) %>% 
#   distinct(Num.mothers) %>% 
#   pull(Num.mothers)

#Number of comparisons for each year gap (mort.yrs)
mom_comps.yr.df <- mom_C.HS %>% group_by(mort.yrs) %>% 
  summarize(comps.yr = sum(all))

#Calculate expected MHSPs 
mom_C.HS.Exp <- mom_C.HS %>% 
  group_by(mort.yrs) %>% 
  summarize(cross.cohort.instances = n()) %>% 
  mutate(cum.surv = adult.survival ^ mort.yrs) %>% 
  mutate(Ct_y = cross.cohort.instances * cum.surv) %>% 
  mutate(C_yr = Ct_y/mom_age.gap.denom) %>% 
  left_join(mom_comps.yr.df, by = "mort.yrs") %>% 
  mutate(Exp.R = (cum.surv*comps.yr)/Nf.truth)

mom_weighted.mean.C <- sum(mom_C.HS.Exp$C_yr)
mom_total.comps <- sum(mom_comps.yr.df$comps.yr)

#(mom.Exp.HS <- round((mom_weighted.mean.C * mom_total.comps)/mom_N.truth, 0))
mom.Exp.HS <- round(sum(mom_C.HS.Exp$Exp.R), 0)


#Calculate Exp(R) for MPOPs
mom_C.PO <- mom_comps.all %>% dplyr::filter(type == "PO")

#Need to use Na to get the correct number of expected POs 
# Nf.truth <- pop.size.tibble %>% filter(year == 85) %>% 
#   distinct(Female.adult.pop) %>% 
#   pull(Female.adult.pop)

mom.Exp.PO <- round(sum(mom_C.PO$all)/Nf.truth, 0)

#---------------------Dad--------------------------------
#Calculate Exp(R) for MHSPs
dad_C.HS <- dad_comps.all %>% dplyr::filter(type == "HS")

#Summarize umber of distinct mortality years
dad_age.gap.denom <- dad_C.HS %>% distinct(mort.yrs) %>% 
  summarize(s = sum(mort.yrs)) %>% 
  pull(s)

#What's the true abundance?
# dad_N.truth <- pop.size.tibble %>% filter(year == 85) %>% 
#   distinct(Male.adult.pop) %>% 
#   pull(Male.adult.pop)

#Number of comparisons for each year gap (mort.yrs)
dad_comps.yr.df <- dad_C.HS %>% group_by(mort.yrs) %>% 
  summarize(comps.yr = sum(all))

#Calculate expected FHSPs 
dad_C.HS.Exp <- dad_C.HS %>% 
  group_by(mort.yrs) %>% 
  summarize(cross.cohort.instances = n()) %>% 
  mutate(cum.surv = adult.survival ^ mort.yrs) %>% 
  mutate(Ct_y = cross.cohort.instances * cum.surv) %>% 
  mutate(C_yr = Ct_y/dad_age.gap.denom) %>% 
  left_join(dad_comps.yr.df, by = "mort.yrs") %>% 
  mutate(Exp.R = (cum.surv*comps.yr)/Nm.truth)

dad_weighted.mean.C <- sum(dad_C.HS.Exp$C_yr)
dad_total.comps <- sum(dad_comps.yr.df$comps.yr)

#(dad.Exp.HS <- round((dad_weighted.mean.C * dad_total.comps)/dad_N.truth, 0))
dad.Exp.HS <- round(sum(dad_C.HS.Exp$Exp.R), 0) #Sum expected HAPs per year


#Calculate Exp(R) for MPOPs
dad_C.PO <- dad_comps.all %>% dplyr::filter(type == "PO")

dad.Exp.PO <- round(sum(dad_C.PO$all)/Nm.truth, 0)

return(list(mom.Exp.HS, mom.Exp.PO, dad.Exp.HS, dad.Exp.PO))
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


misassign.ages <- function(samples){

  #------------Use VonBertalanffy growth function to assign lengths to individuals from ages--------------
  #Set up growth function
  ages <- 0:max.age
  vonBert <- vbFuns(param = "Typical")
  tmp <- growthFunShow("vonBertalanffy","Typical")
  
  #Use parameters from Brown and Gruber 1988
  mean_length.at.age <- vonBert(t = ages, Linf = 317.65, K = 0.057, t0 = -2.302)
  mean_length.subset <- mean_length.at.age[4:length(mean_length.at.age)]
  
  sd.vec <- c(3, 10, 9, age.cv*mean_length.subset) #Try 10 so we don't misassign across too many ages
  
  length.at.age_df <- tibble(mean.length = mean_length.at.age,
                             sd.length = sd.vec) %>% #Name as mean.1 and mean.2 to allow for join later
    mutate(age.x = ages)

  #Save vonBert parameters
  Linf_val <- 317.65
  t0_val <- -2.302
  K_val <- 0.057
  
  #Assign lengths using VonBert function
  samples.int <- samples %>% left_join(length.at.age_df, by = "age.x") %>% 
    mutate(indv.length = rnorm(n = nrow(samples), mean = mean.length, sd = sd.length)) %>% 
    mutate(indv.length = ifelse(indv.length >= Linf_val, Linf_val - 1, indv.length)) %>% 
    dplyr::select(-c(mean.length, sd.length))

  #Specify reverse VonBert function that will assign ages based on length
  rev.VonBert <- function(t0, Lt, Linf, K){
    age <- round(t0 + (log(1 - Lt/Linf)/-K), 0)
    return(age)
  }

  #Misassign ages from function specified above
  samples.miss <- samples.int %>% mutate(age.miss = rev.VonBert(t0 = t0_val, Linf = Linf_val, K = K_val, Lt = indv.length)) %>% 
    mutate(age.miss = ifelse(age.miss < 0, 0, age.miss)) %>% 
    mutate(birth.year.miss = capture.year - age.miss)  
    
    return(samples.miss)
}

calc.truth <- function(truth.df = truth.df){
  #If there's no lambda in the model, then the true value of abundance is the mean over the estimated years. Here, we make this correction for the scenarios that do not include lambda in the model; for those that do, the values in truth.df are correct.
  if(scenario == "scenario_1_model.validation" | 
     scenario == "scenario_1.2.1" | 
     scenario == "scenario_1.2.2" | 
     scenario == "scenario_1.2.3"){
    
    #Make truth equal to mean over years
    #est.year.calibrate is the first year of data aka T0. This is different than the object estimation.year, which is what we loop over and change.
    Nf.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.female.pop = mean(Female.adult.pop)) %>% 
      pull(mean.female.pop)
    
    Nfb.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.mothers = mean(Num.mothers)) %>% 
      pull(mean.mothers)
    
    Nm.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.male.pop = mean(Male.adult.pop)) %>% 
      pull(mean.male.pop)
    
    Nmb.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.fathers = mean(Num.fathers)) %>% 
      pull(mean.fathers)
    
    truth.iter <- truth.iter %>% mutate(all.truth = ifelse(parameter == "Nf", Nf.truth, 
                                                           ifelse(parameter == "Nfb1" | parameter == "Nfb2", Nfb.truth, 
                                                                  ifelse(parameter == "Nm", Nm.truth, 
                                                                         ifelse(parameter == "Nmb1" | parameter == "Nmb2", Nmb.truth, all.truth)))))
    cat("\nUsing average abundance for truth\n")
  } else cat("\nUsing year-specific truth, rather than average.\n")
  
  
  ###ADDED 10/03/2023
  
  #Calculate truth for lambda
  lambda.pop.df <- pop_size.df %>% dplyr::filter(iteration == iter)
  
  lambda.truth <- (lambda.pop.df$Total.adult.pop[n_yrs]/lambda.pop.df$Total.adult.pop[est.year.calibrate])^(1/(n_yrs - est.year.calibrate))
  
  lambda.truth.df <- tibble(estimation.year = estimation.year,
                            iteration = iter,
                            seed = rseed,
                            parameter = "lambda",
                            all.truth = lambda.truth)
  
  
  #Make dataframe of truth
  truth.iter.temp <- pop.size.tibble %>% dplyr::filter(year == estimation.year) %>% 
    mutate(Nfbt = Num.mothers,
           Nft = Female.adult.pop,
           Nmt = Male.adult.pop) %>%
    dplyr::select(Nfbt,
                  Nft,
                  Nmt,
                  estimation.year = year,
                  iteration = iteration,
                  seed = seed) %>% 
    pivot_longer(cols = starts_with("N"),
                 names_to = "parameter",
                 values_to = "all.truth") 
  
  truth.iter <- pop.size.tibble %>% dplyr::filter(year == est.year.calibrate) %>% 
    mutate(Nfb0 = Num.mothers,
           Nf0 = Female.adult.pop,
           Nm0 = Male.adult.pop) %>%
    dplyr::select(Nfb0,
                  Nf0,
                  Nm0,
                  estimation.year = year,
                  iteration = iteration,
                  seed = seed) %>%
    pivot_longer(cols = starts_with("N"),
                 names_to = "parameter",
                 values_to = "all.truth") %>% 
    bind_rows(truth.iter.temp, lambda.truth.df) %>% 
    mutate(T0 = est.year.calibrate,
           sampling.scheme = s.scheme,
           sample.prop = sample.proportion)
  
  #----------END ADDED 10/3/2023
  # N0.truth <- pop_size.df %>% dplyr::filter(iteration == iter,
  #                                           year == est.year.calibrate) %>% 
  #   dplyr::select(Nm0 = Male.adult.pop, 
  #                 Nf0 = Female.adult.pop) %>% 
  #   pivot_longer(cols = starts_with("N"),
  #                names_to = "parameter",
  #                values_to = "all.truth") %>% 
  #   mutate(estimation.year = estimation.year,
  #          sampling.scheme = s.scheme,
  #          sample.prop = sample.proportion,
  #          T0 = est.year.calibrate,
  #          iteration = iter,
  #          seed = rseed)
  
  (truth.iter_all.params <- truth.df %>% dplyr::filter(sampling.scheme == s.scheme,
                                                       iteration == iter,
                                                       sample.prop == sample.proportion) %>% 
      mutate(estimation.year = estimation.year) %>% 
      dplyr::mutate(T0 = ref.yr + 1) %>% #Should be the same as est.year.calibrate
      dplyr::select(-c(surv_min, surv_max, population.growth)) %>% 
      bind_rows(truth.iter) %>% 
      dplyr::select(parameter, all.truth, estimation.year, T0, iteration, sampling.scheme, sample.prop, seed) %>% #Change order of columns
      dplyr::arrange(parameter))
  
  
  results.temp <- model.summary2 %>% left_join(truth.iter_all.params, by = c("parameter", "iteration", "seed")) %>% 
    mutate(estimation.sim = ifelse(est == 1, "T0", 
                                   ifelse(est == 2, "T0-10",
                                          ifelse(est == 3, "present-5",
                                                 ifelse(est == 4, "present", NA))))) %>% 
    mutate(full.siblings = full.sibs,
           diff.cohort_full.siblings = diff.cohort.sibs,
           imposters = as.numeric(NA)) %>% 
    dplyr::select(parameter,
                  estimation.sim, 
                  Q50,
                  all.truth,
                  HPD2.5,
                  HPD97.5,
                  T0,
                  estimation.year,
                  sampling.scheme,
                  iteration,
                  full.siblings,
                  diff.cohort_full.siblings,
                  imposters,
                  mean,
                  sd,
                  Q2.5,
                  Q97.5,
                  Rhat,
                  neff,
                  sample.prop,
                  seed) 
  
  if(test.decoys == "yes"){
    
    charlatan_HSPs <- imposters.df %>%
      dplyr::filter(iteration == iter,
                    sampling.scheme == s.scheme,
                    sample.prop == sample.proportion)
    
    if(filter.decoys == "yes"){
      
      #See what would happen if we filtered them by age difference
      charlatan_HSPs <- charlatan_HSPs %>% dplyr::filter(niece.nephew_birth.year - aunt.unc_birth.year < year.gap.threshold)
      
    }
    
    results.temp <- results.temp %>% mutate(imposters = nrow(charlatan_HSPs))
    
  }
  
  
  
  if(HS.only == "yes"){
    
    results.temp <- results.temp %>% 
      mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi", "Nf0", "Nft"), mom.HSPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2", "Nm0", "Nmt"), dad.HSPs, mom.HSPs + dad.HSPs)),
             HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi", "Nf0", "Nft"), mom.Exp.HS, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2", "Nm0", "Nmt"), dad.Exp.HS, mom.Exp.HS + dad.Exp.HS)),
             POPs_detected = NA,
             POPs_expected = NA,
             scenario = scenario)
    
  } else if(HS.only != "yes"){
    
    results.temp <- results.temp %>% 
      mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi", "Nf0", "Nft"), mom.HSPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2", "Nm0", "Nmt"), dad.HSPs, mom.HSPs + dad.HSPs)),
             HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi", "Nf0", "Nft"), mom.Exp.HS, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2", "Nm0", "Nmt"), dad.Exp.HS, mom.Exp.HS + dad.Exp.HS)),
             POPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi", "Nf0", "Nft"), mom.POPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2", "Nm0", "Nmt"), dad.POPs, mom.POPs + dad.POPs)),
             POPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi", "Nf0", "Nft"), mom.Exp.PO, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2", "Nm0", "Nmt"), dad.Exp.PO, mom.Exp.PO + dad.Exp.PO)),
             scenario = scenario)
    
  }
  
  return(results.temp)
}
