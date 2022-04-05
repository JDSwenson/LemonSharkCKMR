####---------------Filter samples -------------####
#Filters for full sibs, but also identifies appropriate sample sets for each offspring birth year (for PO) and each younger sib birth year (for half-sib)
filter.samples <- function(samples){

  #------------------Filter 1: full siblings----------------------
    NoFullSibs.df <- samples %>% distinct(mother.x, father.x, .keep_all = TRUE) %>%
    as_tibble() #If there is more than one individual with the same mother AND father, then only keep one.
    total.sibs <- nrow(samples) - nrow(NoFullSibs.df)
      
    dif.cohort.sibs <- samples %>% group_by(mother.x, father.x) %>%
      filter(n() > 1) %>% 
      distinct(birth.year) %>%
      ungroup() %>% 
      group_by(mother.x, father.x) %>% 
      summarize(distinct.birth.years = n()) %>% 
      filter(distinct.birth.years > 1) %>% 
      nrow()

    #The vector of years for which we need to split samples into potential parents and offspring i.e. offspring birth years.
  OffBirth.years <- NoFullSibs.df %>% distinct(birth.year) %>% 
    arrange(birth.year) %>%
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
      mutate(age.in.OffBirth.year = age.x + (OffBirth.year - capture.year)) %>% 
      filter(age.in.OffBirth.year >= repro.age) %>% 
      mutate(relation = "parent")
    
    #Dataframe of offspring born in this year
    offspring[[y]] <- NoFullSibs.df %>% 
      mutate(age.in.OffBirth.year = age.x + (OffBirth.year - capture.year)) %>% #unnecessary line, I think
      filter(birth.year == OffBirth.year) %>%
      mutate(relation = "offspring")
      
    #Add dataframe of all potential offspring and parents to a list where each dataframe corresponds to an offspring birth year
    PO.samps.list[[y]] <- rbind(parents[[y]], offspring[[y]]) %>% 
      dplyr::select(indv.name, 
                    birth.year, 
                    capture.year,
                    age.in.OffBirth.year,
                    mother.x, 
                    father.x,
                    sex,
                    relation)
    
    #Not sure if this is needed, but add dataframe of potential half-sib samples for the year to a list where each dataframe corresponds to a specific birth year
    # HS.samps.list[[y]] <- NoFullSibs.df %>%
    #   filter(birth.year <= OffBirth.year) #Filters for offspring that were born before or in Offbirth.year
  }
  
  #Name list elements to correspond to sample year
  names(PO.samps.list) <- paste0("PO.samples.year_", c(OffBirth.years))
  #names(HS.samps.list) <- paste0("HS.samples.year_", c(OffBirth.years))
  
  #Not sure if this is needed, but create dataframes of samples for PO and HS comparisons. Think only need this for half-sibs
#  PO.samps.df <- bind_rows(PO.samps.list) #Confirmed that it's the correct number of rows
  HS.samps.df <- NoFullSibs.df
    
  print(paste0("There are ", total.sibs, " pairs of full siblings in the dataset, ", dif.cohort.sibs, " of these pairs were born in different years."))
  return(list(PO.samps.list, HS.samps.df)) #Return list of possible parents and offspring for each year, and dataframe of potential half-sibs

}



####--------------Build pairwise comparison matrices------------####
#Input is the filtered list from above for POs, and the filtered dataframe for HS
build.pairwise <- function(filtered.samples.PO.list = PO.samps.list, filtered.samples.HS.df = HS.samps.df){

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
  select(offspring.birth, parent.capture.year) %>% 
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
  select(offspring.birth, parent.capture.year) %>% 
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
pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
  select(older.sib = indv.name, 
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
  select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  rename(yes = freq)

dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
  select(older.sib.birth, younger.sib.birth) %>%
  plyr::count() %>% 
  rename(yes = freq)


#Make dataframes for negative comparisons
#Sex-specific
# mom_negatives.HS <- pairwise.df_all.info.HS %>% filter(older.sib.mom != younger.sib.mom & older.sib.birth != younger.sib.birth) %>% #filter for same cohort is repetitive
#   select(older.sib.birth, younger.sib.birth) %>% 
#   plyr::count()
# 
# dad_negatives.HS <- pairwise.df_all.info.HS %>% filter(older.sib.dad != younger.sib.dad & older.sib.birth != younger.sib.birth) %>% 
#   select(older.sib.birth, younger.sib.birth) %>% 
#   plyr::count()

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
  dplyr::filter(mort.yrs < (2*repro.age)) %>% #Added 04/01/2022 to remove comparisons that could be confused for grand-parents
  mutate(type = "HS",
         parent = "mother")

HS.dad.pairwise.df <- dad_comps.HS %>%
  dplyr::rename(ref.year = younger.sib.birth,
                mort.yrs = year_gap) %>% 
  dplyr::select(ref.year, all, yes, mort.yrs) %>% 
  dplyr::filter(mort.yrs < (2*repro.age)) %>% #Added 04/01/2022 to remove comparisons that could be confused for grand-parents
  mutate(type = "HS",
         parent = "father")

#Combine PO and HS dataframes for each sex
mom_comps.all <- bind_rows(HS.mom.pairwise.df, PO.mom.pairwise.df) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  arrange(desc(ref.year), mort.yrs)

dad_comps.all <- bind_rows(HS.dad.pairwise.df, PO.dad.pairwise.df) %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  arrange(desc(ref.year), mort.yrs)

return(list(mom_comps.all, dad_comps.all))
}


filter_fullSibs <- function(samples, positives.HS){
#Remove full sib comparisons -- 10/27/2021: potentially problematic bc we're only removing the comparison, not the individuals. Better to identify full sibs and remove the individuals.
#positives <- positives %>% filter(Ind_1_mom != Ind_2_mom | Ind_1_dad != Ind_2_dad) # Keeps only comparisons where either the mother or the father is different
#nrow(positives)

#Identify full sibs
full_sibs <- positives.HS %>% filter(Ind_1_mom == Ind_2_mom & Ind_1_dad == Ind_2_dad)

print(paste0("There are ", nrow(full_sibs), " full siblings that will be removed from the sample dataframe."))

full_sibs2 <- c(full_sibs$Ind_1, full_sibs$Ind_2) %>% 
  as_tibble() %>% 
  rename(indv.name = value)

#Create new sample dataframe without full sibs for half-sib CKMR
filtered_samples.HS <- samples %>% anti_join(full_sibs2, by = "indv.name")


return(filtered_samples.HS)
}

#build.pairwise.PO <- function(samples){
#  samples %>% filter(age.x > )
#}