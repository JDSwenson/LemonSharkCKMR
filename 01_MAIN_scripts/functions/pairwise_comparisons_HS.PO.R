####-------------Construct pairwise comparison matrix--------------####
filter.samples <- function(samples){
  NoFullSibs.df <- samples %>% distinct(mother.x, father.x, .keep_all = TRUE) %>%
    as_tibble()

  PO.samps.list <- list()
  HS.samps.list <- list()
  parents <- list()
  offspring <- list()
  
  for(y in 1:length(estimated.years)){
    est.year <- estimated.years[y]
    
    parents[[y]] <- NoFullSibs.df %>% mutate(est.year = est.year, age.in.est.year = age.x + (est.year - capture.year)) %>%  #Create column for age in year of estimation to filter samples
      filter(age.in.est.year >= repro.age & capture.year == est.year) 
    
    offspring[[y]] <- NoFullSibs.df %>% mutate(est.year = est.year, age.in.est.year = age.x + (est.year - capture.year)) %>%  #Create column for age in year of estimation to filter samples
      filter(birth.year == est.year)
      
    PO.samps.list[[y]] <- rbind(parents[[y]], offspring[[y]])
    
    HS.samps.list[[y]] <- NoFullSibs.df %>% mutate(est.year = est.year) %>% 
      filter(birth.year <= est.year)
    
  }
  names(PO.samps.list) <- paste0("PO.samples.year_", c(estimated.years))
  names(HS.samps.list) <- paste0("HS.samples.year_", c(estimated.years))
  
  return(list(PO.samps.list, HS.samps.list))

}




build.pairwise <- function(filtered.samples.PO = filtered.samples.PO, filtered.samples.HS = filtered.samples.HS){
  PO.mom.pairwise.list <- list()
  HS.mom.pairwise.list <- list()

  
  PO.dad.pairwise.list <- list()
  HS.dad.pairwise.list <- list()
  
  for(y in 1:length(estimated.years)){
    est.year <- estimated.years[y]
    
    pairwise.df.PO <- data.frame(t(combn(filtered.samples.PO[[y]]$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
    colnames(pairwise.df.PO) <- c("Ind_1", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.PO)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
Ind1_birthyears.PO <- filtered.samples.PO[[y]] %>%
  select(indv.name, birth.year, age.x, mother.x, father.x, parent.capture.year = capture.year) %>% #select relevant columns only
  dplyr::rename("Ind_1" = indv.name, "Ind_1_birth" = birth.year, "Ind_1_age" = age.x, "Ind_1_mom" = mother.x, "Ind_1_dad" = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.PO <- pairwise.df.PO %>% left_join(Ind1_birthyears.PO, by = "Ind_1") %>% 
  left_join(sample.df_all.info, by = "indv.name") %>% 
  dplyr::rename("Ind_2" = indv.name, "Ind_2_birth" = birth.year, "Ind_2_age" = age.x, "Ind_2_mom" = mother.x, "Ind_2_dad" = father.x) %>%
  select(Ind_1, Ind_1_birth, Ind_1_age, Ind_2, Ind_2_birth, Ind_2_age, Ind_1_mom, Ind_2_mom, Ind_1_dad, Ind_2_dad, parent.capture.year) %>% 
  filter(Ind_1_birth != Ind_2_birth) %>% #Filter intra-cohort comparisons (for now)
  filter(Ind_2_birth == est.year)

#head(pairwise.df_all.info)

#Sex-specific PO
#Extract unique parents for sampled individuals
unique_mothers <- filtered.samples.PO[[y]] %>% 
  distinct(mother.x) %>% 
  pull()

unique_fathers <- filtered.samples.PO[[y]] %>% 
  distinct(father.x) %>% 
  pull()

sampled.parents <- filtered.samples.PO[[y]] %>% filter(indv.name %in% unique_mothers | indv.name %in% unique_fathers) %>% 
  pull(indv.name)

#Can whittle down comparisons/samples here, if needed
parents.df <- filtered.samples.PO[[y]] %>% filter(indv.name %in% sampled.parents) %>%
  select(sampled.parent.id = indv.name, parent.birth.year = birth.year, parent.capture.year = capture.year)

# mothers.df <- samples %>% filter(indv.name %in% sampled.parents & sex == "F") %>% 
#   select(sampled.parent.id = indv.name, parent.birth.year = birth.year, parent.capture.year = capture.year)
# 
# fathers.df <- samples %>% filter(indv.name %in% sampled.parents & sex == "M") %>%
#   select(sampled.parent.id = indv.name, parent.birth.year = birth.year, parent.capture.year = capture.year)

positives.PO <- filtered.samples.PO[[y]] %>% filter(mother.x %in% sampled.parents | father.x %in% sampled.parents) %>% #Select only sampled offspring with sampled parent
  mutate(sampled.parent = ifelse(mother.x %in% sampled.parents, "mother", "father")) %>% #determine which parent has been sampled
  mutate(sampled.parent.id = ifelse(sampled.parent == "mother", mother.x, father.x)) %>% #extract id of sampled parent
  select(offspring = indv.name, offspring.birth.year = birth.year, offspring.capture.year = capture.year, sampled.parent, sampled.parent.id) %>% #select columns of interest
  left_join(parents.df, by = "sampled.parent.id") #join with parent info

mom_positives.PO <- positives.PO %>% filter(sampled.parent == "mother") %>% #filter for positives with mother
  select(offspring.birth.year, parent.birth.year, parent.capture.year) %>% 
  plyr::count() 

dad_positives.PO <- positives.PO %>% filter(sampled.parent == "father") %>% #filter for positives with father
  select(offspring.birth.year, parent.birth.year, parent.capture.year) %>% 
  plyr::count() 

#Make dataframes for all comparisons
#Sex-specific
mom_comps.PO <- pairwise.df_all.info.PO %>% #filter for same cohort is repetitive
  select(parent.birth.year = Ind_1_birth, offspring.birth.year = Ind_2_birth, parent.capture.year) %>% 
  plyr::count() %>% 
  rename(all = freq) %>% 
  left_join(mom_positives.PO, by = c("parent.birth.year", "offspring.birth.year", "parent.capture.year")) %>% 
  rename(yes = freq) %>%
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(no = all - yes) %>% 
  mutate(no = replace_na(no, 0)) %>% 
  filter(offspring.birth.year - parent.birth.year >= repro.age)


dad_comps.PO <- pairwise.df_all.info.PO %>% #filter for same cohort is repetitive
  select(parent.birth.year = Ind_1_birth, offspring.birth.year = Ind_2_birth, parent.capture.year) %>% 
  plyr::count() %>% 
  rename(all = freq) %>% 
  left_join(dad_positives.PO, by = c("parent.birth.year", "offspring.birth.year", "parent.capture.year")) %>% 
  rename(yes = freq) %>%
  mutate(yes = replace_na(yes, 0)) %>% 
  mutate(no = all - yes) %>% 
  mutate(no = replace_na(no, 0)) %>% 
  filter(offspring.birth.year - parent.birth.year >= repro.age)


PO.mom.pairwise.list[[y]] <- mom_comps.PO
names(PO.mom.pairwise.list)[[y]] <- paste0("PO.mom.pairwise.year_", est.year)


PO.dad.pairwise.list[[y]] <- dad_comps.PO
names(PO.dad.pairwise.list)[[y]] <- paste0("PO.dad.pairwise.year_", est.year)





#--------------Half sibling--------------------------
pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS[[y]]$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("Ind_1", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
Ind1_birthyears.HS <- filtered.samples.HS[[y]] %>%
  select(indv.name, birth.year, age.x, mother.x, father.x, parent.capture.year = capture.year) %>% #select relevant columns only
  dplyr::rename("Ind_1" = indv.name, "Ind_1_birth" = birth.year, "Ind_1_age" = age.x, "Ind_1_mom" = mother.x, "Ind_1_dad" = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(Ind1_birthyears.HS, by = "Ind_1") %>% 
  left_join(sample.df_all.info, by = "indv.name") %>% 
  dplyr::rename("Ind_2" = indv.name, "Ind_2_birth" = birth.year, "Ind_2_age" = age.x, "Ind_2_mom" = mother.x, "Ind_2_dad" = father.x) %>%
  select(Ind_1, Ind_1_birth, Ind_1_age, Ind_2, Ind_2_birth, Ind_2_age, Ind_1_mom, Ind_2_mom, Ind_1_dad, Ind_2_dad, parent.capture.year) %>% 
  filter(Ind_1_birth != Ind_2_birth) %>% #Filter intra-cohort comparisons (for now)
  filter(Ind_2_birth == est.year)

#Extract positive half-sib comparisons
positives.HS <- pairwise.df_all.info.HS %>% filter(Ind_1_mom == Ind_2_mom | Ind_1_dad == Ind_2_dad) 
#nrow(positives)

####----------------Split dataframes into final form for model----------####
#Sex-specific half-sib
mom_positives.HS <- positives.HS %>% filter(Ind_1_mom == Ind_2_mom) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count() 

dad_positives.HS <- positives.HS %>% filter(Ind_1_dad == Ind_2_dad)  %>%
  select(Ind_1_birth, Ind_2_birth) %>%
  plyr::count()

#Make dataframes for negative comparisons
#Sex-specific
mom_negatives.HS <- pairwise.df_all.info.HS %>% filter(Ind_1_mom != Ind_2_mom & Ind_1_birth != Ind_2_birth) %>% #filter for same cohort is repetitive
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

dad_negatives.HS <- pairwise.df_all.info.HS %>% filter(Ind_1_dad != Ind_2_dad & Ind_1_birth != Ind_2_birth) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

mom_comps.HS <- mom_positives.HS %>% 
  rename(yes = freq) %>% 
  full_join(mom_negatives.HS, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)

dad_comps.HS <- dad_positives.HS %>% 
  rename(yes = freq) %>% 
  full_join(dad_negatives.HS, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)



HS.mom.pairwise.list[[y]] <- mom_comps.HS
names(HS.mom.pairwise.list)[[y]] <- paste0("HS.mom.pairwise.year_", est.year)

HS.dad.pairwise.list[[y]] <- dad_comps.HS
names(HS.dad.pairwise.list)[[y]] <- paste0("HS.dad.pairwise.year_", est.year)
}


return(list(PO.mom.pairwise.list, PO.dad.pairwise.list, HS.mom.pairwise.list, HS.dad.pairwise.list))
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