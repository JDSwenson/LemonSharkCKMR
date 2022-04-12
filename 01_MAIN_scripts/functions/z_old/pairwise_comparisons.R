
####-------------Construct pairwise comparison matrix--------------####

build.pairwise <- function(samples){
#Create dataframe of pairwise comparisons with just individual IDs
pairwise.df <- data.frame(t(combn(samples$indv.name, m=2))) # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df) <- c("Ind_1", "indv.name") #Rename columns so they can easily be joined
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
Ind1_birthyears <- samples %>%
  dplyr::select(indv.name, birth.year, age.x, mother.x, father.x) %>% #select relevant columns only
  dplyr::rename("Ind_1" = indv.name, "Ind_1_birth" = birth.year, "Ind_1_age" = age.x, "Ind_1_mom" = mother.x, "Ind_1_dad" = father.x) 

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info <- pairwise.df %>% left_join(Ind1_birthyears, by = "Ind_1") %>% 
  left_join(sample.df_all.info, by = "indv.name") %>% 
  dplyr::rename("Ind_2" = indv.name, "Ind_2_birth" = birth.year, "Ind_2_age" = age.x, "Ind_2_mom" = mother.x, "Ind_2_dad" = father.x) %>% 
  dplyr::select(Ind_1, Ind_1_birth, Ind_1_age, Ind_2, Ind_2_birth, Ind_2_age, Ind_1_mom, Ind_2_mom, Ind_1_dad, Ind_2_dad) %>% 
  filter(Ind_1_birth != Ind_2_birth) #Filter intra-cohort comparisons (for now)

#head(pairwise.df_all.info)

#Extract positive half-sib comparisons
positives <- pairwise.df_all.info %>% filter(Ind_1_mom == Ind_2_mom | Ind_1_dad == Ind_2_dad) 
#nrow(positives)

####----------------Split dataframes into final form for model----------####
#Sex-specific
mom_positives <- positives %>% filter(Ind_1_mom == Ind_2_mom) %>% 
  dplyr::select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count() 

dad_positives <- positives %>% filter(Ind_1_dad == Ind_2_dad)  %>%
  dplyr::select(Ind_1_birth, Ind_2_birth) %>%
  plyr::count()

#Make dataframes for negative comparisons
#Sex-specific
mom_negatives <- pairwise.df_all.info %>% filter(Ind_1_mom != Ind_2_mom & Ind_1_birth != Ind_2_birth) %>% #filter for same cohort is repetitive
  dplyr::select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

dad_negatives <- pairwise.df_all.info %>% filter(Ind_1_dad != Ind_2_dad & Ind_1_birth != Ind_2_birth) %>% 
  dplyr::select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

mom_comps <- mom_positives %>% 
  rename(yes = freq) %>% 
  full_join(mom_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no) #%>% 
  #filter(yes > 0)

dad_comps <- dad_positives %>% 
  rename(yes = freq) %>% 
  full_join(dad_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no) #%>% 
  #filter(yes > 0)


return(list(pairwise.df_all.info, positives, mom_comps, dad_comps))
}


filter_indvs <- function(samples, positives){
#Remove full sib comparisons -- 10/27/2021: potentially problematic bc we're only removing the comparison, not the individuals. Better to identify full sibs and remove the individuals.
#positives <- positives %>% filter(Ind_1_mom != Ind_2_mom | Ind_1_dad != Ind_2_dad) # Keeps only comparisons where either the mother or the father is different
#nrow(positives)

#Identify full sibs
full_sibs <- positives %>% filter(Ind_1_mom == Ind_2_mom & Ind_1_dad == Ind_2_dad)

print(paste0("There are ", nrow(full_sibs), " full siblings that will be removed from the sample dataframe."))

full_sibs2 <- c(full_sibs$Ind_1, full_sibs$Ind_2) %>% 
  as_tibble() %>% 
  rename(indv.name = value)

#Create new sample dataframe without full sibs
filtered_samples <- samples %>% anti_join(full_sibs2, by = "indv.name")

return(filtered_samples)
}