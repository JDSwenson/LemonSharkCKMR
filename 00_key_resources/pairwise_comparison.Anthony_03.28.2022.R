####-------------ANTHONY SEVEQUE, MISSISSIPPI STATE UNIVERSITY, 2022--------------####
####-------------FORWARD IN TIME POPULATION SIMULATION AND KINSHIP ANALYSIS FOR CKMR USE IN LARGE TERRESTRIAL MAMMALS--------------####

####----------CONSTRUCT PAIRWISE COMPARISON MATRIX----------####

# Anthony: change all filter to subset at some point, to keep the script base R?

build.pairwise <- function(samples){

#Create dataframe of pairwise comparisons with just individual IDs

pairwise.df <- data.frame(t(combn(samples$indv.name, m=2))) # Generates all the possible combinations. m = 2 means pair comparisons. t(x) gives the dimensions of the matrix (dataframe).
colnames(pairwise.df) <- c("Ind_1", "indv.name") #Rename columns so they can easily be joined
head(pairwise.df)

# Ind 2 is called "indv.name" just to call it from the samples dataset at a later stage

#Create dataframe that will be used to extract the birth years for the younger individual from each pairwise comparison using joins.
Ind1_birthyears <- samples %>%
  select(indv.name, birth.year, sampling_year, mother, father, sex) %>% #select relevant columns only
  dplyr::rename("Ind_1" = indv.name, "Ind_1_birth" = birth.year, "Ind_1_sampling" = sampling_year, "Ind_1_mom" = mother, "Ind_1_dad" = father, "Ind_1_sex" = sex) 

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all comparisons (not all relevant for now) and individual data.

pairwise.df_all <- pairwise.df %>%
  left_join(Ind1_birthyears, by = "Ind_1") %>% 
  left_join(samples, by = "indv.name") %>% 
  dplyr::rename("Ind_2" = indv.name, "Ind_2_birth" = birth.year, "Ind_2_sampling" = sampling_year, "Ind_2_mom" = mother, "Ind_2_dad" = father, "Ind_2_sex" = sex) %>% 
  select(Ind_1, Ind_1_birth, Ind_1_sampling, Ind_1_sex, Ind_2, Ind_2_birth, Ind_2_sampling, Ind_2_sex, Ind_1_mom, Ind_1_dad, Ind_2_mom, Ind_2_dad) 

head(pairwise.df_all)

####----------PARENT OFFSPRING PAIRS----------####

# Remove comparisons that are "not biologically possible"

# Important: only adults should be sampled as potential parents (Ind 1). 
# If juveniles are sampled (genotyped), N goes up quasi exponentially (see Waples & Feutry, 2021)

pop.pairwise.df_filt <- pairwise.df_all %>%
  filter ((Ind_1_sampling - Ind_1_birth) >= repro.age) %>% # Remove potential parents that were harvested before reaching age at maturity
  filter ((Ind_2_birth - Ind_1_birth) >= repro.age) %>% # Remove potential parents that did not reach age at maturity when offspring was born
  filter ((Ind_2_birth - Ind_1_birth) <= max.age) %>% # Remove potential parents that were dead when offspring was born (can happen if more than one sampling occasion)
  filter (Ind_1_sampling >= Ind_2_birth) %>% # Remove potential parents that were sampled before offspring was born (because of lethal sampling). For now, let's assume that adults can reproduce and then die in the same year (could make the difference later for males and females).
  filter(Ind_1_birth != Ind_2_birth) #Filter intra-cohort comparisons (redundant with line 2 but keep it for clarity's sake)

# Extract positive POP comparisons
POP_positives_grouped <- pop.pairwise.df_filt %>%
  filter(Ind_1 == Ind_2_mom | Ind_1 == Ind_2_dad)

nrow(POP_positives_grouped)

###--- Sexes grouped ---###

# Extract positive POP comparisons
POP_positives <- pop.pairwise.df_filt %>%
  filter(Ind_1 == Ind_2_mom | Ind_1 == Ind_2_dad) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count() 

sum(POP_positives$freq)

# Make dataframes for negative comparisons
POP_negatives <- pop.pairwise.df_filt %>%
  filter(Ind_1 != Ind_2_mom & Ind_1 != Ind_2_dad) %>%
  select(Ind_1_birth, Ind_2_birth) %>%
  plyr::count()

sum(POP_negatives$freq)

POP_comps <- POP_positives %>% 
  dplyr::rename(yes = freq) %>% 
  full_join(POP_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  dplyr::rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>%
  mutate(all = yes + no) %>%
  dplyr::relocate(yes, .before=no)

sum(POP_comps$yes)

###--- Sex-specific ---###

POP_mom_positives <- pop.pairwise.df_filt %>%
  filter(Ind_1 == Ind_2_mom) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count() 

sum(POP_mom_positives$freq)

POP_dad_positives <- pop.pairwise.df_filt %>%
  filter(Ind_1 == Ind_2_dad)  %>%
  select(Ind_1_birth, Ind_2_birth) %>%
  plyr::count()

sum(POP_dad_positives$freq)

#Make dataframes for negative comparisons
POP_mom_negatives <- pop.pairwise.df_filt %>%
  filter(Ind_1 != Ind_2_mom) %>%
  filter(Ind_1_sex == "F") %>% # One condition for an individual to be a "potential mothers" is to be female. 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

POP_dad_negatives <- pop.pairwise.df_filt %>%
  filter(Ind_1 != Ind_2_dad) %>% 
  filter(Ind_1_sex == "M") %>%
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

POP_mom_comps <- POP_mom_positives %>% 
  dplyr::rename(yes = freq) %>% 
  full_join(POP_mom_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  dplyr::rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>%
  mutate(all = yes + no)

sum(POP_mom_comps$yes) # Just to verify with the positive pairwise comparison below

POP_dad_comps <- POP_dad_positives %>% 
  dplyr::rename(yes = freq) %>% 
  full_join(POP_dad_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  dplyr::rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)

sum(POP_dad_comps$yes) # Just to verify with the positive pairwise comparison below


sum(POP_mom_comps$all) + sum(POP_dad_comps$all) 
nrow (pop.pairwise.df_filt) # Should match

####----------HALF SIBLINGS PAIRS----------####

# Since we are sampling over a long time span, we need to remove comparisons with a number of years between os birth and ys birth
# higher than max.age - repro.age (i.e. by the time Ind2 is born, parent of Ind1 is assured to be dead).

hsp.pairwise.df_filt <- pairwise.df_all %>%
  filter ((Ind_2_birth - Ind_1_birth) <= (max.age - repro.age)) %>%
  filter(Ind_1_birth != Ind_2_birth) #Filter intra-cohort comparisons

# Next, we need to identify FSP, and remove one of the individuals (random) from the dataframe (not just the pairwise comparison, but the individual!) 
# This is to prevent dependence between two FSP and any third animal that is HSP (if HSP of one, then automatically HSP of the other)

full_sibs <- hsp.pairwise.df_filt %>% 
  filter (Ind_1_mom == Ind_2_mom & Ind_1_dad == Ind_2_dad)

print(paste0("There are ", nrow(full_sibs), " full siblings removed from the sample dataframe."))

hsp.pairwise.df_filt <- hsp.pairwise.df_filt[!(hsp.pairwise.df_filt$Ind_1 %in% full_sibs$Ind_1),]
hsp.pairwise.df_filt <- hsp.pairwise.df_filt[!(hsp.pairwise.df_filt$Ind_2 %in% full_sibs$Ind_1),]

###--- Sexes grouped ---###

HSP_positives <- hsp.pairwise.df_filt %>%
  filter(Ind_1_mom == Ind_2_mom | Ind_1_dad == Ind_2_dad) %>% 
  filter(Ind_1_mom != Ind_2_mom | Ind_1_dad != Ind_2_dad) %>% # Keep only comparisons where either the mother or the father is different (to remove full siblings)
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count() 

sum(HSP_positives$freq)

#Make dataframes for negative comparisons
HSP_negatives <- hsp.pairwise.df_filt %>% 
  filter(Ind_1_mom != Ind_2_mom & Ind_1_dad != Ind_2_dad) %>%
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

HSP_comps <- HSP_positives %>% 
  dplyr::rename(yes = freq) %>% 
  full_join(HSP_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  dplyr::rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)


###--- Sex-specific ---###

HSP_mom_positives <- hsp.pairwise.df_filt %>%
  filter(Ind_1_mom == Ind_2_mom) %>% 
  filter(Ind_1_dad != Ind_2_dad) %>% # Remove full siblings
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count() 

sum(HSP_mom_positives$freq)

HSP_dad_positives <- hsp.pairwise.df_filt %>%
  filter(Ind_1_dad == Ind_2_dad)  %>%
  filter(Ind_1_mom != Ind_2_mom) %>% # Remove full siblings
  select(Ind_1_birth, Ind_2_birth) %>%
  plyr::count()

sum(HSP_dad_positives$freq)

#Make dataframes for negative comparisons
HSP_mom_negatives <- hsp.pairwise.df_filt %>%
  filter(Ind_1_mom != Ind_2_mom) %>%
  filter(Ind_1_dad != Ind_2_dad) %>% # Remove paternal half sibling combinations too (because we know they can't be MHSB)?
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

HSP_dad_negatives <- hsp.pairwise.df_filt %>%
  filter(Ind_1_dad != Ind_2_dad) %>% 
  filter(Ind_1_mom != Ind_2_mom) %>% # Remove maternal half sibling combinations too?
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

HSP_mom_comps <- HSP_mom_positives %>% 
  dplyr::rename(yes = freq) %>% 
  full_join(HSP_mom_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  dplyr::rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)

sum(HSP_mom_comps$yes) # Just to verify with the positive pairwise comparison below

HSP_dad_comps <- HSP_dad_positives %>% 
  dplyr::rename(yes = freq) %>% 
  full_join(HSP_dad_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
  dplyr::rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)

sum(HSP_dad_comps$yes) # Just to verify with the positive pairwise comparison below


return(list(pairwise.df_all,
            POP_positives, POP_comps, POP_mom_comps, POP_dad_comps, 
            HSP_positives, HSP_comps, HSP_mom_comps, HSP_dad_comps))
}
