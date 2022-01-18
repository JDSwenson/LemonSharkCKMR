
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
#Sex-specific, focus on gap year for prior
mom_positives.gap <- positives %>% filter(Ind_1_mom == Ind_2_mom) %>% 
  mutate(yr_gap = Ind_2_birth - Ind_1_birth) %>% 
  dplyr::select(yr_gap) %>% 
  plyr::count() 

dad_positives.gap <- positives %>% filter(Ind_1_dad == Ind_2_dad)  %>%
  mutate(yr_gap = Ind_2_birth - Ind_1_birth) %>% 
  dplyr::select(yr_gap) %>% 
  plyr::count() 

#Make dataframes for negative comparisons
#Sex-specific
mom_negatives.gap <- pairwise.df_all.info %>% filter(Ind_1_mom != Ind_2_mom & Ind_1_birth != Ind_2_birth) %>% #filter for same cohort is repetitive
  mutate(yr_gap = Ind_2_birth - Ind_1_birth) %>% 
  dplyr::select(yr_gap) %>% 
  plyr::count() 

dad_negatives.gap <- pairwise.df_all.info %>% filter(Ind_1_dad != Ind_2_dad & Ind_1_birth != Ind_2_birth) %>% 
  mutate(yr_gap = Ind_2_birth - Ind_1_birth) %>% 
  dplyr::select(yr_gap) %>% 
  plyr::count() 

mom_comps.gap <- mom_positives.gap %>% 
  rename(yes = freq) %>% 
  full_join(mom_negatives.gap, by = "yr_gap") %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no) %>% 
  filter(yes > 0)

mc <- mom_comps.gap %>% mutate(all.adj = all * (Adult.survival^yr_gap)) %>% 
  mutate(mom.N = all.adj/yes)

mc.all <- sum(mc$all.adj)
mc.pos <- sum(mc$yes)
mc.mean <- mean(mc$mom.N)
mc.sd <- sd(mc$mom.N)

mc.prior <- mc.all/mc.pos

dad_comps.gap <- dad_positives.gap %>% 
  rename(yes = freq) %>% 
  full_join(dad_negatives.gap, by = "yr_gap") %>% 
  rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no) %>% 
  filter(yes > 0)

dc <- dad_comps.gap %>% mutate(all.adj = all * (Adult.survival^yr_gap)) %>% 
  mutate(dad.N = all.adj/yes)

dc.all <- sum(dc$all.adj)
dc.pos <- sum(dc$yes)
dc.mean <- mean(dc$dad.N)
dc.sd <- sd(dc$dad.N)

dc.prior <- dc.all/dc.pos


#Sex-specific, w/ individual years
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
  mutate(all = yes + no)

dad_comps <- dad_positives %>% 
  rename(yes = freq) %>% 
  full_join(dad_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>%   rename(no = freq) %>% 
  mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
  mutate(all = yes + no)


return(list(pairwise.df_all.info, positives, mom_comps, dad_comps, mc.prior, mc.sd, dc.prior, dc.sd, mom_comps.gap, dad_comps.gap))
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