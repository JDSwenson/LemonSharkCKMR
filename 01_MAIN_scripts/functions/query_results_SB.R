####------------- Examine simulation results ---------------####
plot(pop.size.tibble$population_size, pch=19, ylim=c(0.9*min(pop.size.tibble$population_size), 1.1*max(pop.size.tibble$population_size))) #Plot population size through time

#nrow(YOY.df)/length(mothers) # Average fecundity for last year; remember whether they've skipped breeding

#Make dataframe of abundance for each year
pop.size.tibble <- pop.size.tibble %>% mutate(Total.adult.pop = Male.adult.pop + Female.adult.pop)

#Calculate population growth for whole population
total.lambda <- NULL
for(l in 2:nrow(pop.size.tibble)){ 
  total.lambda.1 <- pop.size.tibble$population_size[l]/pop.size.tibble$population_size[l-1]
  total.lambda <- c(total.lambda, total.lambda.1)
}

#Calculate population growth for adults only
adult.lambda <- NULL
for(l in 2:nrow(pop.size.tibble)){ 
  adult.lambda.1 <- pop.size.tibble$Total.adult.pop[l]/pop.size.tibble$Total.adult.pop[l-1]
  adult.lambda <- c(adult.lambda, adult.lambda.1)
}

#Add NA to first element of population growth vectors
total.lambda <- c(NA, total.lambda)
adult.lambda <- c(NA, adult.lambda)

#Add population growth per year to pop.size dataframe
pop.size.tibble$total.lambda <- total.lambda
pop.size.tibble$adult.lambda <- adult.lambda

#plot(total.lambda[(burn.in+1):n_yrs], pch=19)
#abline(h=1, lty=3)

(mean.total.lam <- mean(pop.size.tibble$total.lambda[(burn.in+1):n_yrs], na.rm=T)) # mean population growth for whole population
sd(pop.size.tibble$total.lambda[(burn.in+1):n_yrs], na.rm=T) # sd Lambda


#Calculate SURVIVAL for each year
sVec <- NULL #Make empty vector to save yearly survival rates

#Store annual survival of adults
for(yr in 2:length(loopy.list)){
  sVec[yr] <- length(which(loopy.list[[yr]]$Survival=='S' & loopy.list[[yr]]$age.x>=repro.age))/length(which(loopy.list[[yr]]$age.x>=repro.age))
}
sVec <- c(NA, sVec) #Add NA to survival vector for first year of simulation

length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age))/length(which(loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age)) #survival of juveniles

length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x==0))/length(which(loopy.list[[n_yrs]]$age.x==0)) #survival of YOY


#Calculate breeding interval
#Initialize sample dataframes
BI.df <- NULL
BI.df_temp <- NULL

for(i in sample.years){ #Extract all YOY for each year sampled
    BI.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
      dplyr::filter(age.x == 0)

        #Combine all
    BI.df <- bind_rows(BI.df, BI.df_temp)
}

#Extract unique mothers for each birth year
BI.df <- BI.df %>% as_tibble() %>% 
  arrange(birth.year) %>% 
  distinct(mother.x, birth.year, .keep_all = TRUE)

#-----------------Identify all relatives in population since estimation year and find breeding interval--------------------#
pairwise.BI.df <- data.frame(t(combn(BI.df$indv.name, m=2)))  # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.BI.df) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(BI.df)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.BI <- BI.df %>%
  select(older.sib = indv.name, 
         older.sib.birth = birth.year, 
         older.sib.age = age.x, 
         older.sib.mom = mother.x,
         older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.BI.all.info <- pairwise.BI.df %>% left_join(OlderSib_birthyears.BI, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(BI.df, by = "indv.name") %>% 
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
                younger.sib.dad) %>% 
  dplyr::filter(older.sib.birth != younger.sib.birth)


#Extract all positive maternal half-sib comparisons from sampled years
positives.mom.BI <- pairwise.BI.all.info %>% filter(older.sib.mom == younger.sib.mom) %>% 
  mutate(yr.gap = younger.sib.birth - older.sib.birth) %>% 
  mutate(BI = ifelse(yr.gap %% 2 == 0, "even", "odd")) %>% #If the year gap between ha;f-sibs is divisible by 2, then call the breeding interval (BI) "even"
  dplyr::distinct(older.sib.birth, younger.sib.birth, younger.sib.mom, .keep_all = TRUE) #Duplicative of earlier; making sure cohort size isn't biasing things here

positives.mom.BI %>% group_by(BI) %>% summarize(n()) #How many positive comparisons for 

skipped.moms <- positives.mom.BI %>% dplyr::filter(BI == "even") %>% 
  pull(younger.sib.mom) #Should be the same as the older sib mom
annual.moms <- positives.mom.BI %>% dplyr::filter(BI == "odd") %>% 
  pull(younger.sib.mom) #Should be the same as the older sib mom

#Which moms ONLY reproduced in even years?
skipped.only.moms <- skipped.moms[which(!skipped.moms %in% annual.moms)]
all.moms <- unique(positives.mom.BI$younger.sib.mom)

#Percent of individuals that only breed when years are evenly spaced
psi.truth <- length(skipped.only.moms)/length(all.moms)