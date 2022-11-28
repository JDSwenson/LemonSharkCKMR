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