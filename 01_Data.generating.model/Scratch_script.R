#Find aunt/niece pairs
#Make dataframe of sampled individuals parents and the sampled indv's birth year (so we can subset from loopy.list). Preserve column names for join
#Would have to insert this into the sampling scheme loop ... 
sampled.indvs.rents.df <- sample.info2 %>% 
  dplyr::select(mother.x, father.x, birth.year) %>% 
  arrange(birth.year)

#Make a vector with all unique parents
rents.vec <- unique(c(sampled.indvs.rents.df$mother.x, sampled.indvs.rents.df$father.x))


#Make a dataframe of all loopy.list dataframes from the oldest sampled indv's birth year until present
b.y <- min(sampled.indvs.rents.df$birth.year)

#Bind rows of loopy.list together so we can grab the info from all parents in the dataset
#Confirmed it's doing what we want it to
loopy.list.merged.df <- bind_rows(loopy.list[b.y:n_yrs]) %>%
  distinct(indv.name, .keep_all = TRUE) %>% #Only keep one copy of each individual
  dplyr::filter(indv.name %in% rents.vec) %>% #Only keep individuals that are parents of sampled individuals
  as_tibble() %>% 
  rename(grandmother = mother.x, grandfather = father.x) #The parents of the sampled indvidual's parents are the grandparents of the sampled individuals

#Split into separate dfs for maternal and paternal lineages. Keep the sampeld indv as the reference point.
mothers.grandparents.df <- loopy.list.merged.df %>% dplyr::filter(sex == "F") %>% 
  dplyr::select(mother = indv.name, 
         mom.grandmother = grandmother, 
         mom.grandfather = grandfather)

fathers.grandparents.df <- loopy.list.merged.df %>% dplyr::filter(sex == "M") %>% 
  dplyr::select(father = indv.name,
         dad.grandmother = grandmother,
         dad.grandfather = grandfather)

#Add grandparents to each sample dataframe
sample.info_w_grandparents <- sample.info2 %>% 
  rename(mother = mother.x, father = father.x) %>% 
  inner_join(mothers.grandparents.df, by = "mother") %>% #add maternal grandmother info to the df, joining by mom
  inner_join(fathers.grandparents.df, by = "father") %>% #add paternal grandmother info to the df, joining by dad
  dplyr::select(indv.name, mother, father, mom.grandmother, mom.grandfather, dad.grandmother, dad.grandfather, birth.year, capture.year, sampling.scheme, iteration, seed, sample.prop, ref.yr) %>% #Keep most columns but remove the x from mother and father
  mutate(both.rents = paste0(mother, father), #Reduce parents to one variable
         both.mom.grandparents = paste0(mom.grandmother, mom.grandfather),
         both.dad.grandparents = paste0(dad.grandmother, dad.grandfather)) %>% 
  group_by(mother, father, birth.year) %>% #Filter for full sibs 
  slice_sample(n = 1) %>% #Keep one sibling. The pairwise comparison script just keeps one indv from each litter, so we can do that here too.
  ungroup()

#Q for Ben or Anthony: is it good that we're only keeping one indv from each litter? Or can we easily integrate full sibs into the likelihood?
#Check how many full sibs there are after filtering
sample.info_w_grandparents %>% 
  group_by(mother, father, birth.year) %>%
  summarize(n=n()) %>%
  dplyr::filter(n>1)

###---- Testing ----###
#Save each combos of parents as vector objects to ease code later
both.parents.vec <- sample.info_w_grandparents$both.rents
both.mom.grandparents.vec <- sample.info_w_grandparents$both.mom.grandparents
both.dad.grandparents.vec<- sample.info_w_grandparents$both.dad.grandparents

#Save dataframe of sampled aunt/uncles
sampled_aunts_or_uncs.df <- sample.info_w_grandparents %>% dplyr::filter(both.rents %in%  both.mom.grandparents.vec | both.rents %in% both.dad.grandparents.vec) %>% #If the sampled individual's parents are the grandparents of another sampled individual, then this is indv is an aunt or uncle.
  mutate(shared.relation = both.rents) %>% #Save the shared relation for joins
  dplyr::select(aunt.unc = indv.name, #Rename columns for joins
                aunt.unc_mother = mother, 
                aunt.unc_father = father, 
                shared.relation,
                aunt.unc_birth.year = birth.year)

#Save dataframe of sampled nieces/nephews. Easier to split maternal and paternal and then join. First focus on maternal line.
sampled_nieces_or_nephews_maternal.df <- sample.info_w_grandparents %>% dplyr::filter(both.mom.grandparents %in% both.parents.vec) %>% #If the indvidual's maternal grandparents are the parents to another sampled individual, then this is a niece/nephew
  mutate(shared.relation = both.mom.grandparents) %>% #Save the shared relation for joins and for double-checking code
  dplyr::select(niece.nephew = indv.name, #Rename columns for joins
                shared.relation, 
                niece.nephew_maternal.grandmother= mom.grandmother, 
                niece.nephew_maternal.grandfather = mom.grandfather,
                niece.nephew_birth.year = birth.year)

#Isolate paternal nieces and nephews, then join with maternal niece/nephews
  all_sampled_nieces_or_nephews.df <- sample.info_w_grandparents %>% 
    dplyr::filter(both.dad.grandparents %in% both.parents.vec) %>%  #If the individual's paternal grandparents are the parents to another sampled individual, then this is a niece/nephew
    mutate(shared.relation = both.dad.grandparents) %>% #Save shared relation for join
    dplyr::select(niece.nephew = indv.name, #Rename for join
                  shared.relation, 
                  niece.nephew_paternal.grandmother = dad.grandmother, 
                  niece.nephew_paternal.grandfather = dad.grandfather,
                  niece.nephew_birth.year = birth.year) %>% 
    full_join(sampled_nieces_or_nephews_maternal.df) #Join with maternal nieces and nephews. If related through grandparents in the maternal line, then the paternal line will be NA and vice versa.

  #Combine niece/nephew dataframe with aunt/uncle dataframe
  aunt.unc_niece.nephew_pw.comps.all <- sampled_aunts_or_uncs.df %>% 
    full_join(all_sampled_nieces_or_nephews.df, by = c("shared.relation")) #This is the key join. It makes a dataframe where we join by the shared relation, and then we should have the aunt/uncle and niece/nephew name and information in one dataframe, with the positives linked up.

#GOT EM
charlatan_HSPs <- aunt.unc_niece.nephew_pw.comps.all %>% 
    dplyr::select(older.sib.birth = aunt.unc_birth.year,
                  younger.sib.birth = niece.nephew_birth.year)
  
#See what would happen if we filtered them by age difference
charlatan_HSPs %>% dplyr::filter(younger.sib.birth - older.sib.birth < 12)

#Seems like we could basically eliminate them. Probably a function of the four year sampling scheme ... we'll see how the full simulations turn out ... 
  
  
  
  
  
    
shit <- tibble(aunt = as.character(NA))

#Now, isolate instances where the mother AND father of the sampled individual (aunt) matches the grandmother and grandfather of any sampled individual (niece)
for(i in 1:nrow(sample.info_w_grandrents)){
  #Try matching at once with which
sample.info_w_grandrents[which(sample.info_w_grandrents$both.rents[i] %in% sample.info_w_grandrents$both.mom.grandrents), ]
  

which(sample.info_w_grandrents$both.rents %in% sample.info_w_grandrents$both.dad.grandrents)
which(sample.info_w_grandrents$both.rents %in% sample.info_w_grandrents$both.mom.grandrents)

#Try a long loop to triple check that this is true
  for(j in 1:nrow(sample.info_w_grandrents)){
    
    if(sample.info_w_grandrents$mother.x[i] == sample.info_w_grandrents$mom.grandmother[j] & sample.info_w_grandrents$father.x[i] == sample.info_w_grandrents$mom.grandfather[j]){
    
       shit.temp <- tibble(aunt = sample.info_w_grandrents$indv.name[i])
       shit <- bind_rows(shit, shit.temp)
    }
  }
}


#Check proportion of conformists and non-conformists that are breeding each year
for(i in 1:length(loopy.list)){
Newmoms.vec <- loopy.list[[i]] %>% dplyr::filter(age.x == 0) %>% pull(mother.x)

Pot.conf.moms <- loopy.list[[i]] %>% dplyr::filter(sex == "F", 
                                                    age.x == repro.age,
                                                    repro.strategy == "conformist")


Newmoms.conf.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec, 
                                                 age.x == repro.age,
                                                 repro.strategy == "conformist")

Pot.Nonconf.moms <- loopy.list[[i]] %>% dplyr::filter(sex == "F", 
                                                   age.x == repro.age,
                                                   repro.strategy == "non-conformist")


Newmoms.Nonconf.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec, 
                                                     age.x == repro.age,
                                                     repro.strategy == "non-conformist")

moms.all.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec)

moms.all <- moms.all.df %>% group_by(repro.strategy) %>% 
  summarize(n = n()) %>% 
  mutate(perc.conf = round(n/sum(n) * 100, 0)) %>% 
  dplyr::filter(repro.strategy == "conformist") %>% 
  pull(perc.conf)


print(paste0("Year ", i, " prop conformists moms: ", round((nrow(Newmoms.conf.df)/nrow(Pot.conf.moms))*100), "%"))
print(paste0("Year ", i, " prop non-conformists moms: ", round((nrow(Newmoms.Nonconf.df)/nrow(Pot.Nonconf.moms))*100), "%"))
print(paste0("Year ", i, " percent conformist all moms: ", moms.all, "%"))
}



Newmoms.df %>% gghistogram(x = "age.x", color = "repro.strategy", fill = "repro.strategy")





data1 %>% dplyr::filter(sex == "F" & age.x >= repro.age)

init.pop2 %>% dplyr::filter(sex == "F" & age.x >= repro.age) %>% 
  group_by(repro.strategy) %>% 
  summarize(num = n())
  

data1.temp %>% dplyr::filter(sex == "F", 
                        age.x == repro.age,
                        is.na(repro.strategy) == TRUE,
                        is.na(repro.cycle) == TRUE)


data1.temp %>% dplyr::filter(sex == "F",
                             age.x == repro.age)


loopy.list[[90]] %>% dplyr::filter(sex == "F") %>% 
  dplyr::select(birth.year, repro.cycle)


parents.tibble_all %>% arrange(num.off) %>% 
  group_by(parent.sex, year) %>% 
  summarize(mean.off = mean(num.off)) %>% 
  dplyr::filter(parent.sex == "mother") %>% 
  summarize(mean(mean.off))

pop.size.tibble_all %>% dplyr::filter(iteration == 2) %>% 
  View()


truth.all %>% View()

pop.size.tibble_all %>% View()


#Check age distribution
g1 <- loopy.list[[10]] %>% gghistogram(x = "age.x", binwidth = 1)
g2 <- loopy.list[[20]] %>% gghistogram(x = "age.x", binwidth = 1)
g3 <- loopy.list[[30]] %>% gghistogram(x = "age.x", binwidth = 1)
g4 <- loopy.list[[70]] %>% gghistogram(x = "age.x", binwidth = 1)
g5 <- loopy.list[[80]] %>% gghistogram(x = "age.x", binwidth = 1)
g6 <- loopy.list[[90]] %>% gghistogram(x = "age.x", binwidth = 1)

ggarrange(g1, g2, g3, g4, g5, g6, common.legend = TRUE)


#Check that females breed, on average, with two males (they do)
loopy.list[[80]] %>% dplyr::filter(age.x == 0) %>% 
  distinct(mother.x, father.x, .keep_all = TRUE) %>% 
  group_by(mother.x) %>% 
  summarize(num.fathers = n()) %>% 
  summarize(mean = mean(num.fathers, na.rm = TRUE))

#Check that all males and females breed for the first time at age 12 (they do)
yr = 80
pt.temp <- parents.tibble_all %>% dplyr::filter(year == yr) %>% 
  pull(parent)

loopy.list[[yr]] %>% dplyr::filter(indv.name %in% pt.temp) %>% 
  arrange(age.x) %>% 
  head()

sample.df_all.info %>% group_by(age.x) %>% 
  summarize(num = n()) %>% 
  arrange(age.x)


length(s1.2)
length(s2.2)
length(s3.2)
length(s4.2)


library(jagsUI)





BI.df.87 <- BI.df %>% dplyr::filter(birth.year == 87)
BI.df.88 <- BI.df %>% dplyr::filter(birth.year == 88)
BI.df.89 <- BI.df %>% dplyr::filter(birth.year == 89)
BI.df.90 <- BI.df %>% dplyr::filter(birth.year == 90)

#How many mothers gave birth in both years?
BI.df.87 %>% dplyr::filter(mother.x %in% BI.df.89$mother.x) %>% 
  nrow()

BI.df.89 %>% dplyr::filter(mother.x %in% BI.df.90$mother.x) %>% 
  nrow()

#How many mothers from 87 were alive in 89?
BI.df.87 %>% dplyr::filter(mother.x %in% loopy.list[[89]]$indv.name) %>% 
  nrow()

BI.df.88 %>% dplyr::filter(mother.x %in% BI.df.90$mother.x) %>% 
  nrow()

#Make sure simulation is doing what I think
head(sims.list.1)
head(sims.list.2)
head(sims.list.3)
head(sims.list.4)
head(parents.tibble)
head(results2)
head(pop.size.tibble)
head(sample.info)
head(mom.comps.tibble)
head(dad.comps.tibble)
tail(mom.comps.tibble)
tail(dad.comps.tibble)


#Check that downsample, HS, PO are working as expected
mom_comps.all %>% group_by(type) %>%
  summarize(sum(yes))

#Should be the same as above
mom.comps.tibble %>% dplyr::filter(iteration == 2, sample.size.juvs == 800) %>% 
  group_by(type) %>%
  summarize(sum(yes))

#Dads
dad_comps.all %>% group_by(type) %>%
  summarize(sum(yes))

#Should be the same as above
dad.comps.tibble %>% dplyr::filter(iteration == 2, sample.size.juvs == 800) %>% 
  group_by(type) %>%
  summarize(sum(yes))

#Check that the correct age classes were sampled
sample.info %>% dplyr::filter(age.x < repro.age) %>% 
  group_by(age.x) %>% 
  summarize(n())





pop.size.1.summary <- pop.size.1 %>% dplyr::filter(year == 85) %>% 
  dplyr::distinct(year, iter, .keep_all = TRUE) %>%
  dplyr::select(!Total.adult.pop) %>% 
 # dplyr::rename(Nm = Male.adult.pop, Nf = Female.adult.pop) %>% 
  pivot_longer(cols = ends_with("adult.pop"), names_to = "col", values_to = "truth.N") %>% 
  mutate(parameter = ifelse(col == "Male.adult.pop", "Nm", "Nf")) %>% 
  dplyr::select(parameter, truth.N, iteration = iter)

pop.size.1.lam <- pop.size.1 %>% dplyr::filter(year >= 85) %>% 
  group_by(iter) %>% 
  summarize(true.lam = mean(adult.lambda)) %>% 
  dplyr::rename(iteration = iter)


results.1 <- results.1 %>% left_join(pop.size.1.summary, by = c("iteration", "parameter")) %>% 
  left_join(pop.size.1.lam, by = "iteration") %>% 
  mutate(truth = ifelse(parameter == "Nf" | parameter == "Nm", truth.N, 
                        ifelse(parameter == "lam", true.lam, truth))) 


pop.size.2.summary <- pop.size.2 %>% dplyr::filter(year == 85) %>% 
  dplyr::distinct(year, iter, .keep_all = TRUE) %>%
  dplyr::select(!Total.adult.pop) %>% 
  # dplyr::rename(Nm = Male.adult.pop, Nf = Female.adult.pop) %>% 
  pivot_longer(cols = ends_with("adult.pop"), names_to = "col", values_to = "truth.N") %>% 
  mutate(parameter = ifelse(col == "Male.adult.pop", "Nm", "Nf")) %>% 
  dplyr::select(parameter, truth.N, iteration = iter)


results.2 <- results.2 %>% left_join(pop.size.2.summary, by = c("iteration", "parameter")) %>% 
  left_join(pop.size.2.lam, by = "iteration") %>% 
  mutate(truth = ifelse(parameter == "Nf" | parameter == "Nm", truth.N, 
                        ifelse(parameter == "lam", true.lam, truth))) 

pop.size.2.lam <- pop.size.2 %>% dplyr::filter(year >= 85) %>% 
  group_by(iter) %>% 
  summarize(true.lam = mean(adult.lambda)) %>% 
  dplyr::rename(iteration = iter)




purpose1 <- "HS.PO_refined.samples_all.comps"
purpose2 <- "HS.PO_refined.samples_downsample"
purpose3 <- "HS.only_refined.samples_all.comps"
purpose4 <- "HS.only_refined.samples_downsample"




#--------------Posterior predictive distribution examples------------------------
#-----------------------Example 1----------------------
#Get data
data(longley)
gnp <- longley$GNP
employed <- longley$Employed
n <- length(employed)
data <- list(gnp=gnp,employed=employed,n=n)

#Identify filepath of model file
modfile <- tempfile()

#Write model
#Note calculation of discrepancy stats fit and fit.new
#(sums of residuals)
writeLines("
model{

  #Likelihood
  for (i in 1:n){ 

    employed[i] ~ dnorm(mu[i], tau)     
    mu[i] <- alpha + beta*gnp[i]
    
    res[i] <- employed[i] - mu[i]   
    emp.new[i] ~ dnorm(mu[i], tau)
    res.new[i] <- emp.new[i] - mu[i]

  }
    
  #Priors
  alpha ~ dnorm(0, 0.00001)
  beta ~ dnorm(0, 0.00001)
  sigma ~ dunif(0,1000)
  tau <- pow(sigma,-2)
  
  #Derived parameters
  fit <- sum(res[])
  fit.new <- sum(res.new[])

}
", con=modfile)

#Set parameters to monitor
params <- c('alpha','beta','sigma','fit','fit.new')

#Run analysis

out <- jags(data = data,
            inits = NULL,
            parameters.to.save = params,
            model.file = modfile,
            n.chains = 3,
            n.adapt = 100,
            n.iter = 1000,
            n.burnin = 500,
            n.thin = 2)

#Examine output summary
out

#Posterior predictive check plot
pp.check(out, observed = 'fit', simulated = 'fit.new')

#-------------------------Example 2------------------------------------
set.seed(02162018)
num.obs <- 50
true.mu <- 0
true.sigmasq <- 1
y <- rnorm(num.obs, mean = true.mu, sd = sqrt(true.sigmasq))

#Specify priors
M <- 0
S <- 100
C <- 100000

dataList = list(y = y, Ntotal = num.obs, M = M, S = S, C = C)

modelString = "model {
for ( i in 1:Ntotal ) {
y[i] ~ dnorm(mu, 1/sigma^2) # sampling model
}
mu ~ dnorm(M,1/S^2)
sigma ~ dunif(0,C)
} "
writeLines( modelString, con='NORMmodel.txt')

initsList <- function(){
  # function for initializing starting place of theta
  # RETURNS: list with random start point for theta
  return(list(mu = rnorm(1, mean = 0, sd = 100), sigma = runif(1,0,1000)))
}

library(rjags)
library(runjags)

jagsModel <- jags.model( file = "NORMmodel.txt", data = dataList,
                         inits =initsList, n.chains = 2, n.adapt = 100)

update(jagsModel, n.iter = 500)

num.mcmc <- 1000
codaSamples <- coda.samples( jagsModel, variable.names = c('mu', 'sigma'), n.iter = num.mcmc)

#Calculate posterior predictive distribution
posterior.mu <- codaSamples[[1]][,'mu']
posterior.sigma <- codaSamples[[1]][,'sigma']
posterior.pred <- rnorm(num.mcmc, mean = posterior.mu, sd = posterior.sigma)
prob.greater <- mean(posterior.pred > -0.2)

#For a posterior predictive interval, I should be able to use a binomial distribution with each set of parameter values from the posterior distribution, and the CKMR equation. This will give the number of predicted kin pairs from the model.