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