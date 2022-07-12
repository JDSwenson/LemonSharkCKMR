#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

#-------------------Set values for Leslie matrix---------------------#
juv.ages <- repro.age - 1 #Years of being a juvenile
ages <- c(0:max.age) #Total ages
adult.ages <- length(ages) - (juv.ages + 1)
vec.mean <- c(YOY.survival, juvenile.survival, adult.survival, leslie.fecundity)   #vector of the means

#Set CV and SD for survival
YOY.survival.sd <- YOY.surv.cv * YOY.survival #survival standard deviation
juv.survival.sd <- juv_adult.surv.cv * juvenile.survival
adult.survival.sd <- juv_adult.surv.cv * adult.survival

#Set CV and SD for fecundity
fec.cv <- fec.cv
fecundity.sd <- fec.cv * leslie.fecundity #Fecundity per breeding cycle i.e. ff

#create matrix of SD values
vec.sd <- c(YOY.survival.sd, juv.survival.sd, adult.survival.sd, fecundity.sd)
mat.sd <- diag(vec.sd)  

#Loop over potential correlations
# for(c in 1:length(corr.vec)){
#   correlation <- corr.vec[c]
#       
#   #create correlation matrix (assuming correlation = -0.5) ====
#   corr.fec.surv <- correlation
#   corr.mat <- matrix(c(1, corr.fec.surv, corr.fec.surv, 1), nrow=2, ncol=2, byrow=T)
#   # get covariance matrix
#   covar.mat <- mat.sd %*% corr.mat %*% mat.sd
      
  # fecSurv.draw <- mvrnorm(n=n.draws, mu=vec.mean, Sigma=covar.mat) %>% 
  #   as_tibble() %>% 
  #   rename(fec = V1, surv = V2) %>% 
  #   mutate(case = correlation) %>% 
  #   mutate(surv = ifelse(surv >=1, 0.99, surv)) #Make sure that survival is never greater than 1

lambda.df <- NULL #Initialize dataframe to store lambda values
lambda.temp.df <- NULL

batchSize <-  leslie.fecundity * mean(num.mates) #Adult fecundity -- average fecundity per breeding (2.91) x average number of mates (2)
  
#female offspring only
leslie.ff <- batchSize/2 #Adult fecundity - female offspring only (assume equal sexes)

#Loop over misspecified values
YOY.miss <- c(-0.10, -0.05, 0.05, 0.10)

for(y in 1:length(YOY.miss)){
  
  YOY.survival.miss <- YOY.survival + (YOY.survival*YOY.miss[y])

  YOY.survival.draw <- rnorm(n = n.draws, mean = YOY.survival.miss, sd = YOY.survival.sd)
  
  for(z in 1:length(YOY.survival.draw)){
    #Input to Leslie matrix
    survival.vec <- c(YOY.survival.draw[z], rep(juvenile.survival, times = juv.ages), rep(adult.survival, times = adult.ages - 1), 0)
    
    fecund.vec <- c(rep(0, times = repro.age), rep(leslie.ff, times = maxAge - juv.ages))

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = survival.vec, #survival
  mx = fecund.vec
)

A1_pre <- make_Leslie_matrix(Leslie_input)
#View(A1)
#A1_post <- pre_to_post(Amat = A1_pre, S0 = YOY.survival)
#View(A1_post)

#Calculate dominant eigenvalue (i.e. lambda) from transition matrix
(lambda1 <- as_tibble(lambda1(A1_pre)) %>% 
  rename(lambda = value))

#lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)
lambda.temp.df <- lambda1 %>% 
  mutate(misspecified.parameter = "YOY.survival",
         misspecified.percent = YOY.miss[y],
         mean.value = YOY.survival.miss,
         CV.value = YOY.surv.cv,
         YOY.survival = YOY.survival.draw[z], 
         juvenile.survival = juvenile.survival, 
         adult.survival = adult.survival, 
         fecundity = leslie.ff,
         iteration = z)

#Store lambda dataframe
lambda.df <- rbind(lambda.df, lambda.temp.df)

} #End loop over survival draws
print(paste0("Finished iteration all iterations with misspecification set to ", YOY.miss[y]*100, "%"))
} #End loop over misspecification percent


#Extract range of lambda values
(low.lambda <- min(lambda.df$lambda))
(high.lambda <- max(lambda.df$lambda))
(mean.lambda <- mean(lambda.df$lambda))
(lambda.sd <- sd(lambda.df$lambda))