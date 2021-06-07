library(rjags)
library(R2jags)

#input dataframes are mom_positives, dad_positives, mom_negatives, and dad_negatives. Also adult_positives and adult_negatives

#Do sex-aggregated first
#Need priors for:
#Number of adults (uninformative)
#survival (beta -- conjugate prior for binomial)

##CKMR code
cat("model{

	# Adult numbers (uninformative prior)
	for(i in 1:nad) {
		Nadult[i] ~ dnorm(0,1.0E-6)
	}

	M ~ dbeta(1,1)

	for(j in 1:ncom) {
	  POPs[j] ~ dbern((OldEno[j] * ((M)^ParYear[j]))/(Nadult[BYear[j]])) # Repro chances reduced by chance that adult died prior to spawning
	}

}", file = "C:\\Users\\Ben\\Google Drive\\Research\\HAL CKMR\\Estimation Model\\JAGS_Test\\CKMRmultiY.jags")


# Data for JAGS
data = list(nad = nad,
            ncom = ncom,
            BYear = BYear,
            POPs = POPs,
            OldEno = ROs,
            ParYear = ParYear)

# Initial values
inits = function() {
  list(Nadult = rep(1000,nad),
       M = 0.5)
}

# Parameters to follow
params = c("Nadult", "M")

Out = jags(data, inits, params, "C:\\Users\\Ben\\Google Drive\\Research\\HAL CKMR\\Estimation Model\\JAGS_Test\\CKMRmultiY.jags", 
           n.burnin = nburn, n.chains = nchains, n.iter = niter, parallel = parallel, n.cores = n.cores, verbose = verbose)








#Other lecture code
cat("
model {
 B0 ~ dnorm(0,1.0E-6)        # uniformative prior
 Binf ~ dnorm(0,1.0E-6)      # uniformative prior
 m ~ dnorm(0,1.0E-6)         # uninformative prior 
 prec ~ dgamma(0.001,0.001)  # prior for the precision
 B[1] <- B0                  # initial year
 Bobs[1] ~ dnorm(B[1],prec)  # draw from normal

 for( i in 1:13)             # for remaining years
 { 	
   B[i+1] <- B[i]+(4*m/Binf)*(1-B[i]/Binf)*B[i]-Cobs[i]
   Bobs[i+1] ~ dnorm(B[i+1], prec) # Bobs, drawn from a normal
 }
     }
",file="LogisticModel.txt")


# Load initial values
#
jags.inits = function()
  list(B0=100,Binf=500,m=50,prec=100)
#
#
# Specify parameter values to track
#
jags.params = c("Binf","m","B0","B")
#
#
# Run MCMC for model
#
LogisticModel.fit = jags(
  data=jags.data,
  inits=jags.inits,
  jags.params,
  n.iter=1000,
  model.file="LogisticModel.txt"
) 
#
# Examine summary output
#
print(LogisticModel.fit)
plot(LogisticModel.fit)
#
jagsfit.mcmc = as.mcmc(LogisticModel.fit)
xx = as.matrix(jagsfit.mcmc)
head(xx)
#
par(mfrow=c(1,2))
hist(xx[,"Binf"],col="dodgerblue",xlab="Binf",main="Binf")
hist(xx[,"m"],col="dodgerblue",xlab="MSY",main="MSY")
par(mfrow=c(1,1))
#
# Plot fit
#
Time = 1982:1995
plot(Time,arrowtooth$Bobs,
     ylim=c(0,max(arrowtooth$Bobs)),
     xlab="Year",
     ylab="Biomass (Thousands of mt)",
     pch=15)
title("Arrowtooth Flounder Biomass")
lines(1982:1995,predict(Bobs.nls),lwd=3,col=3) # From nls above
lines(1982:1995,B.process,lwd=3,col=2)         # From optim ps above
lines(1982:1995,B.obserr,lwd=3,col=4)          # From optim oe above

#
# Add 95% Bayesian Posterior Intervals
#
BiomassMCMC = xx[,c(1,7,8,9,10,11,12,13,14,2,3,4,5,6)]
head(BiomassMCMC)
PI.lower  = apply(BiomassMCMC,2,quantile,0.025)
PI.median = apply(BiomassMCMC,2,quantile,0.5)
PI.upper  = apply(BiomassMCMC,2,quantile,0.975)
#
lines(Time,PI.lower,lwd=3,lty=2,col=5)
lines(Time,PI.median,lwd=3,lty=1,col=5)
lines(Time,PI.upper,lwd=3,lty=2,col=5)
#
#
legend("bottomright",c("nls","optim ps","optim oe","Bayes"),
       lty=1,lwd=3,col=c(3,2,4,5),bty="n")
#
#
# Explore a Range of Fishing Mortality Rates
# for the System with Parameter Estimates
# Following the Bayes Method
#
Fmort = seq(0,.45,by=0.01)
Blast = rep(0,length(Fmort))
for(k in seq(length(Fmort))){
  B0 = 96.7
  Binf = 512.4
  m = 58.3
  B = rep(B0,100)
  for(i in seq(100)){
    dP = (4*m/Binf)*(1-B[i]/Binf)*B[i]	
    B[i+1] = B[i]+dP-Fmort[k]*B[i]
  }
  Blast[k] = B[100]
}
par(mfrow=c(2,1),mar=c(0, 4, 4, 2))
plot(Fmort,Blast,xlab="Fishing Mortality",ylab="Equilibrium Biomass",
     axes=F,
     main="Equilibrium Levels vs. Fishing Mortality",type="l",lwd=2)
box()
axis(side=2)
par(mar=c(4, 4, 0, 2))
plot(Fmort,Fmort*Blast,xlab="Fishing Mortality",ylab="Equilibrium Yield",
     type="l",lwd=2)
par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)