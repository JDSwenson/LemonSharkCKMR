## Diagnostics
library(coda)
rm(list=ls())

#-----------------------------Trace plots---------------------------------------------------
#Specify parameters to plot
jags_params_4plot <- c("Nf", "Nm", "surv") #Specify parameters

#Specify save location for pdf of plots
tracePlot.file <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Diagnostic.plots/TracePlots_", today, "_400samples.pdf")
pdf(file = tracePlot.file)

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(jags.model.400)){
  diag_plots(jags.model.400[[i]], jags_params_4plot, layout = "4x1")
}

dev.off()

#-----------------------------Autocorrelation---------------------------------------------------
#Specify save location for pdf of plots
autocorr.file <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Diagnostic.plots/Autocorrelation.plots_", today, "_400samples_longChain.pdf")
pdf(file = autocorr.file)

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(jags.model.400)){
  autocorr.plot(jags.model.400[[i]])
}

dev.off()

#Visual check for one iteration
autocorr.diag(jags.model.400[[5]], lags = c(0, 1, 5, 10, 15, 20))

#-----------------------------Cross-correlation---------------------------------------------------
#Specify save location for pdf of plots
crosscorr.file <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Diagnostic.plots/Cross_correlation.plots_", today, "_400samples_longChain.pdf")
pdf(file = crosscorr.file)

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(c in 1:length(jags.model.400)){
  crosscorr.plot(jags.model.400[[c]])
}

dev.off()

#Visual check for one iteration
crosscorr(jags.model.400[[5]])

#-----------------------------Effective chain length---------------------------------------------------
#Initialize dataframes
effectiveSize.df = EF.temp <- NULL
for(j in 1:length(jags.model.400)){
  EF.temp <- data.frame(t(effectiveSize(jags.model.400[[j]]))) %>% 
    mutate(iteration = j)
  effectiveSize.df <- rbind(effectiveSize.df, EF.temp)
}
effectiveSize.df

#-----------------------------Geweke---------------------------------------------------
# Tests the means of the first and last parts of the chain and compares their z scores. We expect that the means will be the same; otherwise, our model may not be converging until a later sample is drawn.
# the test statistic is the standard z score for the equality of two means - value greater than 1.64 (p<0.1) or 1.96 (p<0.05) two-tailed
# We do NOT want to reject the null here, as that means the beginning and end of the chains is very different

geweke.df = geweke.temp.c1 = geweke.temp.c2 <- NULL
g.thresh <- 1.65 # Set threshold for geweke failure

# Loop through all chains and iterations and store geweke diagnostic results
for(w in 1:length(jags.model.400)){
  geweke.temp.c1 <- data.frame(t(geweke.diag(jags.model.400[[w]])[[1]][[1]])) %>% 
    mutate(iteration = w, chain = 1)
  geweke.temp.c2 <- data.frame(t(geweke.diag(jags.model.400[[w]])[[2]][[1]])) %>% 
    mutate(iteration = w, chain = 2)
  geweke.df <- rbind(geweke.df, geweke.temp.c1, geweke.temp.c2)
}

# Identify index for each instance where the geweke diagnostic failed
Nf.Gew.vec <- which(abs(geweke.df$Nf) > g.thresh)
Nm.Gew.vec <- which(abs(geweke.df$Nm) > g.thresh)
surv.Gew.vec <- which(abs(geweke.df$surv) > g.thresh)
total.Gew.vec <- unique(c(Nf.Gew.vec, Nm.Gew.vec, surv.Gew.vec))

#Save iteration number of each failure so can analyze results
geweke.failed = G.temp <- NULL
for(f in 1:length(total.Gew.vec)){
  g.index <- total.Gew.vec[f]
  G.temp <- geweke.df$iteration[g.index]
  geweke.failed <- c(geweke.failed, G.temp)
}
geweke.failed

#Summary
#N Females
geweke.df %>% summarize(`0.10` = sum(abs(Nf) >1.64))
geweke.df %>% summarize(`0.05` = sum(abs(Nf) >1.96))

#N Males
geweke.df %>% summarize(`0.10` = sum(abs(Nm) >1.64))
geweke.df %>% summarize(`0.05` = sum(abs(Nm) > 1.96))

#Survival
geweke.df %>% summarize(`0.10` = sum(abs(surv) > 1.64))
geweke.df %>% summarize(`0.05` = sum(abs(surv) > 1.96))

#Save plots of geweke diagnostic
#Specify save location for pdf of plots
geweke.file <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Diagnostic.plots/Geweke_", today, "_400samples_longChain.pdf")
pdf(file = geweke.file) #Open pdf file for plotting

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(jags.model.400)){
  geweke.plot(jags.model.400[[i]])
}

dev.off() #Close pdf file

#-----------------------------Gelman & Rubin---------------------------------------------------
# The Gelman diagnostic calculates the potential scale reduction factor (PSRF) for each variable. The PSRF estimates a factor by which the scale of the distribution might be reduced if the simulations were run for an infinite number of iterations. As the number of iterations approaches infinity, the PSRF should decline to 1.  ... it's kind of like an ANOVA, where it compares the within-chain and between-chain variance. This is essentially the same as the rhat metric.

#Calculate gelman diagnostic for each iteration
gelman.df = gelman.temp <- NULL

for(g in 1:length(jags.model.400)){

  gelman.temp <- data.frame(t(gelman.diag(jags.model.400[[g]])[[1]])) %>% 
    rownames_to_column(var = "type") %>% 
    mutate(iteration = g)
  
  gelman.df <- rbind(gelman.df, gelman.temp)
  }

#Save plots of gelman diagnostic
#Specify save location for pdf of plots
gelman.file <- paste0("G://My Drive/Personal_Drive/R/CKMR/Model.validation/Diagnostic.plots/Gelman_", today, "_400samples_longChain.pdf")
pdf(file = gelman.file) #Open pdf file for plotting

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(jags.model.400)){
  gelman.plot(jags.model.400[[i]])
}

dev.off() #Close pdf file

#-----------------------------Heidelberger & Welch---------------------------------------------------
#Heidelberger and Welch’s convergence diagnostic
# The HW diagnostic examines whether the draws from the posterior come from a stationary distribution. 
# It iteratively removes proportions of the samples, and reports the iteration at which we should start the chain i.e. increase the burn-in period by the iteration reported here.
# We do NOT want to reject the null.

heidel.diag(jags.model.400[[1]])


#-----------------------------Raftery-Lewis---------------------------------------------------
# Raftery-Lewis test
# Will determine the appropriate additional burn-in and thinning rate
# M = number of iterations that should be discarded
# N = total number of iterations that should be run for each variable
# Nmin = the minimum number of iterations that should be run for each variable
# I = the increase in number of iterations needed to reach convergence.
raftery.diag(jags.model.400[[1]])


superdiag(jags.model.400[[1]])
##===============================================================================

#Charlotte's code
# Nyears <- 16
# N.saved <- length(jags.model$sims.list$deviance)/2
# xx <- 1:N.saved
# 
# labs <- c("final year abundance", "long-term trend", "deviance")
# 
# outs <- array(NA, c(N.saved, 2, 3)) #Dimensions are row, column, matrix; here, 
# for(i in 1:N.saved)
# {
#   outs[i,1,1] <- jags.model$sims.list$N[i,Nyears]
#   outs[i,2,1] <- jags.model$sims.list$N[N.saved+i,Nyears]
#   
#   outs[i,1,2] <- jags.model$sims.list$u[i]
#   outs[i,2,2] <- jags.model$sims.list$u[N.saved+i]
#   
#   outs[i,1,3] <- jags.model$sims.list$deviance[i]
#   outs[i,2,3] <- jags.model$sims.list$deviance[N.saved+i]
# }
# 
# par(mfrow=c(2,2))
# for(k in 1:3) 
# {
#   yy <- outs[,1,k]
#   plot(xx, yy, xlab="cycle number", ylab="", main=labs[k], type='b', pch=16, ylim=range(outs[,,k]))
#   yy <- outs[,2,k]
#   lines(xx, yy, type='b',pch=16, col="gray50")
# }
# 