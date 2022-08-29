library(tidyverse)
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

#Initial N
N = 1000 #Has to be something but doesn't really matter what
mat = 3

#Values to run through
psis = seq(0,1,by = 0.05)
phis = seq(0.1,0.9,by = 0.2)
as = seq(5)

#Create data frame to fill
combos = expand.grid(psi = psis, phi = phis, a = as)
combos$Navail.AllAnn = NA
combos$Navail.AllSkip = NA
#combos$Navail.AllPsiRatio = NA
combos$Navail.PsiRatioAnns =NA
combos$Navail.PsiRatioSkips =  NA
combos$Navail.PsiRatioSkipsAll = NA

#Run through all possible combos
for(i in seq(dim(combos)[1]))
{
  #Pull values for this run
  phi = combos$phi[i]
  psi = combos$psi[i]
  a = combos$a[i]
  
  #Calculate the stable age distribution
  maxage = which(phi^seq(1000) <= 0.001)[1] #set the maximum age for calculating init age props to where <= 1% of inds in cohort are left
  if(is.na(maxage)){next} #If phi = 0 | 1 this obvi doesn't work
  age.props = (phi^seq(0,maxage))/sum((phi^seq(0,maxage)))
  if(length(age.props) < mat){next} #Ditto if phi is too low for the set maturity
  
  #The number of spawners that would be available if all were annual is the number that are mature
  Navail.AllAnn = sum(N * age.props[mat:length(age.props)])
  #The number of spawners that would be available if all were skips is the number that are mature AND on-cycle
  Navail.AllSkip = sum(N * age.props[seq(from = mat, to = length(age.props), by = a)])
  #The number of annual spawners that are available for a given psi is the number of mature indivs alive for a given year times (1-psi)
  Navail.PsiRatioAnns = sum(N * age.props[mat:length(age.props)] * (1 - psi))
  #The number of skip spawners that are available for a given psi is the number of mature indivs alive for a given year times psi 
  Navail.PsiRatioSkips =  sum(N * age.props[seq(from = mat, to = length(age.props), by = a)] * psi)
  #The total number of spawners for a given psi is the sum of the annuals and skips
  Navail.PsiRatioSkipsAll = Navail.PsiRatioAnns + Navail.PsiRatioSkips
  
  #Save everything
  combos$Navail.AllAnn[i] = Navail.AllAnn
  combos$Navail.AllSkip[i] = Navail.AllSkip
  #combos$Navail.AllPsiRatio[i] = Navail.PsiRatio
  combos$Navail.PsiRatioAnns[i] =Navail.PsiRatioAnns
  combos$Navail.PsiRatioSkips[i] =  Navail.PsiRatioSkips
  combos$Navail.PsiRatioSkipsAll[i] = Navail.PsiRatioSkipsAll
  
}

#Get the relative proportions of the total spawning pool for boths annuals and skips
combos$SkipProp = combos$Navail.PsiRatioSkips / combos$Navail.PsiRatioSkipsAll
combos$AnnProp = combos$Navail.PsiRatioAnns / combos$Navail.PsiRatioSkipsAll

#Add in mathematical relationships to show that they're equivalent to seq()-based approach
combos$Nsim.skips = (N* combos$phi^(mat-1) * ((1 - combos$phi) / (1 - combos$phi^combos$a))) * combos$psi
combos$Nsim.anns = (N* combos$phi^(mat-1) * (1- combos$psi))

#For two indivs with an age difference equal to a (the lowest on-cycle difference allowed since we're not doing same cohort comparisons),
#the prob of a shared parent is the total number of spawners available reduced by the prob that their shared parent died during that period (phi^a)
combos$HS.Prob.OnCycle.aY = combos$phi^combos$a / combos$Navail.PsiRatioSkipsAll

#For two indivs with an age difference equal to 1 (the lowest off-cycle difference allowed),
#the prob of a shared parent is the number of annual spawners each year reduced by the chance that their shared parent from that pool died (phi^1)
#and further reduced by the prob that the older indiv's parent belonged to the skip-spawning contingent 
#Which is a function of the proportion of annual spawners each year relative to the total spawning pool
combos$HS.Prob.OffCycle.1Y = ((combos$phi * combos$AnnProp) / combos$Navail.PsiRatioAnns)

#Remove off-cycle comparisons for a=1 because they're impossible and therefore will never be encountered
combos[which(combos$a == 1),]$HS.Prob.OffCycle.1Y = NA


#Add labeller functions for faceting
alabeller = function(x){return(paste0("a = ",as.character(x)))}
philabeller = function(x){return(paste0("phi = ",as.character(x)))}

#Plot the proportions of each parental type as a function of a, phi, and psi
plot.Nprops = ggplot(combos)+
  geom_line(aes(x = psi, y = SkipProp), col = '#619CFF', size = 1)+
  geom_line(aes(x = psi, y = AnnProp), col = '#F8766D', size = 1)+
  geom_point(data = data.frame(x = rep(0.5,2), y = rep(0.5,2), type = c("Annual", "Skip")), #Trick the plot into showing a legend
             aes(x = x, y = y, col = type), alpha = 0)+
  scale_color_manual("Spawner type",values = c("Annual" = '#F8766D', "Skip" = '#619CFF'))+
  scale_y_continuous("Proportion of adult population")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  facet_grid(rows = vars(a), cols = vars(phi),
             labeller = labeller(a = alabeller, phi = philabeller))+
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8),
        legend.position = 'bottom')+
  ggtitle("Proportional of each spawner type in annual breeding pool\nby phi, psi, and a")

plot.Nprops


#Plot the relative total breeding pool size as a function of a, phi, and psi
plot.NtotRel = ggplot(combos)+
  geom_line(aes(x = psi, y = Navail.PsiRatioSkipsAll), col = 1, size = 1)+
  scale_y_continuous("Relative adult N (all spawner types")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  facet_grid(rows = vars(a), cols = vars(phi),
             labeller = labeller(a = alabeller, phi = philabeller))+
  theme(strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom')+
  ggtitle("Relative size of adult breeding pool\nby phi, psi, and a")

plot.NtotRel

#Plot the relative probabilities of on- and off-cycle MHSPs (or PHSPs, just one-sided) by a, phi, and psi
#use lowest possible age difference while avoiding same-cohort comparisons
plot.Prel = ggplot(combos)+
  geom_line(aes(x = psi, y = HS.Prob.OnCycle.aY), col = '#39B600', size = 1)+
  geom_line(aes(x = psi, y = combos$HS.Prob.OffCycle.1Y), col = '#C77CFF', size = 1)+
  geom_point(data = data.frame(x = rep(0.0,2), y = rep(0.0,2), type = c("On-Cycle", "Off-Cycle")), #Trick the plot into showing a legend
             aes(x = x, y = y, col = type), alpha = 0)+
  scale_color_manual("Comparison type",
                     values = c("On-Cycle" = '#39B600', "Off-Cycle" = '#C77CFF'),
                     labels = c("On-Cycle" = "On-Cycle\n(agediff = a)\n", "Off-Cycle" = 'On-Cycle\n(agediff = 1)\n'))+
  scale_y_continuous("Relative probabilities of being a MHSP")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  facet_grid(rows = vars(a), cols = vars(phi),
             labeller = labeller(a = alabeller, phi = philabeller))+
  theme(strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom')+
  ggtitle("Relative probabilities of being a MHSP\nby phi, psi, and a")
#Note that left-right difference is merely effect of phi on total parental pop size

plot.Prel

#Can also be plotted as a scaled value (i.e., removing pop. size effect)
combos.test = combos %>% 
  group_by(phi, a) %>%
  mutate(HS.Prob.OnCycle.aY.scaled = HS.Prob.OnCycle.aY / max(c(HS.Prob.OnCycle.aY,HS.Prob.OffCycle.1Y), na.rm = T),
         HS.Prob.OffCycle.1Y.scaled = HS.Prob.OffCycle.1Y / max(c(HS.Prob.OnCycle.aY,HS.Prob.OffCycle.1Y), na.rm = T))

plot.PrelScaled = ggplot(combos.test)+
  geom_line(aes(x = psi, y = HS.Prob.OnCycle.aY.scaled), col = '#39B600', size = 1)+
  geom_line(aes(x = psi, y = HS.Prob.OffCycle.1Y.scaled), col = '#C77CFF', size = 1)+
  geom_point(data = data.frame(x = rep(0.0,2), y = rep(0.0,2), type = c("On-Cycle", "Off-Cycle")), #Trick the plot into showing a legend
             aes(x = x, y = y, col = type), alpha = 0)+
  scale_color_manual("Comparison type",
                     values = c("On-Cycle" = '#39B600', "Off-Cycle" = '#C77CFF'),
                     labels = c("On-Cycle" = "On-Cycle\n(agediff = a)\n", "Off-Cycle" = 'On-Cycle\n(agediff = 1)\n'))+
  scale_y_continuous("Relative probabilities of being a MHSP")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  facet_grid(rows = vars(a), cols = vars(phi),
             labeller = labeller(a = alabeller, phi = philabeller))+
  theme(strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom')+
  ggtitle("Relative probabilities of being a MHSP\nby phi, psi, and a\n(scaled to remove effect of phi on total adult N)")

plot.PrelScaled

# plotdir = "D:\\Google Drive\\Research\\Postdoc\\Meetings\\CKMRgrp_082422\\"
# daydate = Sys.Date() #today's date
# daydate = paste0(substr(daydate,6,7), substr(daydate,9,10), substr(daydate,3,4)) #Chop it all up to get the format right
# 
# png(paste0(plotdir, "\\", "CKMR_SkipSpawn_Sims_Nprops", "_", daydate, ".png"), #Open a .png instead of plot window
#     width = 10.25, height = 7.5, units = "in", res = 600)
# plot.Nprops
# dev.off()
# 
# png(paste0(plotdir, "\\", "CKMR_SkipSpawn_Sims_NtotRel", "_", daydate, ".png"), #Open a .png instead of plot window
#     width = 10.25, height = 7.5, units = "in", res = 600)
# plot.NtotRel
# dev.off()
# 
# png(paste0(plotdir, "\\", "CKMR_SkipSpawn_Sims_Prel", "_", daydate, ".png"), #Open a .png instead of plot window
#     width = 10.25, height = 7.5, units = "in", res = 600)
# plot.Pabs
# dev.off()
# 
# png(paste0(plotdir, "\\", "CKMR_SkipSpawn_Sims_PrelScaled", "_", daydate, ".png"), #Open a .png instead of plot window
#     width = 10.25, height = 7.5, units = "in", res = 600)
# plot.PrelScaled
# dev.off()

