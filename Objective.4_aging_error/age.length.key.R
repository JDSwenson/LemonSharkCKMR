library(tidyverse)
library(FSA) #Mainly for von bertalanffy growth function (vbFuns)
library(RColorBrewer)

rm(list=ls())

#Set path to data file
input.file <- "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/data/Main_lemon_shark.csv"

#Read in lemon shark data and specify column types with short string representation (c is character and d is double; D is date, but is difficult to work with)
lemon_data <- read_csv(input.file, col_types = "ccdccccddccc")

lemon_data %>% separate(`Tube Label`, sep=";", into = c("Juvenile_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4")) -> lemon_data2 #Separate Tube labels into capture history (there is one tube label per sample instance)

colnames(lemon_data2) <- c("PIT_tag","Capture_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4", "Capture_Year", "Capture_Date", "Island", "Site", "Sex", "PCL_cm", "TL_cm", "DOB", "Father", "Mother") #relabel columns
#head(lemon_data2)

lemon_ref <- lemon_data2[which(lemon_data2$DOB!=""),] #Subset for juveniles with known birth dates

#head(lemon_ref)
#length(lemon_data2[,1]) #How many total bimini lemon sharks in the dataset
#length(lemon_ref[,1]) #How many bimini lemon sharks with known birth dates

Juv_ref <- lemon_ref %>% dplyr::select(Capture_ID,
                                       Capture_Year,
                                       Sex,
                                       PCL_cm,
                                       TL_cm,
                                       DOB)

Juv_ref %>% separate(DOB, sep="/", into=c("DOB", "DOB_Yr2")) -> Juv_ref #Separate uncertain birth years (e.g. 1993/1994) into separate columns

#Filter out individuals to retain just those with one assigned birth year and also length information
Age.length.ref <- Juv_ref %>% dplyr::filter(is.na(DOB_Yr2) == TRUE, 
                          is.na(PCL_cm) != TRUE,
                          is.na(Capture_Year) != TRUE) %>% 
  dplyr::select(!DOB_Yr2) %>% 
  mutate(DOB = as.numeric(DOB)) %>% 
  mutate(age.at.capture = Capture_Year - DOB)


#Look at mean and sd of lengths at different ages
Age.length.ref %>% group_by(age.at.capture) %>% 
  summarize(mean.PCL = mean(PCL_cm, na.rm = TRUE),
            sd.PCL = sd(PCL_cm, na.rm = TRUE),
            number = n())

Age.length.ref %>% ggplot(aes(x = age.at.capture, y = PCL_cm)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x)


min(Age.length.ref$PCL_cm)
max(Age.length.ref$PCL_cm)

#Create bins for length
length_bins <- c(0, 
                 seq(from = round(min(Age.length.ref$PCL_cm), 0),
                to = round(max(Age.length.ref$PCL_cm), 0),
                by = 5))

#Make the last bin big enough to cover all the "+" values
length_bins[length(length_bins)] <- max(Age.length.ref$PCL_cm) + 1

#Create labels for the lengths
length_tags_temp <- NULL
for(i in 1:length(length_bins)){
  length_tags_temp[i] <- paste0(length_bins[i], "-", length_bins[i+1])
}

#Need to remove the last one because there should be one less label than bin
#length_tags_temp[length(length_tags_temp)] <- paste0(max(length_bins), "+")
length_tags <- length_tags_temp[1:length(length_tags_temp) - 1]

#Assign a length bin to each individual
Age.length.ref2 <- Age.length.ref %>% mutate(
  length_bin = cut(Age.length.ref$PCL_cm,
      breaks = length_bins,
      labels = length_tags,
      right = FALSE
      ))

#Create an age-length count matrix, with the number of individuals at each age in each length bin
(age.length.count <- Age.length.ref2 %>% dplyr::count(length_bin, age.at.capture) %>% 
  dplyr::arrange(length_bin, age.at.capture) %>% 
  pivot_wider(names_from = age.at.capture,
              values_from = n,
              names_prefix = "age_") %>% 
  mutate_at(vars(contains("age")), ~replace_na(.,0)) %>% 
  column_to_rownames(var = "length_bin"))

#Convert the above matrix into probability densities
age.length.density <- NULL
for(d in 1:nrow(age.length.count)){
  ald.vec <- round(age.length.count[d,]/rowSums(age.length.count[d,]), 3)
  age.length.density <- rbind(age.length.density, ald.vec)
}

#Age-length probability matrix
age.length.density





#--------------------Lemon VonBert--------------------------------#
#Using parameter values from Brown & Gruber 1988; they used the same VonBert parameterization as the one specified by "typical" in the below function
#Save von bertalanffy growth function; NOTE that this uses PCL, not total length
#Double-checked and the values-at-age/length are the same as the paper.
#Average percent error (APE) was 3.4%
ages <- 0:50
vonBert <- vbFuns(param = "Typical")
tmp <- growthFunShow("vonBertalanffy","Typical")
plot(vonBert(t = ages, Linf = 317.65, K = 0.057, t0 = -2.302), pch = 19, type = "b", main = tmp, ylab = "Length", xlab = "Age")

#Copied from above to see sd of length-at-age
#Look at mean and sd of lengths at different ages
Age.length.ref %>% group_by(age.at.capture) %>% 
  summarize(mean.PCL = mean(PCL_cm, na.rm = TRUE),
            sd.PCL = sd(PCL_cm, na.rm = TRUE),
            number = n())

#Vector of mean length at age from VonBert curve
mean_length.at.age <- vonBert(t = ages, Linf = 317.65, K = 0.057, t0 = -2.302)





####--------------------Pick up here 08/23/2022--------------------------####
#Need to make the sd smaller for older ages -- too many overlapping distributions


#Vector of sd length at age, partially from data partially arbitrary
sd_length.at.age <- c(3, 9, 9, rep(10, times = 48)) #Try 10 so we don't misassign across too many ages

#Combine above into tibble
length.at.age_df <- tibble(mean.1 = mean_length.at.age, 
                           mean.2 = mean_length.at.age, 
                           sd = sd_length.at.age) %>% #Name as mean.1 and mean.2 to allow for join later
  mutate(age.1 = ages,
         age.2 = ages)

#Dataframe of means
laf.tib.mean <- t(combn(length.at.age_df$mean.1, m=2)) %>% 
  as_tibble() %>% 
  dplyr::rename(mean.1 = V1, mean.2 = V2)

#Dataframe of sds
laf.tib.sd <- t(combn(length.at.age_df$sd, m=2)) %>% 
  as_tibble() %>% 
  dplyr::rename(sd.1 = V1, sd.2 = V2)

#Pairwise comparison dataframe of means and standard deviations for each age
(laf.tib <- laf.tib.mean %>% bind_cols(laf.tib.sd) %>% 
  left_join(length.at.age_df, by = "mean.1") %>% 
  dplyr::select(mean.1, mean.2 = mean.2.x, sd.1, sd.2, age.1) %>% 
  left_join(length.at.age_df, by = "mean.2") %>% 
  dplyr::select(mean.1 = mean.1.x, mean.2, sd.1, sd.2, age.1 = age.1.x, age.2))


#----------------Calculate overlap among distributions using Monte Carlo--------------------
#Function from: https://rpsychologist.com/calculating-the-overlap-of-two-normal-distributions-using-monte-carlo-integration

ovl.all_df = ovl.temp <- NULL

for(i in 1:nrow(laf.tib)){
  # Numerical integration using monte carlo methods
  set.seed(456456)
  n <- 100000
  mu1 <- laf.tib$mean.1[i]
  sd1 <- laf.tib$sd.1[i]
  mu2 <- laf.tib$mean.2[i]
  sd2 <- laf.tib$sd.2[i]
  
  xs <- seq(min(mu1 - 3*sd1, mu2 - 3*sd2), max(mu1 + 3*sd1, mu2 + 3*sd2), length.out=n)
  f1 <- dnorm(xs, mean=mu1, sd=sd1) # dist1
  f2 <- dnorm(xs, mean=mu2, sd=sd2) # dist2
  
  ps <- matrix(c(runif(n, min(xs), max(xs)), runif(n, min=0, max=max(f1,f2)) ), ncol=2) # sample x,y from uniform dist
  
  z1<- ps[,2] <= dnorm(ps[,1], mu1, sd1) # dist1
  z2<- ps[,2] <= dnorm(ps[,1], mu2, sd2) # dist 2
  z12 <- z1 | z2 # both dists
  z3 <- ps[,2] <= pmin(dnorm(ps[,1], mu1, sd1), dnorm(ps[,1], mu2, sd2)) # overlap
  
  o <- (sum(z3)/sum(z1) + sum(z3)/sum(z2))/2
  ovl.temp <- tibble(age.1 = laf.tib$age.1[i], age.2 = laf.tib$age.2[i], overlap = o)
  ovl.all_df <- bind_rows(ovl.all_df, ovl.temp)
print(paste0("Finished iteration ", i))
}

ovl.all_df

#Make reformatted dataframe where each age is in age.1, regardless of whether larger or smaller than age.2
overlap.ages <- ovl.df %>% dplyr::filter(overlap > 0.001) #Subset for distributions with overlap

age.dist.denom = denom.temp = ovl.switch.temp = ovl.switch <- NULL
for(j in ages){
  
  #Calculate total overlap of distributions with one another to use as denominator for setting probability of age misassignment. Start with probability of getting it wrong on the lower end.
  # ovl.above <- overlap.ages %>% dplyr::filter(age.1 == j) %>%
  #   dplyr::filter(row_number() != 1) %>% 
  #   summarize(ovl.over = sum(overlap)) %>% 
  #   pull(ovl.over)
  # 
  # #What's the probability of misassigning to be too high?  
  # ovl.below <-  overlap.ages %>% dplyr::filter(age.2 == j) %>%
  #   dplyr::filter(row_number() != n()) %>% 
  #   summarize(ovl.under = sum(overlap)) %>% 
  #   pull(ovl.under)
  # 
  # denom.temp <- tibble(age = j, denominator = (ovl.above + ovl.below + 1))
  # 
  # age.dist.denom <- bind_rows(age.dist.denom, denom.temp)
  
  #Make new dataframe where the ages from age.2 are flipped, and add it to our original dataframe of probability overlap
  ovl.switch.temp <- overlap.ages %>% dplyr::filter(age.2 == j) %>% 
    dplyr::rename(age.1 = age.2, age.2 = age.1)

  ovl.switch <- bind_rows(ovl.switch, ovl.switch.temp)
  
  }

#Make new dataframe where age.1 refers to the age being misassigned, and not just the younger individual
overlap.ages_corr <- overlap.ages %>% bind_rows(ovl.switch) %>% 
  dplyr::arrange(age.1, age.2)


#--------------Calculate age misassignment probability------------------#
denom_df <- overlap.ages_corr %>% group_by(age.1) %>% 
  summarize(denominator = n())

age.miss_df <- overlap.ages_corr %>%  
  left_join(denom_df, by = "age.1") %>% 
  mutate(prob.ovlp = overlap/denominator)

age.miss_df %>% group_by(age.1) %>% 
  summarize(sum(prob.ovlp)) %>% 
  View()


#NEXT: check into whether the values make sense
laf.df = laf.temp <- NULL

#Draws from a Normal distribution
for(i in 1:nrow(length.at.age_df)){
  
  laf.temp <- tibble(age = rep(i-1, times = 1000)) %>% mutate(length = rnorm(n = 1000, mean = length.at.age_df$mean[i], sd = length.at.age_df$sd[i]))
  
  laf.df <- bind_rows(laf.df, laf.temp)
}

laf.df

#Check which ages to add
overlap.ages_corr %>% dplyr::filter(age.1 == 40)

laf.df %>% dplyr::filter(age >= 20) %>% 
  mutate(age = factor(age)) %>% 
  ggplot(aes(x = length, colour = age)) +
  geom_density()












#--------------------------DON'T NEED ANYTHING BELOW HERE (I think)----------------------------------#
#Trying to add values progressively to make up a denominator ... prob more complicated than it should be
ovlp.temp = ovl.add.temp = ovl.add <- NULL

#Start For loop
for(k in 2:(length(ages)-1)){
  
  a <- ages[k]
  
  min.ovl <- overlap.ages_corr %>% dplyr::filter(age.1 == a) %>%
    slice_min(age.2) %>% 
    pull(age.2)

  max.ovl <- overlap.ages_corr %>% dplyr::filter(age.1 == a) %>%
    slice_max(age.2) %>% 
    pull(age.2)

for(m in seq(from = min.ovl, to = a-1, by = 1)){
  for(n in seq(from = a+1, to = max.ovl, by = 1)){
    
    ovlp <- overlap.ages_corr %>% dplyr::filter(age.1 == m, age.2 == n)
    ovlp.temp <- bind_rows(ovlp.temp, ovlp)
    
     }
}

ovl.add.temp <- tibble(age = a, add.ovl = sum(ovlp.temp$overlap))
ovl.add <- bind_rows(ovl.add, ovl.add.temp)

ovlp.temp <- NULL

print(paste0("Finished with age ", a))
  }

#Relevant dataframes
ovl.add
ovl.all_df
overlap.ages_corr

denom.final_df <- age.dist.denom %>% left_join(ovl.add, by = "age") %>% 
  replace_na(list(add.ovl = 0)) %>% 
  mutate(total_denom = denominator + add.ovl)

age.miss_df <- overlap.ages_corr %>% dplyr::rename(age = age.1) %>% 
  left_join(denom.final_df, by = "age") %>% 
  dplyr::select(age.1 = age, age.2, overlap, total_denom) %>% 
  mutate(prob.ovlp = overlap/total_denom)


age.miss_df %>% group_by(age.1) %>% 
  summarize(sum(prob.ovlp)) %>% 
  View()






laf.df = laf.temp <- NULL

#Draws from a Normal distribution
for(i in 1:nrow(length.at.age_df)){
  
  laf.temp <- tibble(age = rep(i-1, times = 1000)) %>% mutate(length = rnorm(n = 1000, mean = length.at.age_df$mean[i], sd = length.at.age_df$sd[i]))

  laf.df <- bind_rows(laf.df, laf.temp)
}

laf.df


#visualize
n <- 51
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


laf.df %>% mutate(age = factor(age)) %>%  
  gghistogram(x = "length", color = "age", fill = "age", palette = col_vector, add_density = TRUE) +
  theme(legend.position = "none")

 






#### Age Slicing function
#Use the below function and reference this page: https://stecf.jrc.ec.europa.eu/c/document_library/get_file?uuid=0efd7160-4555-41f5-8234-6e74db26e632&groupId=43805

ageIt=function(len=NULL,n=NULL,Linf=NULL,K=NULL,t0=0.0,timing=0.5,plusGroup=30){
  ## expected age at length adjust to beginning of year
  age=pmax(pmin(floor(t0-log(1-pmin(len/Linf,.9999999))/K+timing),plusGroup),0)
  
  ## calculate frequencies
  res=aggregate(n, list(age=age), sum)
  
  return(res)}                  





#------------------------From Liz--------------------------------#
#Multiplicative lets it change more with age
Linf <- 100
K <- 0.25
t0 <- -0.05

ages <- seq(1,19)
len <- seq(1,100)
# ---
vonbert <-function(x,Linf,K,t0=0){
  
  
  y=Linf*(1.0-exp(-K*(x-t0)))
  return(y)
  
}
# ---

mean.len.age<-vonbert(ages, Linf, K, t0)   #vector of length at age
sd.len.add <- 18
# additive
age.len.mat.add <- sapply(mean.len.age, function(x) qnorm(p=seq(0.05, 0.95, by=0.05), mean=x, sd=sd.len.add))

# multiplicative
sd.len.mult <- 0.2
#Each cell is length; column is age; each row is percentile, going by 5%
age.len.mat.mult <- exp(sapply(log(mean.len.age), function(x) qnorm(p=seq(0.05, 0.95, by=0.05), mean=x, sd=sd.len.mult)) )

#I can loop through the dataframe length.at.age_df and make the matrix of quantiles from Liz
#I'd have the true age for each individual; I'm looking to get a matrix where the diagonal is the highest value and the row sums to 1
#I'd take the true age, then sample from a vector of misassigned ages and pass it the probabilities; it will be flatter with higher ages
#The exact matrix Liz sent is probably early on in the process

int_f <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

m.a1 <- 30
sd.a1 <- 5
m.a2 <- 45
sd.a2 <- 5
m.a3 <- 50
sd.a3 <- 5

integrate(int_f, -Inf, Inf, mu1=m.a1, mu2=m.a2, sd1=sd.a1, sd2=sd.a3)