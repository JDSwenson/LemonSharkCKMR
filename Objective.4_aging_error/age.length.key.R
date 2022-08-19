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

#Vector of sd length at age, partially from data partially arbitrary
sd_length.at.age <- c(3, 9, 9, rep(10, times = 48)) #Try 10 so we don't misassign across too many ages

#Combine above into tibble
length.at.age_df <- tibble(mean = mean_length.at.age, sd = sd_length.at.age)

laf.df = laf.temp <- NULL

for(i in 1:nrow(length.at.age_df)){
  
  laf.temp <- tibble(age = rep(i-1, times = 1000)) %>% mutate(length = rnorm(n = 1000, mean = length.at.age_df$mean[i], sd = length.at.age_df$sd[i]))

  laf.df <- bind_rows(laf.df, laf.temp)
}

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