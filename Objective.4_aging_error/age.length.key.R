library(tidyverse)
library(FSA) #Mainly for von bertalanffy growth function (vbFuns)

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
