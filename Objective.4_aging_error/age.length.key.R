#To quickly navigate:
#Search for DELETE for sections that can likely be deleted
#Search for VARIABLES for variables that can probably be moved to top of script
#Search for Troubleshoot for areas that may have issues
#rm(list=ls())
library(tidyverse)
library(FSA) #Mainly for von bertalanffy growth function (vbFuns)

rm(list=ls())

#Set path to data file
input.file <- "G://My Drive/Personal_Drive/R/CKMR/Objective.5_lemon_shark_data/data/Main_lemon_shark.csv"

#Save von bertalanffy growth function
vonBert <- vbFuns(param = "Typical")

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

Age.length.ref %>% dplyr::filter(age.at.capture == 7)

Age.length.ref %>% group_by(age.at.capture) %>% 
  summarize(mean.PCL = mean(PCL_cm, na.rm = TRUE),
            sd.PCL = sd(PCL_cm, na.rm = TRUE),
            number = n())

Age.length.ref %>% ggplot(aes(x = age.at.capture, y = PCL_cm)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y~x)
