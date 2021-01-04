############################################################
        ####  Matrix population models  ####
############################################################
#Two useful websites: 
#1) https://kevintshoemaker.github.io/NRES-470/LECTURE7.html
#2) https://rstudio-pubs-static.s3.amazonaws.com/100473_88bd83e70c434c67a3c03efa8e6b04f7.html

#NEXT TO DO (9/24/2020):
#Adapt fishSim CKMR script to estimate constant abundance
#Work through fishSim markdown - add all relevant chunks, including model code, import, and vizualization

rm(list=ls())
library(popbio)

# Syntax for projecting abundance using a transition matrix (NOTE: this code won't run until we specify the terms on the right)
#Year1 <- projection_matrix %*% Abundance_year0  # matrix multiplication!


####### First, build a simple transition matrix filled with 0s, one column and one row per stage

#----------------Initialize matrix w/ plus group------------------------
s_YOY = 0.6    # YOY survival
s_adult = 0.85  # Adult (mature) survival 
#Change length.out to ramp up from s_YOY to s_adult, and then edit times so the total number of survival rates is 11
s_juv = c(round(seq(s_YOY, s_adult, length.out = 3), 2), rep(s_adult, times = 8))    # Juvenile survival - ramp up from YOY survival to adult survival; with 13 columns corresponding to years 0 (YOY) - 12 (mature adults), we set 11 different survival rates, starting with the YOY survival and moving to adult. The for loop in the function below starts at index 2, so doesn't repeat YOY survival
a_m = 12         # Age at sexual maturity
b_max = 1.5     #Females recruited into the population. Assume equal sex ratio. Calculated from: average fecundity (6) divided by two (for skipped-breeding) divided by two (for females vs males)
age_x = a_m + 1 # The "plus-group" stage (aka self-loop group) 

init_matrix = function(s_j, s_a, a_mat, age_xx, b_maxx){
  A_matrix = matrix(data = 0, nrow = age_xx, ncol = age_xx) # dimension matrix and fill with zeros
  A_matrix[1, (a_mat:age_xx)] = b_maxx # assign birth rates to mature ages in first row of matrix
  A_matrix[2,1] <- s_YOY
  for(ii in 2:(a_mat-1)) A_matrix[ii+1, ii] = s_j[ii] # assign juvenile survival rates
  A_matrix[age_xx,a_mat] = s_a # adult survival assumed for maturing animals transitioning into plus-group
  A_matrix[age_xx, age_xx] = s_a # adult survival assumed for plus-group
  return(A_matrix)
}

A = init_matrix(s_j = s_juv, s_a = s_adult, a_mat = a_m, 
                age_xx = age_x, b_maxx = b_max) # Call function to initialize matrix A
A

#-----------------Initialize matrix w/ all ages -----------------------
s_YOY = 0.6    # YOY survival
s_adult = 0.85  # Adult (mature) survival 
#Change length.out to ramp up from s_YOY to s_adult, and then edit times so the total number of survival rates is 11
s_juv = c(round(seq(s_YOY, s_adult, length.out = 3), 2), rep(s_adult, times = 8))    # Juvenile survival - ramp up from YOY survival to adult survival; with 13 columns corresponding to years 0 (YOY) - 12 (mature adults), we set 11 different survival rates, starting with the YOY survival and moving to adult. The for loop in the function below starts at index 2, so doesn't repeat YOY survival
a_m = 12         # Age at sexual maturity
b_max = 1.5     #Females recruited into the population. Assume equal sex ratio. Calculated from: average fecundity (6) divided by two (for skipped-breeding) divided by two (for females vs males)
age_x <- max_age


init_matrix = function(s_j, s_a, a_mat, age_xx, b_maxx){
  A_matrix = matrix(data = 0, nrow = age_xx, ncol = age_xx) # dimension matrix and fill with zeros
  A_matrix[1, (a_mat:age_xx)] = b_maxx # assign birth rates to mature ages in first row of matrix
  A_matrix[2,1] <- s_YOY
  for(ii in 2:(a_mat-1)) A_matrix[ii+1, ii] = s_j[ii] # assign juvenile survival rates
  #A_matrix[age_xx,a_mat] = s_a # adult survival assumed for maturing animals transitioning into plus-group
  for(jj in a_mat:(age_xx-1)) A_matrix[jj+1, jj] = s_adult # assign juvenile
  A_matrix[age_xx, age_xx] = s_a # adult survival assumed for plus-group
  return(A_matrix)
}

A = init_matrix(s_j = s_juv, s_a = s_adult, a_mat = a_m, 
                age_xx = age_x, b_maxx = b_max) # Call function to initialize matrix A
View(A)

#Calculate dominant eigenvalue (i.e. population growth rate) from transition matrix
lambda(A)

#Calculate stable age structure of the population
stable_age = stable.stage(A)
stable_age
sum(stable_age) #should equal 1




####### Then we specify initial abundances for the three age classes
init_pop_size <- 1000
aa1 <- round(stable_age[12]*init_pop_size,0)
aa2 <- round(stable_age[13]*init_pop_size,0)
Abundance_year0 <- c(rep(0,times=11), aa1, aa2)
Abundance_year0


######
# Now we can run the code for real

Year1 <- A %*% Abundance_year0  # matrix multiplication!
Year1


######### Project another year
Year2 <- A %*% Year1  # matrix multiplication!
Year2


##########
# Use a FOR loop to project many years into the future

nYears <- 100                 # set the number of years to project
TMat <- A    # define the projection matrix
InitAbund <- Abundance_year0 # define the initial abundance

## NOTE: the code below can be re-used without modification:
allYears <- matrix(0,nrow=nrow(TMat),ncol=nYears+1)     # build a storage array for all abundances!
allYears[,1] <- InitAbund  # set the year 0 abundance                                    
for(t in 2:(nYears+1)){   # loop through all years
  allYears[,t] <-  TMat %*% allYears[,t-1]
}

allYears[,nYears]




#Need to work on this plot ... 9/17/20
plot(1,1,pch="",ylim=c(0,max(allYears)),xlim=c(0,nYears+1),xlab="Years",ylab="Abundance",xaxt="n")  # set up blank plot

cols <- rainbow(ncol(TMat))    # set up colors to use

for(s in 1:ncol(TMat)){
  points(allYears[s,],col=cols[s],type="l",lwd=2)     # plot out each life stage abundance, one at a time
}

axis(1,at=seq(1,nYears+1),labels = seq(0,nYears))   # label the axis
legend("top",col=cols,lwd=rep(2,ncol(TMat)),legend=paste("Age ",seq(0:ncol(TMat)-1)))  # put a legend on the plot