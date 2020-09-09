### Practicing Matrix Population Modeling
#### A.K. Teffer

setwd("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/")

#Load packaages and set directory
library(expm)    
library(ggplot2)
library(tidyverse)
library(reshape2)
library(plotly)

## Define parameters
p = 0.5  # prevalence or other pathogen metric, 0-1 (current "optimal" pathogen metric = 0.775, moderately prevalent agent but affect pop at at very high prev)
t = 0.00  # absolute value of degrees away from optimal. Current optimal temp anomaly = 0
p <- seq(0.0, 1.0, 0.2)  # pathogen prevalence or likelihood of presence or virulence degree - j as each level
t <- seq(0.0, 2.0, 0.2)  # temperature anomaly - k as each level

### Fecundity within and between populations (a,b) 
Fa <- 0.7
Fb <- 0.7
Fab <- 0.7
Fba <- 0.7

### Survival by population/stage (0, 1+)
Sa0 <- function(p,t){0.8 - p * 0.7 - t * 0.02}
Sa1 <- function(p,t){0.8 - p * 0.7 - t * 0.02}
Sb0 <- function(p,t){0.8 - p * 0.7 - t * 0.02}
Sb1 <- function(p,t){0.8 - p * 0.7 - t * 0.02}

### Movement by population/stage
Mab0 <- function(p,t){0.5 - p * 0.4 - t * 0.05}
Mab1 <- function(p,t){0.7 - p * 0.4 - t * 0.05}
Mba0 <- function(p,t){0.5 - p * 0.4 - t * 0.05}
Mba1 <- function(p,t){0.7 - p * 0.4 - t * 0.05}

## Create our matrix (take the blue pill)
mtx <- function(p,t){
  out <- matrix(
        c(0, Sa0(p,t), 0, Mab0(p,t),
        Fa, Sa1(p,t), Fab, Mab1(p,t),
        0, Mba0(p,t), 0, Sb0(p,t),
        Fba, Mba1(p,t), Fb, Sb1(p,t)),
      nrow=4)
}

mtx1 <- mtx(p,t) #test it
mtx1

## Characterize population attributes (r, stable stage structure)
eigen(mtx1)
eigen(mtx(p,t))$values[1]  #extract value

### Create empty matrix to store values per year for 100 years
#### Initial population structure: 20,50,40,30
#### Loop over 100 years per parameter

Q <- 100  # years
Z <- 0  # starting row for building matrix
W <- length(p) * length(t)*Q  # total rows in output matrix is number of values in P and T

#### Create matrix to hold loop output - name columns
dat <- tibble(j=NA, k=NA, i=NA, eig = NA, mat1 = NA, mat2 = NA, mat3 = NA, mat4 = NA, .rows=W)

for (j in p) {                            # for each level of the pathogen metric set
  for (k in t) {                          # for each level of temperature anomaly set
    mat <- matrix(NA, nrow=4, ncol=100)   # 4x100 matrix of population structure per year
    mat[,1] <- c(20,50,40,30)             # initial population structure
    for(i in 2:Q){                        # beginning with 2nd year up to 100
      mat[,i] <- mtx(j,k) %*% mat[,i-1]   # apply mtx fxn to the previous col in mat (pop str) and fill in current col with output
      Z = Z+1                             # Increase Z by 1 (next step for outter matrix, e.g. 0+1)
      dat$j[Z] <- j                       # put current "j" in the Z row of dat in j col
      dat$k[Z] <- k                       # put current "k" in the Z row of dat in k col
      dat$i[Z] <- i                       # put current "i" in the Z row of dat in i col
      dat$eig[Z] <- eigen(mtx(j,k))$values[1]
      dat$mat1[Z] <- mat[1,i]             # put the new pop str for iteration i in the final 4 col
      dat$mat2[Z] <- mat[2,i]
      dat$mat3[Z] <- mat[3,i]
      dat$mat4[Z] <- mat[4,i]
    } 
  }
}

### r ~ 0.98-1.02 for sweet spot of population stability
### Plot heatmap to identify optimal parameter values
dat2<-data.frame(dat)
names(dat2)
ggplot(dat2) +
  geom_tile(aes(j, k, fill=eig))

### Eaxamine relationship of eigenvalue (r) with pathogen and temperature metrics
ggplot(dat2) +
  geom_point(aes(j, eig))
ggplot(dat2) +
  geom_point(aes(k, eig))
plot_ly(dat2, x=~k, y=~j, z=~eig, type="scatter3d", mode="markers", color=~eig)








### Bonus plots - ignore
### Plot age distributions in 2 populations over time
pop.100 <- as.data.frame(mat)
rownames(pop.100) <- c("PopA0", "PopA1", "PopB0", "PopB1")
colnames(pop.100) <- c(1:100)

### Add row names and melt into long form to plot
pop.df <- tibble::rownames_to_column(pop.100, "pop")
pop.df2 <- melt(data = pop.df,
                id.vars = "pop")

### Plot population size change over 100 cycles
ggplot(pop.df2, aes(x=variable, y=value, fill=pop)) + 
  geom_bar(stat = "identity") +
  labs(y="Population size", x="Cycle") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Next - identify parameter values that give r ~ 0.98-1.02 - heatmap to find sweet spot
  
  

