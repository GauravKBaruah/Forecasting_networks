rm(list=ls())
source("01_functions_cluster.R")
library(statmod)
require(deSolve)    ## for integrating ordinary differential equations
require(stringr)
require(tidyr)      ## for efficient data manipulation & plotting
library(cowplot)    ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)

# change if necessary
mydir =   '../networks_sizes_8_34'    
myfiles = list.files(path = mydir, pattern = "*.csv", full.names = TRUE)
webfiles <- myfiles

fact_2 <- expand.grid(`web` = webfiles) %>% 
  as_tibble %>%
  mutate(`Nestedness`= 0,
         `Connectance` = 0,
         `network_size`= 0)

# a loop to go over the rows of the dataframe to calculate nestedness and network size for each network
for (r in 1:nrow(fact_2)){
  g <- adj.mat(myfiles[which(myfiles == fact_2$web[r])]) #network web names
  Aspecies <- dim(g)[2] # no of animal species
  Pspecies <- dim(g)[1] # no of animal species
  
  fact_2$network_size[r] <- Aspecies+Pspecies
  fact_2$Nestedness[r] <- round(nestedness_NODF(g),2)
  fact_2$Connectance[r] <- Connectance(g)
}

webfiles <- as.character(fact_2$web)

fact <- expand.grid(Temperature=20,
                  `web` = webfiles,
                  `replicates`=1+(1:1)*100) 


for(r in 1:nrow(fact)){
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[r])])  

  Aspecies <- dim(g)[2]        
  Pspecies <- dim(g)[1]
  
  degree.plants<-degree.animals<-numeric()
  for(i in 1:Pspecies){
    degree.plants[i]<-sum(g[i,])
  } 
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j])
  }
  
  degree <- c(degree.animals, degree.plants)
  nestedness <- nestedness_NODF(g)
  c <- Connectance(g)
  network_size <- Aspecies + Pspecies
  
  #parameters for modelling: intialisation of the model

  sigma<-runif((Aspecies+Pspecies), 0.05,0.09) 

  na <- rep(1, Aspecies) 
  np <- rep(1, Pspecies) 
  muA <- runif(Aspecies, 16, 24)   #initial mean phenotypic optimum trait values
  muP <- runif(Pspecies, 16, 24)   #initial mean phenotypic optimum trait values
  
  h2 <- runif((Pspecies + Aspecies), 0.1, 0.25)
  #parameters below are taken from another paper Akesson et al 2021 Nat Comms.
  bw  <- 2
  aw  <- 0
  gi <- 1
  ki <- 0.1 # mortality rate
  w <- 7    # mutualism interaction width
  Temp <- fact$Temperature[1]  
  
  Amatrix <- mat.comp(g)$Amatrix  # competition matrix , aij, for animals
  Pmatrix <- mat.comp(g)$Pmatrix  # competition matrix , aij, for plants
  mut.strength = 1.5 # average mutualistic strength
  
  dt <- 0.1
  initial_temperature <- fact$Temperature[1]
  env_var <- 2 # temperature fluctuation range 
  rate <- 0 # rate of temperature increase --> constant temperature
  
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,web=fact$web[r],
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,A=Aspecies,P=Pspecies,degree=degree, env_var=env_var,
               var="high",gi=gi,bw=bw,dt=dt,rate=rate,initial_temperature=initial_temperature,
               mut.strength=mut.strength,nestedness=nestedness,C=C,network_size=network_size)
  
  ic <- c(na,np,muA,muP) # initial conditions coerced into a vector
  
  tmax <- 1e4 # time to integrate equations for
  
  out <- eqs_euler(time = tmax, y = ic, pars = params)
  
  networkString = str_extract(fact$web[r], regex("M_PL_\\d+(_\\d+)?"))
  saveRDS(out, paste(networkString, "_constant_high_h2_010_025.rds", sep = ""))

}


#-------------------------------------------------------------------------------
#
#                  plotting time series
#
#-------------------------------------------------------------------------------

out <- readRDS("M_PL_061_33_constant_high_h2_010_025.rds")                      # time series to plot 

#------------------------------------------------------------
#------------------------------------------------------------
#                  color palettes                                               
#------------------------------------------------------------
#------------------------------------------------------------

diverseColors = c("#FF0000", "#00FF00", "#0000FF",
                  "#FFFF00", "#FF00FF", "#00FFFF",
                  "#333333", "#777777", "#BFFFBB",
                  "#FFCCCC", "#CCFFCC", "#CCCCFF",
                  "#FFFFCC", "#FFCCFF", "#CCFFFF",
                  "plum", "orange", "wheat4", "brown")

animalColors = c("#FF0000", "#FF8B00", "#BA5430", "#FFD100", "#FF5D00",
                 "#FF7400", "#AF1400", "#FFA200", "#FFB900", "#FFFF80",
                 "#FFE800", "#FFFF00", "#FFFF2A", "#FA4680", "#FFFFD5")

plantColors = c("#00A600", "#2DB600", "#63C600", "#A0D600", "#E6E600",
                "#E8C32E", "#EBB25E", "#EDB48E", "#F0C9C0", "#52F252")

ts.plot(cbind(out$Na[1:10000,], out$Np[1:10000,]), gpars = list(xlab = "time", ylab = "abundance", col = c(animalColors[1:6], plantColors[1:2])), lwd = 2)

ts.plot(cbind(out$muA, out$muP), gpars = list(xlab = "time", ylab = "species optimum temperature", col = c(animalColors[1:6], plantColors[1:2])), lwd = 2)

