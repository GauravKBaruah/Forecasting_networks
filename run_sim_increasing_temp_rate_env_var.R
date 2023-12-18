rm(list=ls())
source("01_functions_cluster.R")
library(statmod)
require(deSolve)                    # for integrating ordinary differential equations
require(tidyr)                      # for efficient data manipulation & plotting
require(tidyverse)
library(cowplot)                    # for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(wesanderson)

# change if necessary
network_dir = '../networks_sizes_8_34'
#-------------------------------------------------------------------------------
webfiles = list.files(path = network_dir, pattern = "*.csv", full.names = TRUE)

fact <- expand.grid(Temperature=20,
                    `web` = webfiles,
                    network_size = 0,
                    connectance = 0,
                    nestedness = 0,
                    h2 = c("none", "rand"),
                    rate = seq(0.075, 0.2, 0.025),
                    env_var = seq(0, 7.5, 1.5),                                 
                    `replicates`=1+(1:10)*100
                    ) %>% arrange_all()

for (r in 1:nrow(fact)){
  g <- adj.mat(webfiles[which(webfiles == fact$web[r])]) 
  Aspecies <- dim(g)[2]
  Pspecies <- dim(g)[1]
  
  fact$network_size[r] <- Aspecies+Pspecies
  fact$nestedness[r] <- round(nestedness_NODF(g),2)
  fact$connectance[r] <- Connectance(g)
}

# for storage of extinction time points of every species in respective network
extinction_tp_list <- list(rep(0, (fact$network_size)[1]))

for (i in 1:10) {
  for(j in 1:720) {
    network_size <- (fact$network_size)[i*720]
    extinction_tp_list[[(i-1)*720 + j]] <- rep(0, network_size)
  }
}

fact <- fact %>% mutate(extinction_tps = extinction_tp_list)

# simulation of time series
system.time(for(r in 1:nrow(fact)){
  print(r)
  g<-adj.mat(webfiles[which(webfiles == fact$web[r])])  
  
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
  
  # read in data from pre-run with constant temperature
  networkString = str_extract(fact$web[r], regex("M_PL_\\d+(_\\d+)?")) 
  nw_data_base <- readRDS(paste(networkString, "_constant_high_h2_010_025.rds", sep = ""))
  time_points_base <- length(nw_data_base$Temp)
  animal_abundances_base_final <- nw_data_base$Na[time_points_base,]
  plant_abundances_base_final <- nw_data_base$Np[time_points_base,]
  animal_traits_base_final <- nw_data_base$muA[time_points_base,]
  plant_traits_base_final <- nw_data_base$muP[time_points_base,]
  
  na <- animal_abundances_base_final 
  np <- plant_abundances_base_final  
  muA <- animal_traits_base_final    
  muP <- plant_traits_base_final     
  
  if (fact$h2[r] == "none") {
    h2 <- rep(0, Aspecies + Pspecies)
  } else {
      h2 <- runif(Aspecies + Pspecies, 0.2, 0.4)
  }
  
  #parameters below are taken from another paper Akesson et al 2021 Nat Comms.
  bw  <- 2
  aw  <- 0
  gi <- 1
  ki <- 0.1                                                                     # mortality rate
  w <- 7                                                                        # mutualism interaction width
  Temp <- fact$Temperature[1]                                               
                                                
  Amatrix <- mat.comp(g)$Amatrix                                                # competition matrix , aij, for animals
  Pmatrix <- mat.comp(g)$Pmatrix                                                # competition matrix , aij, for plants
  mut.strength = 1.5                                                            # average mutualistic strength
                                                
  tmax <- 5e3                                                                   # time to integrate equations for
  dt <- 0.1
  initial_temperature <- fact$Temperature[1]
  
  env_var <- fact$env_var[r]                                                    # temperature fluctuation range 
  rate <- fact$rate[r]                                                          # temperature increase
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,web=fact$web[r],
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,A=Aspecies,P=Pspecies,degree=degree, env_var=env_var,
               var="high",gi=gi,bw=bw,dt=dt,rate=rate,initial_temperature=initial_temperature,
               mut.strength=mut.strength,nestedness=nestedness,C=C,network_size=network_size)
  
  ic <- c(na,np,muA,muP)                                                        # initial conditions coerced into a vector
  
  out <- eqs_euler(time = tmax, y = ic, pars = params)                        # generation of time series
  fact <- calculate_extinction_tps(out, network_size, fact, r)                  
  
  #------------     Quick testing    -------------------------------------------
  #
  # ts.plot(out$Na)
  # ts.plot(out$Np)
  # ts.plot(out$Temp)
  #
  #-----------------------------------------------------------------------------
  
})

# adjust date 
saveRDS(fact, "increasing_fluctuating_231204.rds")

#-------------------------------------------------------------------------------
#
#                  plotting time series
#
#-------------------------------------------------------------------------------

out <- readRDS(file = "M_PL_061_33_rand_0.075_1.5_101_time_series_231205.rds" ) # time series to plot 

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

ts.plot(cbind(out$Na[0:900,], out$Np[0:900, ]), gpars = list(xlab = "time", ylab = "abundance", col = c(animalColors[1:6], plantColors[1:2])), lwd = 2)

