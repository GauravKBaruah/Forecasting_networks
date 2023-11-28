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
                    rate = seq(0.05, 0.2, 0.025),
                    env_var = c(0.25, 2.5, 5),                                 
                    `replicates`=1+(1:3)*100
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
  for(j in 1:126) {
    network_size <- (fact$network_size)[i*126]
    extinction_tp_list[[(i-1)*126 + j]] <- rep(0, network_size)
  }
}

fact <- fact %>% mutate(extinction_tps = extinction_tp_list)

# simulation of time series
system.time(for(r in 1:nrow(fact)){
  
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
      h2 <- runif(Aspecies + Pspecies, 0.25, 0.4)
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
                                                
  tmax <- 5e3                                                                    # time to integrate equations for
  dt <- 0.1
  initial_temperature <- fact$Temperature[1]
  
  env_var <- fact$env_var[r]                                                    # temperature fluctuation range 
  rate <- fact$rate[r]                                                          # temperature increase
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,web=fact$web[r],
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,A=Aspecies,P=Pspecies,degree=degree, env_var=env_var,
               var="high",gi=gi,bw=bw,dt=dt,rate=rate,initial_temperature=initial_temperature,
               mut.strength=mut.strength,nestedness=nestedness,C=C,network_size=network_size)
  
  ic <- c(na,np,muA,muP)                                                        # initial conditions coerced into a vector
  
  out <- eqs_euler_2(time = tmax, y = ic, pars = params)                          # generation of time series
  fact <- calculate_extinction_tps(out, network_size, fact, r)                  
  
  #------------     Quick testing    -------------------------------------------
  #
  # ts.plot(out$Na)
  # ts.plot(out$Np)
  # ts.plot(out$Temp)
  #
  #-----------------------------------------------------------------------------
  
})

saveRDS(fact, "increasing_fluctuating_231128.rds")

#-------------------------------------------------------------------------------

#     plotting

#-------------------------------------------------------------------------------

myColors = c("#FF0000", "#00FF00", "#0000FF",
             "#FFFF00", "#FF00FF", "#00FFFF",
             "#333333", "#777777", "#BFFFBB",
             "#FFCCCC", "#CCFFCC", "#CCCCFF",
             "#FFFFCC", "#FFCCFF", "#CCFFFF",
             "plum", "orange", "wheat4", "brown") 


ts.plot(out$Na, gpars = list(xlab = paste("increasing", networkString, "high", fact$h2[r],  "rate: ", rate, "env_var: ", env_var, sep = ", "), ylab = "abundance, animals", col = myColors)) 
ts.plot(out$Np, gpars = list(xlab = paste("increasing", networkString, "high", fact$h2[r],  "rate: ", rate, "env_var: ", env_var, sep = ", "), ylab = "abundance, plants", col = myColors))
ts.plot(out$Temp, gpars = list(xlab = paste("increasing", networkString, "high", fact$h2[r],  "rate: ", rate, "env_var: ", env_var, sep = ", "), ylab = "abundance, animals", col = myColors)) 


#------------------------------------------------------------------------------
#
#         miscellaneous
#
#------------------------------------------------------------------------------

# system.time()

# ggplot(data = as_tibble(out$Na)) 
#   geom_line(mapping = aes(x = out$Temp, y = , color = factor(h2)))
# 
# ggplot(data = dat_01) +
#   geom_point(mapping = aes(x = temperature, y = plant_richness, color = factor(h2)))
# 
# ggplot(data = dat_01) +
#   geom_point(mapping = aes(x = temperature, y = animal_richness, color = factor(h2)))

# plot.ts(out$Na, ylab = "abundance, animals")
# plot.ts(out$Np, ylab = "abundance, plants")
# plot.ts(out$muA, ylab = "trait, animals")
# plot.ts(out$muP, ylab = "trait, plants")
# plot.ts(out$Temp, ylab = "Temperature[°C]")

#------------------------------------------------------------------------

# ggplot(networkData$Na, aes(x = Time, y = )) + 
#   geom_line(aes(color = variable), size = 1) +
#   scale_color_manual(values = myColors) +
#   theme_minimal()

#save(out, file ="../output/warming_increasing_temp.RData")

# https://www.climate.gov/news-features/understanding-climate/climate-change-global-temperature
#--> 5/9 * 6 projected increase in temperature between 2021 and 101 

#------------------------------------------------------------------------

# fact_2<- expand.grid(`web` =webfiles[1:153]) %>%
#   as_tibble %>%
#   mutate(`Nestedness`= 0,
#          `Connectance` = 0,
#          `network_size`=0)
# # a loop to go over the rows of the dataframe to calculate nestedness and network size for each network
# for (r in 1:nrow(fact_2)){
#   g<-adj.mat(myfiles[which(myfiles == fact_2$web[r])]) #network web names
#   Aspecies<- dim(g)[2] # no of animal species
#   Plantspecies<- dim(g)[1]
#   
#   fact_2$network_size[r]<-Aspecies+Plantspecies
#   fact_2$Nestedness[r]<-  round(nestedness_NODF(g),2)
#   fact_2$Connectance[r] <- Connectance(g)
# }- 

# a loop to go over the rows of the dataframe to calculate nestedness and network size for each network
# fact_2 <- expand.grid(`web` = webfiles) %>% #153]) %>% rate, variance
#   as_tibble %>%
#   mutate(`Nestedness`= 0,
#          `Connectance` = 0,
#          `network_size`= 0)

# ts.plot(networkData$muA, gpars = list(xlab = paste("increasing", network, var_string, h2_string, sep = ", "), ylab = "trait, animals", col = myColors))
# ts.plot(networkData$muP, gpars = list(xlab = paste("increasing", network, var_string, h2_string, sep = ", "), ylab = "trait, plants", col = myColors))
# ts.plot(networkData$Temp, gpars = list(xlab = paste("increasing", network, var_string, h2_string, sep = ", "), ylab = "temperature", col = myColors))

