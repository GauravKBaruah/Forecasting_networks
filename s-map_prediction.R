rm(list=ls())
source("../src/01_functions_cluster.R")
library(statmod)
require(deSolve)                    # for integrating ordinary differential equations
require(tidyr)                      # for efficient data manipulation & plotting
require(tidyverse)
library(cowplot)                    # for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)

#------------------------------------------------------------------------------------------------------------------------

#     generation of time series

#------------------------------------------------------------------------------------------------------------------------

# change if necessary
network_dir = 'small_networks'
#-------------------------------------------------------------------------------
webfiles = list.files(path = network_dir, pattern = "*.csv", full.names = TRUE)

fact <- expand.grid(Temperature=20,
                    `web` = webfiles,
                    network_size = 0,
                    connectance = 0,
                    nestedness = 0,
                    h2 = c("none", "rand"),
                    rate = c(0.05, 0.2),
                    env_var = c(1, 7)                              
) %>% arrange_all()

for (r in 1:nrow(fact)){
  g <- adj.mat(webfiles[which(webfiles == fact$web[r])]) 
  Aspecies <- dim(g)[2]
  Pspecies <- dim(g)[1]
  
  fact$network_size[r] <- Aspecies+Pspecies
  fact$nestedness[r] <- round(nestedness_NODF(g),2)
  fact$connectance[r] <- Connectance(g)
}

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
  
  tmax <- 3e3                                                                    # time to integrate equations for
  dt <- 0.1
  initial_temperature <- fact$Temperature[1]
  
  env_var <- fact$env_var[r]                                                    # temperature fluctuation range 
  rate <- fact$rate[r]                                                          # temperature increase
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,web=fact$web[r],
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,A=Aspecies,P=Pspecies,degree=degree, env_var=env_var,
               var="high",gi=gi,bw=bw,dt=dt,rate=rate,initial_temperature=initial_temperature,
               mut.strength=mut.strength,nestedness=nestedness,C=C,network_size=network_size)
  
  ic <- c(na,np,muA,muP)                                                        # initial conditions coerced into a vector
  
  out <- eqs_euler_2(time = tmax, y = ic, pars = params)                        # generation of time series
  
  saveRDS(out, paste(networkString, "_", fact$h2[r], "_", rate, "_", env_var, "_s-map_231130.rds", sep = ""))
  
  #------------     Quick testing    -------------------------------------------
  #
  # ts.plot(out$Na)
  # ts.plot(out$Np)
  # ts.plot(out$Temp)
  #
  #-----------------------------------------------------------------------------
  
})

# output files, contain time series data
network_output = list.files(path = "network_output", pattern = regex("M_PL_\\d+(_\\d+)?_(rand)?(none)?_0.\\d+_\\d_s-map_231130.rds"), full.names = T) 

#------------------------------------------------------------------------------------------------------------------------

#     s-map prediction on time series

#------------------------------------------------------------------------------------------------------------------------

ts <- 1                                                                         # time_series index, adjust accordingly
ts_data_raw <- readRDS(network_output[ts])
ts_data <- cbind(ts_data_raw$Na, ts_data_raw$Np)

colnames <- c("A1")

for(i in 2:ncol(ts_data_raw$Na)) {
  colnames <- c(colnames, paste("A", i, sep = ""))
}

for(i in 1:ncol(ts_data_raw$Np)) {
  colnames <- c(colnames, paste("P", i, sep = ""))
}

colnames(ts_data) <- colnames

targ_col <- 8                                                                   # insert column of species, animals before plants in ts_data

Embedding <- colnames                                                           # change 'colnames' to specific columns if necessary
Edim <- length(Embedding)

ts_data <- ts_data[, Embedding]

coeff_names <- sapply(colnames(ts_data), function(x) paste("d", 
                                                     colnames(ts_data)[targ_col], "d", x, sep = ""))

block <- cbind(ts_data[2:dim(ts_data)[1], targ_col], ts_data[1:(dim(ts_data)[1] - 1),])
norm_consts <- apply(block, 2, function(x) sd(x))
block <- as.data.frame(apply(block, 2, function(x) (x - mean(x)) / sd(x)))

lib <- 1:dim(block)[1]
pred <- 1:dim(block)[1]
theta <- 8

coeff <- array(0, dim = c(length(pred), Edim))
colnames(coeff) <- coeff_names
coeff <- as.data.frame(coeff)

lm_svdsolve <- function(y, x, ws, subset = seq_along(y)) {
  x <- x[subset, ]
  y <- y[subset]
  ws <- ws[subset]
  
  A <- cbind(1, x) * ws
  A_svd <- svd(A)
  s <- A_svd$d
  s_inv <- matrix(0, nrow = dim(x)[2] + 1, ncol = dim(x)[2] + 1)
  for(i in seq_along(s)) {
    if(s[i] >= max(s) * 1e-5)
      s_inv[i,i] <- 1/s[i]
  }
  
  coeff <- A_svd$v %*% s_inv %*% t(A_svd$u) %*% (ws * y)
  coeff <- t(coeff)
  colnames(coeff) <- c("const", colnames(x))
  return (coeff)
}

for (ipred in 1:length(pred)) {
  libs = lib[-pred[ipred]]
  q <- matrix(as.numeric(block[pred[ipred], 2:dim(block)[2]]),
              ncol = Edim, nrow = length(libs), byrow = TRUE)
  distances <- sqrt(rowSums((block[libs, 2:dim(block)[2]] - q)^2))
  dbar <- mean(distances)
  Ws <- exp(-theta*distances/dbar)
  svd_fit <- lm_svdsolve(block[libs, 1], block[libs, 2:dim(block)[2]], Ws)
  coeff[ipred, ] <- svd_fit[-1]
}

coeff <- cbind(pred, coeff)
colnames(coeff)[1] <- "t"

for (i in 1:length(coeff_names)) {
  cd = str_extract(coeff_names[i], regex("(?<=d)(A|P)(?=\\dd)"))                # clade of dependent species
  ci = str_extract(coeff_names[i], regex("(?<=d(A|P)\\dd)(A|P)"))               # clade of independent species
  
  snum_d = str_extract(coeff_names[i], regex("(?<=d(A|P))\\d(?=d)"))            # which dependent species
  snum_i = str_extract(coeff_names[i], regex("(?<=d(A|P)\\dd(A|P))\\d"))        # which independent species
  
  col_comp = "#cb1249"
  col_mut = "#48ab49"
  
  par(mfrow = c(1,1))
  coefflims = c(-0.5, 0.5)
  trange <- 1:1000
  plot(coeff[trange,"t"],coeff[trange,i + 1],type="l",col=if(cd == ci) col_comp else col_mut ,xlab="time",
       ylab=bquote(partialdiff*.(cd)[.(snum_d)] / partialdiff*.(ci)[.(snum_i)]),
       ylim=coefflims,xlim=range(trange),lwd=2)
  abline(a=0 ,b=0 , lty="dashed", col="black", lwd=.5)
}

assay <- "Some Assay"
ylab <-bquote(.(cd) ~ AC50 ~ (mu*M))
plot(0, xlab = xlab)

par(mfrow = c(1,1))
coefflims = c(-0.5, 0.5)
trange <- 1:1000
plot(coeff[trange,"t"],coeff[trange,8],type="l",col="blue",xlab="time",
     ylab=expression(partialdiff*A[2] / partialdiff*P[1]),
     ylim=coefflims,xlim=range(trange),lwd=2)
abline(a=0 ,b=0 , lty="dashed", col="black", lwd=.5)

#------------------------------------------------------------------------------------------------------------------------

#     plotting of time series

#------------------------------------------------------------------------------------------------------------------------

myColors = c("#FF0000", "#00FF00", "#0000FF",
             "#FFFF00", "#FF00FF", "#00FFFF",
             "#333333", "#777777", "#BFFFBB",
             "#FFCCCC", "#CCFFCC", "#CCCCFF",
             "#FFFFCC", "#FFCCFF", "#CCFFFF",
             "plum", "orange", "wheat4", "brown") 

for (i in 1:length(network_output)) {
  out <- readRDS(network_output[i])
  network_string <- str_extract(network_output[i], regex("M_PL_\\d+(_\\d+)?"))   
  h2 <- str_extract(network_output[i], regex("(rand)|(none)"))
  rate <- str_extract(network_output[i], regex("0\\.\\d+"))
  env_var <- str_extract(network_output[i], regex("(?<=_)(1|7)(?=_)"))
  
  
  ts.plot(out$Na, gpars = list(xlab = paste("increasing", network_string, "high", h2, rate, env_var, sep = "_"), ylab = "abundance, animals", col = myColors)) 
  ts.plot(out$Np, gpars = list(xlab = paste("increasing", network_string, "high", h2, rate, env_var, sep = "_"), ylab = "abundance, plants", col = myColors))
  ts.plot(out$Temp, gpars = list(xlab = paste("increasing", network_string, "high", h2, rate, env_var, sep = "_"), ylab = "temperature", col = myColors)) 
}


#------------------------------------------------------------------------------------------------------------------------

#     miscellaneous

#------------------------------------------------------------------------------------------------------------------------


#str_view_all('M_PL_061_34_none_0.05_7_s-map_231130', regex("M_PL_\\d+(_\\d+)?_(rand)?(none)?_0.\\d+_\\d_s-map_231130.rds")) # _(rand | none)"))

             