rm(list=ls())
source("../src/01_functions_cluster.R")
source("function_predictions.R")
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
network_dir = 'smallestAndLargest/' #'small_networks'
#-------------------------------------------------------------------------------
webfiles = list.files(path = network_dir, pattern = "*.csv", full.names = TRUE)

fact <- expand.grid(Temperature=20,
                    `web` = webfiles,
                    network_size = 0,
                    connectance = 0,
                    nestedness = 0,
                    h2 = c("none", "rand"),
                    rate = c(0.075, 0.2),                                       # lowest rate --> look at effect of different fluctuations
                    env_var = seq(0, 7.5, 1.5)                              
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

for(i in 1:48) {
    network_size <- (fact$network_size)[i]
    extinction_tp_list[[i]] <- rep(0, network_size)
}


fact <- fact %>% mutate(extinction_tps = extinction_tp_list)

system.time(for(r in 43:nrow(fact)){
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
  network_string = str_extract(fact$web[r], regex("M_PL_\\d+(_\\d+)?")) 
  nw_data_base <- readRDS(paste(network_string, "_constant_high_h2_010_025.rds", sep = ""))
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
  fact <- calculate_extinction_tps(out, network_size, fact, r)                  
  
  saveRDS(out, paste(network_string, "_", fact$h2[r], "_", rate, "_", env_var, "_time_series_231204.rds", sep = ""))
  #------------     Quick testing    -------------------------------------------
  #
  # ts.plot(out$Na)
  # ts.plot(out$Np)
  # ts.plot(out$Temp)
  #
  #-----------------------------------------------------------------------------
})

saveRDS(fact, "network_output/231204_prediction/fact_predictions_231204.rds")

#------------------------------------------------------------------------------------------------------------------------

#     preparation

#------------------------------------------------------------------------------------------------------------------------

networks_time_series = list.files(path = "network_output/231204_prediction/", pattern = regex("M_PL_\\d+(_\\d+)?_(rand)?(none)?_0\\.\\d+_\\d(\\.\\d+)?_time_series_231204.rds"), full.names = T) 

fact <- readRDS("network_output/231204_prediction/fact_predictions_231204.rds")

fact <- fact %>% mutate(
  spearman_correlation_degree = 0,
  spearman_correlation_rate = 0,
  spearman_correlation_eigen = 0
  
)

#------------------------------------------------------------------------------------------------------------------------

#     degree, prediction of extinction sensitivity

#------------------------------------------------------------------------------------------------------------------------

large_60_11 <- adj.mat("smallestAndLargest/M_PL_060_11.csv")
small_61_33 <- adj.mat("smallestAndLargest/M_PL_061_33.csv")  

degrees_large_60_11 <- c(colSums(large_60_11), rowSums(large_60_11))
degrees_small_61_33 <- c(colSums(small_61_33), rowSums(small_61_33))

# ---------------------- large network, degree ranks ---------------------------

ranks_large_60_11_final <- rank(degrees_large_60_11)
ranks_large_60_11_temp <- unique(sort(ranks_large_60_11_final))
ranks_large_60_11_replace <- list(list(0,0))

for(i in 1:length(ranks_large_60_11_final)) {
  ranks_large_60_11_replace[[i]] <- list(ranks_large_60_11_temp[i], i)
}

for(i in 1:length(degrees_large_60_11)) {
  for(j in 1:length(ranks_large_60_11_replace)) {
    if(ranks_large_60_11_final[i] == ranks_large_60_11_replace[[j]][[1]]) {
      ranks_large_60_11_final[i] <- ranks_large_60_11_replace[[j]][[2]]
      break
    }
  }
}

# ---------------------- small network, degree ranks ---------------------------

ranks_small_61_33_final <- rank(degrees_small_61_33)
ranks_small_61_33_temp <- unique(sort(ranks_small_61_33_final))
ranks_small_61_33_replace <- list(list(0,0))

for(i in 1:length(ranks_small_61_33_final)) {
  ranks_small_61_33_replace[[i]] <- list(ranks_small_61_33_temp[i], i)
}

for(i in 1:length(degrees_small_61_33)) {
  for(j in 1:length(ranks_small_61_33_replace)) {
    if(ranks_small_61_33_final[i] == ranks_small_61_33_replace[[j]][[1]]) {
      ranks_small_61_33_final[i] <- ranks_small_61_33_replace[[j]][[2]]
      break
    }
  }
}
# ------------------------------------------------------------------------------
# --------------------------------- predictions --------------------------------
# ------------------------------------------------------------------------------


# ---------------------- large network, spearman -------------------------------

for(ts in 1:nrow(fact) / 2) {
  ranks_extinctions <- ranking_empirical(fact, ts)
  fact$spearman_correlation_degree[ts] <- cor(ranks_extinctions, ranks_large_60_11_final, method = "spearman")
}

# ---------------------- small network, spearman -------------------------------

for(ts in ((nrow(fact)/2) + 1):nrow(fact)) {
  ranks_extinctions <- ranking_empirical(fact, ts)
  fact$spearman_correlation_degree[ts] <- cor(ranks_extinctions, ranks_small_61_33_final, method = "spearman")
}

#------------------------------------------------------------------------------------------------------------------------

#     rate of change, prediction of extinction sensitivity

#------------------------------------------------------------------------------------------------------------------------

for (ts in 1:length(networks_time_series)) {
  ts_data_raw <- readRDS(networks_time_series[ts])

  ranks_extinctions <- ranking_empirical(fact, ts)
  
  rates_of_change_animals <- matrix(0, nrow = nrow(ts_data_raw$Na), ncol = ncol(ts_data_raw$Na))  
  rates_of_change_plants <- matrix(0, nrow = nrow(ts_data_raw$Np), ncol = ncol(ts_data_raw$Np))  
  
  for(j in 2:nrow(ts_data_raw$Na)) {
    rates_of_change_animals[j,] <- abs((ts_data_raw$Na)[j,] - (ts_data_raw$Na)[j-1,]) / (ts_data_raw$Na)[j-1,]
    rates_of_change_plants[j,] <- abs((ts_data_raw$Np)[j,] - (ts_data_raw$Np)[j-1,]) / (ts_data_raw$Np)[j-1,]
  }
  
  rates_of_change <- cbind(rates_of_change_animals, rates_of_change_plants)
  
  tp_begin_collapse <- min(fact$extinction_tps[[ts]])
  
  average_rates_of_change <- rep(0, ncol(rates_of_change))
  
  for(j in 1:length(average_rates_of_change)) {
    average_rates_of_change[j] <- sum(rates_of_change[(tp_begin_collapse - 80):(tp_begin_collapse-1),j]) / 80
  }
  
  fact$spearman_correlation_rate[ts] <- cor(ranks_extinctions, rank(-average_rates_of_change), method = "spearman")
}


#------------------------------------------------------------------------------------------------------------------------

#     s-map, on time series

#------------------------------------------------------------------------------------------------------------------------

for (ts in 1:length(networks_time_series)) { #ts <- 1                                 # time_series index, adjust accordingly
  print(ts)
  
  ranks_extinctions <- ranking_empirical(fact, ts)
  
  ts_data_raw <- readRDS(networks_time_series[ts])
  network_identifier = str_extract(networks_time_series[ts], regex("M_PL_\\d+(_\\d+)?_(rand)?(none)?_0.\\d+_\\d")) 
  
  ts_data <- cbind(ts_data_raw$Na, ts_data_raw$Np)
  
  tp_begin_collapse <- min(fact$extinction_tps[[ts]])
  
  colnames <- c("A1")
  
  for(i in 2:ncol(ts_data_raw$Na)) {
    colnames <- c(colnames, paste("A", i, sep = ""))
  }
  
  for(i in 1:ncol(ts_data_raw$Np)) {
    colnames <- c(colnames, paste("P", i, sep = ""))
  }
  
  colnames(ts_data) <- colnames
  
  jacobians <- list(matrix(0, length(colnames), length(colnames)))                # to store jacobians for every timepoint
  
  # for(tp in 2:dim(ts_data_raw$Na)[1] - 1) {                                       # prepare list of jacobians for storage, dim(ts_data_raw$Na)[1] - 1 different timepoints
  #   jacobians[[tp]] <- matrix(0, length(colnames), length(colnames))
  # }
  
  for(tp in 2:80) {
    jacobians[[tp]] <- matrix(0, length(colnames), length(colnames))
  }
  
  for (targ_col in 1:length(colnames)) {                                          # targ_col =  column of dependent species, animals before plants in ts_data
    print(targ_col)
    Embedding <- colnames                                                         # change 'colnames' to specific columns if necessary
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
    
    for (ipred in (tp_begin_collapse - 80):(tp_begin_collapse - 1)) { # 1:length(pred)) {
      libs = lib[-pred[ipred]]
      q <- matrix(as.numeric(block[pred[ipred], 2:dim(block)[2]]),
                  ncol = Edim, nrow = length(libs), byrow = TRUE)
      distances <- sqrt(rowSums((block[libs, 2:dim(block)[2]] - q)^2))
      dbar <- mean(distances)
      Ws <- exp(-theta*distances/dbar)
      svd_fit <- lm_svdsolve(block[libs, 1], block[libs, 2:dim(block)[2]], Ws)
      jacobians[[ipred - (tp_begin_collapse - 80) + 1]][targ_col,] <- coeff[ipred, ] <- svd_fit[-1]              # calculate partial derivatives and store in jacobian matrix for time point 'ipred'
    }                                                                                                            # partial derivatives of targ_col with respect to other cols in coeff[ipred, ]
    
    average_eigen_sensitivities <- rep(0, ncol(ts_data))
    
    for (i in 1:length(average_eigen_sensitivities)) {
      sum_eigen_entries <- 0
      for (j in 1:80) {
        sum_eigen_entries <- sum_eigen_entries + abs(Re(eigen(jacobians[[j]])$vectors[,1]))[i]
      }
      average_eigen_sensitivities[i] <- sum_eigen_entries / 80
    }
    
    fact$spearman_correlation_eigen[ts] <- cor(ranks_extinctions, rank(-average_eigen_sensitivities), method = "spearman")
    
    # coeff <- cbind(pred, coeff)
    # colnames(coeff)[1] <- "t"
    # 
    # for (i in 1:length(coeff_names)) {
    #   cd = str_extract(coeff_names[i], regex("(?<=d)(A|P)(?=\\dd)"))              # clade of dependent species
    #   ci = str_extract(coeff_names[i], regex("(?<=d(A|P)\\dd)(A|P)"))             # clade of independent species
    #   
    #   snum_d = str_extract(coeff_names[i], regex("(?<=d(A|P))\\d+(?=d)"))         # which dependent species
    #   snum_i = str_extract(coeff_names[i], regex("(?<=d(A|P)\\dd(A|P))\\d+"))     # which independent species
    #   
    #   col_comp = "#cb1249"
    #   col_mut = "#48ab49"
    #   
    #   par(mfrow = c(1,1))
    #   coefflims = c(-1, 1)
    #   trange <- 1:1000
    #   plot(coeff[trange,"t"],coeff[trange,i + 1],type="l",col=if(cd == ci) col_comp else col_mut ,
    #        main = network_identifier,
    #        xlab="time",
    #        ylab=bquote(partialdiff*.(cd)[.(snum_d)] / partialdiff*.(ci)[.(snum_i)]),
    #        ylim=coefflims,xlim=range(trange),lwd=2)
    #   abline(a=0 ,b=0 , lty="dashed", col="black", lwd=.5)
    # }
  }
}  

saveRDS(fact, "network_output/231205_prediction/fact_predictions_231205.rds")

#------------------------------------------------------------------------------------------------------------------------

#     analysis of jacobians

#------------------------------------------------------------------------------------------------------------------------

myColors = c("#FF0000", "#00FF00", "#0000FF",
             "#FFFF00", "#FF00FF", "#00FFFF",
             "#333333", "#777777", "#BFFFBB",
             "#FFCCCC", "#CCFFCC", "#CCCCFF",
             "#FFFFCC", "#FFCCFF", "#CCFFFF",
             "plum", "orange", "wheat4", "brown") 

Re(eigen(jacobians[[10]])$values[1])
tempr<-matrix(0,nrow = 80, ncol = 8)
for(i in 1:80) {
tempr[i,]<-  abs(Re(eigen(jacobians[[i]])$vectors[,1]))
}
#------------------------------------------------------------------------------------------------------------------------
ts.plot(tempr, col = myColors, lwd=3)
ts.plot(ts_data[(tp_begin_collapse-80):(tp_begin_collapse-1),], col = myColors)
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

#     plotting of prediction correlation

#------------------------------------------------------------------------------------------------------------------------

data <- readRDS("network_output/231205_prediction/fact_predictions_231205.rds")

# degree
data %>%
  ggplot(aes(factor(rate), spearman_correlation_degree,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(size = "network size", colour = "heritability") +
  xlab("rate of temperature increase") + 
  ylab("spearman correlation")

# rate
data %>%
  ggplot(aes(factor(rate), spearman_correlation_rate,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(size = "network size", colour = "heritability") +
  xlab("rate of temperature increase") + 
  ylab("spearman correlation")

# eigen
data %>%
  ggplot(aes(factor(rate), spearman_correlation_eigen,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(size = "network size", colour = "heritability") +
  xlab("rate of temperature increase") + 
  ylab("spearman correlation")

#------------------------------------------------------------------------------------------------------------------------

#     miscellaneous

#------------------------------------------------------------------------------------------------------------------------


#str_view_all('M_PL_061_34_none_0.05_7_s-map_231130', regex("M_PL_\\d+(_\\d+)?_(rand)?(none)?_0.\\d+_\\d_s-map_231130.rds")) # _(rand | none)"))

             