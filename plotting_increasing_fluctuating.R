library(tidyverse)

#data from the simulations, extinction tp of every species was saved, insert appropriate file
data_raw <- readRDS(file = "increasing_fluctuating_231204.rds")   

data <- expand.grid(Temperature=20,
                    `web` = unique(data_raw$web),
                    network_size = 0,
                    connectance = 0,
                    nestedness = 0,
                    h2 = c("none", "rand"),
                    rate = seq(0.075, 0.2, 0.025),
                    env_var = seq(0, 7.5, 1.5)  
            ) %>% arrange_all()

#preparation for *
extinction_tp_list <- list(rep(0, (data_raw$network_size)[1]))
for (i in 1:10) {                                                               # i --> no. networks
  for(j in 1:72) {                                                              # j --> no. parameter combinations for each single network
    network_size <- (data_raw$network_size)[i * 720]                            #        == ((size of data_raw) / i) / no. replicates  
    extinction_tp_list[[(i-1)*72 + j]] <- rep(0, network_size)
  }
}

# *
# 1. transfer network_size etc. to data --> no changes needed
# 2. calculate average extinction tp over 'by' replicates for each run  
pointer_data <- 1
for (i in seq(1, nrow(data_raw), by=10)) {
  
  data[pointer_data, "network_size"] <- data_raw[i,"network_size"] 
  data[pointer_data, "connectance"] <- data_raw[i,"connectance"] 
  data[pointer_data, "nestedness"] <- data_raw[i,"nestedness"] 
  
  network_size <- (data_raw$network_size)[i]
  for(pointer_species in 1:network_size) {
    sum_extinction_tps <- 0
    for(replicate in 1:10) {
      sum_extinction_tps <- sum_extinction_tps + (data_raw$extinction_tps)[[(i-1) + replicate]][pointer_species]
    }
    extinction_tp_list[[pointer_data]][pointer_species] <- sum_extinction_tps / 10
  }
  pointer_data <- pointer_data + 1
}

data <- data %>% mutate(
  extinction_tps = extinction_tp_list
)

data <- data %>% mutate(
  tp_first_extinction = 0,
  tp_last_extinction = 0,
  tp_avg_extinction = 0
)

for (i in 1:nrow(data)) {
  data[i, "tp_first_extinction"] <- unlist(data[i, "extinction_tps"]) %>% min
  data[i, "tp_last_extinction"] <- unlist(data[i, "extinction_tps"]) %>% max
  data[i, "tp_avg_extinction"] <- unlist(data[i, "extinction_tps"]) %>% mean
}

#------------------------------------------------------------
#
#                  Plotting
#
#------------------------------------------------------------

#first extinction
data %>%
  ggplot(aes(rate, tp_first_extinction,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(size = "network size", colour = "heritability") +
  xlab("rate") + 
  ylab("timepoint of first extinction") 


#last extinction
data %>%
  ggplot(aes(rate, tp_last_extinction,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(size = "network size", colour = "heritability") +
  xlab("rate") + 
  ylab("timepoint of last extinction")


#average extinction, for thesis
data %>%
  ggplot(aes(rate, tp_avg_extinction,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(size = "network size", colour = "heritability") +
  xlab("rate") + 
  ylab("average timepoint of extinction")

#first extinction, rate vs. fluctuation
data %>% 
  filter(rate > 0.05) %>% 
  ggplot(aes(rate, env_var)) +
  geom_raster(aes(fill = tp_first_extinction)) +
  scale_fill_gradientn(colours = (heat.colors(10))) +
  labs(fill = "timepoint of\nfirst extinction") +
  xlab("rate") + 
  ylab("fluctuation") + 
  theme_classic() 
  
#last extinction, rate vs. fluctuation
data %>% 
  filter(rate > 0.05) %>% 
  ggplot(aes(rate, env_var)) +
  geom_raster(aes(fill = tp_last_extinction)) +
  scale_fill_gradientn(colours = (heat.colors(10))) +
  labs(fill = "timepoint of\nlast extinction") +
  xlab("rate") + 
  ylab("fluctuation") + 
  theme_classic()

#average extinction tp, rate vs. fluctuation
data %>% 
  filter(rate > 0.05) %>% 
  ggplot(aes(rate, env_var)) +
  geom_raster(aes(fill = tp_avg_extinction)) +
  scale_fill_gradientn(colours = (heat.colors(10))) +
  labs(fill = "average timepoint\nof extinction") +
  xlab("rate") + 
  ylab("fluctuation") +
  theme_classic() 
