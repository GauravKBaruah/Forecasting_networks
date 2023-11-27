library(tidyverse)

#data from the simulations, extinction tp of every species was saved
data_raw <- readRDS(file = "increasing_fluctuating_231123.rds")   

data <- expand.grid(Temperature=20,
                    `web` = unique(data_raw$web),
                    network_size = 0,
                    connectance = 0,
                    nestedness = 0,
                    h2 = c("none", "rand"),
                    rate = seq(0.05, 0.2, 0.025),
                    env_var = seq(0.5, 3, 0.5)
            ) %>% arrange_all()

#preparation for *
extinction_tp_list <- list(rep(0, (data_raw$network_size)[1]))
for (i in 1:10) {
  for(j in 1:84) {
    network_size <- (data_raw$network_size)[i*504]
    extinction_tp_list[[(i-1)*84 + j]] <- rep(0, network_size)
  }
}

# *
# 1. transfer network_size etc. to data --> no changes needed
# 2. calculate average extinction tp over 6 replicates for each run  
pointer_data <- 1
for (i in seq(1, nrow(data_raw), by=6)) {
  
  data[pointer_data, "network_size"] <- data_raw[i,"network_size"] 
  data[pointer_data, "connectance"] <- data_raw[i,"connectance"] 
  data[pointer_data, "nestedness"] <- data_raw[i,"nestedness"] 
  
  network_size <- (data_raw$network_size)[i]
  for(pointer_species in 1:network_size) {
    sum_extinction_tps <- 0
    for(replicate in 1:6) {
      sum_extinction_tps <- sum_extinction_tps + (data_raw$extinction_tps)[[(i-1) + replicate]][pointer_species]
    }
    extinction_tp_list[[pointer_data]][pointer_species] <- sum_extinction_tps / 6
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
  geom_point(alpha = 0.25,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(title = "Time point of first extinction vs. rate of temperature increase, facets show temperature variation") +
  xlab("rate") + 
  ylab("first extinction")


#last extinction
data %>%
  ggplot(aes(rate, tp_last_extinction,
             colour = h2)) +
  geom_point(alpha = 0.25,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(title = "Time point of last extinction vs. rate of temperature increase, facets show temperature variation") +
  xlab("rate") + 
  ylab("last extinction")


#average extinction
data %>%
  filter(network_size < 20) %>% 
  ggplot(aes(rate, tp_avg_extinction,
             colour = h2)) +
  geom_point(alpha = 0.35,
             aes(size = network_size)) +
  facet_wrap(~env_var) + 
  theme_bw() +
  labs(title = "Time point of average extinction vs. rate of temperature increase, facets show temperature variation") +
  xlab("rate") + 
  ylab("average extinction")

data %>% 
  #filter(network_size > 20) %>% 
  ggplot(aes(h2, tp_avg_extinction)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(aes(ymin = min(tp_avg_extinction,
                              ymax = max(tp_avg_extinction)))) +
  geom_jitter(alpha = 0.25,
              width = 0.15,
              height = 0,
              aes(colour = factor(network_size))) +
  facet_wrap(~rate, nrow = 2) +
  labs(title = "Time point of average extinction vs. heritability, facets show rate of temperature increase") +
  xlab("heritability") + 
  ylab("time point average extinction")
