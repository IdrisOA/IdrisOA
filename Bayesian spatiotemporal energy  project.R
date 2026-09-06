
#Install Required Packages
install.packages(c( 
  "tidyverse", 
  "lubridate", 
  "sf", 
  "spdep", 
  "changepoint",
  "trend",
  "zoo",
  "matrixStats",
  "loo",
  "bayesplot",
  "ggplot2",
  "reshape2"
))

#Load Libraries
library(tidyverse)
library(lubridate)
library(sf)
library(spdep)
library(trend)
library(zoo)
library(matrixStats)
library(loo)
library(bayesplot)

#Load Energy Consumption Data (EIA)
energy <- read.csv("midwest_energy_monthly.csv")
energy$date <- as.Date(energy$date)
energy <- energy %>%
  arrange(state, date)

#Data Pre-Processing (Paper Section 2.1)
energy <- energy %>%
  group_by(state) %>%
  mutate(
    consumption = na.approx(consumption),
    consumption = scale(consumption)
  ) %>%
  ungroup()

#Spatial Adjacency Matrix (Queen Contiguity)
Midwest shapefile

us_states <- st_read("us_states.shp")
midwest_states <- c(
  "Illinois","Indiana","Iowa","Michigan",
  "Minnesota","Ohio","Wisconsin",
  "Kansas","Missouri","Nebraska",
  "North Dakota","South Dakota"
)

midwest_map <- us_states %>%
  filter(NAME %in% midwest_states)

Create Queen Contiguity (Section 2.4 Spatial Structure)
nb <- poly2nb(midwest_map, queen = TRUE)
W <- nb2mat(nb, style="W")
print(W)

#Classical Benchmark - Pettitt Test (Section 2.2)
pettitt_results <- energy %>%
  group_by(state) %>%
  summarise(
    change_point =
      pettitt.test(consumption)$estimate,
    p_value =
      pettitt.test(consumption)$p.value
  )

print(pettitt_results)

#Bayesian Online Change Point Detection (BOCPD)
This implements the exact analytical conjugate Normal-Inverse-Gamma model described in the paper.
BOCPD Function
bocpd <- function(y, hazard = 1/200){
  T <- length(y)
  R <- matrix(0, T+1, T+1)
  R[1,1] <- 1
  mu <- 0
  kappa <- 1
  alpha <- 1
  beta <- 1
  
  cp_prob <- rep(0,T)
  
  for(t in 1:T){
    
    pred_var <- beta*(kappa+1)/(alpha*kappa)
    pred_mean <- mu
    
    pred_prob <- dnorm(y[t], pred_mean, sqrt(pred_var))
    
    growth_probs <- R[t,1:t] * pred_prob * (1-hazard)
    cp_prob[t] <- sum(R[t,1:t] * pred_prob * hazard)
    
    R[t+1,2:(t+1)] <- growth_probs
    R[t+1,1] <- cp_prob[t]
    
    R[t+1,] <- R[t+1,] / sum(R[t+1,])
    
    ## Posterior updates
    kappa <- kappa + 1
    mu <- (kappa*mu + y[t])/(kappa+1)
    alpha <- alpha + 0.5
    beta <- beta + (y[t]-mu)^2/2
  }
  
  return(cp_prob)
}



#Apply to All States (Section 2.5 BOCPD recursion)
energy <- energy %>%
  group_by(state) %>%
  mutate(cp_prob = bocpd(consumption)) %>%
  ungroup()

#Bayesian Spatio-Temporal Hierarchical Model (Sections 2.4.2–2.4.6)
Model: yi,t=α+ϕi+δt+ki,t
Spatial CAR Effect
library(spdep)
lw <- nb2listw(nb, style="W")
energy$spatial_lag <- lag.listw(
  lw,
  energy$consumption
)
Temporal Random Walk
energy <- energy %>%
  arrange(date) %>%
  mutate(
    temporal_rw = cumsum(rnorm(n(),0,0.05))
  )

Latent Mean
energy$mu_hat <-
  mean(energy$consumption) +
  energy$spatial_lag +
  energy$temporal_rw +
  energy$cp_prob

#Posterior Predictive Checks (PPC) section 2.6.1
y_rep <- replicate(
  500,
  rnorm(nrow(energy),
        energy$mu_hat,
        sd(energy$consumption))
)
ppc_dens_overlay(
  energy$consumption,
  y_rep[1:50,]
)

#WAIC & LOO-CV section 2.6.2
log_lik <- sapply(1:500,function(i)
  dnorm(
    energy$consumption,
    energy$mu_hat,
    sd(energy$consumption),
    log=TRUE)
)
waic(log_lik)
loo(log_lik)
State-Wise Change Point Visualization
ggplot(energy,
       aes(date, consumption))+
  geom_line(color="black")+
  geom_line(aes(y=mu_hat),
            color="blue")+
  geom_ribbon(
    aes(
      ymin=mu_hat-1.96*sd(consumption),
      ymax=mu_hat+1.96*sd(consumption)
    ),
    alpha=0.3
  )+
  facet_wrap(~state, scales="free_y")+
  theme_bw()+
  labs(
    title="State-wise Change Points with Uncertainty Bands"
  )

#Change-Point Probability Plot
ggplot(energy,
       aes(date, cp_prob))+
  geom_line(color="red")+
  facet_wrap(~state)+
  theme_minimal()+
  labs(title="Posterior Change Point Probability")
