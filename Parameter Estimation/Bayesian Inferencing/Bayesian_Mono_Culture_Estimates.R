# Connor Tynan - 2023

# Bayesian methods account for the uncertainty that is present in our observed data.
# Stan enables one to call powerful sampling methods. Create Stan script -> Call in R/Python
# Key ingredients for Bayesian models:
# Process Model:     Modelling the underlying process 
# Measurement Model: Want to model/measure how error in our model is distributed
# Parameter Model:   The prior distribution. Specifying our beliefs about the variables. 
#                    Beliefs are updated from sampling the posterior

#### Preamble ####

rm(list = ls())  # Creating fresh workspace (Deletes current workspace)

# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Initializing Stan packages
library(rstan)
options(mc.cores = parallel::detectCores()) # Options for RStan
rstan_options(auto_write = TRUE)            # Options for RStan
library(bayesplot)
library("shinystan")
  
# Miscellaneous / QOL
library(formatR) # Formatting code

# Importing data
library(readxl)

#### Data ####

# Importing Excel Data
# This script is only interested in the mono-culture data

# Staphylococcus Aureus Mono-Culture
sa_mono_df <-
  read_excel(
    "221108-Staph-Pseudomonas interaction data.xlsx",
    sheet = "S. aureus mono-culture",
    skip = 1
  )

# Pseudomonas Aeruginosa Mono-Culture
pa_mono_df <-
  read_excel(
    "221108-Staph-Pseudomonas interaction data.xlsx",
    sheet = "P. aeruginosa mono-culture",
    skip = 1
  )

# Parsing the data
# Renaming columns for clarity
col_names = c('time', 'rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6')
colnames(sa_mono_df)[which(names(sa_mono_df) == "Time (h)")] <-
  "time"
names(sa_mono_df)[1:7] <- col_names
colnames(pa_mono_df)[which(names(pa_mono_df) == "Time (h)")] <-
  "time"
names(pa_mono_df)[1:7] <- col_names

# Using averaged replicate data (across all 6 replicates)
# sa_mono_df$mean <- rowMeans(sa_mono_df[,2:7], na.rm=TRUE)
# sa_mono_df[,2:7] <- NULL

# Using averaged replicate data (across all 6 replicates)
# pa_mono_df$mean <- rowMeans(pa_mono_df[,2:7], na.rm=TRUE)
# pa_mono_df[,2:7] <- NULL

# Using individual replicate data - Comment out to measure averaged data
sa_mono_df$mean <- sa_mono_df[, 2]
sa_mono_df[, 2:7] <- NULL
sa_mono_df$mean <- rowMeans(sa_mono_df[, 2], na.rm = TRUE)
sa_mono_df <- na.omit(sa_mono_df)

# Using individual replicate data - Comment out to measure averaged data
pa_mono_df$mean <- pa_mono_df[, 2]
pa_mono_df[, 2:7] <- NULL
pa_mono_df$mean <- rowMeans(pa_mono_df[, 2], na.rm = TRUE)
pa_mono_df <- na.omit(pa_mono_df)

####

# Normalising the data - Dividing through by initial time population
sa_mono_df[, 2] = sa_mono_df[, 2] / (filter(sa_mono_df, time == 0) %>% select(-time) %>% unlist)
pa_mono_df[, 2] = pa_mono_df[, 2] / (filter(pa_mono_df, time == 0) %>% select(-time) %>% unlist)

# Please interchange between one of the two species

species_df = sa_mono_df  # Species of interest for estimation

#### Parameter Estimation ####

# Preparing information for parameter estimation

# Initial time
t0 = 0.0  
# Remaining times
ts = filter(species_df, time > 0) %>% select(time) %>% unlist  
# Growth data
z = filter(species_df, time > 0) %>% select(-time)             
# Initial population size (should be == 1, if the data is normalised properly)
y0 = filter(species_df, time == 0) %>% select(-time) %>% unlist 

# Number of samples
nSamples = nrow(species_df) - 1
# Number of replicates
nReplicates = 1 

# Initialising Stan Model
gLV_mono <-
  stan_model("Case_1_gLV.stan", model_name = "ConnorModel")

# MCMC estimation
estimates <- sampling(
  object = gLV_mono,
  data = list (
    T  = nSamples,
    nReplicates = nReplicates,
    y0 = array(y0, dim = 1),
    z  = z[, 1:nReplicates],
    t0 = t0,
    ts = ts
  ),
  seed = 13371,
  chains = 20,
  iter = 30000,
  warmup = 15000,
  init_r = 2,
  show_messages = TRUE
)

#### Results ####

# Specifying type of and printing results
parametersToPlot = c("theta", "sigma", "lp__")
print(estimates, pars = parametersToPlot)

# Initialising ShinyStan visualisation library
launch_shinystan(estimates)

# Comparing fitted data to experimental data

# Preparing results for visualisation
xdata <- data.frame(size = unlist(z[, 1:1]),
                    well = as.vector(matrix(
                      rep(1:1, nSamples),
                      nrow = nSamples, byrow = TRUE
                    )),
                    time = rep(ts, 1))

# Predictions
pred <- as.data.frame(estimates, pars = "z_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(
    lb = quantile(value, probs = 0.25),
    median = quantile(value, probs = 0.5),
    ub = quantile(value, probs = 0.75)
  ) %>%
  bind_cols(xdata)

# Aesthetics
theme_set(theme_bw())

# Line plot of results with error quantile error bars
p1 <- ggplot(pred, aes(x = time, y = size))
p1 <- p1 + geom_point() +
  labs(x = "Time (hrs)", y = "Normalised Population Size") +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.15) +
  theme(panel.grid = element_blank())
