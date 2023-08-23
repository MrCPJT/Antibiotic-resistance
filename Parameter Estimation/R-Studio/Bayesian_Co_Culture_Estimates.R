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
library(ggmcmc)

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
# This script is only interested in the co-culture data

# Staphylococcus Aureus Co-Culture
sa_co_df <-
  read_excel(
    "Staph-Pseudomonas-Data.xlsx",
    sheet = "S. aureus co-culture",
    skip = 1
  )

# Pseudomonas Aeruginosa Co-Culture
pa_co_df <-
  read_excel(
    "Staph-Pseudomonas-Data.xlsx",
    sheet = "P. aeruginosa co-culture",
    skip = 1
  )

# Parsing the data
# Renaming columns for clarity
col_names = c('time', 'rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6')
colnames(sa_co_df)[which(names(sa_co_df) == "Time (h)")] <-
  "time"
names(sa_co_df)[1:7] <- col_names
colnames(pa_co_df)[which(names(pa_co_df) == "Time (h)")] <-
  "time"
names(pa_co_df)[1:7] <- col_names

## Using averaged replicate data (across all 6 replicates)
sa_co_df$mean <- rowMeans(sa_co_df[,2:7], na.rm=TRUE)
sa_co_df[,2:7] <- NULL

## Using averaged replicate data (across all 6 replicates)
pa_co_df$mean <- rowMeans(pa_co_df[,2:7], na.rm=TRUE)
pa_co_df[,2:7] <- NULL

####

# Normalising the data - Dividing through by initial time population
pa_co_df[, 2] = pa_co_df[, 2] / (filter(pa_co_df, time == 0) %>% select(-time) %>% unlist)
sa_co_df[, 2] = sa_co_df[, 2] / (filter(sa_co_df, time == 0) %>% select(-time) %>% unlist)

#### Parameter Estimation ####

# Preparing information for parameter estimation

# Initial time
t0 = 0.0  
# Remaining times
ts = filter(sa_co_df, time > 0) %>% select(time) %>% unlist  
# Growth data
z = filter(sa_co_df,time>0) %>% select(-time)
z = cbind(z,filter(pa_co_df,time>0) %>% select(-time))             
# Initial population size
y0 = c(sa_co_df[1,2],pa_co_df[1,2])
y0 <- as.numeric(unlist(y0))

# Number of samples
nSamples = nrow(sa_co_df) - 1
# Number of replicates
nReplicates = 2

# Initialising Stan Model
gLV_co <-
  stan_model("gLV_Case_1c.stan", model_name = "ConnorModel")

# MCMC estimation
# 50/50 Train/test split
estimates <- sampling(
  object = gLV_co,
  data = list (
    T  = nSamples,
    nReplicates = nReplicates,
    y0 = y0,
    z  = z[, 1:nReplicates],
    t0 = t0,
    ts = ts
  ),
  seed = 13371,
  chains = 15,
  iter = 25000,
  warmup = 12500,
  init_r = 2,
  show_messages = TRUE
)

#### Results ####

# Specifying type of and printing results
parametersToPlot = c("theta", "sigma", "lp__")
print(estimates, pars = parametersToPlot)

# Initialising ShinyStan visualisation library
launch_shinystan(estimates)
