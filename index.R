# Connor Tynan - November 2022 - % Bayesian Inference

# Bayesian methods account for the uncertainty that is present in our observed data.
# Stan enables one to call powerful sampling methods. Create Stan script -> Call in R/Python
# Key ingredients for Bayesian models:
# Process Model: Modelling the underlying process (differential equations, linear regression models, etc)
# Measurement Model: Want to model/measure how error in our model is distributed
# Parameter Model: The prior distribution. Specifying our beliefs about the variables. Beliefs are updated from sampling the posterior

rm(list = ls())  # Clearing workspace

# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(formatR) # Formatting code

# Initializing Stan package
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Importing data
library(readxl)

sp_mono_df <-
  read_excel(
    "221108-Staph-Pseudomonas interaction data.xlsx",
    sheet = "S. aureus mono-culture",
    skip = 1
  )

pa_mono_df <-
  read_excel(
    "221108-Staph-Pseudomonas interaction data.xlsx",
    sheet = "P. aeruginosa mono-culture",
    skip = 1
  )

# Parsing the data - Adjusting column names, setting to long/melt form, etc
col_names = c('time', 'rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6')

colnames(sp_mono_df)[which(names(sp_mono_df) == "Time (h)")] <- "time"
names(sp_mono_df)[1:7] <- col_names

colnames(pa_mono_df)[which(names(pa_mono_df) == "Time (h)")] <- "time"
names(pa_mono_df)[1:7] <- col_names

sp_mono_df$mean <- rowMeans(sp_mono_df[,2:7], na.rm=TRUE)
sp_mono_df[,2:7] <- NULL

sp_mono_df[,2] = sp_mono_df[,2]/(filter(sp_mono_df,time==0) %>% select(-time) %>% unlist)

# Putting data into long/melt format for manipulation
long_sp_mono <- sp_mono_df %>% gather(mean,size,-time)
glimpse(long_sp_mono)

# Variables / Parameters for estimation
t0 = 0.0  # Initial time
ts = filter(sp_mono_df,time>0) %>% select(time) %>% unlist  # Remaining times
z = filter(sp_mono_df,time>0) %>% select(-time)             # Growth data and respective rep

y0 = filter(sp_mono_df,time==0) %>% select(-time) %>% unlist # Initial population size of reps

nSamples = nrow(sp_mono_df) - 1 # Number of samples
nReplicates = 1 # Number of replicates

# Calling of stan model (gLV.stan)
gLV_mono <- stan_model("gLV.stan", model_name = "ConnorModel")
  
# MCMC estimation of theta[1:2], growth and inter-interactivity parameters respectively
estimates <- sampling(object = gLV_mono,
                      data = list (
                        T  = nSamples,
                        nReplicates = nReplicates,
                        y0 = array(y0,dim=1),
                        z  = z[,1:nReplicates],
                        t0 = t0,
                        ts = ts
                      ),
                      seed = 1241928,
                      chains = 4,
                      iter = 1500,
                      warmup = 750,
                      init_r=1.5,
                      show_messages = TRUE
)

# Collecting fit data

parametersToPlot = c("theta","sigma","lp__")
print(estimates, pars = parametersToPlot)

# Visualisation of fitting

draws <- as.array(estimates, pars=parametersToPlot)
mcmc_trace(draws)

color_scheme_set("brightblue")
mcmc_scatter(draws,pars=c('theta[1]','theta[2]'))

xdata <- data.frame(size = unlist(z[,1:1]),
                    well = as.vector(matrix(rep(1:1,nSamples),
                                            nrow=nSamples,byrow=TRUE)),time = rep(ts,1))

pred <- as.data.frame(estimates, pars = "z_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.25),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.75)) %>%
  bind_cols(xdata)

theme_set(theme_bw())

p1 <- ggplot(pred, aes(x = time, y = size))
p1 <- p1 + geom_point() +
  labs(x = "Time (hrs)", y = "Normalised Population Size") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none")
p1 + geom_line(aes(x = time, y = median)) +
geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.15) +
  theme(panel.grid = element_blank()) + scale_x_continuous(breaks=seq(0,15,1), limits = c(0,14)) + 
  scale_y_continuous(breaks=seq(0,30000,5000))

# Acknowledgements
# Code inspired by: https://shug3502.github.io/blog/DifferentialEqnsStan
