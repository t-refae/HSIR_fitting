#### user defined params (how to solicit input?) ####

rm(list=ls())
library(dplyr)

## set population size
N <- 1e5

## set (rate) parameters
params <- c(beta = 1.2,
            gamma = 0.4,
            cv = 1)

## set times
tf  <- 50
tps <- seq(1, tf, by=0.1)

# set ICs
init_state <- c(S = 1-1/N, I = 1/N, R=0)


#### max time to capture all dynamics ####

beta_range <- data.frame(beta = c(seq(0.3, 0.6, by=0.15), 
                                  seq(0.6, 1.2, by=0.3), 
                                  seq(0.9, 1.8, by=0.45))
              )

gamma_range <- data.frame(gamma = c(0.2, 0.3, 0.4))

R0_range <- c(1.5, 3, 4.5)

beta_gamma_join <- cross_join(beta_range, gamma_range) %>%
  mutate(beta = round(beta, digits=2),
         gamma= round(gamma, digits=2),
         R0 = round(beta / gamma, digits=2)) %>%
  filter(R0 %in% R0_range) %>%
  arrange(R0, beta, gamma) %>%
  distinct()

cv_range <- data.frame(cv = c(0, 0.5, 1, 2))

# full range of parameter combinations explored
bgcv_join <- cross_join(beta_gamma_join, cv_range)


# different types of fitting to perform

# to quantify wrt n_runs: VOI? NPI (NOT NPI_release) ?
fitting_types <- c("full", "partial", "NPI_release")
n_fitting_types <- c(1, 3, 4)

n_to_runs <- nrow(bgcv_join)*sum(n_fitting_types)

n_to_runs * 5 / 24
