
source("code/functions.R")

# Data ----
pairs <- data.table::fread("data/transmission_pairs.csv")[, -1]

# Informative priors ----
mu_mean_priors <- c(4.848, 5.0)
mu_sigma_priors <- c(0.610, 0.9)
par1_mean_priors <- c(2.305, 1.548)
par1_sigma_priors <- c(0.439, 0.178)

types <- unique(pairs$type)

stan_data <- lapply(types, function(i) {
  dat <- pairs |> filter(type==i)
  data_list <- list(
    N = nrow(dat),
    EL = dat[, "EL"],
    ER = dat[, "ER"] + 1,
    CL = dat[, "CL"],
    CR = dat[, "CR"] + 1,
    S1L = dat[, "S1"],
    S1R = dat[, "S1"] + 1,
    mu_mean_prior = mu_mean_priors,
    mu_sigma_prior = mu_sigma_priors,
    par1_mean_prior = par1_mean_priors,
    par1_sigma_prior = par1_sigma_priors
  )
  init_list = list(
    e_raw = runif(nrow(dat), .1, .5),
    c_raw = runif(nrow(dat), .5, .95),
    s1_raw = runif(nrow(dat), .5, .95)
  )
  list(data=dat, data_list=data_list, init_list=init_list)
})
names(stan_data) <- types




