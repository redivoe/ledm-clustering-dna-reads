pkgs <- c("stringdist", "matrixStats", "MCMCpack")
is_available <- sapply(pkgs, requireNamespace, quietly = TRUE)
if(!all(is_available)){
  stop("Missing packages: ", paste0(pkgs[!is_available], collapse = " "))
}
library(Rcpp)
library(purrr)
library(dplyr)

source("../functions/ledm-em.R")


nk <- 15
p <- c(50, 100, 150)
K <- c(10, 50)
theta <- list(c(0.6, 0.11, 0.18, 0.11),
              c(0.8, 0.06, 0.08, 0.06),
              c(0.9, 0.03, 0.04, 0.03),
              c(0.95, 0.015, 0.02, 0.015))

n_reps <- 100

settings <- cross(list("nk" = nk,
                       "p" = p,
                       "K" = K,
                       "theta" = theta,
                       "rep" = 1:n_reps))

set.seed(1)
data <- map(settings, \(x) sample_ledm(nk = rep(x$nk, x$K), p = x$p, theta = x$theta))

settings_t <- transpose(settings)
settings_t$data <- data
settings_df <- map(settings_t[c("nk", "p", "K", "rep")], list_simplify) |>
  as_tibble() |> 
  mutate(theta = match(x = settings_t$theta, table = unique(settings_t$theta)),
         across(.cols = everything(), .fns = \(x) factor(x, ordered = TRUE)))
data_df <- as_tibble(transpose(settings_t$data)) |> 
  rename(Z_true = Z,
         cl_true = cl)
out_df <- bind_cols(settings_df, data_df)

save(list = c("data", "settings", "out_df"), file = "sim-data.RData")

