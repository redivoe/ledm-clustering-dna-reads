pkgs <- c("stringdist", "matrixStats", "MCMCpack")
is_available <- sapply(pkgs, requireNamespace, quietly = TRUE)
if(!all(is_available)){
  stop("Missing packages: ", paste0(pkgs[!is_available], collapse = " "))
}
library(Rcpp)
library(purrr)
library(furrr)
sourceCpp("../functions/ledm-functions-band-patched.cpp")
source("../functions/ledm-em.R")

nk <- 15
p <- 100
K <- 10
theta <- c(0.9, 0.03, 0.04, 0.03)

n_reps <- 100

set.seed(14)
data <- map(1:n_reps, \(i) sample_ledm(nk = rep(nk, K), p = p, theta = theta))

plan(multisession(workers = 16))

out <- future_map(data, \(x) {
  if (is.null(getOption("ledm_cpp_loaded"))) {
    Rcpp::sourceCpp("../functions/ledm-functions-band-patched.cpp")
    options(ledm_cpp_loaded = TRUE)
  }
  ledm_em_sel(
    X = x$X,
    p = p,
    K_seq = 1:20,
    selection_criterion = "bic",
    B = 10,
    cl_init_method = "stringdist",
    Z_init_method = "consensus"
  )
}, .options = furrr_options(seed = 1))

save(list = c("data", "out"), file = "sim-model-selection-k10.RData")

