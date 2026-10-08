pkgs <- c("stringdist", "matrixStats", "MCMCpack")
is_available <- sapply(pkgs, requireNamespace, quietly = TRUE)
if(!all(is_available)){
  stop("Missing packages: ", paste0(pkgs[!is_available], collapse = " "))
}
library(Rcpp)
library(purrr)
library(dplyr)
library(furrr)
sourceCpp("../functions/ledm-functions-band-patched.cpp")
source("../functions/ledm-em.R")


ledm_em_band <- function(X,
                         p,
                         K,
                         B_seq,
                         Z_init_method = "consensus",
                         cl_init_method = "stringdist",
                         stringdist_method = "lv",
                         init = NULL,
                         seed = NULL,
                         max_em_iter = 100,
                         eps = 1e-3) {
  
  out_band <- vector(mode = "list", length = length(B_seq))
  for(i in 1:length(B_seq)){
    start_time <- Sys.time()
    out_band[[i]] <- ledm_em(
      X = X,
      p = p,
      K = K,
      B = B_seq[i],
      Z_init_method = Z_init_method,
      cl_init_method = cl_init_method,
      stringdist_method = stringdist_method,
      init = init,
      seed = seed,
      max_em_iter = max_em_iter,
      eps = eps
    )
    end_time <- Sys.time()
    out_band[[i]]$elapsed_time <- difftime(end_time, start_time, units = "secs") |>
      as.numeric()
  }

  return(out_band)
}


nk <- 15
p <- 100
K <- 10
theta <- list(c(0.6, 0.11, 0.18, 0.11),
              c(0.9, 0.03, 0.04, 0.03))

n_reps <- 100

settings <- cross(list("nk" = nk,
                       "p" = p,
                       "K" = K,
                       "theta" = theta,
                       "rep" = 1:n_reps))

set.seed(3)
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

# minimum possible B for every dataset
# map_dbl(data, \(x) max(abs(sapply(x$X, length) - p)))

B_seq <- seq(10, 100, by = 10)

plan(multisession(workers = 16))

out_band <- future_map2(data,
                        settings,
                        \(x, y) {
                          if (is.null(getOption("ledm_cpp_loaded"))) {
                            Rcpp::sourceCpp("../functions/ledm-functions-band-patched.cpp")
                            options(ledm_cpp_loaded = TRUE)
                          }
                          ledm_em_band(X = x$X,
                                       p = y$p,
                                       K = y$K,
                                       B_seq = B_seq,
                                       cl_init_method = "stringdist",
                                       Z_init_method = "consensus",
                                       seed = 123)
                        },
                        .options = furrr_options(seed = 1))

save(list = c("out_df", "out_band", "B_seq"),
     file = "sim-band.RData")

