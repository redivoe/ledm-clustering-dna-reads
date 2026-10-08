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

load("sim-data.RData")

plan(multisession(workers = 16))

out_ledm <- future_map2(data,
                        settings,
                        \(x, y) {
                             if (is.null(getOption("ledm_cpp_loaded"))) {
                               Rcpp::sourceCpp("../functions/ledm-functions-band-patched.cpp")
                               options(ledm_cpp_loaded = TRUE)
                               }
                             ledm_em(X = x$X,
                                     p = y$p,
                                     K = y$K,
                                     B = 10,
                                     cl_init_method = "stringdist",
                                     Z_init_method = "consensus")
                             },
                           .options = furrr_options(seed = 321))

save(list = c("out_ledm"), file = "sim-ledm-reps.RData")
