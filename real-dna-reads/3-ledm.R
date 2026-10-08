library(Rcpp)
library(purrr)
library(stringr)
library(furrr)

sourceCpp("../functions/ledm-functions-band-patched.cpp")
source("../functions/ledm-em.R")

load("output/subset-perturbed-K50.RData")


dnastring_to_list <- function(dnastring, nucleotides = c("A", "C", "G", "T")){
  map(dnastring, str_split_1, pattern = "") |> 
    map(\(x) match(x, table = nucleotides))
}

X <- map(list(data$read, data$read_perturbed_5, data$read_perturbed_10, data$read_perturbed_15),
         dnastring_to_list)

K_seq <- 40:60
out_grid <- expand.grid(dataset_id = seq_along(X),
                        K = K_seq)

plan(multisession(workers = 16))

out <- future_map(1:nrow(out_grid),
                  \(j){
                    if (is.null(getOption("ledm_cpp_loaded"))) {
                      Rcpp::sourceCpp("../functions/ledm-functions-band-patched.cpp")
                      options(ledm_cpp_loaded = TRUE)
                    }
                    ledm_em(X = X[[out_grid$dataset_id[j]]],
                            p = 110,
                            K = out_grid$K[j],
                            B = 10,
                            cl_init_method = "stringdist",
                            Z_init_method = "consensus",
                            seed = 11)
                  },
                  .options = furrr_options(seed = 1))

save(list = c("out", "out_grid"), file = "out-perturbed-ledm.RData")
