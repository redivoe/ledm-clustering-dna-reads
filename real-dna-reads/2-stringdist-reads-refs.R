library(stringdist)
library(furrr)
library(tidyverse)

refs <- readLines(con = "data-nbt17/nbt17-refs-preprocessed.txt")
reads <- readLines(con = "data-nbt17/nbt17-reads-preprocessed.txt")

find_closest_reads <- function(ref, dist_threshold = 10){
  d <- stringdist(ref, reads, method = "lv", nthread = 1)
  return(which(d <= dist_threshold))
}

refs_subset <- refs[1001:1500]

options(future.globals.maxSize = 2000 * 1024 ^ 2)
plan(multicore, workers = 16)
idx_reads <- future_map(.x = refs_subset,
                        .f = find_closest_reads,
                        .options = furrr_options(scheduling = 2, seed = NULL))

data <- tibble("idx_reads" = idx_reads) |>
  mutate(cl = 1:n(),
         ref = refs_subset) |>
  unnest(idx_reads) |>
  mutate(read = reads[idx_reads])

save(list = c("data"), file = "subset-stringdist.RData")
load("output/subset-stringdist.RData")

set.seed(321)
data <- data |>
  nest_by(cl) |>
  mutate(cluster_size = nrow(data)) |>
  ungroup() |>
  slice_sample(n = 50) |>
  unnest(data) |> 
  filter(!str_detect(read, "N"),
         !str_detect(ref, "N"))

perturbe_read <- function(read, idx){
  nchar_read <- nchar(read)
  operation <- sample.int(n = 3, size = 1)
  if(operation == 1){
    return(paste0(str_sub(string = read, start = 1, end = idx-1),
                  str_sub(string = read, start = idx+1, end = nchar_read),
                  collapse = ""))
  }else if(operation == 2){
    str_sub(string = read, start = idx, end = idx) <- sample(x = setdiff(c("A", "C", "G", "T"), str_sub(read, idx, idx)), size = 1)
    return(read)
  }else if(operation == 3){
    return(paste0(str_sub(string = read, start = 1, end = idx),
           sample(x = c("A", "C", "G", "T"), size = 1),
           str_sub(string = read, start = idx+1, end = nchar_read),
           collapse = ""))
  }
}

set.seed(10)
data <- data |> 
  rowwise() |> 
  mutate(perturbed_idx = list(sample.int(n = nchar(read), size = 15, replace = FALSE)),
         read_perturbed_5 = reduce(.x = perturbed_idx[1:5], .f = perturbe_read, .init = read),
         read_perturbed_10 = reduce(.x = perturbed_idx[1:10], .f = perturbe_read, .init = read),
         read_perturbed_15 = reduce(.x = perturbed_idx, .f = perturbe_read, .init = read),
         dist_lv_5 = stringdist(a = read_perturbed_5, b = read, method = "lv"),
         dist_lv_10 = stringdist(a = read_perturbed_10, b = read, method = "lv"),
         dist_lv_15 = stringdist(a = read_perturbed_15, b = read, method = "lv")) |> 
  ungroup()

save(list = c("data"), file = "output/subset-perturbed-K50.RData")

