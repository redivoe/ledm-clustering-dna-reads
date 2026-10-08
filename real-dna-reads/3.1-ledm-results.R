library(tidyverse)
library(ggh4x)
library(gt)
source("../functions/align-labels.R")

load("output/subset-perturbed-K50.RData")
load("output/out-perturbed-ledm.RData")

out_grid <- out_grid |> 
  mutate(out = out) |> 
  rowwise() |> 
  mutate(bic = out$bic,
         cl = list(out$cl),
         cl_init = list(out$cl_init),
         Z = list(out$Z),
         Z_init = list(out$Z_init),
         theta = list(out$theta)) |> 
  group_by(dataset_id) |>
  mutate(bic_seq = list(bic), K_seq = list(K)) |> 
  slice_min(order_by = bic, n = 1)


(p_sel <- out_grid |> 
  select(dataset_id, bic_seq, K_seq) |> 
  unnest(c(bic_seq, K_seq)) |> 
  ggplot(aes(x = K_seq, y = bic_seq))+
  geom_line()+
  geom_point()+
  facet_grid2(cols = vars(dataset_id),
              scales = "free_y",
              independent = "y",
              labeller = labeller(dataset_id = c("1" = "Original", "2" = "+5 edits", "3" = "+10 edits", "4" = "+15 edits")))+
  labs(x = "K", y = "BIC")+
  theme_bw()+
  theme(strip.background = element_rect(fill = NA)))

ggsave(filename = "output/real-ledm-selection.pdf", plot = p_sel, width = 9, height = 2.5)

nucleotides <- c("A", "C", "G", "T")
Z_true <- map(unique(data$ref[order(data$cl)]), str_split_1, pattern = "") |> 
  map(\(x) match(x, table = nucleotides)) |> 
  do.call(args = _, what = rbind)

out_tab <- out_grid |> 
  rowwise() |> 
  mutate(K_effective = length(unique(cl)),
         ari = mclust::adjustedRandIndex(x = cl, y = data$cl),
         out_align_init = list(align_labels_int(reference = as.integer(factor(data$cl)), labels = cl_init)), 
         out_align = list(align_labels_int(reference = as.integer(factor(data$cl)), labels = cl)),
         Z = list(Z[out_align$permutation[1:nrow(Z_true)], ]),
         Z_init = list(Z_init[out_align_init$permutation[1:nrow(Z_true)], ]),
         accuracy_z = mean(Z == Z_true),
         accuracy_z_init = mean(Z_init == Z_true)) |>
  select(dataset_id, K, K_effective, ari, accuracy_z, accuracy_z_init)


out_tab |> 
  select(dataset_id, accuracy_z, accuracy_z_init) |> 
  pivot_longer(cols = starts_with("accuracy"),
               names_to = "metric",
               values_to = "value") |> 
  pivot_wider(names_from = dataset_id,
              values_from = value,
              names_prefix = "dataset_") |> 
  gt() |> 
  fmt_number(columns = where(is.numeric),
             decimals = 2) |> 
  cols_label(`metric` = "",
             `dataset_1` = "Original",
             `dataset_2` = "+5 edits",
             `dataset_3` = "+10 edits",
             `dataset_4` = "+15 edits")

# editing probabilities
map(out_grid$theta, \(x) sum(x[-1]))
