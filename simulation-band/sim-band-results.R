library(tidyverse)
library(patchwork)
source("../functions/align-labels.R")

load("sim-band.RData")

time <- map_depth(out_band, 2, \(x) x$elapsed_time) |> 
  purrr::map(list_simplify)
cl <- map_depth(out_band, 2, \(x) x$cl)
Z <- map_depth(out_band, 2, \(x) x$Z)
loglik <- map_depth(out_band, 2, \(x) tail(x$loglik, 1)) |> 
  purrr::map(list_simplify)
out_df$out <- pmap(.l = list(time, cl, Z, loglik),
                   .f = ~ tibble("B" = B_seq, "time" = ..1, "cl" = ..2, "Z" = ..3, "loglik" = ..4))
out_df <- out_df |> 
  unnest(out) |> 
  rowwise() |> 
  mutate(ari = mclust::adjustedRandIndex(x = cl, y = cl_true),
         out_align = list(align_labels(reference = cl_true, labels = cl)),
         Z = list(Z[out_align$permutation, ]),
         accuracy_z = mean(Z == Z_true)) |> 
  select(nk, p, K, rep, theta, ari, accuracy_z, time, loglik, B) |> 
  group_by(rep, theta) |> 
  arrange(rep, theta, B) |> 
  mutate(diff_loglik = (last(loglik) - loglik) / last(loglik)) |> 
  ungroup()


theta_labs <- list(c(0.6, 0.11, 0.18, 0.11), c(0.9, 0.03, 0.04, 0.03)) |> 
  map_chr(\(x) paste0("(", paste0(x, collapse = ", "), ")"))

plottable <- out_df |> 
  group_by(theta, B) |> 
  summarise(across(.cols = c(time, ari, accuracy_z, diff_loglik),
                   .fns = list("mean" = median, "min" = min, "max" = max),
                   .names = "{.col}_{.fn}"),
            .groups = "drop") |> 
  pivot_longer(cols = -c(theta, B), names_to = c("stat", ".value"), names_pattern = "(.+)_(.+)$") |> 
  mutate(stat = factor(stat, levels = c("time", "ari", "accuracy_z", "diff_loglik")),
         theta = factor(theta, levels = c(1, 2), labels = theta_labs))

p_band <- map2(c("time", "ari"), c("Time (seconds)", "ARI"),
               \(x, y) plottable |> 
                 filter(stat == x) |> 
                 ggplot(aes(x = B, y = mean, ymin = min, ymax = max, col = theta, lty = theta))+
                 geom_linerange(alpha = 0.5)+
                 geom_line()+
                 geom_point()+
                 scale_x_continuous(breaks = B_seq)+
                 scale_color_viridis_d(end = 0.7)+
                 labs(y = y, col = expression(theta), lty = expression(theta))+
                 theme_bw())
               
(p_band_layout <- p_band[[1]] + p_band[[2]] + plot_layout(guides = "collect") &
  theme(legend.position = 'top'))
ggsave(filename = "sim-band.pdf", plot = p_band_layout, width = 6, height = 3)
