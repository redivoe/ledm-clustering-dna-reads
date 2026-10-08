library(tidyverse)
library(latex2exp)
source("../functions/align-labels.R")

load("sim-data.RData")
load("sim-ledm-reps.RData")

out_ledm <- as_tibble(transpose(out_ledm)) |> 
  select(-c(theta))

out_df <- bind_cols(out_df, out_ledm) |> 
  rowwise() |> 
  mutate(ari_cl = mclust::adjustedRandIndex(x = cl_true, y = cl),
         ari_cl_init = mclust::adjustedRandIndex(x = cl_true, y = cl_init),
         out_align = list(align_labels_int(reference = cl_true, labels = cl)),
         out_align_init = list(align_labels_int(reference = cl_true, labels = cl_init)),
         Z = list(Z[out_align$permutation, ]),
         Z_init = list(Z_init[out_align_init$permutation, ]),
         accuracy_z = mean(Z == Z_true),
         accuracy_z_init = mean(Z_init == Z_true)) |> 
  ungroup()


theta_labs <- unique(transpose(settings)$theta) |> 
  map_chr(\(x) paste0("(", paste0(x, collapse = ", "), ")"))

theta_labs_parsed <- paste0("theta == '", theta_labs, "'")


(p_z <- out_df |> 
    mutate(theta = factor(theta, levels = 1:4, labels = theta_labs_parsed)) |> 
    ggplot(aes(x = p, y = accuracy_z, color = factor(K)))+
    geom_boxplot(outlier.shape = 1,
                 outlier.alpha = 0.7,
                 staplewidth = 0.5,
                 width = 0.5,
                 position = position_dodge(width = 0.65))+
    facet_grid(cols = vars(theta),
               labeller = label_parsed)+
               # labeller = labeller(theta = function(x) paste0(expression(theta)," = ", x)))+
    scale_color_viridis_d(end = 0.7)+
    labs(x = "p",
         y = "Accuracy",
         color = "K") +
    theme_bw()+
    theme(strip.background = element_rect(fill = NA)))

ggsave(plot = p_z, filename = "sim-ledm-z-accuracy.pdf",
       width = 8.5, height = 2.5)

(p_accuracy_improvement <- out_df |> 
  mutate(theta = factor(theta, levels = 1:4, labels = theta_labs)) |> 
  ggplot(aes(x = accuracy_z_init, y = accuracy_z, col = theta))+
  geom_abline(slope = 1, intercept = 0, lty = 3)+
  geom_point(pch = 1, alpha = 0.7)+
  labs(x = "Accuracy initialization", y = "Accuracy final", col = expression(theta))+
  lims(x = c(0.2, 1), y = c(0.2, 1))+
  coord_equal()+
  theme_bw())
ggsave(plot = p_accuracy_improvement,
       filename = "accuracy-improvement.pdf",
       width = 5, height = 3)


p_labs_parsed <- paste0("p == '", levels(out_df$p), "'")

(p_accuracy_improvement_hist <- out_df |> 
  mutate(theta = factor(theta, levels = 1:4, labels = theta_labs_parsed),
         p = factor(p, labels = p_labs_parsed),
         accuracy_improvement = accuracy_z - accuracy_z_init) |> 
  ggplot(aes(x = accuracy_improvement, y = after_stat(density), fill = theta))+
  geom_vline(xintercept = 0, lty = 3)+
  geom_histogram(alpha = 0.7, position="identity")+
  facet_grid(rows = vars(p), cols = vars(theta),
             labeller = label_parsed)+
  labs(y = NULL, x = "Accuracy improvement (final - initial)")+
  theme_bw()+
  theme(strip.background = element_rect(fill = NA),
        legend.position = "none"))
ggsave(plot = p_accuracy_improvement_hist, filename = "accuracy-improvement-hist.pdf",
       width = 8, height = 5.5)



out_ledm <- out_df |> 
  select(nk, p, K, theta, rep, ari_cl, ari_cl_init, accuracy_z, accuracy_z_init)
save(list = "out_ledm", file = "sim-ledm-reps-post-processed.RData")


# 
# out_df |> 
#   mutate(acc_improv = accuracy_z - accuracy_z_init) |> 
#   slice_min(order_by = acc_improv, n = 10) |> 
#   select(p, K, rep, theta, accuracy_z, accuracy_z_init, ari_cl, ari_cl_init)
# 
# ex <- out_df |> 
#   filter(K == 50, rep == 50, theta == 2, p == 50)
# table(as.integer(factor(ex$cl_init[[1]], ex$cl[[1]])
# length(unique(ex$cl_init[[1]]))
# length(unique(ex$cl[[1]]))
# tabulate(ex$cl_init[[1]])
# tabulate(ex$cl[[1]])
# ?tabulate
