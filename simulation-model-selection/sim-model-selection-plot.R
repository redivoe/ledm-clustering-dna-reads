library(tidyverse)
library(latex2exp)

load("sim-model-selection-k10.RData")

# Solo BIC:
K_bic <- tibble(reps = seq_along(out),
                K_selected_bic = map_int(out, "K_selected_bic"),
                K_true = 10)

load("sim-model-selection-k50.RData")

K_bic <- bind_rows(K_bic,
                   tibble(reps = seq_along(out),
                          K_selected_bic = map_int(out, "K_selected_bic"),
                          K_true = 50)) |> 
  mutate(correct = factor(K_selected_bic == K_true, levels = c("TRUE", "FALSE")))

(p_bic_selection <- ggplot(K_bic, aes(x = K_selected_bic, fill = correct)) +
  geom_bar(width = 0.5) +
  facet_grid(cols = vars(K_true),
             scales = "free_x", space = "free_x",
             labeller = labeller(K_true = function(x) paste0("K = ", x)))+
  scale_x_continuous(breaks = c(1:20, 41:60))+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_manual(values = c("Light Sky Blue", "Light Salmon 1 "))+
  labs(x = TeX("Selected number of groups $(\\hat{K})$"),
       y = "Frequency",
       fill = TeX("$\\hat{K} = K$"))+
  theme_bw()+
  theme(strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 10),
        legend.position = "none"))

ggsave(plot = p_bic_selection,
       filename = "sim-model-selection-bic.pdf",
       width = 8.5,
       height = 3.5)

