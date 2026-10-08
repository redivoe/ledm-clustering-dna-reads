library(tidyverse)

load("output/out-perturbed-ledm.RData")


dataset_names <- c("Original", "+5 edits", "+10 edits", "+15 edits")
editing_labels <- c(K = "Keep",
                    D = "Deletion",
                    S = "Substitution",
                    I = "Insertion")


selected_fit_indices <- integer(length(dataset_names))
selected_K <- integer(length(dataset_names))

for (dataset_id in seq_along(dataset_names)) {
  candidate_indices <- which(out_grid$dataset_id == dataset_id)
  if (length(candidate_indices) == 0L) {
    stop("No fitted models found for dataset_id = ", dataset_id, ".")
  }

  candidate_bic <- vapply(
    out[candidate_indices],
    function(fit) fit$bic,
    FUN.VALUE = numeric(1)
  )
  selected_index <- candidate_indices[which.min(candidate_bic)]
  selected_fit_indices[dataset_id] <- selected_index
  selected_K[dataset_id] <- out_grid$K[selected_index]
}

cat("BIC-selected number of components:\n")
print(
  data.frame(
    dataset = dataset_names,
    K = selected_K
  )
)

theta_rows <- vector("list", length(dataset_names) * 4L)
theta_row_index <- 1L
diagnostic_rows <- vector("list", length(dataset_names))
read_uncertainty_rows <- vector("list", length(dataset_names))

for (dataset_id in seq_along(dataset_names)) {
  fit <- out[[selected_fit_indices[dataset_id]]]
  K <- selected_K[dataset_id]

  theta_values <- fit$theta[c("K", "D", "S", "I")]
  for (editing_state in names(theta_values)) {
    theta_rows[[theta_row_index]] <- data.frame(
      dataset_id = dataset_id,
      dataset = dataset_names[dataset_id],
      editing_state = editing_state,
      editing_operation = editing_labels[editing_state],
      estimate = unname(theta_values[editing_state]),
      stringsAsFactors = FALSE
    )
    theta_row_index <- theta_row_index + 1L
  }

  
}

theta_results <- do.call(rbind, theta_rows)
rownames(theta_results) <- NULL
theta_results$dataset <- factor(
  theta_results$dataset,
  levels = dataset_names
)
theta_results$editing_operation <- factor(
  theta_results$editing_operation,
  levels = unname(editing_labels)
)

theta_results$panel <- ifelse(
  theta_results$editing_state == "K",
  "Keep probability",
  "Editing probabilities"
)

paper_theta_plot <- ggplot2::ggplot(
  theta_results,
  ggplot2::aes(
    x = dataset,
    y = estimate,
    color = editing_operation,
    group = editing_operation
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.7) +
  ggplot2::facet_wrap(
    ggplot2::vars(panel),
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::labs(
    x = "Number of additional edits per read",
    y = "Estimated probability",
    color = NULL
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 11),
    legend.position = "top",
    legend.spacing.x = grid::unit(0.15, "cm"),
    panel.spacing = grid::unit(0.6, "cm")
  )

ggplot2::ggsave(
  "output/real-ledm-editing-probabilities-paper.pdf",
  plot = paper_theta_plot,
  width = 8,
  height = 3.6,
  device = cairo_pdf
)

