
# assumes both reference and labels are integer vectors
align_labels_int <- function(reference, labels) {
  K <- max(c(reference, labels))
  reference <- factor(reference, levels = seq_len(K))
  labels <- factor(labels, levels = seq_len(K))
  
  tab <- unclass(table(reference, labels))

  permutation <- clue::solve_LSAP(tab, maximum = TRUE)
  
  return(list("permutation" = permutation,
              "labels_permuted" = match(x = labels, table = permutation)))
}


reference <- c(1, 1, 2, 2, 3, 3)
labels <- c(1, 1, 3, 3, 4, 5)

out <- align_labels_int(reference, labels)
out$labels_permuted
out$permutation |> unclass()

