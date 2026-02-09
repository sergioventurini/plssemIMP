make_conditional_boot_indices <- function(data, R) {
  row_miss <- apply(data, 1, function(x) any(is.na(x)))

  idx_complete   <- which(!row_miss)
  idx_incomplete <- which(row_miss)

  n_complete   <- length(idx_complete)
  n_incomplete <- length(idx_incomplete)

  indices_list <- vector("list", R + 1)
  indices_list[[1]] <- seq_len(nrow(data))

  for (b in seq_len(R)) {
    indices_list[[b + 1]] <- c(
      sample(idx_complete,   n_complete,   replace = TRUE),
      sample(idx_incomplete, n_incomplete, replace = TRUE)
    )
  }

  indices_list
}
