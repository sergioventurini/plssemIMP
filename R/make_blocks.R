make_blocks <- function(model) {
  model_csem <- cSEM::parseModel(model)
  meas_matrix <- model_csem$measurement
  indicators <- model_csem$indicators

  blocks <- list()
  for (i in 1:nrow(meas_matrix)) {
    blocks[[i]] <- indicators[meas_matrix[i, ] == 1]
  }
  names(blocks) <- paste0("B", 1:nrow(meas_matrix))

  blocks
}
