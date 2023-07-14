slice_samp <- function(.data, n, replace = TRUE) {
  dplyr::slice_sample(.data = .data, n = n, replace = replace)
}