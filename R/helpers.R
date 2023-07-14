slice_samp <- function(.data, n, replace = TRUE) {
  dplyr::slice_sample(.data = .data, n = n, replace = replace)
}

`%dopar%` <- foreach::`%dopar%`

str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}
