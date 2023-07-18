slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}

# `%dopar%` <- foreach::`%dopar%`

str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}
