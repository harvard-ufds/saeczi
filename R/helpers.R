slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}

# `%dopar%` <- foreach::`%dopar%`

str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}


samp_by_domain <- function(data_list, freq) {
  
  raw <- mapply(slice_samp,
         .data = data_list,
         n = freq,
         SIMPLIFY = F) 
  
  bound <- do.call(rbind, raw)
  
  return(bound)
  
} 
