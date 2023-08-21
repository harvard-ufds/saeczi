boot_rep <- function(pop_boot, samp_dat, domain_level, boot_lin_formula, boot_log_formula) {

  boot_truth <- stats::setNames(stats::aggregate(response ~ domain, data = pop_boot,
                                      FUN = mean), c("domain", "domain_est"))

  by_domains <- split(pop_boot, f = pop_boot$domain)

  num_domains <- length(by_domains)
  num_plots <- data.frame(table(samp_dat[ ,domain_level]))
  boot_data <- purrr::map2_dfr(.x = by_domains, .y = num_plots$Freq, slice_samp)
  #boot_data <- samp_by_domain(by_domains, num_plots$Freq)
  
  # make this tryCatch recursive
  fit_zi_tc <- function(samp_dat, pop_dat, lin_formula, log_formula , domain_level) {
    return(
      tryCatch(
        {
          fit_zi(boot_data, pop_boot, boot_lin_formula, boot_log_formula, domain_level)
        },
        error = function(cond) {
          boot_data <- purrr::map2_dfr(.x = by_domains, .y = num_plots$Freq, slice_samp)
          #boot_data <- samp_by_domain(by_domains, num_plots$Freq)
          fit_zi(boot_data, pop_boot, boot_lin_formula, boot_log_formula, domain_level)
        }
      )
    )
  }
  
  boot_samp_fit  <- tryCatch(
    {
      fit_zi_tc(boot_data, pop_boot, boot_lin_formula, boot_log_formula, domain_level)
    },
    error = function(cond) {
      boot_data <- purrr::map2_dfr(.x = by_domains, .y = num_plots$Freq, slice_samp)
      #boot_data <- samp_by_domain(by_domains, num_plots$Freq)
      fit_zi_tc(boot_data, pop_boot, boot_lin_formula, boot_log_formula, domain_level)
    }
  )

  squared_error <- merge(x = boot_samp_fit$pred, y = boot_truth, by = "domain", all.x = TRUE) |>
    transform(sq_error = (Y_hat_j - domain_est)^2)

  squared_error <- squared_error[ , c("domain", "sq_error")]

  return(squared_error)
}
