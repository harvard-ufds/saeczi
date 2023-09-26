# fit_zi function
fit_zi <- function(samp_dat,
                   pop_dat,
                   lin_formula,
                   log_formula,
                   domain_level) {
  
  Y <- deparse(lin_formula[[2]])
  lin_X <- unlist(str_extract_all_base(deparse(lin_formula[[3]]), "\\w+"))
  log_X <- unlist(str_extract_all_base(deparse(log_formula[[3]]), "\\w+"))
  
  # function will always treat domain_level as the random intercept
  rand_intercept <- paste0("( 1 | ", domain_level, " )")
  
  # of form y ~ x_1 + ... + x_n + (1 | domain_level)
  lin_reg_formula <- stats::as.formula(
    paste0(deparse(lin_formula[[2]]), " ~ ",
           paste(lin_X, collapse = " + "), " + ",
           rand_intercept)
  )
  
  # of form y != 0 ~ x_1 + ... + x_n + (1 | domain_level)
  log_reg_formula <- stats::as.formula(
    paste0(deparse(log_formula[[2]]), " != 0 ~ ",
           paste(log_X, collapse = " + "), " + ",
           rand_intercept)
  )
  
  # creating nonzero version of our sample data set
  nz <- samp_dat[samp_dat[ , Y] > 0, ]
  
  # fit linear mixed model on nonzero data
  lmer_nz <- suppressMessages(lme4::lmer(lin_reg_formula, data = nz))
  
  # Fit logistic mixed effects on ALL data
  glmer_z <- suppressMessages(
    lme4::glmer(log_reg_formula, data = samp_dat, family = "binomial")
  )
  
  lin_pred <- stats::predict(lmer_nz, pop_dat, allow.new.levels = TRUE)
  log_pred <- stats::predict(glmer_z, pop_dat, type = "response")
  
  unit_level_preds <- lin_pred*log_pred
  
  # d x 2 dataframe
  # where d = # of domains
  zi_domain_preds <- data.frame(
    domain = pop_dat[ , domain_level, drop = T],
    unit_level_preds = unit_level_preds)
  
  zi_domain_preds <- stats::setNames(stats::aggregate(unit_level_preds ~ domain, data = zi_domain_preds,
                                                      FUN = mean), c("domain", "Y_hat_j"))
  
  return(list(lmer = lmer_nz, glmer = glmer_z, pred = zi_domain_preds))
  
}

# base version of dplyr::slice_sample
slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}


# base version of stringr::str_extract_all
str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}


# bootstrap rep helper
boot_rep <- function(pop_boot,
                     samp_dat,
                     domain_level,
                     boot_lin_formula,
                     boot_log_formula) {
  
  boot_truth <- stats::setNames(stats::aggregate(response ~ domain, data = pop_boot,
                                                 FUN = mean), c("domain", "domain_est"))
  
  by_domains <- split(pop_boot, f = pop_boot$domain)
  
  num_domains <- length(by_domains)
  num_plots <- data.frame(table(samp_dat[ , domain_level]))
  boot_data_ls <- purrr::map2(.x = by_domains, .y = num_plots$Freq, slice_samp)
  boot_data <- do.call("rbind", boot_data_ls)
  
  # make this tryCatch recursive
  fit_zi_tc <- function(samp_dat, pop_dat, lin_formula, log_formula , domain_level) {
    return(
      tryCatch(
        {
          fit_zi(boot_data, pop_boot, boot_lin_formula, boot_log_formula, domain_level)
        },
        error = function(cond) {
          boot_data_ls <- purrr::map2(.x = by_domains, .y = num_plots$Freq, slice_samp)
          boot_data <- do.call("rbind", boot_data_ls)
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
      boot_data_ls <- purrr::map2(.x = by_domains, .y = num_plots$Freq, slice_samp)
      boot_data <- do.call("rbind", boot_data_ls)
      fit_zi_tc(boot_data, pop_boot, boot_lin_formula, boot_log_formula, domain_level)
    }
  )
  
  squared_error <- merge(x = boot_samp_fit$pred, y = boot_truth, by = "domain", all.x = TRUE) |>
    transform(sq_error = (Y_hat_j - domain_est)^2)
  
  squared_error <- squared_error[ , c("domain", "sq_error")]
  
  return(squared_error)
  
}


# helper function to extract model coefs for bootstrap data generation
mse_coefs <- function(lmer_model, glmer_model) {
  
  # from lmer model
  beta_hat <- lmer_model@beta # linear model coefficients
  model_summary_df <- data.frame(summary(lmer_model)$varcor)
  
  sig2_mu_hat <- model_summary_df[1, ]$vcov
  sig2_eps_hat <- subset(model_summary_df, grp == "Residual")$vcov
  
  # from glmer model
  alpha_1 <- glmer_model@beta
  
  b_i <- lme4::ranef(glmer_model)[[1]][,1]
  b_domain_levels <- rownames(lme4::ranef(glmer_model)[[1]])
  
  return(list(
    beta_hat = beta_hat, sig2_mu_hat = sig2_mu_hat,
    sig2_eps_hat = sig2_eps_hat, alpha_1 = alpha_1,
    b_i = b_i, domain_levels = b_domain_levels))
  
}



