fit_zi <- function(samp_dat, pop_dat, lin_formula, log_formula , domain_level) {

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
    unit_level_preds = unit_level_preds) |>
    dplyr::group_by(domain) |>
    dplyr::summarise(Y_hat_j = mean(unit_level_preds)) |>
    dplyr::ungroup()

  return(list(lmer = lmer_nz, glmer = glmer_z, pred = zi_domain_preds))

}
