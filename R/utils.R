# fast samp-by-grp
samp_by_grp <- function(samp, pop, dom_nm, B) {
  
  num_plots <- dplyr::count(samp, !!rlang::sym(dom_nm))
  # our boot_pop_data has column name domain as its group variable
  setup <- dplyr::count(pop, domain) |> 
    dplyr::left_join(num_plots, by = c(domain = dom_nm)) |>
    dplyr::mutate(add_to = dplyr::lag(cumsum(n.x), default = 0)) |>
    dplyr::rowwise() |> 
    dplyr::mutate(map_args = list(list(n.x, n.y, add_to))) 
  
  all_samps <- vector("list", length = B)
  
  for (i in 1:B) {
    ids <- setup |>
      dplyr::mutate(samps = purrr::pmap(.l = map_args, .f = \ (x, y, z) {
        sample(1:x, size = y, replace = TRUE) + z
      })) |>
      dplyr::pull(samps) |>
      unlist()
    
    out <- pop[ids, ]
    all_samps[[i]] <- out
  }
  
  return(all_samps)
}



# fit_zi function

# don't do prediction here
# predict_zi

# take the mean of the pixels in that county
# then predict on those means

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
  glmer_z <- suppressMessages(lme4::glmer(log_reg_formula, data = samp_dat, family = "binomial"))
  
  # dont do this
  unit_level_preds <- setNames(
    stats::predict(lmer_nz, pop_dat, allow.new.levels = TRUE) * stats::predict(glmer_z, pop_dat, type = "response"),
    as.character(pop_dat[ , domain_level, drop = T])
  ) 
  
  # idea: just return model params and fit later

  
  zi_domain_preds <- aggregate(unit_level_preds, by = list(names(unit_level_preds)), FUN = mean)
  
  names(zi_domain_preds) <- c("domain", "Y_hat_j")
  
  return(list(lmer = lmer_nz, glmer = glmer_z, pred = zi_domain_preds))
  
}


# predict_zi <- function(mod1, mod2, data) {
#   
# }

# base version of dplyr::slice_sample
slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}


# base version of stringr::str_extract_all
str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}


# bootstrap rep helper
boot_rep <- function(boot_samp,
                     pop_boot,
                     domain_level,
                     boot_lin_formula,
                     boot_log_formula,
                     boot_truth) {
  
  # capture warnings and messages silently when bootstrapping
  fit_zi_capture <- capture_all(fit_zi)
  
  boot_samp_fit  <- tryCatch(
    {
      fit_zi_capture(boot_samp,
                     pop_boot,
                     boot_lin_formula,
                     boot_log_formula,
                     domain_level)
    },
    error = function(cond) {
      zi_domain_preds <- boot_truth
      zi_domain_preds$domain_est <- NA
      names(zi_domain_preds) <- c("domain", "Y_hat_j")
      list(result = list(lmer = NA,
                         glmer = NA,
                         pred = zi_domain_preds),
           log = cond)
    
    }
  )
  
  squared_error <- merge(x = boot_samp_fit$result$pred,
                         y = boot_truth,
                         by = "domain",
                         all.x = TRUE) |>
    transform(sq_error = (Y_hat_j - domain_est)^2)
  
  squared_error <- squared_error[ , c("domain", "sq_error")]
  
  return(list(sqerr = squared_error, log = boot_samp_fit$log))
  
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


capture_all <- function(.f){
  
  .f <- purrr::as_mapper(.f)
  
  function(...){
    
    try_out <- suppressMessages(suppressWarnings(
      try(.f(...), silent = TRUE)
    ))
    
    res <- rlang::try_fetch(
      .f(...),
      error = function(err) rlang::abort("Failed.", parent = err),
      warning = function(warn) warn,
      message = function(message) message,
    )
    
    out <- list(
      result = NULL,
      log = NULL
    )
    
    if("error" %in% class(res)) {
      stop(res$message)
    } else if (!any(c("warning", "message") %in% class(res))){
      out$result <- try_out
      out$log <- NA
    } else {
      out$result <- try_out
      out$log <- res$message
    }

    return(out)
    
  }
  
}



