#' Fit a zero-inflation estimator to a dataset
#'
#' Calculate the domain predictions using the zero-inflation estimator, and outputs
#' those domain-level predictions, in the form of a dataframe. It contains
#' the estiamtes for each domain, as well as the mean squared error estimates should
#' the user choose. The output of the function is a list, with the first item being
#' said dataframe, and the second being the R squared value of the model.
#'
#' @param samp_dat A dataframe with domains, predictor variables, and the response variable of a sample
#' @param pop_dat A dataframe with domains and predictor variables of a population
#' @param lin_formula A string of the formula for the linear regression model
#' @param log_formula A string of the formula for the logistic regression model
#' @param domain_level A string of the column name in the dataframes that reflect the domain level
#' @param B An integer of the number of reps desired for the bootstrap
#' @param mse_est A boolean that specifies if the user
#' @param boot_type IDK
#'
#' @details The arguments `lin_formula`, and `log_formula`
#' can be unquoted or quoted. The function can handle both forms.
#'
#' The two datasets (pop_dat and samp_dat) must have the same column names for the domain level,
#' as well as the predictor variables for the function to work.
#' @export
unit_zi <- function(samp_dat, pop_dat, lin_formula, log_formula = lin_formula, domain_level = "COUNTYFIPS",
                    B = 100, mse_est = F, boot_type = "parametric"){

  if(!("formula" %in% class(lin_formula))) {
    lin_formula <- stats::as.formula(lin_formula)
    message("lin_formula was converted to class 'formula'")
  }

  if(!("formula" %in% class(log_formula))) {
    log_formula <- stats::as.formula(log_formula)
    message("log_formula was converted to class 'formula'")
  }

  # creating strings of original X, Y names
  Y <- deparse(lin_formula[[2]])
  lin_X <- unlist(str_extract_all_base(deparse(lin_formula[[3]]), "\\w+"))
  log_X <- unlist(str_extract_all_base(deparse(log_formula[[3]]), "\\w+"))

  original_pred <- fit_zi(samp_dat, pop_dat, lin_formula, log_formula, domain_level)

  #R2 <- r.squaredGLMM(original_pred$lmer)[1]
  R2 <- 0

  if (mse_est == T & boot_type == "parametric") {

    # MSE estimation -------------------------------------------------------------
    zi_model_coefs <- mse_coefs(original_pred$lmer, original_pred$glmer)

    params_and_domain <- data.frame(dom = zi_model_coefs$domain_levels,
                                    b_i = zi_model_coefs$b_i)

    joined_pop_bi <- merge(x = data.frame(dom = pop_dat[ , domain_level, drop = T]),
                           y = params_and_domain, by = "dom", all.x = TRUE)


    x_log_matrix <- stats::model.matrix(
      stats::as.formula(paste0(" ~ ", paste(log_X, collapse = " + "))),
      data = pop_dat[ , log_X, drop = F]
    )

    x_lin_matrix <- stats::model.matrix(
      stats::as.formula(paste0(" ~ ", paste(lin_X, collapse = " + "))),
      data = pop_dat[ , lin_X, drop = F]
    )

    # generating random errors
    individual_random_errors <- data.frame(
      dom = pop_dat[ , domain_level, drop = T],
      individual_random_errors = stats::rnorm(
        n = length(pop_dat[ , domain_level, drop = T]),
        mean = 0,
        sd = sqrt(zi_model_coefs$sig2_eps_hat)
      )
    )

    area_random_errors <- data.frame(
      dom = zi_model_coefs$domain_levels,
      area_random_errors = stats::rnorm(
        length(zi_model_coefs$domain_levels),
        mean = 0,
        sd = sqrt(zi_model_coefs$sig2_mu_hat)
      )
    )

    random_effects <- merge(x = individual_random_errors, y = area_random_errors, by = "dom", all.x = TRUE)

    # predict probability of non-zeros
    p_hat_i <- exp(x_log_matrix %*% zi_model_coefs$alpha_1 + joined_pop_bi$b_i)/
      (1 + exp(x_log_matrix %*% zi_model_coefs$alpha_1 + joined_pop_bi$b_i))

    # Generate corresponding deltas
    delta_i_star <- stats::rbinom(length(p_hat_i), 1, p_hat_i)

    boot_data_generation_params <- list(random_effects = random_effects,
                                        p_hat_i = p_hat_i,
                                        delta_i_star = delta_i_star)

    linear_preds <- x_lin_matrix %*% zi_model_coefs$beta_hat +
      boot_data_generation_params$random_effects$area_random_errors +
      boot_data_generation_params$random_effects$individual_random_errors

    boot_pop_response <- as.vector(
      linear_preds * boot_data_generation_params$delta_i_star
    )

    # bootstrap population data
    boot_pop_data <- data.frame(
      domain = pop_dat[ , domain_level, drop = T],
      response = boot_pop_response
    )

    # domain level estimates for bootstrap population data
    boot_pop_param <- stats::setNames(stats::aggregate(response ~ domain,
              data = boot_pop_data,
              FUN = mean), c("domain", "domain_est"))

    #diff <- boot_pop_param$domain_est - truth$biomass


    ## bootstrapping -------------------------------------------------------------

    boot_pop_data <- cbind(pop_dat, boot_pop_data)

    # creating bootstrap formula to be used to fit zi-model to bootstrap samples
    boot_lin_formula <- stats::as.formula(paste0("response ~ ", paste(lin_X, collapse = " + ")))
    boot_log_formula <- stats::as.formula(paste0("response ~ ", paste(log_X, collapse = " + ")))


    # furrr with progress bar
    boot_rep_with_progress_bar <- function(x) {
      p <- progressr::progressor(steps = length(x))

      x |> furrr::future_map_dfr( ~{
        p()
        boot_rep(boot_pop_data, samp_dat, domain_level,
                 boot_lin_formula, boot_log_formula)
      },
      .options = furrr::furrr_options(seed = TRUE))
    }

    progressr::with_progress({
      mse_df <- boot_rep_with_progress_bar(1:B)
    })

    # furrr without progress bar
    # mse_df <- 1:B |>
    #   furrr::future_map_dfr(~ boot_rep(boot_pop_data, samp_dat, domain_level,
    #                         boot_lin_formula, boot_log_formula))

    ########foreach and dopar
    # mse_df <-  foreach::foreach(i = 1:B,
    #                   .combine = 'rbind') %dopar% {
    #                     boot_rep(boot_pop_data, samp_dat, domain_level, boot_lin_formula, boot_log_formula)
    #                   }


    final_df <- stats::setNames(stats::aggregate(sq_error ~ domain,
                                          data = mse_df,
                                          FUN = function(x) sum(x)/B), c("domain", "mse"))

    final_df <- merge(x = final_df, y = original_pred$pred, by = "domain", all.x = TRUE)

    final_df <- stats::setNames(final_df[ ,c("domain", "mse", "Y_hat_j")], c("domain", "mse", "est"))

  } else if (mse_est == T & boot_type == "vanilla") {

    ############# to do

  }
  else {
    final_df <- original_pred$pred
  }

  return(list(final_df, original_pred$lmer, original_pred$glmer))

}
