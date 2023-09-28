#' Fit a zero-inflation estimator to a dataset
#'
#' Calculate the domain predictions using the zero-inflation estimator, and outputs
#' those domain-level predictions, in the form of a dataframe. It contains
#' the estimates for each domain, as well as the mean squared error estimates should
#' the user choose. The output of the function is a list, with the first item being
#' said dataframe, and the second being the R squared value of the model.
#'
#' @param samp_dat A dataframe with domains, predictor variables, and the response variable of a sample
#' @param pop_dat A dataframe with domains and predictor variables of a population
#' @param lin_formula model formula for the linear regression model
#' @param log_formula model formula for the logistic regression model
#' @param domain_level A string of the column name in the dataframes that reflect the domain level
#' @param B An integer of the number of reps desired for the bootstrap
#' @param mse_est A boolean that specifies if the user
#' @param parallel Compute MSE estimation in parallel
#'
#' @details The arguments `lin_formula`, and `log_formula`
#' can be unquoted or quoted. The function can handle both forms.
#'
#' The two datasets (pop_dat and samp_dat) must have the same column names for the domain level,
#' as well as the predictor variables for the function to work.
#' @export unit_zi
#' @import stats
#' @importFrom progressr progressor with_progress
#' @importFrom furrr future_map furrr_options 
#' @importFrom purrr map
#' @importFrom methods is

unit_zi <- function(samp_dat,
                    pop_dat,
                    lin_formula,
                    log_formula = lin_formula,
                    domain_level = "COUNTYFIPS",
                    B = 100,
                    mse_est = FALSE,
                    parallel = FALSE) {
  
  funcCall <- match.call() 

  if(!("formula" %in% class(lin_formula))) {
    lin_formula <- as.formula(lin_formula)
    message("lin_formula was converted to class 'formula'")
  }

  if(!("formula" %in% class(log_formula))) {
    log_formula <- as.formula(log_formula)
    message("log_formula was converted to class 'formula'")
  }
  
  if (parallel && is(future::plan(), "sequential")) {
    warning("In order for the internal processes to be run in parallel a `future::plan()` must be specified by the user")
    message("See <https://future.futureverse.org/reference/plan.html> for reference on how to use `future::plan()`")
  }

  # creating strings of original X, Y names
  Y <- deparse(lin_formula[[2]])
  
  lin_X <- unlist(str_extract_all_base(
    deparse(lin_formula[[3]]),
    "\\w+"
  ))
  
  log_X <- unlist(str_extract_all_base(
    deparse(log_formula[[3]]),
    "\\w+"
  ))

  original_pred <- fit_zi(
    samp_dat,
    pop_dat,
    lin_formula,
    log_formula,
    domain_level
  )

  if (mse_est == T) {
    
    zi_model_coefs <- mse_coefs(
      original_pred$lmer,
      original_pred$glmer
    )

    params_and_domain <- data.frame(
      dom = zi_model_coefs$domain_levels,
      b_i = zi_model_coefs$b_i
    )

    joined_pop_bi <- merge(
      x = data.frame(dom = pop_dat[ , domain_level, drop = T]),
      y = params_and_domain,
      by = "dom",
      all.x = TRUE
    )

    x_log_matrix <- model.matrix(
      as.formula(paste0(" ~ ", paste(log_X, collapse = " + "))),
      data = pop_dat[ , log_X, drop = F]
    )

    x_lin_matrix <- model.matrix(
      as.formula(paste0(" ~ ", paste(lin_X, collapse = " + "))),
      data = pop_dat[ , lin_X, drop = F]
    )

    # generating random errors
    individual_random_errors <- data.frame(
      dom = pop_dat[ , domain_level, drop = T],
      individual_random_errors = rnorm(
        n = length(pop_dat[ , domain_level, drop = T]),
        mean = 0,
        sd = sqrt(zi_model_coefs$sig2_eps_hat)
      )
    )

    area_random_errors <- data.frame(
      dom = zi_model_coefs$domain_levels,
      area_random_errors = rnorm(
        length(zi_model_coefs$domain_levels),
        mean = 0,
        sd = sqrt(zi_model_coefs$sig2_mu_hat)
      )
    )

    random_effects <- merge(
      x = individual_random_errors,
      y = area_random_errors,
      by = "dom",
      all.x = TRUE
    )

    # predict probability of non-zeros
    p_hat_i <- exp(x_log_matrix %*% zi_model_coefs$alpha_1 + joined_pop_bi$b_i)/
      (1 + exp(x_log_matrix %*% zi_model_coefs$alpha_1 + joined_pop_bi$b_i))

    # Generate corresponding deltas
    delta_i_star <- rbinom(length(p_hat_i), 1, p_hat_i)

    boot_data_generation_params <- list(
      random_effects = random_effects,
      p_hat_i = p_hat_i,
      delta_i_star = delta_i_star
    )
    
    # remove columns if lme4 removes them in the initial fitting process
    x_lin_matrix <- x_lin_matrix[, colnames(model.matrix(original_pred$lmer))]

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
    boot_pop_param <- setNames(
      aggregate(
        response ~ domain,
        data = boot_pop_data,
        FUN = mean
      ),
      c("domain", "domain_est")
    )


    ## bootstrapping -------------------------------------------------------------

    boot_pop_data <- cbind(pop_dat, boot_pop_data)

    # creating bootstrap formula to be used to fit zi-model to bootstrap samples
    boot_lin_formula <- as.formula(
      paste0(
        "response ~ ",
        paste(lin_X, collapse = " + ")
      )
    )
    
    boot_log_formula <- as.formula(
      paste0(
        "response ~ ",
        paste(log_X, collapse = " + ")
      )
    )
    # define these before bootstrap
    boot_truth <- stats::setNames(stats::aggregate(response ~ domain, data = boot_pop_data,
                                                   FUN = mean), c("domain", "domain_est"))
    
    by_domains <- split(boot_pop_data, f = boot_pop_data$domain)

    # furrr with progress bar
    boot_rep_with_progress_bar <- function(x) {

      p <- progressor(steps = length(x))

      res <- x |> future_map( ~{
        p()
        boot_rep(
          boot_pop_data,
          samp_dat,
          domain_level,
          boot_lin_formula,
          boot_log_formula,
          boot_truth,
          by_domains
        )
        },
        .options = furrr_options(seed = TRUE))
    
      res_df <- do.call("rbind", res)
      
    }

    with_progress({
      mse_df <- boot_rep_with_progress_bar(1:B)
    })

    final_df <- setNames(
      aggregate(sq_error ~ domain,
                data = mse_df,
                FUN = function(x) sum(x)/B),
      c("domain", "mse")
      )

    final_df <- merge(
      x = final_df,
      y = original_pred$pred,
      by = "domain",
      all.x = TRUE
    )

    final_df <- setNames(
      final_df[ ,c("domain", "mse", "Y_hat_j")],
      c("domain", "mse", "est")
    )

  } else {

    final_df <- original_pred$pred

  }

  out <- list(
    call = funcCall,
    res = final_df,
    lin_mod = original_pred$lmer,
    log_mod = original_pred$glmer
  )
  
  structure(out, class = "zi_mod")

}

#' @export
print.zi_mod <- function(x, ...) {
  
  cat("\nCall:\n")
  console_cat(deparse(x$call))
  cat("\n\n")
  
  cat("Linear Model: \n")
  cat("- Fixed effects: \n")
  print(summary(x$lin_mod)$coefficients[ ,1])
  cat("\n")
  cat("- Random effects: \n")
  print(summary(x$lin_mod)$varcor)
  cat("\n")
  
  cat("Logistic Model: \n")
  cat("- Fixed effects: \n")
  print(summary(x$log_mod)$coefficients[ ,1])
  cat("\n")
  cat("- Random effects: \n")
  print(summary(x$log_mod)$varcor)
  cat("\n")
  
}

#' #' @export
#' summary.zi_mod <- function(obj, ...) {
#'   
#'   
#'   
#' }
#' 
#' #' @export
#' print.summary.zi_mod <- function(x, ...) {
#'   
#'   
#'   
#' }


