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
#' 
#' @returns 
#' An object of class `zi_mod` with defined `print()` and `summary()` methods. 
#' The object is structured like a list and contains the following elements:
#' 
#' * call: The original function call
#' 
#' * res: A data.frame containing the estimates and mse estimates
#' 
#' * bootstrap_log: A list containing any messages and warnings that occur during the bootstrap process
#' 
#' * lin_mod: The modeling object used to fit the original linear model
#' 
#' * log_mod: The modeling object used to fit the original logistic model
#' 
#' @examples 
#' data(pop)
#' data(samp)
#' 
#' lin_formula <- DRYBIO_AG_TPA_live_ADJ ~ tcc16 + elev
#' 
#' result <- saeczi(samp,
#'                  pop, 
#'                  lin_formula,
#'                  log_formula = lin_formula,
#'                  domain_level = "COUNTYFIPS",
#'                  mse_est = FALSE)
#'
#' @export saeczi
#' @import stats
#' @importFrom progressr progressor with_progress
#' @importFrom furrr future_map furrr_options future_map2
#' @importFrom purrr map map2 map_dfr
#' @importFrom methods is

saeczi <- function(samp_dat,
                   pop_dat,
                   lin_formula,
                   log_formula = lin_formula,
                   domain_level,
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
    message("In order for the internal processes to be run in parallel a `future::plan()` must be specified by the user")
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
  
  all_preds <- unique(lin_X, log_X)
  
  original_out <- fit_zi(
    samp_dat,
    lin_formula,
    log_formula,
    domain_level
  )
  
  mod1 <- original_out$lmer
  mod2 <- original_out$glmer
  .data <- pop_dat[ ,c(all_preds, domain_level)]
  
  unit_level_preds <- setNames(
    stats::predict(mod1, newdata = .data, allow.new.levels = TRUE) * stats::predict(mod2, newdata = .data, allow.new.levels = TRUE, type = "response"),
    as.character(.data[ , domain_level, drop = T])
  )
  
  zi_domain_preds <- aggregate(unit_level_preds, by = list(names(unit_level_preds)), FUN = mean)
  names(zi_domain_preds) <- c("domain", "Y_hat_j")
  
  original_pred <- zi_domain_preds
  
  if (mse_est) {
    
    zi_model_coefs <- mse_coefs(
      original_out$lmer,
      original_out$glmer
    )
    
    params_and_domain <- data.frame(
      dom = zi_model_coefs$domain_levels,
      b_i = zi_model_coefs$b_i
    )
    
    # new levels need to be set to zero here
    joined_pop_bi <- merge(
      x = data.frame(dom = pop_dat[ , domain_level, drop = T]),
      y = params_and_domain,
      by = "dom",
      all.x = TRUE
    ) 
    
    joined_pop_bi[is.na(joined_pop_bi$b_i), "b_i"] <- 0
    
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
    # tweak for allowing new levels
    pop_doms <- unique(pop_dat[[domain_level]])
    all_doms <- unique(pop_doms, zi_model_coefs$domain_levels)
    area_random_errors <- data.frame(
      dom = all_doms,
      area_random_errors = rnorm(
        length(all_doms),
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
    p_hat_i <- 1/(1 + exp(-(x_log_matrix %*% zi_model_coefs$alpha_1 + joined_pop_bi$b_i)))
    
    # Generate corresponding deltas
    delta_i_star <- rbinom(length(p_hat_i), 1, p_hat_i)
    
    boot_data_generation_params <- list(
      random_effects = random_effects,
      p_hat_i = p_hat_i,
      delta_i_star = delta_i_star
    )
    
    # remove columns if lme4 removes them in the initial fitting process
    x_lin_matrix <- x_lin_matrix[, colnames(model.matrix(original_out$lmer))]
    
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
    
    # create bootstrap samples 
    boot_samp_ls <- samp_by_grp(samp_dat, boot_pop_data, domain_level, B) 
    
    # goal is to not pass boot_pop_data to the map at all
    
    # furrr with progress bar
    boot_rep_with_progress_bar <- function(x, boot_lst) {
      
      p <- progressor(steps = length(x))
      
      res <- 
        furrr::future_map(.x = boot_lst,
                          .f = \(.x) {
                            p()
                            boot_rep(boot_samp = .x,
                                     domain_level,
                                     boot_lin_formula,
                                     boot_log_formula)
                          },
                          .options = furrr_options(seed = TRUE))
      
      beta_lm_mat <- res |>
        map_dfr(.f = ~ .x$params$beta_lm) |>
        as.matrix()
      
      beta_glm_mat <- res |>
        map_dfr(.f = ~ .x$params$beta_glm) |>
        as.matrix()
      
      u_lm <- res |> 
        map_dfr(.f = ~ .x$params$u_lm) |> 
        as.matrix()
      
      u_glm <- res |> 
        map_dfr(.f = ~ .x$params$u_glm) |> 
        as.matrix()
      
      # sometimes u_lm will have fewer domains once it is filtered
      # down to positive response values
      u_lm[is.na(u_lm)] <- 0
      
      preds_full <- generate_mse(.data = boot_pop_data,
                                 truth = boot_truth,
                                 domain_level = domain_level,
                                 beta_lm_mat = beta_lm_mat,
                                 beta_glm_mat = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 lin_X = lin_X,
                                 log_X = log_X)
      
      log_lst <- res |>
        map(.f = ~ .x$log)
      
      list(preds_full, log_lst)
      
    }
    
    if (parallel) {
      
      with_progress({
        boot_res <- boot_rep_with_progress_bar(x = 1:B,
                                               boot_lst = boot_samp_ls)
      }) 
      
    } else {
      
      res <- 
        purrr::map(.x = boot_samp_ls,
                   .f = \(.x) { 
                     boot_rep(boot_samp = .x,
                              domain_level,
                              boot_lin_formula,
                              boot_log_formula)
                   },
                   .progress = list(
                     type = "iterator",
                     clear = TRUE
                   ))
      
      beta_lm_mat <- res |>
        map_dfr(.f = ~ .x$params$beta_lm) |>
        as.matrix()
      
      beta_glm_mat <- res |>
        map_dfr(.f = ~ .x$params$beta_glm) |>
        as.matrix()
      
      u_lm <- res |> 
        map_dfr(.f = ~ .x$params$u_lm) |> 
        as.matrix()
      
      u_glm <- res |> 
        map_dfr(.f = ~ .x$params$u_glm) |> 
        as.matrix()
      
      # sometimes u_lm will have fewer domains once it is filtered
      # down to positive response values
      u_lm[is.na(u_lm)] <- 0
      
      preds_full <- generate_mse(.data = boot_pop_data,
                                 truth = boot_truth,
                                 domain_level = domain_level,
                                 beta_lm_mat = beta_lm_mat,
                                 beta_glm_mat = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 lin_X = lin_X,
                                 log_X = log_X)
      
      log_lst <- res |>
        map(.f = ~ .x$log)
      
      boot_res <- list(preds_full, log_lst)
      
    }
    
    final_df <- setNames(
      boot_res[[1]],
      c("domain", "mse")
    )
    
    final_df <- merge(
      x = final_df,
      y = original_pred,
      by = "domain",
      all.x = TRUE
    )
    
    final_df <- setNames(
      final_df[ ,c("domain", "mse", "Y_hat_j")],
      c("domain", "mse", "est")
    )
    
    bootstrap_log <- boot_res[[2]]
    
  } else {
    
    final_df <- original_pred
    bootstrap_log <- NA
    
  }
  
  
  out <- list(
    call = funcCall,
    res = final_df,
    bootstrap_log = bootstrap_log,
    lin_mod = original_out$lmer,
    log_mod = original_out$glmer
  )
  
  structure(out, class = "zi_mod")
  
}

#' @export
print.zi_mod <- function(x, ...) {
  
  cat("\nCall:\n")
  cat(deparse(x$call))
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

#' @export
summary.zi_mod <- function(object, ...) {
  
  print(summary(object$lin_mod))
  cat("\n")
  print(summary(object$log_mod))
  
}



