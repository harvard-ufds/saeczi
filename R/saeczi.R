#' Fit a zero-inflation estimator.
#'
#' @param samp_dat A data.frame with domains, auxiliary variables, and the response variable of a sample
#' @param pop_dat A data.frame with domains and auxiliary variables of a population.
#' @param lin_formula Formula. Specification of the response and fixed effects of the linear regression model
#' @param log_formula Formula. Specification of the response and fixed effects of the logistic regression model
#' @param domain_level String. The column name in samp_dat and pop_dat that encodes the domain level
#' @param B Integer. The number of bootstraps to be used in MSE estimation.
#' @param mse_est Logical. Whether or not MSE estimation should happen.
#' @param estimand String. Whether the estimates should be 'totals' or 'means'.
#' @param parallel Logical. Should the MSE estimation be computed in parallel
#'
#' @returns
#' An object of class `zi_mod` with defined `print()` and `summary()` methods.
#' The object is structured like a list and contains the following elements:
#'
#' * call: The original function call
#'
#' * res: A data.frame containing the estimates and mse estimates
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
#' result <- saeczi(samp_dat = samp,
#'                  pop_dat = pop,
#'                  lin_formula = lin_formula,
#'                  log_formula = lin_formula,
#'                  domain_level = "COUNTYFIPS",
#'                  mse_est = FALSE)
#'
#' @export saeczi
#' @import stats
#' @importFrom rlang sym
#' @importFrom dplyr summarise group_by mutate left_join
#' @importFrom progressr progressor with_progress
#' @importFrom furrr future_map furrr_options future_map2
#' @importFrom purrr map map2 map_dfr

saeczi <- function(samp_dat,
                   pop_dat,
                   lin_formula,
                   log_formula = lin_formula,
                   domain_level,
                   B = 100L,
                   mse_est = FALSE,
                   estimand = "means",
                   parallel = FALSE) {

  funcCall <- match.call()

  check_inherits("data.frame", samp_dat, pop_dat)
  check_inherits("formula", lin_formula, log_formula)
  check_inherits("character", domain_level, estimand)
  check_inherits("integer", B)
  check_inherits("logical", mse_est, parallel)
  
  check_parallel(parallel)
  check_re(pop_dat, samp_dat, domain_level)

  if(!(estimand %in% c("means", "totals"))) {
    stop("Invalid estimand, must be either 'means' or 'totals'")
  }

  Y <- toString(lin_formula[[2]])

  lin_X <- unlist(str_extract_all_base(deparse(lin_formula[[3]]), "\\w+"))
  log_X <- unlist(str_extract_all_base(deparse(log_formula[[3]]), "\\w+"))
  rand_intercept <- paste0("( 1 | ", domain_level, " )")
  lin_formula <- reformulate(c(lin_X, rand_intercept), response = Y)
  log_formula <- reformulate(c(log_X, rand_intercept), response = paste0(Y, "!= 0"))

  all_preds <- unique(c(lin_X, log_X))

  original_out <- fit_zi(samp_dat,
                         lin_formula,
                         log_formula,
                         domain_level)

  mod1 <- original_out$lmer
  mod2 <- original_out$glmer

  .data <- pop_dat[, c(all_preds, domain_level)]

  original_pred <- collect_preds(mod1, mod2, estimand, .data, domain_level)

  if (mse_est) {

    boot_pop_data <- generate_boot_pop(original_out,
                                       pop_dat,
                                       domain_level,
                                       log_X,
                                       all_preds)

    boot_lin_formula <- reformulate(c(lin_X, rand_intercept), "response")
    boot_log_formula <- reformulate(c(log_X, rand_intercept), "response != 0")

    if (estimand == "means") {
      boot_truth <- boot_pop_data |>
        group_by(!!rlang::sym(domain_level)) |>
        summarise(domain_est = mean(response))
    } else {
      boot_truth <- boot_pop_data |>
        group_by(!!rlang::sym(domain_level)) |>
        summarise(domain_est = sum(response))
    }

    boot_samp_ls <- samp_by_grp(samp_dat, boot_pop_data, domain_level, B)

    if (parallel) {
      with_progress({
        boot_res <- boot_rep_par(x = 1:B,
                                 boot_lst = boot_samp_ls,
                                 domain_level,
                                 boot_lin_formula,
                                 boot_log_formula,
                                 boot_pop_data,
                                 boot_truth,
                                 estimand,
                                 lin_X,
                                 log_X)
        })

      names(boot_res) <- c("preds")

    } else {

      res <-
        map(.x = boot_samp_ls,
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
        map_dfr(.f = ~ .x$beta_lm) |>
        as.matrix()

      beta_glm_mat <- res |>
        map_dfr(.f = ~ .x$beta_glm) |>
        as.matrix()

      u_lm <- res |>
        map_dfr(.f = ~ .x$u_lm) |>
        as.matrix()

      u_glm <- res |>
        map_dfr(.f = ~ .x$u_glm) |>
        as.matrix()

      u_lm[is.na(u_lm)] <- 0

      preds_full <- generate_mse(.data = boot_pop_data,
                                 truth = boot_truth,
                                 domain_level = domain_level,
                                 beta_lm_mat = beta_lm_mat,
                                 beta_glm_mat = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 lin_X = lin_X,
                                 log_X = log_X,
                                 estimand = estimand)
    

      boot_res <- list(preds = preds_full)

    }

    mse_df <- setNames(boot_res$preds,
                       c(domain_level, "mse"))

    final_df <- mse_df |>
      left_join(original_pred, by = domain_level)

  } else {

    final_df <- original_pred

  }

  out <- list(
    call = funcCall,
    res = final_df,
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
  out <- list(
    lin_mod = summary(object$lin_mod),
    log_mod = summary(object$log_mod)
  )

  class(out) <- "summary.zinf_bayes"
  out
}

#' @export
print.summary.zi_mod <- function(x, ...) {
  print(x$lin_mod)
  cat("\n")
  print(x$log_mod)
}
