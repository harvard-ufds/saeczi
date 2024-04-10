#' Generates B-many bootstrap samples
#' 
#' Ensures that each group in a bootstrap sample has the same number of rows as 
#' that group did in the original sample data. 
#' 
#' @param samp The sample data
#' @param pop The bootstrap population data
#' @param dom_nm Character string of the domain identifier name as it appears in samp and pop
#' @param B Integer. The number of bootstrap samples to produce
#' 
#' @return A list of length B where each item is a unique bootstrap sample
#' @noRd
#' 
samp_by_grp <- function(samp, pop, dom_nm, B) {
  
  num_plots <- dplyr::count(samp, !!rlang::sym(dom_nm))
  diff <- unique(pop$domain)[!unique(pop$domain) %in% unique(samp[[dom_nm]])]
  if (length(diff) != 0) {
    to_append <- data.frame(
      doms = diff,
      n = rep(0, length(diff))
    )
    colnames(to_append) <- c(dom_nm, "n")
    num_plots <- rbind(num_plots, to_append)
  }
  
  # our boot_pop_data has column name domain as its group variable
  setup <- dplyr::count(pop, !!rlang::sym(dom_nm)) |> 
    dplyr::left_join(num_plots, by = dom_nm) |>
    dplyr::mutate(add_to = dplyr::lag(cumsum(n.x), default = 0)) |>
    dplyr::rowwise() |> 
    dplyr::mutate(map_args = list(list(n.x, n.y, add_to))) 
  
  all_samps <- vector("list", length = B)
  pop_ordered <- pop[order(pop[[dom_nm]]), ]
  
  for (i in 1:B) {
    ids <- setup |>
      dplyr::mutate(samps = purrr::pmap(.l = map_args, .f = \ (x, y, z) {
        sample(1:x, size = y, replace = TRUE) + z
      })) |>
      dplyr::pull(samps) |>
      unlist()
    
    out <- pop_ordered[ids, ]
    all_samps[[i]] <- out
  }
  
  return(all_samps)
  
}


#' Fits the two models in a zi-estimator
#' 
#' @param samp_dat Sample data
#' @param lin_formula Formula to be used in the linear regression model
#' @param log_formula Formula to be used in the logistic regression model
#' @param domain_level Character. Domain identifier name.
#' 
#' @return A list containing the two model objects
#' @noRd
#' 
fit_zi <- function(samp_dat,
                   lin_formula,
                   log_formula,
                   domain_level) {
  
  Y <- deparse(lin_formula[[2]])
  
  # creating nonzero version of our sample data set
  nz <- samp_dat[samp_dat[[Y]] > 0, ]
  
  # fit linear mixed model on nonzero data
  lmer_nz <- suppressMessages(
    lme4::lmer(lin_formula, data = nz)
  )
  
  # Fit logistic mixed effects on ALL data
  glmer_z <- suppressMessages(
    lme4::glmer(log_formula, data = samp_dat, family = 'binomial')
  )
  
  return(list(lmer = lmer_nz, glmer = glmer_z))
  
}

#' Generates mse estimates
#' 
#' @param .data The bootstrap population data
#' @param truth A data.frame containing the true domain level values from the bootstrap population data
#' @param domain_level Character. Domain identifier name
#' @param beta_lm_mat A matrix containing the fixed effects coefficients resulting from fitting the linear model to each bootstrap sample
#' @param beta_glm_mat A matrix containing the fixed effects coefficients resulting from fitting the logistic model to each bootstrap sample
#' @param u_lm A matrix containing the random effects values for each domain resulting from fitting the linear model to each bootstrap sample
#' @param u_glm A matrix containing the random effects values for each domain resulting from fitting the linear model to each bootstrap sample
#' @param lin_X A character vector with the names of the predictor variables used in the linear model
#' @param log_X A character vector with the names of the predictor variables used in the logistic model
#' @param estimand A string specifying whether the estimates should be 'totals' or 'means'
#' 
#' @return A data.frame with the mse estimate for every domain
#' @noRd
generate_mse <- function(.data,
                         truth,
                         domain_level,
                         beta_lm_mat,
                         beta_glm_mat,
                         u_lm,
                         u_glm,
                         lin_X,
                         log_X,
                         estimand) {
  
  boot_pop_by_dom <- split(.data, f = .data[[domain_level]])
  
  design_mat_ls <-  boot_pop_by_dom |> 
    map(.f = function(.x) {
      dmat_lm <- model.matrix(~., .x[ ,lin_X, drop = FALSE])
      dmat_glm <- model.matrix(~., .x[ ,log_X, drop = FALSE])
      return(list(design_mat_lm = dmat_lm,
                  design_mat_glm = dmat_glm))
    })
  
  n_doms <- length(colnames(u_glm))
  dom_order <- names(boot_pop_by_dom)
  u_lm <- u_lm[, order(match(colnames(u_lm), dom_order))]
  u_glm <- u_glm[, order(match(colnames(u_glm), dom_order))]
  
  if (ncol(u_lm) != ncol(u_glm)) {
    diff <- base::setdiff(colnames(u_glm), colnames(u_lm))
    to_append <- matrix(rep(NA, times = nrow(u_lm) * length(diff)), ncol = length(diff))
    colnames(to_append) <- diff
    u_lm <- cbind(u_lm, to_append)
  }
  
  dom_res_wide <- generate_preds(beta_lm = beta_lm_mat,
                                 beta_glm = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 design_mats = design_mat_ls,
                                 J = n_doms,
                                 estimand = estimand)
  
  truth_ordered <- truth[order(match(truth[[domain_level]], dom_order)), ]
  truth_vec <- truth_ordered$domain_est
  
  mse <- (dom_res_wide - truth_vec)^2 |>
    rowMeans(na.rm = TRUE)
  
  res_doms <- data.frame(
    domain = truth_ordered[[domain_level]],
    mse = mse
  )
  
  return(res_doms)
  
}

#' Generate the Bootstrap population data
#' 
#' @param original_out List containing original model objects
#' @param pop_dat The population data frame
#' @param domain_level Character. The domain column names in pop_dat
#' @param log_X Vector of characters containing logistic model predictor names
#' @param all_preds Vector of characters containing all predictor names
#' 
#' @return The population bootstrap data
#' @noRd

generate_boot_pop <- function(original_out, 
                              pop_dat,
                              domain_level,
                              log_X,
                              all_preds) {
  
  zi_mod_coefs <- mse_coefs(original_out$lmer, original_out$glmer)
  
  params_and_domain <- setNames(
    zi_mod_coefs$b_i,
    zi_mod_coefs$domain_levels
  )
  
  pop_b_i <- data.frame(
    dom = pop_dat[ , domain_level, drop = TRUE],
    b_i = params_and_domain[pop_dat[ , domain_level, drop = TRUE]]
  )
  
  pop_b_i[is.na(pop_b_i$b_i), "b_i"] <- 0
  
  x_matrix <- model.matrix(
    as.formula(paste0(" ~ ", paste(all_preds, collapse = " + "))),
    data = pop_dat[ , all_preds, drop = FALSE]
  )
  
  indv_re <- data.frame(
    dom = pop_dat[ , domain_level, drop = TRUE],
    eps_ij = rnorm(nrow(pop_dat), 0, sqrt(zi_mod_coefs$sig2_eps_hat))
  )
  
  # tweak for allowing new levels
  pop_doms <- unique(pop_dat[[domain_level]])
  all_doms <- unique(pop_doms, zi_mod_coefs$domain_levels)
  
  area_re_lkp <- setNames(
    rnorm(length(all_doms), 0, sqrt(zi_mod_coefs$sig2_mu_hat)),
    all_doms
  )
  
  rand_effs <- data.frame(
    indv_re,
    u_j = area_re_lkp[pop_dat[ , domain_level, drop = TRUE]]
  )
  
  x <- x_matrix[ , c("(Intercept)", log_X)] %*% zi_mod_coefs$alpha_1 + pop_b_i$b_i
  p_hat_i <- binomial()$linkinv(x)
  
  delta_i_star <- rbinom(length(p_hat_i), 1, p_hat_i)
  
  boot_dat_params <- list(
    random_effects = rand_effs,
    p_hat_i = p_hat_i,
    delta_i_star = delta_i_star
  )
  
  linear_preds <- (x_matrix[, colnames(model.matrix(original_out$lmer))] %*% zi_mod_coefs$beta_hat) +
    boot_dat_params$random_effects$u_j + boot_dat_params$random_effects$eps_ij
  
  boot_pop_data <- data.frame(
    pop_dat[ , c(domain_level, all_preds)],
    response = linear_preds * boot_dat_params$delta_i_star
  ) 
  
  return(boot_pop_data)
  
}

#' Bootstrap procedure for the parallel option
#' 
#' @param x The vector 1:B where B is the number of total bootstraps
#' @param boot_lst A list where each element contains a bootstrap sample
#' @param domain_level Character. Domain identifier name
#' @param boot_lin_formula The formula to be used for the linear model
#' @param boot_log_formula The formula to be used for the logistic model
#' @param boot_pop_data The bootstrap population data
#' @param boot_truth A data.frame containing the true domain level values from the bootstrap population data
#' @param estimand A string specifying whether the estimates should be 'totals' or 'means'
#' @param lin_X A character vector with the names of the predictor variables used in the linear model
#' @param log_X A character vector with the names of the predictor variables used in the logistic model
#' 
#' @return A list containing the mse estimates data.frame.
#' @noRd
#' 
boot_rep_par <- function(x,
                         boot_lst,
                         domain_level,
                         boot_lin_formula,
                         boot_log_formula,
                         boot_pop_data,
                         boot_truth,
                         estimand,
                         lin_X,
                         log_X) {
  
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
                             log_X = log_X,
                             estimand = estimand)
  
  return(list(preds = preds_full))
  
}




#' Sample n rows from a data.frame
#' 
#' @param .data The data.frame to sample from
#' @param n The number of rows to sample
#' @param replace Logical. Should selected rows be added back to the data.frame for the next row selection?
#' 
#' @return A data.frame with n randomly chosen rows
#' @noRd
#' 
slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}


#' Extract all matches from a string
#' 
#' @param string String to extract matches from
#' @param pattern Regex pattern to use for extractions
#' 
#' @return A vector of the matches
#' @noRd
str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}


#' Format model parameters for mse estimation
#' 
#' @param .fit A custom list containing model fits
#' @param ref  A list containing names to be used to fill in when a model fit failed
#' 
#' @return A list containing all of the properly formated parameters
#' @noRd
#' 
mod_param_fmt <- function(.fit, ref = NULL) {
  
  if (!is.null(.fit)) {
    .lmer <- .fit$lmer
    .glmer <- .fit$glmer
    
    beta_lm <- lme4::fixef(.lmer)
    beta_glm <- lme4::fixef(.glmer)
    
    ref_lm <- lme4::ranef(.lmer)[[1]]
    ref_glm <- lme4::ranef(.glmer)[[1]]
    
    u_lm <- setNames(
      ref_lm[ ,1],
      rownames(ref_lm)
    )
    u_glm <- setNames(
      ref_glm[ ,1],
      rownames(ref_glm)
    ) 
  } else {
    lm_terms <- ref$.lm[!grepl("\\|", ref$.lm)]
    glm_terms <- ref$.glm[!grepl("\\|", ref$.glm)]
    beta_lm <- setNames(
      rep(NA, times = length(lm_terms) + 1),
      c('(Intercept)', lm_terms)
    )
    beta_glm <- setNames(
      rep(NA, times = length(glm_terms) + 1),
      c('(Intercept)', glm_terms)
    )
    u_lm <- setNames(
      rep(NA, times = length(ref$d)),
      ref$d
    )
    u_glm <- u_lm
  }
  
  list(beta_lm = beta_lm,
       beta_glm = beta_glm,
       u_lm = u_lm,
       u_glm = u_glm)
  
}


#' Perform a single bootstrap repetition
#' 
#' @param boot_samp data.frame, An individual bootstrap sample.
#' @param domain_level Character. Domain identifier name
#' @param boot_lin_formula The formula to be used for the linear model
#' @param boot_log_formula The formula to be used for the logistic model
#' 
#' @return A list containing the properly formated model parameters from fitting the two models to the sample data.
#' @noRd
#' 
boot_rep <- function(boot_samp,
                     domain_level,
                     boot_lin_formula,
                     boot_log_formula) {
  
  boot_samp_fit  <- tryCatch(
    {
      out <- fit_zi(boot_samp,
                    boot_lin_formula,
                    boot_log_formula,
                    domain_level)
      
      ps <- mod_param_fmt(out)
      return(ps)
      
    },
    error = function(cond) {
      
      doms <- unique(boot_samp[[domain_level]])
      lm_cfs <- labels(terms(boot_lin_formula))
      glm_cfs <- labels(terms(boot_lin_formula))
      
      ps <- mod_param_fmt(.fit = NULL,
                          ref = list(d = doms,
                                     .lm = lm_cfs,
                                     .glm = glm_cfs))
      
      return(ps)
      
    }
  )
  
  return(boot_samp_fit)
  
}


#' Format model coefficients
#' 
#' @param lmer_model lme4::lmer model object
#' @param glmer_model lme4::glmer model object
#' 
#' @return A list with all of the necessary model parameters
#' @noRd
#' 
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

#' A function factory that causes a function to capture messages, warnings, or errors 
#' 
#' @param .f The function 
#' 
#' @return A function that works just like .f but also returns any messages, warnings, or errors that result from using the function
#' @noRd
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

#' Checking if a param inherits a class
#' 
#' @param x The parameter input(s) to check
#' @param what What class to check if the parameter input inherits
#' 
#' @return Nothing if the check is passed, but an error if the check fails
#' @noRd
check_inherits <- function(x, what) {
  for (i in seq_along(x)) {
    if (!inherits(x[[i]], what)) {
      stop(paste0(x[[i]], " needs to be of class ", what))
    }
  }
  invisible(x)
}

#' Checking if parallel functionality is properly set up
#' 
#' @param x The parameter input to check
#' @param call The caller environment to check in
#' 
#' @return Nothing if the check is passed, but an error if the check fails
#' @noRd
check_parallel <- function(x, call = rlang::caller_env()) {
  
  if (x) {
    if (eval(inherits(future::plan(), "sequential"), envir = call)) {
      message("In order for the internal processes to be run in parallel a `future::plan()` must be specified by the user")
      message("See <https://future.futureverse.org/reference/plan.html> for reference on how to use `future::plan()`")
    }
  }
  
  if (x && future::nbrOfWorkers() == 1) {
    warning("Argument `parallel` is set to true, but only one core is being used")
  }
  
  invisible(x)
}

#' Checking random effect column
#' 
#' @param pop_dat The population dataset to check
#' @param domain_level Character. The domain level identifier.
#' 
#' @return Nothing if the check is passed, but and error if it fails
#' @noRd
check_re <- function(pop_dat, samp_dat, domain_level) {
  if (!(domain_level %in% names(pop_dat))) {
    stop(paste0("Column ", domain_level, " does not exist in pop_dat"))
  }
  if (!(domain_level %in% names(samp_dat))) {
    stop(paste0("Column ", domain_level, " does not exist in samp_dat"))
  }
  if (!inherits(pop_dat[[domain_level]], "character") || !inherits(samp_dat[[domain_level]], "character")) {
    stop(paste0("Column ", domain_level, " must be of type `character` in both pop_dat and samp_dat"))
  }
}

#' Fast aggregation
#' 
#' @noRd
agg_stat <- function(vals, nms, .f) {
  agg <- tapply(vals, nms, .f)
  out <- data.frame(nms = names(agg), vals = agg)
  out
}

#' Predict with both models and return result
#' 
#' @param mod1 Linear model object
#' @param mod2 Logistic model object
#' @param estimand Character, either "means" or "totals"
#' @param .data The data to predict on
#' @param domain_level Character name of domain variable in .data
#' 
#' @returns A data frame of results
#' @noRd

collect_preds <- function(mod1, mod2, estimand, .data, domain_level) {
  
  lin_pred <- predict(mod1, newdata = .data, allow.new.levels = TRUE)
  log_pred <- predict(mod2, newdata = .data, type = "response", allow.new.levels = TRUE)
    
  unit_preds <- lin_pred * log_pred
  
  out <- switch(estimand,
                "means" = agg_stat(unit_preds, .data[[domain_level]], mean),
                "totals" = agg_stat(unit_preds, .data[[domain_level]], sum),
                stop())
  
  
  names(out) <- c(domain_level, "est")
  out
  
}