# fast samp-by-grp
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
fit_zi <- function(samp_dat,
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
  lmer_nz <- suppressMessages(
    lme4::lmer(lin_reg_formula, data = nz)
  )
  
  # Fit logistic mixed effects on ALL data
  glmer_z <- suppressMessages(
    lme4::glmer(log_reg_formula, data = samp_dat, family = 'binomial')
  )
  
  return(list(lmer = lmer_nz, glmer = glmer_z))
  
}

# furrr with progress bar
boot_rep_with_progress_bar <- function(x,
                                       boot_lst,
                                       domain_level,
                                       boot_lin_formula,
                                       boot_log_formula,
                                       boot_pop_data,
                                       boot_truth, 
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
                             log_X = log_X,
                             estimand =)
  
  log_lst <- res |>
    map(.f = ~ .x$log)
  
  list(preds_full, log_lst)
  
}


# generate boot_pop
generate_boot_pop <- function(original_out,
                              pop_dat,
                              domain_level,
                              lin_X,
                              log_X) {
  
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
  
  # bootstrap population data
  boot_pop_data <- data.frame(
    domain = pop_dat[ , domain_level, drop = T],
    response = linear_preds * boot_data_generation_params$delta_i_star
  )
  
  boot_pop_data <- cbind(pop_dat, boot_pop_data)
  
  return(boot_pop_data)
  
}



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
  
  boot_pop_by_dom <- split(.data, f = .data$domain)
  
  design_mat_ls <-  boot_pop_by_dom |> 
    map(.f = function(.x) {
      dmat_lm <- model.matrix(~., .x[ ,lin_X, drop = FALSE])
      dmat_glm <- model.matrix(~., .x[ ,log_X, drop = FALSE])
      return(list(design_mat_lm = dmat_lm,
                  design_mat_glm = dmat_glm))
    })
  
  n_doms <- length(colnames(u_lm))
  dom_order <- names(boot_pop_by_dom)
  u_lm <- u_lm[, order(match(colnames(u_lm), dom_order))]
  u_glm <- u_glm[, order(match(colnames(u_glm), dom_order))]
  
  dom_res_wide <- generate_preds(beta_lm = beta_lm_mat,
                                 beta_glm = beta_glm_mat,
                                 u_lm = u_lm,
                                 u_glm = u_glm,
                                 design_mats = design_mat_ls,
                                 J = n_doms,
                                 estimand = estimand)
  
  
  
  
  truth_ordered <- truth[order(match(truth$domain, dom_order)), ]
  truth_vec <- truth_ordered$domain_est
  
  mse <- (dom_res_wide - truth_vec)^2 |>
    rowMeans(na.rm = TRUE)
  

  res_doms <- data.frame(
    domain = truth_ordered$domain,
    mse = mse
  )
  
  return(res_doms)
  
  
}



# base version of dplyr::slice_sample
slice_samp <- function(.data, n, replace = TRUE) {
  .data[sample(nrow(.data), n, replace = replace),]
}


# base version of stringr::str_extract_all
str_extract_all_base <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string))
}

mod_param_fmt <- function(.fit, ref = NULL) {
  
  if (!is.null(.fit)) {
    .lmer <- .fit$result$lmer
    .glmer <- .fit$result$glmer
    
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
    lm_terms <- ref$.lm[!stringr::str_detect(ref$.lm, "\\|")]
    glm_terms <- ref$.glm[!stringr::str_detect(ref$.glm, "\\|")]
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


# bootstrap rep helper
boot_rep <- function(boot_samp,
                     domain_level,
                     boot_lin_formula,
                     boot_log_formula) {
  
  # capture warnings and messages silently when bootstrapping
  fit_zi_capture <- capture_all(fit_zi)
  
  boot_samp_fit  <- tryCatch(
    {
      out <- fit_zi_capture(boot_samp,
                            boot_lin_formula,
                            boot_log_formula,
                            domain_level)
      
      ps <- mod_param_fmt(out)
      return(list(params = ps, log = out$log))
      
    },
    error = function(cond) {
      
      doms <- unique(boot_samp[[domain_level]])
      lm_cfs <- labels(terms(boot_lin_formula))
      glm_cfs <- labels(terms(boot_lin_formula))
      
      ps <- mod_param_fmt(.fit = NULL,
                          ref = list(d = doms,
                                     .lm = lm_cfs,
                                     .glm = glm_cfs))
      
      return(list(params = ps, log = cond))
      
    }
  )
  
  return(list(params_ls = boot_samp_fit$params, log = boot_samp_fit$log))
  
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