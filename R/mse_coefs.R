# helper function use to extract model coefs for bootstrap data generation
mse_coefs <- function(lmer_model, glmer_model) {

  # from lmer model
  beta_hat <- lmer_model@beta # linear model coefficients
  model_summary_df <- data.frame(summary(lmer_model)$varcor)

  sig2_mu_hat <- model_summary_df[1, ]$vcov
  sig2_eps_hat <- dplyr::filter(model_summary_df, grp == "Residual")$vcov

  # from glmer model
  alpha_1 <- glmer_model@beta

  b_i <- plm::ranef(glmer_model)[[1]][,1]
  b_domain_levels <- rownames(ranef(glmer_model)[[1]])

  return(list(
    beta_hat = beta_hat, sig2_mu_hat = sig2_mu_hat,
    sig2_eps_hat = sig2_eps_hat, alpha_1 = alpha_1,
    b_i = b_i, domain_levels = b_domain_levels))
}
