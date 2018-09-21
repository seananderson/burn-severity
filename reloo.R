reloo <- function(m, d, loo_x, obs, ...) {
  J <- length(obs)
  lls <- vector("list", J)
  message(
    J, " problematic observation(s) found.",
    "\nModel will be refit ", J, " times."
  )

  for (j in 1:J) {
    message(
      "\nFitting model ", j, " out of ", J,
      " (leaving out observation ", obs[j], ")"
    )
    omitted <- obs[j]

    fit_j <- fit_model(d = d[d$group != x, , drop = FALSE], model = sm,
         predictors = "xscaled", iter = 500L, chains = 1L, log_lik = FALSE)


      update(
        x,
        data = d[-omitted, , drop = FALSE],
        subset = rep(TRUE, nrow(d) - length(omitted)),
        evaluate = FALSE,
        refresh = 0,
        open_progress = FALSE
      )
    fit_j_call$subset <- eval(fit_j_call$subset)
    fit_j_call$data <- eval(fit_j_call$data)
    if (!is.null(getCall(x)$offset)) {
      fit_j_call$offset <- x$offset[-omitted]
    }
    capture.output(
      fit_j <- suppressWarnings(eval(fit_j_call))
    )

    lls[[j]] <-
      log_lik.stanreg(
        fit_j,
        newdata = d[omitted, , drop = FALSE],
        offset = x$offset[omitted],
        newx = get_x(x)[omitted, , drop = FALSE],
        newz = x$z[omitted, , drop = FALSE], # NULL other than for some stan_betareg models
        stanmat = as.matrix.stanreg(fit_j)
      )
  }

  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- unlist(lapply(lls, log_mean_exp))

  # compute \hat{lpd}_j for each of the held out observations (using log-lik
  # matrix from full posterior, not the leave-one-out posteriors)
  ll_x <- log_lik(
    object = x,
    newdata = d[obs, , drop = FALSE],
    offset = x$offset[obs]
  )
  hat_lpd <- apply(ll_x, 2, log_mean_exp)

  # compute effective number of parameters
  p_loo <- hat_lpd - elpd_loo

  # replace parts of the loo object with these computed quantities
  sel <- c("elpd_loo", "p_loo", "looic")
  loo_x$pointwise[obs, sel] <- cbind(elpd_loo, p_loo, -2 * elpd_loo)
  loo_x$estimates[sel, "Estimate"] <- with(loo_x, colSums(pointwise[, sel]))
  loo_x$estimates[sel, "SE"] <- with(loo_x, {
    N <- nrow(pointwise)
    sqrt(N * apply(pointwise[, sel], 2, var))
  })
  loo_x$diagnostics$pareto_k[obs] <- NA

  return(loo_x)
}
