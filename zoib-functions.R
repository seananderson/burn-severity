get_dat <- function(response, main_predictor,
  file = "ALL_SHOBS_2010_2014_PLOTS_MASTER_w_CO_FINAL_20171029.csv",
  max_predictor = 2500) {

  response <- enquo(response)
  main_predictor <- enquo(main_predictor)

  d <- suppressMessages(read_csv(file))
  d <- select(d, !!response, !!main_predictor, FIRE,
    lat, HeatLoad, slope, preNBR, # Q2
    slope, HeatLoad, # Q3
    QMD_LAF, QMD_all, ba_ha, stems_ha, RedGrayGreenStage.BA.p) %>%
    filter(!is.na(!!response), !is.na(!!main_predictor)) %>%
    rename(y = !!response, x = !!main_predictor, group = FIRE) %>%
    filter(x <= max_predictor) %>%
    mutate(y0 = ifelse(y == 0, 1, 0),
      y1 = ifelse(y == 1, 1, 0),
      yp = ifelse(y == 0 | y == 1, NA, y),
      xscaled = arm::rescale(x),
      group = as.factor(tolower(group)),
      group_id = as.integer(group)) %>%
    arrange(group_id, x)
  names(d) <- tolower(names(d))
  names(d) <- gsub("\\.", "_", names(d))
  d
}

prep_stan_dat <- function(d, predictors = "xscaled") {

  # in case some missing now (e.g. cross-validation):
  d <- d %>% mutate(
      group = as.factor(as.character(group)),
      group_id = as.integer(group)) %>%
    arrange(group_id, x)
  dp <- filter(d, !is.na(yp))

  N0 <- nrow(d)
  N1 <- nrow(d)
  Np <- nrow(dp)

  f0 <- as.formula(paste("y0 ~", paste(predictors, collapse = " + ")))
  f1 <- as.formula(paste("y1 ~", paste(predictors, collapse = " + ")))
  fp <- as.formula(paste("y ~", paste(predictors, collapse = " + ")))
  mm0 <- as.matrix(model.matrix(f0, data = d))
  mm1 <- as.matrix(model.matrix(f1, data = d))
  mmp <- as.matrix(model.matrix(fp, data = dp))

  sd <- list(
    N0 = N0,
    N1 = N1,
    Np = Np,
    J01 = ncol(mm0),
    Jp = ncol(mmp),
    X0_ij = mm0,
    X1_ij = mm1,
    Xp_ij = mmp,
    y0_i = d$y0,
    y1_i = d$y1,
    yp_i = dp$y,

    Ng = length(unique(as.character(d$group))),
    group_id0 = d$group_id,
    group_id1 = d$group_id,
    group_idp = dp$group_id)

  list(stan_dat = sd, f = fp)
}

fit_model <- function(d, model, predictors = "xscaled", iter = 800L, chains = 4L,
  adapt_delta = 0.8, log_lik = TRUE, ...) {
  library(rstan)

  pars <- c("b0_j", "b1_j", "bp_j", "phi", "sigma_z", "z0_g", "z1_g", "zp_g")
  if (log_lik) pars <- c(pars, "log_lik")
  prep <- prep_stan_dat(d, predictors = predictors)
  m <- sampling(model,
    data = prep$stan_dat,
    iter = iter, chains = chains,
    pars = pars,
    control = list(adapt_delta = adapt_delta, max_treedepth = 20), ...
  )
  list(model = m, data = d, f = prep$f, stan_dat = prep$stan_dat)
}

make_predictions <- function(d, f, model, Npred = 200L, use_new_data = TRUE,
  re = TRUE, sd_x = NULL, mean_x = NULL) {

  if (use_new_data) {
    d_pred <- expand.grid(
      xscaled = seq(min(d$xscaled), max(d$xscaled), length.out = Npred),
      group_id = unique(d$group_id), y = 1)
  } else {
    d_pred <- d
    Npred
  }
  e <- rstan::extract(model)

  mm_pred <- as.matrix(model.matrix(f, data = d_pred))

  pred_p <- plogis(mm_pred %*% t(e$bp_j))
  # plot(d_pred$xscaled, apply(pred_p, 1, median), ylim = c(0, 1))
  pred_0 <- plogis(mm_pred %*% t(e$b0_j))
  # plot(d_pred$xscaled, apply(pred_0, 1, median), ylim = c(0, 1))
  pred_1 <- plogis(mm_pred %*% t(e$b1_j))
  # plot(d_pred$xscaled, apply(pred_1, 1, median), ylim = c(0, 1))
  y_pred <- (1 - pred_0) * (pred_1 + (1 - pred_1) * pred_p)
  d_pred$y_pred <- y_pred
  if (re) {
    pred_p_re <- mm_pred %*% t(e$bp_j)
    pred_p_re <- apply(pred_p_re, 1, median)
    zp_g <- apply(e$zp_g, 2, median)
    # zp_g <- tibble(zp_g = zp_g, group_id = seq_along(zp_g))

    pred_0_re <- mm_pred %*% t(e$b0_j)
    pred_0_re <- apply(pred_0_re, 1, median)
    z0_g <- apply(e$z0_g, 2, median)
    # z0_g <- tibble(z0_g = z0_g, group_id = seq_along(z0_g))

    pred_1_re <- mm_pred %*% t(e$b1_j)
    pred_1_re <- apply(pred_1_re, 1, median)
    z1_g <- apply(e$z1_g, 2, median)
    # z1_g <- tibble(z1_g = z1_g, group_id = seq_along(z1_g))
    #
    # d_pred <- inner_join(d_pred, zp_g, by = "group_id")
    # d_pred <- inner_join(d_pred, z0_g, by = "group_id")
    # d_pred <- inner_join(d_pred, z1_g, by = "group_id")
    #
    # d_pred <- mutate(d_pred,
    #   pred_p_re = pred_p_re + zp_g,
    #   pred_0_re = pred_0_re + z0_g,
    #   pred_1_re = pred_1_re + z1_g
    # )

    d_pred$pred_p_re <- plogis(pred_p_re + rep(zp_g, each = Npred))
    d_pred$pred_0_re <- plogis(pred_0_re + rep(z0_g, each = Npred))
    d_pred$pred_1_re <- plogis(pred_1_re + rep(z1_g, each = Npred))
    d_pred$y_pred_re <- (1 - d_pred$pred_0_re) *
      (d_pred$pred_1_re + (1 - d_pred$pred_1_re) * d_pred$pred_p_re)
  }

  pred_df <- data.frame(
    xscaled = d_pred$xscaled,
    group = d_pred$group_id,
    est = apply(y_pred, 1, quantile, probs = 0.5),
    lwr = apply(y_pred, 1, quantile, probs = 0.025),
    upr = apply(y_pred, 1, quantile, probs = 0.975))

  if (re)
    pred_df$est_re = d_pred$y_pred_re

  if (is.null(sd_x))
    pred_df$x <- d_pred$xscaled * 2 * sd(d$x) + mean(d$x)

  pred_df
}

make_pred_plot <- function(pred_df, raw_data, re = TRUE) {
  g <- ggplot(pred_df, aes(x, est, ymin = lwr, ymax = upr)) +
    geom_point(data = d, aes(x, y, colour = group),
      inherit.aes = FALSE, alpha = 0.6) +
    geom_line(lwd = 1.5) +
    geom_ribbon(alpha = 0.5, fill = "grey20") +
    ylim(0, 1) +
    ggsidekick::theme_sleek() +
    ylab("Proportion burned") +
    xlab("Predictor value") +
    labs(colour = "Fire location")
  if (re)
    g <- g + geom_line(aes(y = est_re, group = as.factor(group)),
      col = "grey50", lty = 1, lwd = 0.2)
  g
}

# https://gist.github.com/traversc/1446ebe1dcc2d84bccdca781d4f1fa2a
fastROC <- function(probs, class) {
  class_sorted <- class[order(probs, decreasing = TRUE)]
  TPR <- cumsum(class_sorted) / sum(class)
  FPR <- cumsum(class_sorted == 0) / sum(class == 0)
  data.frame(tpr = TPR, fpr = FPR)
}

# https://gist.github.com/traversc/1446ebe1dcc2d84bccdca781d4f1fa2a
fastAUC <- function(probs, class) {
  x <- probs
  y <- class
  x1 <- x[y == 1]
  n1 <- length(x1)
  x2 <- x[y == 0]
  n2 <- length(x2)
  r <- rank(c(x1, x2))
  (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / n1 / n2
}

make_predictions_dat <- function(model, d, f) {
  e <- rstan::extract(model)
  mm <- as.matrix(model.matrix(f, data = d))

  pred_p <- plogis(mm %*% t(e$bp_j))
  pred_0 <- plogis(mm %*% t(e$b0_j))
  pred_1 <- plogis(mm %*% t(e$b1_j))
  y_pred <- (1 - pred_0) * (pred_1 + (1 - pred_1) * pred_p)
  apply(y_pred, 1, median)
}

make_roc <- function(d, predictions, return_plot = TRUE) {
  thresh <- seq(0.2, 0.8, length.out = 8)
  auc_thresh <- function(thresh = 0.5) {
    d$outcome <- ifelse(d$y > thresh, 1, 0)
    fastAUC(class = d$outcome, probs = predictions)
  }
  a <- sapply(thresh, auc_thresh)

  roc_thresh <- function(thresh = 0.5) {
    d$outcome <- ifelse(d$y > thresh, 1, 0)
    fastROC(class = d$outcome, probs = predictions)
  }

  lab <- paste0("AUC range = ", round(min(a), 2), " - ", round(max(a), 2))
  lab1 <- paste0("mean AUC = ", round(mean(a), 2))
  r <- plyr::ldply(thresh, roc_thresh)
  r$thresh <- rep(thresh, each = length(predictions))
  g <- ggplot(r, aes(fpr, tpr)) +
    geom_line(aes(group = as.factor(thresh), colour = thresh), alpha = 0.9) +
    scale_color_gradient2(mid = "grey75", midpoint = 0.5) +
    geom_text(data = data.frame(fpr = 0.75, tpr = 0.1), label = lab) +
    geom_text(data = data.frame(fpr = 0.75, tpr = 0.15), label = lab1) +
    ggsidekick::theme_sleek()
  if (return_plot)
    return(g)
  else
    return(a)
}
