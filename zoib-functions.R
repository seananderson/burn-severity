get_dat <- function(response, main_predictor,
  file = "RockiesBurnSev_clean_20180216.csv",
  max_predictor = 2500) {

  response <- enquo(response)
  main_predictor <- enquo(main_predictor)

  d <- suppressMessages(read_csv(file))

  if ("pre_NBR" %in% names(d)) # col name was changed
    d <- rename(d, preNBR = pre_NBR)

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

  d <- dplyr::mutate_at(d,
    c("lat",  "heatload",  "slope", "prenbr", "qmd_laf", "qmd_all",
      "ba_ha", "stems_ha", "redgraygreenstage_ba_p"),
    function(x) arm::rescale(x))
  d
}

prep_stan_dat_oib <- function(d, predictors = "xscaled") {

  # in case some missing now (e.g. cross-validation):
  d <- d %>% mutate(
    group = as.factor(as.character(group)),
    group_id = as.integer(group)) %>%
    arrange(group_id, x)
  dp <- filter(d, !is.na(yp))

  # N0 <- nrow(d)
  N1 <- nrow(d)
  Np <- nrow(dp)

  # f0 <- as.formula(paste("y0 ~", paste(predictors, collapse = " + ")))
  f1 <- as.formula(paste("y1 ~", paste(predictors, collapse = " + ")))
  fp <- as.formula(paste("y ~", paste(predictors, collapse = " + ")))
  # mm0 <- as.matrix(model.matrix(f0, data = d))
  mm1 <- as.matrix(model.matrix(f1, data = d))
  mmp <- as.matrix(model.matrix(fp, data = dp))

  sd <- list(
    # N0 = N0,
    N1 = N1,
    Np = Np,
    J01 = ncol(mm1),
    Jp = ncol(mmp),
    # X0_ij = mm0,
    X1_ij = mm1,
    Xp_ij = mmp,
    # y0_i = d$y0,
    y1_i = d$y1,
    yp_i = dp$y,

    Ng = length(unique(as.character(d$group))),
    # group_id0 = d$group_id,
    group_id1 = d$group_id,
    group_idp = dp$group_id)

  list(stan_dat = sd, f = fp)
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

fit_model_oib <- function(d, model, predictors = "xscaled", iter = 800L, chains = 4L,
  adapt_delta = 0.8, log_lik = TRUE, ...) {
  library(rstan)
  pars <- c("b1_j", "bp_j", "phi", "sigma_z", "z1_g", "zp_g")
  if (log_lik) pars <- c(pars, "log_lik")
  prep <- prep_stan_dat_oib(d, predictors = predictors)
  m <- sampling(sm_oib,
    data = prep$stan_dat,
    iter = iter, chains = chains,
    pars = pars,
    control = list(adapt_delta = adapt_delta, max_treedepth = 20), ...
  )
  list(model = m, data = d, f = prep$f, stan_dat = prep$stan_dat)
}

make_predictions_oib <- function(d, f, model, Npred = 200L, use_new_data = TRUE,
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

  pred_p <- plogis(mm_pred %*% t(e$bp_j[, , drop = FALSE]))
  pred_1 <- plogis(mm_pred %*% t(e$b1_j[, , drop = FALSE]))
  y_pred <- pred_1 + (1 - pred_1) * pred_p
  d_pred$y_pred <- y_pred

  pred_df <- data.frame(
    xscaled = d_pred$xscaled,
    group = d_pred$group_id,
    est = apply(y_pred, 1, quantile, probs = 0.5),
    lwr = apply(y_pred, 1, quantile, probs = 0.05),
    upr = apply(y_pred, 1, quantile, probs = 0.95))

  if (is.null(sd_x))
    pred_df$x <- d_pred$xscaled * 2 * sd(d$x) + mean(d$x)

  pred_df
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

plot_zoib_coefs <- function(obj, oib = FALSE) {

  m <- obj$m
  x_name <- paste0(" ", obj$x) # " " for fct order

  terms <- tibble(
    name = colnames(m$stan_dat$Xp_ij),
    num = as.character(seq_along(colnames(m$stan_dat$Xp_ij)))
  )

  models <- tibble(model = c("0", "1", "p"),
    model_full = factor(
      c("Pr(not 0)", "Pr(1)", "Proportion"),
      levels =
        c("Pr(not 0)", "Proportion", "Pr(1)")
    ))
  if (oib)
    models <- filter(models, model != 0) # oib vs zoib

  b <- broom::tidyMCMC(m$model, estimate.method = "median", conf.int = TRUE)
  b_inner <- broom::tidyMCMC(m$model, estimate.method = "median", conf.int = TRUE,
    conf.level = 0.5) %>%
    rename(conf.low.25 = conf.low, conf.high.75 = conf.high) %>%
    select(-estimate, -std.error)

  b <- b %>% inner_join(b_inner, by = "term") %>%
    filter(grepl("b[01p]+_j", term)) %>%
    filter(!grepl("\\[1\\]", term)) %>%
    mutate(num = gsub("b[01p]+_j\\[([0-9]+)\\]", "\\1", term)) %>%
    mutate(model = gsub("b([01p]+)_j\\[([0-9]+)\\]", "\\1", term)) %>%
    left_join(terms, by = "num") %>%
    mutate(name = gsub("xscaled", x_name, name))

  pal <- RColorBrewer::brewer.pal(3, "Set2")
  # pal <- c(pal[4], "grey50", pal[5])

  b %>%
    mutate(estimate = ifelse(model == 0,     -estimate,     estimate)) %>%
    mutate(conf.low = ifelse(model == 0,     -conf.low,     conf.low)) %>%
    mutate(conf.high = ifelse(model == 0,    -conf.high,    conf.high)) %>%
    mutate(conf.high.75 = ifelse(model == 0, -conf.high.75, conf.high.75)) %>%
    mutate(conf.low.25 = ifelse(model == 0,  -conf.low.25,  conf.low.25)) %>%
    mutate(interaction = grepl(":", name)) %>%
    inner_join(models, by = "model") %>%
    mutate(xvar = fct_relevel(name, x_name, after = Inf)) %>%
    ggplot(aes(y = estimate, x = xvar,
      colour = model_full, shape = interaction)) +
    geom_hline(yintercept = 0, lty = 2, col = "grey55") +
    geom_linerange(aes(ymin = conf.low.25, ymax = conf.high.75),
      position = position_dodge(width = 0.4), size = 1.2) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = 0.4)) +
    xlab("") +
    ggsidekick::theme_sleek() +
    coord_flip() +
    scale_color_manual(values = pal, guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c(19, 21)) +
    labs(colour = "Model", shape = "Interaction", y = "Coefficient value") +
    ggtitle(obj$y) +
    guides(shape = FALSE)
}

plot_interaction <- function(obj, int_var, int_lab, xlab, quant = c(0.1, 0.9),
  cols = c("black", "red"), title = "") {

  newdata <- obj$d
  cn <- colnames(obj$m$stan_dat$Xp_ij)
  cn <- cn[!grepl("\\(", cn)]
  cn <- cn[!grepl("\\:", cn)]
  cn <- cn[!grepl("xscaled", cn)]
  newdata[,cn] <- 0
  newdata[, "xscaled"] <- seq(min(newdata[, "xscaled"]),
    max(newdata[, "xscaled"]),
    length.out = nrow(newdata))

  newdata[,int_var] <- quantile(obj$d[, int_var, drop = TRUE], quant[[1]])
  newdata$level <- as.character(quant[[1]])
  newdata1 <- newdata

  for (i in seq(2, length(quant))) {
    newdata2 <- newdata1
    newdata2[, int_var] <- quantile(obj$d[, int_var, drop = TRUE], quant[[i]])
    newdata2$level <- as.character(quant[[i]])
    newdata <- bind_rows(newdata, newdata2)
  }

  p <- make_predictions_oib(d = newdata, f = obj$m$f,
    model = obj$m$model, re = FALSE, use_new_data = FALSE)
  p$level <- newdata$level

  names(cols) <- quant

  g <- ggplot(p, aes(x, est, ymin = lwr, ymax = upr, group = level,
    fill = level)) +
    geom_line(lwd = 1.5, aes(colour = level)) +
    geom_ribbon(alpha = 0.3) +
    ylim(0, 1) +
    ggsidekick::theme_sleek() +
    ylab("Proportion burned") +
    xlab(xlab) +
    labs(colour = paste("Quantile", int_lab),
      fill = paste("Quantile", int_lab)) +
    scale_colour_manual(values = cols, guide = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = cols, guide = guide_legend(reverse = TRUE)) +
    coord_cartesian(expand = FALSE) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0)) +
    ggtitle(title)
  g
}
