
```
  mu_temp = inv_logit(mup);
  A = mu_temp * phi;
  B = (1.0 - mu_temp) * phi;

  for (i in 2:J01) {
    b0_j[i] ~ normal(0, 2);
    b1_j[i] ~ normal(0, 2);
  }
  for (i in 2:Jp)
    bp_j[i] ~ normal(0, 2);

  b0_j[1] ~ normal(0, 5);
  b1_j[1] ~ normal(0, 5);
  bp_j[1] ~ normal(0, 5);

  phi ~ student_t(3, 0, 25);

  y0_i ~ bernoulli_logit(mu0);
  y1_i ~ bernoulli_logit(mu1);
  yp_i ~ beta(A, B);
  
  sigma_z ~ student_t(3, 0, 2);
```

```
ZOI:
y_pred <- (1 - pred_0) * (pred_1 + (1 - pred_1) * pred_p)

OIB:
y_pred <- pred_1 + (1 - pred_1) * pred_p
```

```
 2 * sd(d$x) + mean(d$x)
```

```
redo:
iter = 1000
chains = 4
```
max_predictor = 2500
quants <- c(0.05, 0.5, 0.95)

AUC thresholds...
