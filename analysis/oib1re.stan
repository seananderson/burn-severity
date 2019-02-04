data {
  // int<lower=0> N0; // rows of data
  int<lower=0> N1; // rows of data
  int<lower=0> Np; // rows of data
  int<lower=1> J01; // # fixed covariates disc.
  int<lower=1> Jp; // # fixed covariates cont.
  // matrix[N0,J01] X0_ij; // fixed covariate matrix for mean
  matrix[N1,J01] X1_ij; // fixed covariate matrix for mean
  matrix[Np,Jp] Xp_ij; // fixed covariate matrix for mean
  // int<lower=0,upper=1> y0_i[N0]; // vector to hold 0 vs other observations
  int<lower=0,upper=1> y1_i[N1]; // vector to hold 1 vs other observations
  vector[Np] yp_i; // vector to hold (0, 1) observations

  int<lower=0> Ng; // N groups
  // int<lower=1,upper=Ng> group_id0[N0];
  int<lower=1,upper=Ng> group_id1[N1];
  int<lower=1,upper=Ng> group_idp[Np];
}
parameters {
  // vector[J01] b0_j;
  vector[J01] b1_j;
  vector[Jp] bp_j;
  real<lower=0> phi; // dispersion parameter

  // vector[Ng] z0_g_cent; // group-level deviates
  vector[Ng] z1_g_cent; // group-level deviates
  vector[Ng] zp_g_cent; // group-level deviates
  real<lower=0> sigma_z[2]; // RE sigma
}
transformed parameters {
  // vector[N0] mu0;
  vector[N1] mu1;
  vector[Np] mup;
  vector[Np] mu_temp;

  vector<lower=0>[Np] A; // beta dist. parameter
  vector<lower=0>[Np] B; // beta dist. parameter

  // vector[Ng] z0_g; // group-level deviates
  vector[Ng] z1_g; // group-level deviates
  vector[Ng] zp_g; // group-level deviates

  // mu0 = X0_ij * b0_j;
  mu1 = X1_ij * b1_j;
  mup = Xp_ij * bp_j;

  // z0_g = sigma_z[1] * z0_g_cent;
  z1_g = sigma_z[1] * z1_g_cent;
  zp_g = sigma_z[2] * zp_g_cent;

  // for (i in 1:N0)
    // mu0[i] = mu0[i] + z0_g[group_id0[i]];
  for (i in 1:N1)
    mu1[i] = mu1[i] + z1_g[group_id1[i]];
  for (i in 1:Np)
    mup[i] = mup[i] + zp_g[group_idp[i]];

  mu_temp = inv_logit(mup);
  A = mu_temp * phi;
  B = (1.0 - mu_temp) * phi;
}
model {

  for (i in 2:J01) {
    // b0_j[i] ~ normal(0, 2);
    b1_j[i] ~ normal(0, 2);
  }
  for (i in 2:Jp)
    bp_j[i] ~ normal(0, 2);

  // b0_j[1] ~ normal(0, 5);
  b1_j[1] ~ normal(0, 5);
  bp_j[1] ~ normal(0, 5);

  phi ~ student_t(3, 0, 25);

  // z0_g_cent ~ normal(0, 1);
  z1_g_cent ~ normal(0, 1);
  zp_g_cent ~ normal(0, 1);
  sigma_z ~ student_t(3, 0, 2);

  // y0_i ~ bernoulli_logit(mu0);
  y1_i ~ bernoulli_logit(mu1);
  yp_i ~ beta(A, B);
}
generated quantities {
  // log_lik is for use with the loo package
  // vector[N0] log_lik0;
  vector[N1] log_lik1;
  vector[Np] log_likp;
  vector[N1] log_lik_temp;
  vector[Np+N1] log_lik;

  // for (i in 1:N0)
    // log_lik0[i] = bernoulli_logit_lpmf(y0_i[i] | mu0[i]);
  for (i in 1:N1)
    log_lik1[i] = bernoulli_logit_lpmf(y1_i[i] | mu1[i]);
  for (i in 1:Np)
    log_likp[i] = beta_lpdf(yp_i[i] | A[i], B[i]);

  // log_lik_temp = append_row(log_lik0, log_lik1);
  log_lik = append_row(log_lik1, log_likp);
}
