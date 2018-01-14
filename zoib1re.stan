data {
  int<lower=0> N0; // rows of data
  int<lower=0> N1; // rows of data
  int<lower=0> Np; // rows of data
  int<lower=1> J01; // # fixed covariates disc.
  int<lower=1> Jp; // # fixed covariates cont.
  matrix[N0,J01] X0_ij; // fixed covariate matrix for mean
  matrix[N1,J01] X1_ij; // fixed covariate matrix for mean
  matrix[Np,Jp] Xp_ij; // fixed covariate matrix for mean
  int<lower=0,upper=1> y0_i[N0]; // vector to hold 0 vs other observations
  int<lower=0,upper=1> y1_i[N1]; // vector to hold 1 vs other observations
  vector[Np] yp_i; // vector to hold (0, 1) observations

  int<lower=0> Ng; // N groups
  int<lower=1,upper=Ng> group_id0[N0];
  int<lower=1,upper=Ng> group_id1[N1];
  int<lower=1,upper=Ng> group_idp[Np];
}
parameters {
  vector[J01] b0_j;
  vector[J01] b1_j;
  vector[Jp] bp_j;
  real<lower=0> phi; // dispersion parameter

  vector[Ng] z0_g; // group-level deviates
  vector[Ng] z1_g; // group-level deviates
  vector[Ng] zp_g; // group-level deviates
  real<lower=0> sigma_z[3]; // RE sigma
}
transformed parameters {
  vector[N0] mu0;
  vector[N1] mu1;
  vector[Np] mup;
  vector[Np] mu_temp;

  vector<lower=0>[Np] A; // beta dist. parameter
  vector<lower=0>[Np] B; // beta dist. parameter

  mu0 = X0_ij * b0_j;
  mu1 = X1_ij * b1_j;
  mup = Xp_ij * bp_j;

  for (i in 1:N0)
    mu0[i] = mu0[i] + z0_g[group_id0[i]];
  for (i in 1:N1)
    mu1[i] = mu1[i] + z1_g[group_id1[i]];
  for (i in 1:Np)
    mup[i] = mup[i] + zp_g[group_idp[i]];

  mu_temp = inv_logit(mup);
  A = mu_temp * phi;
  B = (1.0 - mu_temp) * phi;
}
model {
  b0_j ~ normal(0, 5);
  b1_j ~ normal(0, 5);
  bp_j ~ normal(0, 5);
  phi ~ student_t(3, 0, 25);

  z0_g ~ normal(0, 3);
  z1_g ~ normal(0, 3);
  zp_g ~ normal(0, 3);
  sigma_z ~ student_t(3, 0, 3);

  z0_g ~ normal(0, sigma_z[1]);
  z1_g ~ normal(0, sigma_z[2]);
  zp_g ~ normal(0, sigma_z[3]);

  y0_i ~ bernoulli_logit(mu0);
  y1_i ~ bernoulli_logit(mu1);
  yp_i ~ beta(A, B);
}
