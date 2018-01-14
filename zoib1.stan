data {
  int<lower=0> N0; // rows of data
  int<lower=0> N1; // rows of data
  int<lower=0> Np; // rows of data
  int<lower=1> J; // # fixed covariates
  matrix[N0,J] X0_ij; // fixed covariate matrix for mean
  matrix[N1,J] X1_ij; // fixed covariate matrix for mean
  matrix[Np,J] Xp_ij; // fixed covariate matrix for mean
  int<lower=0,upper=1> y0_i[N0]; // vector to hold 0 vs other observations
  int<lower=0,upper=1> y1_i[N1]; // vector to hold 1 vs other observations
  vector[Np] yp_i; // vector to hold (0, 1) observations
}
parameters {
  vector[J] b0_j;
  vector[J] b1_j;
  vector[J] bp_j;
  real<lower=0> phi; // dispersion parameter
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

  mu_temp = inv_logit(mup);
  A = mu_temp * phi;
  B = (1.0 - mu_temp) * phi;
}
model {
  b0_j ~ normal(0, 5);
  b1_j ~ normal(0, 5);
  bp_j ~ normal(0, 5);
  phi ~ student_t(7, 0, 25);

  y0_i ~ bernoulli_logit(mu0);
  y1_i ~ bernoulli_logit(mu1);
  yp_i ~ beta(A, B);
}
