functions {
  vector sir(real t, vector y, vector theta) {
    real S = y[1];
    real I = y[2];
    real beta = theta[1];
    real gamma = 1 / theta[2];
    real cv = theta[3];
    real foi = beta * I * pow(fmax(S, 1e-12), 1 + cv^2);
    return to_vector([-foi, foi - gamma * I, gamma * I, foi]);
  }
}
data {
  int<lower=1> n_days;
  vector[4] y0;
  real t0;
  array[n_days-1] real ts;
  int N;
  array[n_days - 1] int<lower=0> cases;
  int n_fit;
}
parameters {
  real<lower=2, upper=6> D;
  real<lower=0> beta;
  real<lower=0, upper=4> cv;
}
model {
  vector[3] theta = to_vector({beta, D, cv});
  array[n_fit] vector[4] yf
    = ode_rk45_tol(sir, y0, t0, ts[1:n_fit], 1e-8, 1e-10, 100000, theta);
  vector[n_fit] inc;

  inc[1] = fmax(yf[1, 4] - y0[4], 1e-12);
  for (i in 2:n_fit) {
    inc[i] = fmax(yf[i, 4] - yf[i-1, 4], 1e-12);
  }

  cases[1:n_fit] ~ poisson(inc * N);

  beta ~ uniform(0, 2);
  D ~ uniform(2, 6);
  cv ~ uniform(0, 4);
}
generated quantities {
  real R0 = beta * D;
  real gamma = 1 / D;

  array[n_days] vector[4] y;
  vector<lower=0>[n_days - 1] incidence;
  array[n_days-1] int pred_cases;

  y[1] = y0;
  y[2:n_days] = ode_rk45_tol(sir, y0, t0, ts, 1e-8, 1e-10, 100000,
                             to_vector({beta, D, cv}));

  for (i in 1:(n_days-1)) {
    incidence[i] = fmax(y[i+1, 4] - y[i, 4], 1e-12);
  }
  pred_cases = poisson_rng(incidence * N);
}
