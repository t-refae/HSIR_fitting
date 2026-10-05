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
  int<lower=2> n_days;
  vector[4] y0;             // Initial state vector [S, I, R, C] at day 1
  real t0;
  array[n_days-1] real ts;  // Time steps (days 2..n_days)
  int N;
  array[n_days - 1] int<lower=0> cases;
}
transformed data {
  array[0] real x_r;
  array[0] int x_i;
}
parameters {
  // real<lower=0.05, upper=0.8> gamma;
  real<lower=2, upper=6> D;
  real<lower=0> beta;
  real<lower=0, upper=4> cv;
  // real<lower=0> phi; // overdispersion parameter
}
transformed parameters{
  array[n_days-1] vector[4] y;
  vector[3] theta;
  theta[1] = beta;
  theta[2] = D;
  theta[3] = cv;
  
  y = ode_rk45_tol(sir, y0, t0, ts, 1e-8, 1e-10, 100000, theta);  
  
  vector<lower=0>[n_days - 1] incidence;
  
  incidence[1] = fmax(y[1, 4] - y0[4], 1e-12);
  for (i in 2:(n_days-1)){
    incidence[i] = fmax(y[i, 4] - y[i-1, 4], 1e-12); 
  }
}

model {
  // LL
  cases ~ poisson(incidence * N);
  // cases ~ neg_binomial_2(incidence*N, phi);
  
  // priors
  beta ~ normal(1, 1); // tighter prior on beta compared to other fitting procedures
  
  // gamma ~ uniform(0.05, 0.8);
  // gamma ~ normal(0.4, 0.5);
  D ~ uniform(2, 6);
  // cv ~ lognormal(log(0.5-(0.5^2)/2), 0.5); // mean = 0.5, skewed toward 0
  // cv ~ exponential(2); // mean=1/2, skewed toward 0
  cv ~ uniform(0,4);
  // phi ~ exponential(1);
}

generated quantities {
  real R0 = beta*D;
  // real recovery_time = 1/gamma;
  real gamma = 1/D;

  array[n_days-1] int pred_cases;
  pred_cases = poisson_rng(incidence * N);
  // for (i in 1:(n_days-1)) {
  //     pred_cases[i] = neg_binomial_2_rng(incidence[i]*N, phi);
  // }
}
