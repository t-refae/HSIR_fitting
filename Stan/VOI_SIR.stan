functions {  
  vector sir(real t, vector y, vector theta) {
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real C = y[4];
    
    real beta = theta[1];
    real gamma = 1/theta[2];
    
    real dS_dt = -beta * I * S;
    real dI_dt = beta * I * S - gamma * I;
    real dR_dt = gamma * I;
    real dC_dt = beta * I * S;
    
    return to_vector([dS_dt, dI_dt, dR_dt, dC_dt]);
  }
}

data {
  int<lower=2> n_days;
  vector[4] y0;   // Initial state vector [S, I, R, C] at day 1
  real t0;
  array[n_days-1] real ts;  // Time steps (days 2..n_days)
  int N;
  array[n_days-1] int<lower=0> cases;  // Observed cases per day
}

parameters {
  // real<lower=0.05, upper=0.8> gamma;
  real<lower=2, upper=6> D;
  real<lower=0> beta;
}

transformed parameters {
  array[n_days-1] vector[4] y;  // Store the solution of the ODE system
  vector[2] theta;
  theta[1] = beta;
  theta[2] = D;

  // Solve the ODE using ode_bdf with the sir function
  y = ode_bdf(sir, y0, t0, ts, theta);
  
  // Calculate incidence (new cases per day); y0 anchors the day-1 cumulative
  vector<lower=0>[n_days - 1] incidence;
  
  incidence[1] = fmax(y[1, 4] - y0[4], 1e-12);
  for (i in 2:(n_days - 1)) {
    incidence[i] = fmax(y[i, 4] - y[i-1, 4], 1e-12);
  }
}

model {
  // Likelihood: observed cases follow a Poisson distribution
  cases ~ poisson(incidence * N);
  
  // beta ~ uniform(0, 2);
  beta ~ normal(1,1);
  // gamma ~ uniform(0.05, 0.8);
  // gamma ~ normal(0.4, 0.5);
  D ~ uniform(2,6);
}

generated quantities {
  // Calculate R0 and recovery time
  real R0 = beta*D;
  // real recovery_time = 1 / gamma;
  real gamma = 1/D;

  // Predicted cases based on the incidence and Poisson distribution
  array[n_days - 1] int pred_cases;
  pred_cases = poisson_rng(incidence * N);
}
