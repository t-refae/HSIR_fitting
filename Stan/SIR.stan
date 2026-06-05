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
  int<lower=1> n_days;
  vector[4] y0;
  real t0;
  array[n_days] real ts;
  int N;
  array[n_days-1] int<lower=0> cases;
}

parameters {
  real<lower=2, upper=6> D;
  real<lower=0> beta;
}

transformed parameters {
  array[n_days] vector[4] y;
  vector[2] theta;
  theta[1] = beta;
  theta[2] = D;

  y = ode_bdf(sir, y0, t0, ts, theta);
  
  vector<lower=0>[n_days - 1] incidence;
  
  for (i in 1:(n_days - 1)) {
    incidence[i] = y[i+1, 4] - y[i, 4];
  }
}

model {
  // likelihood
  cases ~ poisson(incidence * N);
  
  // priors
  beta ~ uniform(0, 2);
  D ~ uniform(2,6);
}

generated quantities {
  real R0 = beta*D;
  real gamma = 1 / D;

  array[n_days - 1] int pred_cases;
  pred_cases = poisson_rng(incidence * N);
}
