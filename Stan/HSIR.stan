functions {
  vector sir(real t, vector y, vector theta) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real C = y[4];
      
      real beta = theta[1];
      real gamma = 1/theta[2];
      real cv = theta[3];
      
      real dS_dt = -beta * I * S^(1+cv^2);
      real dI_dt =  beta * I * S^(1+cv^2) - gamma * I;
      real dR_dt =  gamma * I;
      real dC_dt =  beta * I * S^(1+cv^2);
      
      return to_vector([dS_dt, dI_dt, dR_dt, dC_dt]);
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
transformed data {
  array[0] real x_r;
  array[0] int x_i;
}
parameters {
  real<lower=2, upper=6> D;
  real<lower=0> beta;
  real<lower=0, upper=4> cv;
}

transformed parameters{
  array[n_days] vector[4] y;
  vector[3] theta;
  theta[1] = beta;
  theta[2] = D;
  theta[3] = cv;
  
  y[1] = y0;
  y[2:n_days] = ode_bdf(sir, y0, t0, ts, theta);
  
  vector<lower=0>[n_days - 1] incidence;
  
  for (i in 1:(n_days-1)){
    incidence[i] = fmax(y[i+1, 4] - y[i, 4], 1e-12);
  }
}

model {
  // LLS
  cases[1:n_fit] ~ poisson(incidence[1:n_fit] * N);

  // priors
  beta ~ uniform(0,2);
  D ~ uniform(2, 6);
  cv ~ uniform(0,4);
}

generated quantities {
  real R0 = beta*D;
  real gamma = 1/D;

  array[n_days-1] int pred_cases;
  pred_cases = poisson_rng(incidence * N);
}
