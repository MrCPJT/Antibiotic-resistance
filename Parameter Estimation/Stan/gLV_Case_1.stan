// Stan file generating gLV Stan model
// Case 1 - For generating growth rate and inter-interactivity paramaters using Mono-Culture data

// Function block:
functions {
  real[] logisticgrowth(real t,
                        real[] y,
                        real[] theta,
                        real[] x_r,
                        int[] x_i
  ) {
    real dydt[x_i[1]];
    for (i in 1:x_i[1]){
      dydt[i] = theta[1] * y[i] * (1-y[i]/theta[2]);
    }
    return dydt;
  }
}

// The input data is a vector 'y' of length 'N'.
// Data block: Information conditioned upon. Lower = Lower bound.

data {
  int<lower=1> T;           // Number of total replicates
  int<lower=0> nReplicates; // Number of replicates to sample
  real y0[nReplicates];     // Initial Population size for replicates
  real z[T,nReplicates];    // Averaging across replicates?
  real t0;                  // Initial time
  real ts[T];               // Time step
}

// Transformed data
transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1] = nReplicates;
}

// Accepted parameters to estimate (Interactivity terms [theta] and error est. [sigma])
parameters {
  real<lower=0> theta[2];
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  real y_hat[T,nReplicates];
  theta[1] ~ cauchy(1,0.01);   // Vague choice of priors
  theta[2] ~ cauchy(100000,1000); // Vague choice of priors
  sigma ~ normal(0,10000);    // Vague choice of priors
  y_hat = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    for (i in 1:nReplicates) {
      z[t,i] ~ normal(y_hat[t,i], sigma);
    }
  }
}

// Generated Quantities
generated quantities{
  real y_pred[T,nReplicates];
  real z_pred[T,nReplicates];
  y_pred = integrate_ode_rk45(logisticgrowth, y0, t0, ts, theta, x_r, x_i );
  for (t in 1:T) {
    for(i in 1:nReplicates){
      z_pred[t,i] = y_pred[t,i] + normal_rng(0,sigma);
    }
  }
}
