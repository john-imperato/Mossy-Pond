
functions {
  // Compute the (T, 3) matrix of forward probabilities for one individual
  matrix forward_prob(int i, // individual index
                      matrix X_surv, vector beta_phi, vector eps_phi,
                      vector logit_detect, vector lambda, vector gam_init,
                      array[] int introduced, array[] int t_intro,
                      array[] int removed, array[] int t_remove,
                      array[] int prim_idx, array[] int any_surveys,
                      array[] int J, array[,] int j_idx, array[,,] int Y,
                      int Jtot, int T) {
    // Define temporary scope vars for function
    array[3] real acc;
    array[T] vector[3] gam;
    vector[T] phi;
    matrix[T, 3] forward;
    array[3, T, 3] real ps;
    array[3, Jtot, 3] real po;
    real p;
    // s = 1 :: not recruited
    // s = 2 :: alive
    // s = 3 :: dead
    
    // Define probs of state S(t+1) | S(t)
    // First index: S(t)
    // Second index: individual
    // Third index: t
    // Fourth index: S(t + 1)
    for (t in 1 : T) {
      phi[t] = inv_logit(X_surv[i,  : ] * beta_phi + eps_phi[t]);
      ps[1, t, 3] = 0; // can't die before being alive
      ps[2, t, 1] = 0; // can't unenter population
      ps[3, t, 1] = 0;
      ps[2, t, 2] = phi[t]; // survive
      ps[2, t, 3] = 1 - phi[t]; // death
      ps[3, t, 2] = 0; // cannot un-die
      ps[3, t, 3] = 1; // dead stay dead
    }
    
    if (introduced[i]) {
      // Individual has been introduced
      // Zero probability of recruiting through t_intro - 1
      for (t in 1 : (t_intro[i] - 1)) {
        ps[1, t, 1] = 1;
        ps[1, t, 2] = 0;
        ps[1, t, 3] = 0;
      }
      
      // Timestep of introduction has Pr(recruiting) = 1
      ps[1, t_intro[i], 1] = 0;
      ps[1, t_intro[i], 2] = 1;
      ps[1, t_intro[i], 3] = 0;
      
      // To avoid NA values, fill in remaining recruitment probs (though they
      // are irrelevant for the likelihood)
      for (t in (t_intro[i] + 1) : T) {
        ps[1, t, 1] = 1;
        ps[1, t, 2] = 0;
        ps[1, t, 3] = 0;
      }
    } else {
      for (t in 1 : T) {
        ps[1, t, 1] = 1 - lambda[t];
        ps[1, t, 2] = lambda[t];
        ps[1, t, 3] = 0;
      }
    }
    
    if (removed[i]) {
      if (t_remove[i] < T) {
        ps[2, t_remove[i] + 1, 2] = 0;
        ps[2, t_remove[i] + 1, 3] = 1;
      }
    }
    
    // Observation probabilities
    for (j in 1 : Jtot) {
      p = inv_logit(logit_detect[j]);
      po[1, j, 1] = 1;
      po[1, j, 2] = 0;
      
      if (prim_idx[j] == t_intro[i]) {
        // Introductions always happen after surveys, so if an individual is
        // released on primary period t, it has a zero probability of
        // detection
        po[2, j, 1] = 1;
        po[2, j, 2] = 0;
      } else {
        po[2, j, 1] = 1 - p;
        po[2, j, 2] = p;
      }
      
      po[3, j, 1] = 1;
      po[3, j, 2] = 0;
    }
    
    for (t in 1 : T) {
      for (k in 1 : 3) {
        for (kk in 1 : 3) {
          if (t == 1) {
            acc[kk] = gam_init[kk];
          } else {
            acc[kk] = gam[t - 1, kk];
          }
          acc[kk] *= ps[kk, t, k];
          if (any_surveys[t]) {
            for (j in 1 : J[t]) {
              // Only increment the log probability with the likelihood
              // if surveys happened
              acc[kk] *= po[k, j_idx[t, j], Y[i, t, j]];
            }
          }
        }
        gam[t, k] = sum(acc);
      }
      forward[t,  : ] = gam[t,  : ]';
    }
    return forward;
  }
  
  real partial_sum_lpmf(array[] int slice_individuals, int start, int end,
                        matrix X_surv, vector beta_phi, vector eps_phi,
                        vector logit_detect, vector lambda, vector gam_init,
                        array[] int introduced, array[] int t_intro,
                        array[] int removed, array[] int t_remove,
                        array[] int prim_idx, array[] int any_surveys,
                        array[] int J, array[,] int j_idx, array[,,] int Y,
                        int Jtot, int T) {
    real loglik = 0;
    matrix[T, 3] forward;
    
    for (i in start : end) {
      forward = forward_prob(i, X_surv, beta_phi, eps_phi, logit_detect,
                             lambda, gam_init, introduced, t_intro, removed,
                             t_remove, prim_idx, any_surveys, J, j_idx, Y,
                             Jtot, T);
      
      loglik += log(sum(forward[T,  : ]));
    }
    
    return loglik;
  }
}
data {
  int<lower=1> M; // Augmented sample size
  int<lower=1> T; // # primary periods
  int<lower=1> maxJ; // Max # of secondary periods
  array[T] int<lower=0, upper=maxJ> J; // # 2ndary periods for each prim. period
  int<lower=1> Jtot; // Total number of surveys
  array[Jtot] int<lower=1, upper=T> prim_idx;
  
  // Observations
  // 0=NA, 1=not detected, 2=detected
  array[M, T, maxJ] int<lower=0, upper=2> Y;
  array[M] int<lower=0, upper=1> introduced; // Indicator for whether introduced
  array[M] int<lower=0, upper=T> t_intro; // When individuals introduced
  array[M] int<lower=0, upper=1> removed;
  array[M] int<lower=0, upper=T> t_remove;
  
  // Index order of surveys (0: NA)
  array[T, maxJ] int<lower=0, upper=Jtot> j_idx;
  array[T] int<lower=0, upper=1> any_surveys;
  
  // Fixed effects design matrices
  int<lower=1> m_detect;
  matrix[Jtot, m_detect] X_detect;
  int<lower=1> m_surv;
  matrix[M, m_surv] X_surv;
  
  int<lower=0, upper=1> any_recruitment;
  int grainsize;
}
transformed data {
  array[M] int Mseq;
  int Tm1 = T - 1; // # primary periods - 1
  vector[3] gam_init = [1, 0, 0]';
  
  for (i in 1 : M) {
    Mseq[i] = i;
  }
}
parameters {
  // Recruitment
  real alpha_lambda;
  real<lower=0> sigma_lambda;
  vector[T] eps_lambda;
  
  // Survival
  vector[m_surv] beta_phi;
  real<lower=0> sigma_phi;
  vector<multiplier=sigma_phi>[T] eps_phi;
  
  // Detection params
  vector[m_detect] beta_detect;
}
transformed parameters {
  vector[Jtot] logit_detect;
  vector<lower=0, upper=1>[T] lambda;
  
  // Probability of entering population
  lambda = any_recruitment
           * inv_logit(alpha_lambda + eps_lambda * sigma_lambda);
  
  // Probability of detection
  logit_detect = X_detect * beta_detect;
}
model {
  // Priors
  alpha_lambda ~ normal(0, 10);
  sigma_lambda ~ normal(0, 1);
  eps_lambda ~ normal(0, sigma_lambda);
  
  beta_phi ~ normal(0, 10);
  sigma_phi ~ normal(0, 1);
  eps_phi ~ normal(0, sigma_phi);
  
  beta_detect ~ normal(0, 10);
  
  // Likelihood
  if (any_recruitment) {
    for (i in 1 : M) {
      // Only include likelihood for individuals that are in population
      if (introduced[i] || removed[i] || any_surveys[t_intro[i]]) {
        target += partial_sum_lpmf(Mseq, i, i, X_surv, beta_phi, eps_phi,
                                   logit_detect, lambda, gam_init, introduced, t_intro,
                                   removed, t_remove, prim_idx, any_surveys, J, j_idx, Y,
                                   Jtot, T);
      }
    }
  }
}
generated quantities {
  array[M, T] al;
  array[M, T] d;
  array[M] alive;
  array[M] w;
  matrix[T, T] survival_probs; // Matrix to store survival probabilities

  // Calculate the survival probabilities between primary periods
  for (t in 1 : T) {
    survival_probs[t, t] = inv_logit(X_surv[t, :] * beta_phi + eps_phi[t]);
    if (t < T) {
      survival_probs[t, t + 1] = inv_logit(X_surv[t + 1, :] * beta_phi + eps_phi[t + 1]);
    }
  }

  for (i in 1 : M) {
    for (t in 2 : T) {
      al[i, t - 1] = s[i, t] == 2;
    }
    for (t in 1 : T) {
      d[i, t] = s[i, t] == al[i, t];
    }
    alive[i] = sum(al[i]);
  }

  for (t in 1 : T) {
    N[t] = sum(al[:, t]);
    B[t] = sum(d[:, t]);
  }
  for (i in 1 : M) {
    w[i] = 1 - !alive[i];
  }
  Nsuper = sum(w);
}

