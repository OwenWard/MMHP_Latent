data{
  int<lower=1> max_Nm; //maximum of number of events for each pair
  int<lower=1> N_til;//number of pairs who have more than "cut_off" interactions
  int<lower=1> N;
  int<lower=1,upper=N> I_fit[N_til];
  int<lower=1,upper=N> J_fit[N_til];
  int<lower=0,upper=max_Nm> Nm[N_til]; //number of events for each pair
  // int<lower=0,upper=N> alpha_id;
  vector[max_Nm] event_matrix[N_til];
  vector[max_Nm] interevent_time_matrix[N_til];
}
parameters{
  //real<lower=0> lambda0; //baseline rate for each pair
  real<lower=0> r_lambda1;
  vector<lower=0,upper=1>[N] f;
  vector<lower=0>[N] gamma;
  vector<lower=0>[N] zeta;
  real<lower=0> eta_1;
  real<lower=0> eta_2;
  real<lower=0> eta_3;
  real<lower=0> beta_delta;
}
transformed parameters{
  vector<lower=0>[N_til] lambda0;
  vector<lower=0>[N_til] lambda1;
  row_vector[2] log_delta;
  vector[N_til] q1; // P(initial state = 1)
  vector[N_til] q2; // P(initial state = 1)
  vector[N_til] alpha; // P(initial state = 1)
  real alpha_max;
  real beta;

  // lambda1 = lambda0*(1+r_lambda1);
  log_delta[1] = log(0.5);
  log_delta[2] = log(0.5);

  for(i in 1:N_til){
    lambda0[i] = gamma[I_fit[i]]+zeta[J_fit[i]];
    lambda1[i] = lambda0[i]*(1+r_lambda1);
    alpha[i] = exp(-eta_2*fabs(f[I_fit[i]]-f[J_fit[i]]))*f[I_fit[i]]*f[J_fit[i]]*eta_1;
    q1[i] = exp(-eta_3*f[I_fit[i]]);
    q2[i] = exp(-eta_3*f[J_fit[i]]);
  }
  alpha_max = max(alpha);
  beta = alpha_max*(1+beta_delta);
}
model{
  row_vector[2] forward[max_Nm]; // Forward variables from forward-backward algorithm
  row_vector[2] probs_1[max_Nm]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm]; // Probability vector for transition to state 2 (inactive state)
  row_vector[2] int_1[max_Nm]; // Integration of lambda when state transit to 1 (active state)
  row_vector[2] int_2[max_Nm]; // Integration of lambda when state transit to 2 (inactive state)
  real R[max_Nm]; // record variable for Hawkes process
  vector[max_Nm] interevent;
  real K0;
  real K1;
  real K2;
  real K3;
  real K4;
  real K5;
  real temp_lambda0;
  real temp_lambda1;
  real temp_alpha;
  real temp_beta;
  real temp_q1;
  real temp_q2;
  row_vector[2] temp_log_delta;
  
  //priors
  //lambda0 ~ lognormal(0,2);
  gamma ~ double_exponential(0, 1);//inv_gamma(3,0.5);
  zeta ~ double_exponential(0, 1);//inv_gamma(3,0.5);
  r_lambda1 ~ lognormal(0,2);
  beta_delta ~ lognormal(0,2);
  eta_1 ~ lognormal(0,1);
  eta_2 ~ lognormal(0,1);
  eta_3 ~ lognormal(0,1);//normal(5,1);
  f[alpha_id] ~ normal(1,0.1);
  
  
  for(i in 1:N_til){ //for each pair
    if(I_fit[i]!=J_fit[i]){
      interevent = interevent_time_matrix[i];

      temp_lambda0 = lambda0[i];
      temp_lambda1 = lambda1[i];
      temp_alpha = alpha[i];
      temp_beta = beta;
      temp_q1 = q1[i];
      temp_q2 = q2[i];
      temp_log_delta = log_delta;
      
      // --- log probability of Markov transition logP_ij(t)
      for(n in 1:Nm[i]){
         probs_1[n][1] = log(temp_q2/(temp_q1+temp_q2)+temp_q1/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[n])); //1->1
         probs_2[n][2] = log(temp_q1/(temp_q1+temp_q2)+temp_q2/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[n])); //2->2
         probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
         probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
      }
      
      // --- R for Hawkes
      R[1] = 0;
      for(n in 2:Nm[i]){
    	  R[n] = exp(-temp_beta*interevent[n])*(R[n-1]+1);
      }
      //intensity_end_window[i,j] = tilde_lambda1+temp_alpha*R[Nm[i,j]+1];
      // Integration of lambda
      for(n in 1:Nm[i]){
        K0 = exp(-(temp_q1+temp_q2)*interevent[n]);
        K1 = (1-exp(-(temp_q1+temp_q2)*interevent[n]))/(temp_q1+temp_q2);
        K2 = (1-exp(-(temp_q1+temp_q2)*interevent[n]))/(temp_q1+temp_q2);
        K3 = R[n]*(exp(temp_beta*interevent[n])-1)/temp_beta;
        K4 = R[n]*(1-exp(-(temp_beta+temp_q1+temp_q2)*interevent[n]))*exp(temp_beta*interevent[n])/(temp_beta+temp_q1+temp_q2);
        K5 = R[n]*(1-exp(-(temp_q1+temp_q2-temp_beta)*interevent[n]))/(temp_q1+temp_q2-temp_beta);
        int_1[n][1] = ((temp_q2^2*temp_lambda1+temp_q2*temp_q1*temp_lambda0)*interevent[n] +
                       (temp_q1^2*temp_lambda1+temp_q2*temp_q1*temp_lambda0)*K0*interevent[n] +
                       (temp_lambda1-temp_lambda0)*temp_q2*temp_q1*K1 + (temp_lambda1-temp_lambda0)*temp_q2*temp_q1*K2 +
                       temp_alpha*K3*(temp_q2^2+temp_q1^2*K0) +
                       temp_alpha*temp_q1*temp_q2*K4 + temp_alpha*temp_q1*temp_q2*K5)/(temp_q1+temp_q2)^2/exp(probs_1[n][1]); //1->1
        int_1[n][2] = ((temp_q2^2*temp_lambda1+temp_lambda0*temp_q1*temp_q2)*interevent[n] -
                       (temp_lambda1*temp_q1*temp_q2+temp_lambda0*temp_q2^2)*K0*interevent[n] +
                       (temp_lambda0-temp_lambda1)*temp_q2^2*K1 + (temp_lambda1-temp_lambda0)*temp_q1*temp_q2*K2 +
                       temp_alpha*temp_q2*K3*(temp_q2-temp_q1*K0) -
                       temp_alpha*temp_q2^2*K4 + temp_alpha*temp_q1*temp_q2*K5)/(temp_q1+temp_q2)^2/exp(probs_1[n][2]); //2->1
        int_2[n][1] = ((temp_q1*temp_q2*temp_lambda1+temp_q1^2*temp_lambda0)*interevent[n] -
                       (temp_q1^2*temp_lambda1+temp_q1*temp_q2*temp_lambda0)*K0*interevent[n] +
                       (temp_lambda1-temp_lambda0)*temp_q1^2*K1 + temp_q1*temp_q2*(temp_lambda0-temp_lambda1)*K2 +
                       temp_alpha*temp_q1*K3*(temp_q2-temp_q1*K0) +
                       temp_alpha*temp_q1^2*K4 - temp_alpha*temp_q2*temp_q1*K5)/(temp_q1+temp_q2)^2/exp(probs_2[n][1]); //1->2
        int_2[n][2] = ((temp_q1*temp_q2*temp_lambda1+temp_lambda0*temp_q1^2)*interevent[n] +
                       (temp_q1*temp_q2*temp_lambda1+temp_lambda0*temp_q2^2)*K0*interevent[n] +
                       (temp_lambda0-temp_lambda1)*temp_q1*temp_q2*K1 + (temp_lambda0-temp_lambda1)*temp_q1*temp_q2*K2 +
                       temp_alpha*temp_q1*temp_q2*K3*(1+K0) -
                       temp_alpha*temp_q1*temp_q2*K4 - temp_alpha*temp_q1*temp_q2*K5)/(temp_q1+temp_q2)^2/exp(probs_2[n][2]); //2->2
      }
  
      //consider n = 1
      forward[1][1] = log(temp_lambda1) + log_sum_exp(probs_1[1]-int_1[1]+temp_log_delta); 
      forward[1][2] = log(temp_lambda0) + log_sum_exp(probs_2[1]-int_2[1]+temp_log_delta); 
      
      if(Nm[i]>1){
        for(n in 2:Nm[i]){
          forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n] - int_1[n]) + log(temp_lambda1+temp_alpha*R[n]);
          forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n] - int_2[n]) + log(temp_lambda0);
        }
      }
      
      target += log_sum_exp(forward[Nm[i]]);
    }
  }
}
