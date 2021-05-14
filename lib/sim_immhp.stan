data{
  // int<lower=1> max_Nm; //maximum of number of events for each pair each window => max(unlist(lapply(return_df$event.times,length))))
  // int<lower=1> N_til; // number of pairs with interactions
  // int<lower=0,upper=max_Nm> Nm[N_til,no_observations]; //number of events for each pair => count_matrix
  // vector[max_Nm+1] time_matrix[N_til,no_observations]; // include termination time in the last entry
  // real<lower=0> max_interevent[N_til];
  int<lower=1> max_Nm; //maximum of number of events for each pair
  int<lower=1> N_til;//number of pairs who have more than "cut_off" interactions
  int<lower=1> N;
  int<lower=1,upper=N> I_fit[N_til];
  int<lower=1,upper=N> J_fit[N_til];
  int<lower=0,upper=max_Nm> Nm[N_til]; //number of events for each pair
  vector[max_Nm] event_matrix[N_til];
  vector[max_Nm] interevent_time_matrix[N_til];
}

parameters{
  vector<lower=0>[N_til] lambda0; //baseline rate for each pair
  vector<lower=0>[N_til] w_lambda;
  vector<lower=0, upper=1>[N_til] w_q1; //CTMC transition rate
  vector<lower=0, upper=1>[N_til] w_q2; //
  vector<lower=0>[N_til] alpha;
  vector<lower=0>[N_til] beta_delta;
  vector<lower=0,upper=1>[N_til] delta_1; // P(initial state = 1)
}
transformed parameters{
  vector<lower=0>[N_til] lambda1;
  vector<lower=0>[N_til] q1;
  vector<lower=0>[N_til] q2;
  vector[N_til] beta;
  
  row_vector[2] log_delta[N_til];
  
  lambda1 = (lambda0).*(1+w_lambda); 
  q2 = (lambda0).*w_q2;
  q1 = (lambda0).*w_q1;
  
  for(i in 1:N_til){
    //delta_1[i] = q2[i]/(q1[i]+q2[i]);
    log_delta[i][1] = log(delta_1[i]);
    log_delta[i][2] = log(1-delta_1[i]);
    beta[i] = alpha[i]*(1+beta_delta[i]);
  }
}
model{
  real integ; // Placeholder variable for calculating integrals
  row_vector[2] forward[max_Nm]; // Forward variables from forward-backward algorithm
  row_vector[2] forward_termination; // Forward variables at termination time
  row_vector[2] probs_1[max_Nm+1]; // Probability vector for transition to state 1 (active state)
  row_vector[2] probs_2[max_Nm+1]; // Probability vector for transition to state 2 (inactive state)
  row_vector[2] int_1[max_Nm+1]; // Integration of lambda when state transit to 1 (active state)
  row_vector[2] int_2[max_Nm+1]; // Integration of lambda when state transit to 2 (inactive state)
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
  real temp_delta_1;
  
  //priors
  w_lambda ~ lognormal(0,2);
  alpha ~ lognormal(0,1);
  beta_delta ~ gamma(1,1);
  //delta_1 ~ beta(2,2);
  w_q1 ~ beta(2,2);
  w_q2 ~ beta(2,2);
  
   for(i in 1:N_til){ //for each pair
    if(I_fit[i]!=J_fit[i]){
      interevent = interevent_time_matrix[i];

      temp_lambda0 = lambda0[i];
      temp_lambda1 = lambda1[i];
      temp_alpha = alpha[i];
      temp_beta = beta[i];
      temp_q1 = q1[i];
      temp_q2 = q2[i];
      temp_log_delta = log_delta[i];
      
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
