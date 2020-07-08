data{
  int<lower=1> N_til;//number of pairs => nrow(unique_pairs_df)
  int<lower=1> no_observations;//total number of observation windows -> max(return_df$observe.id)
  int<lower=1,upper=12> I_fit[N_til];
  int<lower=1,upper=12> J_fit[N_til];
  int<lower=1> max_Nm; //maximum of number of events for each pair each window => max(unlist(lapply(return_df$event.times,length))))
int<lower=0,upper=max_Nm> Nm[N_til,no_observations]; //number of events for each pair => count_matrix
vector[max_Nm+1] interevent_time_matrix[N_til,no_observations]; // include termination time difference in the last entry
vector[max_Nm] event_matrix[N_til,no_observations]; // event times in each observation window
//real<lower=0> delta_window[no_observations]; //length of non-observation period 
real<lower=0> finishing_time[no_observations]; //for each pair, each observation window, what is the finishing time
int<lower=0,upper=12> alpha_id;
vector<lower=0,upper=1>[N_til] delta_1;
}
parameters{
  //real<lower=0.01> lambda0; //baseline rate for each pair
  vector<lower=0>[12] gamma;
  vector<lower=0>[12] zeta;
  //real<lower = 0> scale;
  real<lower=0> w_lambda;
  vector<lower=0,upper=1>[12] f;
  real<lower=0> eta_1;
  real<lower=0> eta_2;
  real<lower=0> eta_3;
  real<lower=0> beta;
  //real<lower=0,upper=1> delta_1; // P(initial state = 1)
  //real<lower=0> tilde_beta;
  //vector<lower=0,upper=1>[N_til] delta_1; // P(initial state = 1)
}
transformed parameters{
  vector<lower=0>[N_til] lambda1;
  vector<lower=0.01>[N_til] lambda0;
  //vector<lower=0,upper=1>[N_til] delta_1; // P(initial state = 1)
  vector[N_til] q1; // P(initial state = 1)
  vector[N_til] q2; // P(initial state = 1)
  vector[N_til] alpha; // P(initial state = 1)
  row_vector[2] log_delta[N_til];
  
  
  for(i in 1:N_til){
    lambda0[i] = gamma[I_fit[i]]+zeta[J_fit[i]];  
    lambda1[i] = lambda0[i]*(1+w_lambda);
    alpha[i] = exp(-eta_2*fabs(f[I_fit[i]]-f[J_fit[i]]))*f[I_fit[i]]*f[J_fit[i]]*eta_1;///(1+exp(-eta_3*(f[I_fit[i]]-f[J_fit[i]])));//eta_1*f[I_fit[i]]*f[J_fit[i]];
    //f_max = (f[I_fit[i]]+f[J_fit[i]]+fabs(f[I_fit[i]]-f[J_fit[i]]))/2;
    //f_min = (f[I_fit[i]]+f[J_fit[i]]-fabs(f[I_fit[i]]-f[J_fit[i]]))/2;
    q1[i] = exp(-eta_3*f[I_fit[i]]);//exp(-(f_max+eta_2*(f[I_fit[i]]-f[J_fit[i]])));
    q2[i] = exp(-eta_3*f[J_fit[i]]);//exp((f_min-eta_3*(f[I_fit[i]]+f[J_fit[i]])));
    //delta_1[i] = q2[i]/(q1[i]+q2[i]);
    log_delta[i][1] = log(delta_1[i]);
    log_delta[i][2] = log(1-delta_1[i]);
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
  real R[max_Nm+1]; // record variable for Hawkes process
  vector[max_Nm+1] interevent;
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
  //lambda0 ~ gamma(0.1,1); // prior has mean 0.1, var 0.1
  gamma ~ double_exponential(0,0.1); // sigma smaller, sparser
  zeta ~ double_exponential(0,0.1);
  w_lambda ~ lognormal(0,2);
  eta_1 ~ lognormal(0,1);
  eta_2 ~ lognormal(0,1);
  eta_3 ~ lognormal(0,1);
  beta ~ lognormal(0,2);
  f[alpha_id] ~ normal(1,0.05);
  
  for(i in 1:N_til){ //for each pair
    temp_lambda0 = lambda0[i];
    temp_lambda1 = lambda1[i];
    temp_alpha = alpha[i];
    temp_beta = beta;
    temp_q1 = q1[i];
    temp_q2 = q2[i];
    temp_log_delta = log_delta[i];
    
    for(j in 1:no_observations){ // for each observation period
      interevent = interevent_time_matrix[i,j];
      if(Nm[i,j]==0){ // there is no event occured in this period
        //--- prepare for forward calculation
        probs_1[1][1] = log(temp_q2/(temp_q1+temp_q2)+temp_q1/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[1])); //1->1
        probs_2[1][2] = log(temp_q1/(temp_q1+temp_q2)+temp_q2/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[1])); //2->2
        probs_1[1][2] = log1m_exp(probs_2[1][2]); //2->1
        probs_2[1][1] = log1m_exp(probs_1[1][1]); //1->2
        R[1] = 0;
        K0 = exp(-(temp_q1+temp_q2)*interevent[1]);
        K1 = (1-exp(-(temp_q1+temp_q2)*interevent[1]))/(temp_q1+temp_q2);
        K2 = (1-exp(-(temp_q1+temp_q2)*interevent[1]))/(temp_q1+temp_q2);
        int_1[1][1] = ((temp_q2^2*temp_lambda1+temp_q2*temp_q1*temp_lambda0)*interevent[1] +
                         (temp_q1^2*temp_lambda1+temp_q2*temp_q1*temp_lambda0)*K0*interevent[1] +
                         (temp_lambda1-temp_lambda0)*temp_q2*temp_q1*K1 + (temp_lambda1-temp_lambda0)*temp_q2*temp_q1*K2)/(temp_q1+temp_q2)^2/exp(probs_1[1][1]); //1->1
        int_1[1][2] = ((temp_q2^2*temp_lambda1+temp_lambda0*temp_q1*temp_q2)*interevent[1] -
                         (temp_lambda1*temp_q1*temp_q2+temp_lambda0*temp_q2^2)*K0*interevent[1] +
                         (temp_lambda0-temp_lambda1)*temp_q2^2*K1 + (temp_lambda1-temp_lambda0)*temp_q1*temp_q2*K2)/(temp_q1+temp_q2)^2/exp(probs_1[1][2]); //2->1
        int_2[1][1] = ((temp_q1*temp_q2*temp_lambda1+temp_q1^2*temp_lambda0)*interevent[1] -
                         (temp_q1^2*temp_lambda1+temp_q1*temp_q2*temp_lambda0)*K0*interevent[1] +
                         (temp_lambda1-temp_lambda0)*temp_q1^2*K1 + temp_q1*temp_q2*(temp_lambda0-temp_lambda1)*K2)/(temp_q1+temp_q2)^2/exp(probs_2[1][1]); //1->2
        int_2[1][2] = ((temp_q1*temp_q2*temp_lambda1+temp_lambda0*temp_q1^2)*interevent[1] +
                         (temp_q1*temp_q2*temp_lambda1+temp_lambda0*temp_q2^2)*K0*interevent[1] +
                         (temp_lambda0-temp_lambda1)*temp_q1*temp_q2*K1 + (temp_lambda0-temp_lambda1)*temp_q1*temp_q2*K2)/(temp_q1+temp_q2)^2/exp(probs_2[1][2]); //2->2
        
        forward_termination[1] = log_sum_exp(temp_log_delta + probs_1[1] - int_1[1]);
        forward_termination[2] = log_sum_exp(temp_log_delta + probs_2[1] - int_2[1]);
        target += log_sum_exp(forward_termination);
        //target += -temp_lambda0*interevent[1]*temp_delta_1-temp_lambda1*interevent[1]*(1-temp_delta_1);
      }else{ 
        // --- log probability of Markov transition logP_ij(t)
        for(n in 1:(Nm[i,j]+1)){
          probs_1[n][1] = log(temp_q2/(temp_q1+temp_q2)+temp_q1/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[n])); //1->1
          probs_2[n][2] = log(temp_q1/(temp_q1+temp_q2)+temp_q2/(temp_q1+temp_q2)*exp(-(temp_q1+temp_q2)*interevent[n])); //2->2
          probs_1[n][2] = log1m_exp(probs_2[n][2]); //2->1
          probs_2[n][1] = log1m_exp(probs_1[n][1]); //1->2
        }
        
        // --- R for Hawkes
        R[1] = 0;
        for(n in 2:(Nm[i,j]+1)){
          R[n] = exp(-temp_beta*interevent[n])*(R[n-1]+1);
        }
        //intensity_end_window[i,j] = tilde_lambda1+temp_alpha*R[Nm[i,j]+1];
        // Integration of lambda
        for(n in 1:(Nm[i,j]+1)){
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
        
        if(Nm[i,j]>1){
          for(n in 2:Nm[i,j]){
            forward[n][1] = log_sum_exp(forward[n-1] + probs_1[n] - int_1[n]) + log(temp_lambda1+temp_alpha*R[n]);
            forward[n][2] = log_sum_exp(forward[n-1] + probs_2[n] - int_2[n]) + log(temp_lambda0);
          }
        }
        
        forward_termination[1] = log_sum_exp(forward[Nm[i,j]] + probs_1[Nm[i,j]+1] - int_1[Nm[i,j]+1]);
        forward_termination[2] = log_sum_exp(forward[Nm[i,j]] + probs_2[Nm[i,j]+1] - int_2[Nm[i,j]+1]);
        
        target += log_sum_exp(forward_termination);
      }
    }
  }
}
