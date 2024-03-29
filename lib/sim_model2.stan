data{
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
  vector<lower=0>[N] gamma;
  vector<lower=0>[N] zeta;
  vector<lower=0,upper=1>[N] f;
  real<lower=0> eta_1;
  real<lower=0> eta_2;
  real<lower=0> eta_3;
  real<lower=0> beta_delta;
}

transformed parameters{
  vector<lower=0>[N_til] alpha;
  real alpha_max;
  real beta;
  for(i in 1:N_til){
    alpha[i] = exp(-eta_2*fabs(f[I_fit[i]]-f[J_fit[i]]))*f[I_fit[i]]*f[J_fit[i]]*eta_1/(1+exp(-eta_3*(f[I_fit[i]]-f[J_fit[i]])));
  }
  alpha_max = max(alpha);
  beta = alpha_max*(1 + beta_delta);
}

model{
  real r[max_Nm]; // record variable for Hawkes process
  //real alpha;
  real lambda_current;
  vector[max_Nm] interevent;
  vector[max_Nm] event;
  real termination;
  
  beta_delta ~ normal(0,1);
  eta_1 ~ normal(0,1);
  eta_2 ~ normal(0,1);
  eta_3 ~ normal(0,1);
  gamma ~ double_exponential(0, 0.1);//inv_gamma(3, 0.5);
  zeta ~ double_exponential(0, 0.1);//inv_gamma(3, 0.5);
  
  for(i in 1:N_til){
    lambda_current = gamma[I_fit[i]]+zeta[J_fit[i]];
    //alpha = exp(-eta_2*fabs(f[I_fit[i]]-f[J_fit[i]]))*f[I_fit[i]]*f[J_fit[i]]*eta_1/(1+exp(-eta_3*(f[I_fit[i]]-f[J_fit[i]])));
    interevent = interevent_time_matrix[i];
    termination = event_matrix[i,Nm[i]];

    target += -lambda_current*termination - alpha[i]/beta*sum((1-exp(-beta*(termination-segment(event_matrix[i],1,Nm[i])))));
    r[1] = 0; 
    target += log(lambda_current+alpha[i]*r[1]);
    if(Nm[i]>1){
      for(n in 2:Nm[i]){
        r[n] = exp(-beta*interevent[n])*(r[n-1]+1); 
        target += log(lambda_current+alpha[i]*r[n]);
      }
    }
  }
}
