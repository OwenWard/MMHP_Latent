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
}
parameters{
  real<lower=0> lambda0; 
  vector<lower=0,upper=1>[12] f;
  real<lower=0> eta_1;
  real<lower=0> eta_2;
  real<lower=0> eta_3;
  real<lower=0> beta;
  //real<lower=0> tilde_beta;
}
model{
  real r[max_Nm+1]; // record variable for Hawkes process
  real alpha;
  real tilde_lambda0;
  //real intensity_end_window[N_til,no_observations]; //intensity by and of observation window
  vector[max_Nm+1] interevent;
  vector[max_Nm] event;
  
  lambda0 ~ lognormal(0,0.5);
  eta_1 ~ lognormal(0,1);
  eta_2 ~ lognormal(0,1);
  eta_3 ~ lognormal(0,1);
  beta ~ normal(0,10);
  //tilde_beta ~ normal(0,20);

  for(i in 1:N_til){
    alpha = exp(-eta_2*fabs(f[I_fit[i]]-f[J_fit[i]]))*f[I_fit[i]]*f[J_fit[i]]*eta_1/(1+exp(-eta_3*(f[I_fit[i]]-f[J_fit[i]])));
    for(j in 1:no_observations){
      if(Nm[i,j]==0){ // there is no event occured in this period
        // if(j==1){
        //   tilde_lambda0 = lambda0;
        // }else{
        //   tilde_lambda0 = lambda0 + (intensity_end_window[i,j-1]-lambda0)*exp(-tilde_beta*(delta_window[j]));
        // }
        tilde_lambda0 = lambda0;
        target += -tilde_lambda0*finishing_time[j];
        //intensity_end_window[i,j] = tilde_lambda0;
      }else{ // there is event occured
        interevent = interevent_time_matrix[i,j];
        // if(j==1){
        //   tilde_lambda0 = lambda0;
        // }else{
        //   tilde_lambda0 = lambda0 + (intensity_end_window[i,j-1]-lambda0)*exp(-tilde_beta*(delta_window[j]));
        // }
        tilde_lambda0 = lambda0;
        target += -tilde_lambda0*finishing_time[j] - alpha/beta*sum((1-exp(-beta*(finishing_time[j]-segment(event_matrix[i,j],1,Nm[i,j])))));
        r[1] = 0; 
        target += log(tilde_lambda0+alpha*r[1]);
        if(Nm[i,j]>1){
          for(n in 2:Nm[i,j]){
            r[n] = exp(-beta*interevent[n])*(r[n-1]+1); 
            target += log(tilde_lambda0+alpha*r[n]);
          }
        }
        //r[Nm[i,j]+1] = exp(-beta*interevent_time_matrix[i,j][Nm[i,j]+1])*(r[Nm[i,j]]+1); 
        //intensity_end_window[i,j] = tilde_lambda0+alpha*r[Nm[i,j]+1];
      }
    }
  }
}
