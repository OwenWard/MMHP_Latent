data{
  int<lower=0> day;
  int<lower=0> Gt[day, 12, 12];
  //real<lower=0> rho; // 1
  real<lower=0> c; // 1
  real<lower=0> sigma; // transition variance = 1
}
parameters{
  vector[12] x[day]; // latent space
  vector<lower=0>[12] delta; // radius
  real<lower=0> rho;
}
transformed parameters{
  real r[12,12];
  for(i in 1:12){
    for(j in 1:12){
      if(delta[i]>=delta[j]){
        r[i,j] = c*(delta[i]+1);
      }else{
        r[i,j] = c*(delta[j]+1);
      }
    }
  }
}
model{
  real d_ij;
  real p_l;
  rho ~ normal(0,1);
  delta ~ lognormal(0,5);
  for(t in 1:day){
    if(t==1){
      x[t] ~ normal(0,sigma);
    }else{
      // latent space transition
      for(i in 1:12){
        x[t][i] ~ normal(x[t-1][i],sigma);
      }
    }
      
    for(i in 1:12){
      for(j in 1:12){
        if(i != j){
          #d_ij, r_ij
          d_ij = fabs(x[t][i]-x[t][j]);
          if(d_ij<=r[i][j]){
            p_l = (1/(1+exp(d_ij-r[i][j]))-rho)*(1-(d_ij/r[i][j])^2)^2+rho;
          }else{
            p_l = rho;
          }
        Gt[t, i, j] ~ poisson(p_l);
        }
      }
    }
  }
}
generated quantities {
  real lambda_d[6, 12, 12];
  vector[12] x_til[6]; // latent space
  real d_ij_til;
  
  for(i in 1:12){
    x_til[1][i] = normal_rng(x[day][i], sigma);
  }
  for(k in 2:6){
    for(i in 1:12){
      x_til[k][i] = normal_rng(x_til[k-1][i], sigma);
    }
  }
  
  for(t in 1:6){
    for(i in 1:12){
      for(j in 1:12){
        if(i != j){
          #d_ij_til, r_ij
          d_ij_til = fabs(x_til[t][i]-x_til[t][j]);
          if(d_ij_til<=r[i][j]){
            lambda_d[t,i,j] = (1/(1+exp(d_ij_til-r[i][j]))-rho)*(1-(d_ij_til/r[i][j])^2)^2+rho;
          }else{
            lambda_d[t,i,j]  = rho;
          }
        }else{
          lambda_d[t,i,j] = 0;
        }
      }
    }
  }
}
