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
  // real<lower=0> sigma;
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
  rho ~ inv_gamma(3, 1);
  // sigma ~ inv_gamma(2, 1);
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
          //d_ij, r_ij
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
