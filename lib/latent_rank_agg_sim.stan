data{
  int<lower=1> n;
  int<lower=0> n_matrix[n,n];
}
parameters{
  vector<lower=0,upper=1>[n] x;
  real<lower=0,upper=1> c;
  real<lower=0,upper=1> beta;
}
model{
  
  c ~ lognormal(0,1);
  beta ~ lognormal(0,2);
  
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        n_matrix[i][j] ~ poisson(c*exp(-beta*(x[i]-x[j]-1)^2));
      }
    }
  }
}
