data{
  int<lower=0> n_matrix[12,12];
}
parameters{
  vector<lower=0,upper=1>[12] x;
  real<lower=0,upper=1> c;
  real<lower=0,upper=1> beta;
}
model{
  
  c ~ lognormal(0,1);
  beta ~ lognormal(0,2);
  
  for(i in 1:12){
    for(j in 1:12){
      if(i != j){
        n_matrix[i][j] ~ poisson(c*exp(-beta*(x[i]-x[j]-1)^2));
      }
    }
  }
}
