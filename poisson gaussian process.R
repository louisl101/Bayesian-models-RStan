stanlist<-list(N=nrow(d),
               M=10,
               TN=d$population,
               TP=d$logpop,
               cyano=d$total_tools,
               site=d$society,
               Dmat=Dmat
)

code<-'
data {
  int<lower=1> N;
  int<lower=1> M;
  int cyano[N];
  real TN[N];
  real TP[N];
  int site[N];
  matrix[M,M] Dmat;
}
parameters{
  real alpha;
  vector[M] gamma_gp;
  real beta;
  real<lower=0> eta2_gp;
  real<lower=0> rho2_gp;
}
model {
   matrix[M,M] sigma_gp;
   vector[N] lambda;
// prior
  alpha ~ normal(0,10);
  beta ~ normal(0,1);
  eta2_gp ~ cauchy(0,1);
  rho2_gp ~ cauchy(0,1);
//hyper-parameters
  for (i in 1:(M-1))
    for (j in (i+1):M) {
      sigma_gp[i,j]=eta2_gp*exp(-rho2_gp*pow(Dmat[i,j],2));
      sigma_gp[j,i]=sigma_gp[i,j];
    }
  for (k in 1:M)
  sigma_gp[k,k]=eta2_gp*exp(-rho2_gp*pow(Dmat[k,k],2))+0.01;
  gamma_gp ~ multi_normal(rep_vector(0,M), sigma_gp);
//parameters
  for (m in 1:N) 
  lambda[m] =log(alpha+ gamma_gp[site[m]]+beta*TP[m]);
// likelihood
  cyano ~ poisson(lambda);
}
generated quantities {
  matrix[M,M] sigma_gp;
  vector[N] lambda;
  vector[N] log_lik;
  for (m in 1:N) 
  lambda[m] =log(alpha+ gamma_gp[site[m]]+beta*TP[m]);
  // log-likelihood
  log_lik += poisson_lpmf(cyano|lambda);
}
'

pois_GP<-stan(model_name = "poi-gp",
              model_code = code,
              data=stanlist,
              iter = 2000,
              warmup = 500,
              algorithm = "NUTS",
              cores = 4,
              control = list(
                stepsize=0.9,
                max_treedepth=25)
)
