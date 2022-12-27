SMSNcode<-
"
data{
int<lower=1> N;
int<lower=1> M;
int<lower=1> L;
int<lower=1,upper=M> Site[N];
int<lower=1,upper=L> Season[N];

real y[N];
real x1[N];
real x2[N];
real x3[N];
real x4[N];
//
int<lower=1> N_p;
int<lower=1> M_p;
int<lower=1> L_p;
int<lower=1,upper=M_p> Site_p[N_p];
int<lower=1,upper=L_p> Season_p[N_p];

real y_p[N_p];
real x1_p[N_p];
real x2_p[N_p];
real x3_p[N_p];
real x4_p[N_p];
}

parameters{
// coef
real b0;
real b1;
real b2;
real b3;
real b4;

// scales
real<lower=0> tau_0;
real<lower=0> tau_1;
real<lower=0> tau_2;
real<lower=0> tau_3;
real<lower=0> tau_4;
real<lower=0> u;// scale factor
real<lower=0> alpha;
real<lower=0> sigma;// scale SN
real<lower=0> sigma1; //scale normal random effect gamma
real<lower=0> sigma2; //scale normal random effect delta
real<lower=0> nu; //scale freedom
real<lower=0> nu_1; //scale freedom
real<lower=0> lambda; //scale freedom
real<lower=0> lambda_1; //scale freedom

// random parameter
vector<offset = 0, multiplier = sigma1>[M] gamma;
vector<offset = 0, multiplier = sigma2>[L] delta;
}

model{
//AuMuliary paramter
vector[N] Mu;

// ridge coef-- prior
//b0 ~ normal(0,1/tau_0);
//b1 ~ normal(0,1/tau_1);
//b2 ~ normal(0,1/tau_2);
//b3 ~ normal(0,1/tau_3);
//b4 ~ normal(0,1/tau_4);
//tau_0 ~ gamma(0.01,0.01);
//tau_1 ~ gamma(0.01,0.01);
//tau_2 ~ gamma(0.01,0.01);
//tau_3 ~ gamma(0.01,0.01);
//tau_4 ~ gamma(0.01,0.01);
b0 ~ normal(0,100);
b1 ~ normal(0,100);
b2 ~ normal(0,100);
b3 ~ normal(0,100);
b4 ~ normal(0,100);

// random effect
for (j in 1:M)
gamma[j]~ normal(0, sigma1);
for (k in 1:L)
delta[L]~ normal(0,sigma2);
sigma1 ~ normal(0, 10);
sigma2 ~ normal(0, 10);

// fixed effect
//u ~ gamma(nu/2,nu/2);  //SMSN-t
//u ~ beta(nu,1);         //SMSN-beta/slash
//u ~ inv_gamma(nu/2,nu/2);         //SMSN-VG
u ~ gamma(nu/2,nu_1/2);  //SMSN-Pearson type VII distribution

//nu ~ gamma(5,1);
nu ~ exponential(lambda);
lambda~uniform(0.1,10);
//nu_1 ~ gamma(5,1);
nu_1 ~ exponential(lambda_1);
lambda_1~uniform(0.1,10);

sigma ~ normal(0,10);
alpha ~ normal(0,1);

for (i in 1:N)
Mu[i]=
b0
+ b1*x1[i]
+ b2*x2[i]
+ b3*x3[i]
+ b4*x4[i]
+gamma[Site[i]]
+delta[Season[i]]
;
//likelihood
for (n in 1:N)
  target+= skew_normal_lpdf(y[n]|Mu[n],sigma/u,alpha);
}

generated quantities{
real y_rep[N];
real y_pred[N_p];
vector[N] log_lik;
vector[N_p] log_lik_pred;
real mu[N];
real mu_p[N_p];
//real gamma_pred[M] = normal_rng(rep_vector(0,M), sigma1);
real delta_pred[L_p] = normal_rng(rep_vector(0,L_p), sigma2);

for (i in 1:N)
  y_rep[i] = skew_normal_rng(
  b0 
  + b1*x1[i]
  + b2*x2[i]
  + b3*x3[i]
  + b4*x4[i]
  +gamma[Site[i]]
  +delta[Season[i]]
  , sigma/u,alpha);
      
for (i in 1:N)
  mu[i]=
  b0 
  + b1*x1[i]
  + b2*x2[i]
  + b3*x3[i]
  + b4*x4[i]
  +gamma[Site[i]]
  +delta[Season[i]]
  ;

for (i in 1:N)
  log_lik[i]= skew_normal_lpdf(
  y[i]|
  b0 
  + b1*x1[i]
  + b2*x2[i]
  + b3*x3[i]
  + b4*x4[i]
  +gamma[Site[i]]
  +delta[Season[i]]
  , sigma/u,alpha); 

for (i in 1:N_p)
  y_pred[i]= skew_normal_rng(
  b0
  + b1*x1_p[i]
  + b2*x2_p[i]
  + b3*x3_p[i]
  + b4*x4_p[i]
  +gamma[Site_p[i]]
  //+delta[Season_p[i]]
  //+gamma_pred[Site_p[i]]
  +delta_pred[Season_p[i]]
  ,sigma/u,alpha);      

//for (i in 1:N_p)
  //log_lik_pred[i]= skew_normal_lpdf(
  //y_p[i]|
  //b0
  //+ b1*x1_p[i]
  //+ b2*x2_p[i]
  //+ b3*x3_p[i]
  //+ b4*x4_p[i]
  //+gamma[Site_p[i]]
  //+delta[Season_p[i]]
  //+gamma_pred[Site_p[i]]
  //+delta_pred[Season_p[i]]
  //,sigma/u,alpha);      

for (i in 1:N_p)
  mu_p[i]=
  b0 
  + b1*x1_p[i]
  + b2*x2_p[i]
  + b3*x3_p[i]
  + b4*x4_p[i]
  +gamma[Site_p[i]]
  //+delta[Season_p[i]]
  //+gamma_pred[Site_p[i]]
  +delta_pred[Season_p[i]]
  ;      
}
"