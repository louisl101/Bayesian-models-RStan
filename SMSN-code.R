##-------------------------my librarys------------------------------
library(rstan)
library(brms)
library(bayesplot)
library(ggplot2)
library(openxlsx)
library(rstanarm)
library(factoextra)
library(FactoMineR)
library(readxl)
library(loo)
library(dplyr)
library(plyr)
library(rethinking)
library(corrplot)
library(vegan)
library(coda)
library(psych)
library(Hmisc)
##-------------------------my functions------------------------------

r_square<-function(pred, obs){
  return (1 - (sum((pred - obs)^2)/sum((obs - mean(obs))^2)) )
}

RMSE <- function(pred, obs) {
  sqrt(mean((pred - obs)^2))
}
NRMSE <- function(pred, obs) {
  sqrt(mean((pred - obs)^2))/(max(obs)-min(obs))
  # sqrt(mean((pred - obs)^2))/mean(obs)
  
}

probs<-function(ppd,lim1,lim2){
  prob1=sum(ppd>=lim1)/length(ppd)
  prob2=sum(ppd>lim2)/length(ppd)
  return(cbind(prob1,prob2))
}

##------------------------- SMSN stancode------------------------------
stancode<-
  '
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
b0 ~ normal(0,1/tau_0);
b1 ~ normal(0,1/tau_1);
b2 ~ normal(0,1/tau_2);
b3 ~ normal(0,1/tau_3);
b4 ~ normal(0,1/tau_4);
tau_0 ~ gamma(0.01,0.01);
tau_1 ~ gamma(0.01,0.01);
tau_2 ~ gamma(0.01,0.01);
tau_3 ~ gamma(0.01,0.01);
tau_4 ~ gamma(0.01,0.01);
//b0 ~ normal(0,100);
//b1 ~ normal(0,100);
//b2 ~ normal(0,100);
//b3 ~ normal(0,100);
//b4 ~ normal(0,100);

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
real mu[N];
//real gamma_pred[M] = normal_rng(rep_vector(0,M), sigma1);
//real delta_pred[L_p] = normal_rng(rep_vector(0,L_p), sigma2);

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
  +delta[Season_p[i]]
  //+gamma_pred[Site_p[i]]
  //+delta_pred[Season_p[i]]
  ,sigma/u,alpha);      
}
'
##------------------------- NTL Trout lakes-SMSN multilevel regression------------------------------
new_d<-read.csv("~/new_d.csv") ## chl-a seed=6 trout
new_d$biomass<-log(new_d$chla)
#
##
new_d$point_num<-as.numeric(as.factor(new_d$lakeid))
new_d$season_num<-as.numeric(as.factor(new_d$date))

train<-subset(new_d,year!=2018)
test<-subset(new_d,year==2018)

train[,c('pH','WT','DO','SDD')]<-decostand(train[,c('pH','WT','DO','SDD')],"standardize",MARGIN = 2)
test[,c('pH','WT','DO','SDD')]<-decostand(test[,c('pH','WT','DO','SDD')],"standardize",MARGIN = 2)
##
stanlist<-list(N=nrow(train),
               M=max(train$point_num),
               L=max(train$season_num),
               y=train$biomass,
               x1=train$WT,
               x2=train$SDD,
               x3=train$DO,
               x4=train$pH,
               
               Site=train$point_num,
               Season=train$season_num,
               
               N_p=nrow(test),
               M_p=max(test$point_num),
               L_p=max(test$season_num),
               y_p=test$biomass,
               x1_p=test$WT,
               x2_p=test$SDD,
               x3_p=test$DO,
               x4_p=test$pH,
               
               Site_p=test$point_num,
               Season_p=test$season_num
)

model_smsn_x <-stan(model_name = "SMSN",
                    model_code = stancode,
                    data=stanlist,
                    chains = 4,
                    iter = 20000,
                    warmup = 5000,
                    algorithm = "NUTS",
                    cores = 6,
                    control = list(
                      stepsize=0.95,
                      max_treedepth=25,
                      adapt_delta=0.99),
                    init = "0",
                    seed =6
)

mc<-extract(model_smsn_x)

plot(x=apply(mc$y_rep,2,median),y=stanlist$y,main=round(r_square(apply(mc$y_rep,2,median),stanlist$y),digits = 4))
# ,xlab=waic(mc$log_lik)[['estimates']][3,1])
abline(a=0,b=1)
plot(x=apply(mc$y_pred,2,median),y=stanlist$y_p,main=round(r_square(apply(mc$y_pred,2,median),stanlist$y_p),digits = 4))
# ,xlab=waic(mc$log_lik_pred)[['estimates']][3,1])
abline(a=0,b=1)

## evaluation-index
r_square(apply(mc$y_rep,2,median),stanlist$y)
RMSE(apply(mc$y_rep,2,median),stanlist$y)

r_square(apply(mc$y_pred,2,median),stanlist$y_p)
RMSE(apply(mc$y_pred,2,median),stanlist$y_p)
##------------------------- regression fig 2400*750------------------------------
## train-data-ready
ppd_trout<-mc$y_rep

ci_trout<-cbind(posterior_interval(ppd_trout,0.95),rep('trout',length(stanlist$y)),apply(ppd_trout,2,median),stanlist$y)
ci<-rbind(ci_trout)
stat_train<-data.frame(down=as.numeric(ci[,1]),
                       up=as.numeric(ci[,2]),
                       mean=as.numeric(ci[,4]),
                       obs=as.numeric(ci[,5]),
                       model=ci[,3],
                       dataset=rep('train',nrow(ci))
)

## test-data-ready
ppd_trout_p<-mc$y_pred
ci_trout_p<-cbind(posterior_interval(ppd_trout_p,0.95),rep('trout',length(stanlist$y_p)),apply(ppd_trout_p,2,median),stanlist$y_p)
ci_p<-rbind(ci_trout_p)
stat_p<-data.frame(down=as.numeric(ci_p[,1]),
                   up=as.numeric(ci_p[,2]),
                   mean=as.numeric(ci_p[,4]),
                   obs=as.numeric(ci_p[,5]),
                   model=ci_p[,3],
                   dataset=rep('test',nrow(ci_p))
)

stat<-rbind(stat_train,stat_p)

m='test'
c="#3C5488FF"

m='train'
c="#B62C19"

p1<-
  ggplot(data=subset(stat,model=='trout'))+
  # geom_point(size=2.5)+
  # geom_errorbar(size=0.5,alpha=0.4)+
  geom_pointrange(aes(x=obs,y=mean,ymin=down, ymax=up,fill=dataset),col='darkgray',size=1.5,fatten =4,shape=21)+
  geom_abline(intercept =0,slope=1,col="black",alpha=0.4,size=1.5)+
  geom_segment(aes(x=log(1),xend=log(1),yend=-4.2,y=log(1)),col="#B62C19",alpha=0.4,size=1.1,linetype=6)+
  geom_segment(aes(x=log(12),xend=log(12),yend=-4.2,y=log(12)),col="#B62C19",alpha=0.4,size=1.1,linetype=6)+
  geom_segment(aes(y=log(1),yend=log(1),xend=-4.2,x=log(1)),col="#B62C19",alpha=0.4,size=1.1,linetype=6)+
  geom_segment(aes(y=log(12),yend=log(12),xend=-4.2,x=log(12)),col="#B62C19",alpha=0.4,size=1.1,linetype=6)+
  # scale_color_manual(values = c("#3C5488FF","#3C884A"))+
  scale_fill_manual(values = c("#3C5488FF","#FF9900"))+
  theme_bw()+
  theme(
    line=element_line(size = 0.5),
    plot.title = element_text(size = 22,face="bold"),
    axis.title=element_text(size = 26),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length=unit(.25, "cm"),
    axis.text = element_text(size=25,face="italic",color="black"),
    axis.title.x =element_blank(),
    axis.title.y =element_blank(),
    # legend.title=element_text(size=14,face="plain",color="black")
    legend.position='none'
  )+
  scale_y_continuous(limits = c(-4.2,6.5),breaks = seq(-4,6,1),label = c('-4','','-2','','0','','2','','4','','6'),expand = c(0,0))+
  scale_x_continuous(limits = c(-4.2,6.5),breaks = seq(-4,6,1),label = c('-4','','-2','','0','','2','','4','','6'),expand = c(0,0))+
  labs(x = "Observed", 
       y = "Predicted",title = "(b) Trout Lake Region lakes")
p1

##------------------------- exceedance fig 1200*700------------------------------
# optimal probabilistic threshold for trout:  0.61 0.56

t1<-0.61
t2<-0.56
#--------exceedance level 1--------
##calibration-probability
par(mfrow=c(1,2))
ppd<-ppd_trout
lim1=log(1)
lim2=log(12)
#-----data ready-------
#level 1
p<-matrix(ncol = 2,nrow = ncol(ppd))
for (i in 1:ncol(ppd)){
  p[i,]<-probs(ppd[,i],lim1,lim2)
}
p<-data.frame(p,y=stanlist$y,p_condition=p[,2]/p[,1])
#
ppd_test<-mc$y_pred
p_test<-matrix(ncol = 2,nrow = ncol(ppd_test))
for (i in 1:ncol(ppd_test)){
  p_test[i,]<-probs(ppd_test[,i],lim1,lim2)
}
p_test<-data.frame(p_test,y=stanlist$y_p,p_condition=p_test[,2]/p_test[,1])

#conditional level 2 --define after calculating t1
p_level2_test<-subset(p_test,X1>t1&y>lim1)
p_level2<-subset(p,X1>t1&y>lim1)
###----- plot data------level 1-----
plot(y=p[,1],x=stanlist$y,
     ylim=c(-0.1,1.05),
     # xlim=c(-max(stanlist$y)-.2,max(stanlist$y)+.2),
     xlim=c(-4.2,10.6),
     ylab="",
     xlab="",
     main="Model Calibration (2018 Data)",
     pch=16,
     cex=2,
     cex.lab=1.2,
     cex.main=1.5,
     xaxt="n",
     yaxt="n",
     xaxs="i",
     yaxs="i",
     bg='black'
)
axis(1,-12:12,
     labels = c("","","-10","","8","","-6","","-4","","-2","","0","","2","","4","","6","","8","","10","",""),
     cex.axis = 1.5,font.axis = 3,gap.axis = 1)
axis(2,seq(0,1,by=0.1), 
     labels = c("0","","0.2","","0.4","","0.6","","0.8","","1"),
     cex.axis = 1.5,font.axis = 3)
abline(v=lim1,lwd=2,col="#DC143C")
abline(h=t1,lwd=2,lty=4,col="#DC143C")

##test-probability
#
plot(y=p_test[,1],x=stanlist$y_p,
     ylim=c(-0.1,1.05),
     # xlim=c(-max(stanlist$y)-.2,max(stanlist$y)+.2),
     xlim=c(-4.2,10.6),
     ylab="",
     xlab='',
     main="Model Validation (2019 Data)",
     pch=1,
     cex=2,
     cex.lab=1.2,
     cex.main=1.5,
     xaxt="n",
     yaxt="n",
     xaxs="i",
     yaxs="i"
)
axis(1,-12:12,
     labels = c("","","-10","","8","","-6","","-4","","-2","","0","","2","","4","","6","","8","","10","",""),
     cex.axis = 1.5,font.axis = 3,gap.axis = 1)
axis(2,seq(0,1,by=0.1), 
     labels = c("0","","0.2","","0.4","","0.6","","0.8","","1"),
     cex.axis = 1.5,font.axis = 3)
abline(v=lim1,lwd=2,col="#DC143C")
abline(h=t1,lwd=2,lty=4,col="#DC143C")

##--------exceedance level 2--------
##calibration-probability
plot(y=p_level2$p_condition,x=p_level2$y,
     ylim=c(-0.1,1.05),
     # xlim=c(-max(stanlist$y)-.2,max(stanlist$y)+.2),
     xlim=c(-4.2,10.6),
     ylab="",
     xlab="",
     main="Model Calibration (2018 Data)",
     pch=16,
     # pch=23,
     cex.lab=1.2,
     cex.main=1.5,
     xaxt="n",
     yaxt="n",
     xaxs="i",
     yaxs="i",
     bg='black'
)
axis(1,-12:12,
     labels = c("","","-10","","8","","-6","","-4","","-2","","0","","2","","4","","6","","8","","10","",""),
     cex.axis = 1.5,font.axis = 3,gap.axis = 1)
axis(2,seq(0,1,by=0.1), 
     labels = c("0","","0.2","","0.4","","0.6","","0.8","","1"),
     cex.axis = 1.5,font.axis = 3)
abline(v=lim2,lwd=2,col="#DC143C")
abline(h=t2,lwd=2,lty=4,col="#DC143C")

##test-probability
plot(y=p_level2_test$p_condition,x=p_level2_test$y,
     ylim=c(-0.1,1.05),
     # xlim=c(-max(stanlist$y)-.2,max(stanlist$y)+.2),
     xlim=c(-4.2,10.6),
     ylab="",
     xlab="",
     main="Model Validation (2019 Data)",
     pch=1,
     cex=2,
     cex.lab=1.2,
     cex.main=1.5,
     xaxt="n",
     yaxt="n",
     xaxs="i",
     yaxs="i"
)
axis(1,-12:12,
     labels = c("","","-10","","8","","-6","","-4","","-2","","0","","2","","4","","6","","8","","10","",""),
     cex.axis = 1.5,font.axis = 3,gap.axis = 1)
axis(2,seq(0,1,by=0.1), 
     labels = c("0","","0.2","","0.4","","0.6","","0.8","","1"),
     cex.axis = 1.5,font.axis = 3)
abline(v=lim2,lwd=2,col="#DC143C")
abline(h=t2,lwd=2,lty=4,col="#DC143C")

