# dynamic-asset-allocation
setwd("C:/Users/32388/Desktop/sxr/ReadMe")
Assets=read.table('asset.txt',header = T,sep = '\t')
head(Assets)
Assets$date=as.Date(Assets$date,"%Y-%m-%d")
#View(Assets)
win=72 
library(data.table)
library(PerformanceAnalytics)
library(CVXR)
Assets=as.data.table(Assets)
#Assets=Assets[1:(nrow(Assets)-win),1:ncol(Assets),with=FALSE] #to do monthly rebalance
#View(Assets)
#Assets=Assets[,-1]
Assets=Assets[(nrow(Assets)-win+1):nrow(Assets),-1,with=FALSE] #do not do monthly reblance
AssetCov=cov(Assets)
AssetRtn=apply(Assets,2,mean)

#View(AssetCov)
#View(AssetRtn)
#########################################

#MVO
Vol_tgt=0.025
ObjFun=function(w){
  -sum(AssetRtn*w)
}

ObjFunGradient=function(w){
  -AssetRtn
}
Equal=function(w){
  sum(w)-1
}
EqualGradient=function(w){
  rep(1,length(w))
}
Inequal=function(w){
  t(w) %*% AssetCov %*% w-Vol_tgt^2
}
InequalGradient=function(w){
  t(w) %*% (AssetCov+t(AssetCov))
}
w0=rep(1/14,14)
library(nloptr)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=Equal,eval_jac_g_eq=EqualGradient,
           eval_g_ineq=Inequal,eval_jac_g_ineq=InequalGradient,
           lb=rep(-1,14),ub=rep(1,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-8,'maxeval'=100000000))
w_MVO=opt$solution
#View(w_MVO)
sum(AssetRtn*w_MVO)
sqrt(t(w_MVO) %*% AssetCov %*% w_MVO)
MCR=AssetCov %*% w_MVO/sqrt(t(w_MVO) %*% AssetCov %*% w_MVO)[1,1]
MVO_202004=data.frame(w_MVO,w_MVO*MCR)
#write.csv(MVO_202004,file = "C:/Users/32388/Desktop/indiproj/MVO_202004.csv")

#risk parity
Assets1=Assets[,-'cash']
AssetCov1=cov(Assets1)
ObjFun=function(w){
  -sum(log(abs(w)))
}

ObjFunGradient=function(w){
  -(1/w)
}

Inequal=function(w){
  t(w) %*% AssetCov1 %*% w-Vol_tgt^2
}
InequalGradient=function(w){
  t(w) %*% (AssetCov1+t(AssetCov1))
}

w0=rep(1/13,13)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_ineq=Inequal,eval_jac_g_ineq=InequalGradient,
           lb=rep(-10,13),ub=rep(10,13),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-8,'maxeval'=100000000))
w_riskpa=opt$solution
#View(w_riskpa)
sum(w_riskpa)
sum(AssetRtn*c(w_riskpa,1-sum(w_riskpa)))
sqrt(t(w_riskpa) %*% AssetCov1 %*% w_riskpa)
MCR=AssetCov1 %*% w_riskpa/sqrt(t(w_riskpa) %*% AssetCov1 %*% w_riskpa)[1,1]
data.frame(w_riskpa,w_riskpa*MCR)
parity_202004=data.frame(w_riskpa,w_riskpa*MCR)
#write.csv(parity_202004,file = "C:/Users/32388/Desktop/indiproj/parity_202004.csv")

#compare return of MVO to risk parity
z=read.table("perf.txt",header = T,sep="\t")
head(z)
return_MVO=ts(z$MVO*100,start=c(2016,12),end=c(2020,4),frequency=12)
return_riskparity=ts(z$riskparity*100,start=c(2016,12),end=c(2020,4),frequency=12)
plot(return_riskparity,main = "portfolio return of risk parity method(%)",col="red",type = 'l')
plot(return_MVO,main = "portfolio return ",col="blue",type = 'l',ylim=c(0,1.8))
lines(return_riskparity,col="red",type = 'l')

################################################################
#quintile portfolio
quintile_portfolio <- function(dataset) {
  X=Assets
  N = ncol(X)
  ranking = sort(colMeans(X), decreasing = TRUE, index.return = TRUE)$ix
  w = rep(0, N)
  w[ranking[1:round(N/5)]] <- 1/round(N/5)
  return(w)
}
w_quintile=quintile_portfolio(Assets)
sum(w_quintile)
sum(AssetRtn*w_quintile)
sqrt(t(w_quintile) %*% AssetCov %*% w_quintile)
MCR=AssetCov %*% w_quintile/sqrt(t(w_quintile) %*% AssetCov %*% w_quintile)[1,1]
quintile=data.frame(w_quintile,w_quintile*MCR)
quintile
#write.csv(quintile,file = "C:/Users/32388/Desktop/indiproj/quintile_202004.csv")


#GMVP (remove cash and long only)
GMVP_portfolio <- function(dataset) {
  X=Assets1
  Sigma = cov(X)  
  w = solve(Sigma, rep(1, nrow(Sigma)))
  w = abs(w)/sum(abs(w))
  return(w)
}
w_GMVP=GMVP_portfolio(Assets)
sum(w_GMVP)
sum(AssetRtn*c(w_GMVP,1-sum(w_GMVP)))
sqrt(t(w_GMVP) %*% AssetCov1 %*% w_GMVP)
MCR=AssetCov1 %*% w_GMVP/sqrt(t(w_GMVP) %*% AssetCov1 %*% w_GMVP)[1,1]
data.frame(w_GMVP,w_GMVP*MCR)
GMVP=data.frame(w_GMVP,w_GMVP*MCR)
GMVP
#write.csv(GMVP,file = "C:/Users/32388/Desktop/indiproj/GMVP_202004.csv")

################################################################
#BL model
Sigma=AssetCov
mu_bl=colMeans(Assets)
mu_bl
T=nrow(Assets)
T_trn<-round(0.7*T) #split data into test part and training part
X_trn<-Assets[1:T_trn,]
X_tst<-Assets[(T_trn+1):T,]

compute_BL=function(pi,c_roll){
  A=1/T_trn
  P=rbind(c(0,0,0,1,0,0,0,0,0,0,0,0,0,0),
          c(0,0,0,0,0,-1,1,0,0,0,0,0,0,0),
          c(0,0,0,0,0,0,0,0,0,0,0,0,1,0))
  v=c(0.01,0.0015,0.01)
  mu_BL=NULL
  Sigma_BL=list()
  for (c in c_roll+1e-6) {
    Omega <- (1/c) * P %*% Sigma %*% t(P)
    mu_ <- pi + A * Sigma %*% t(P) %*% solve(A * P %*% Sigma %*% t(P) + Omega) %*% (v - P %*% pi)
    Sigma_ <- (1+A)*Sigma - 
      (A^2) * Sigma %*% t(P) %*% solve(A * P %*% Sigma %*% t(P) + Omega) %*% P %*% Sigma
    mu_BL <- cbind(mu_BL, mu_)
    Sigma_BL <- c(Sigma_BL, list(Sigma_))
  }
  colnames(mu_BL) <- paste("c =", round(c_roll))
  names(Sigma_BL) <- paste("c =", round(c_roll))
  return(list(mu = mu_BL, Sigma = Sigma_BL))
}
c_roll <- seq(0, T_trn, length.out = 10)
moment_BL <- compute_BL(pi = mu_bl, c_roll)
names(moment_BL)
print(moment_BL$mu)
names(moment_BL$Sigma)
print(moment_BL$Sigma[[1]])


mu_true <- colMeans(X_tst)
error_sm_BL <- colSums(abs((moment_BL$mu - mu_true)^2))
print(error_sm_BL)
plot(c_roll, error_sm_BL, xlab = "c", type = "b", 
     main = "Estimation error in mu based on BL with sample mean")

#compute weights of different assets
portolioBL <- function(mu, Sigma, lmd = 0.5) {
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lmd*quad_form(w, Sigma)),
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

w_BL <- NULL
for (i in 1:length(c_roll)){
  w_BL <- cbind(w_BL, portolioBL(moment_BL$mu[, i], moment_BL$Sigma[[i]]))}
colnames(w_BL) <- colnames(moment_BL$mu)
rownames(w_BL) <- rownames(moment_BL$mu)
w_BL
chart.StackedBar(t(w_BL), ylab = "w", space=0, border = NA,
                 main = "Portfolio (based on BL with sample mean) as a function of uncertainty") #need about 10s to compute

# compute returns of all portfolios
X=as.matrix(Assets)
return_BL <- ts(X %*% w_BL)
return_BL_TR <- return_BL[1:T_trn, ]
return_BL_TEST <- return_BL[-c(1:T_trn), ]
plot(return_BL[,1],col=1,main='return of portfolio by BL method') #whole data set
for (i in 2:10){
  lines(return_BL[,i],col=i)
}
plot(return_BL_TR[,1],col=1,main='return given by BL method-training data',type='l') #training part
for (i in 2:10){
  lines(return_BL_TR[,i],col=i)
}
plot(return_BL_TEST[,1],col=1,main='return given by BL method-test data',type='l',ylim=c(-0.1,0.15)) #testing part
for (i in 2:10){
  lines(return_BL_TEST[,i],col=i)
}

return_BL_TR=ts(return_BL_TR)
return_BL_TEST=ts(return_BL_TEST)
table.AnnualizedReturns(return_BL)
table.AnnualizedReturns(return_BL_TR)
table.AnnualizedReturns(return_BL_TEST)

#############################################################################################
#factor based asset allocation
#case 1 target factor exposure
Factors=read.table('riskfactor.txt',header = T,sep='\t')
Factors=as.data.table(Factors)
head(Factors)
win=72
##Factors=Factors[1:(nrow(Factors)-win),1:ncol(Factors),with=FALSE] #to do monthly rebalance
#View(Factors)
#Factors=Factors[,-1]
Factors=Factors[(nrow(Factors)-win+1):nrow(Factors),-1,with=FALSE] #do not do monthly reblance
AssetExpo=matrix(0,ncol(Assets),ncol(Factors))
rownames(AssetExpo)=colnames(Assets)
colnames(AssetExpo)=colnames(Factors)
#View(AssetExpo)
Spe=matrix(0,win,ncol(Assets))
#View(Spe)

Mapping=function(y){
  dat=data.table(Y=y,Factors)
  full=lm(Y~.,data=dat)
  null=lm(Y~1,data=dat)
  best=step(null,scope=list(upper=full,lower=null),direction = 'both')
  list(coeff=best$coefficients[-1],res=summary(best)$residuals)
}
for(i in 1:ncol(Assets)){
  a=Mapping(data.frame(Assets[,i,with=FALSE])[,1])
  AssetExpo[i,names(a$coeff)]=a$coeff
  Spe[,i]=a$res
}
#View(AssetExpo)
#View(Spe)
library(ggplot2)
risk=read.table("risk.txt",header=T,sep="\t")
ggplot(risk,aes(x=asset,y=V1,fill =risk )) + geom_bar(stat="identity")

FacCov=cov(Factors)
SpeCov=cov(Spe)
TargetExpo=c(0.45,0.12,0.07,0.1,0.29,0.2)
names(TargetExpo)=colnames(Factors)
lambda=0.99
ObjFun=function(w){
  M=w %*% AssetExpo-TargetExpo
  ((1-lambda) * M %*% t(M) + lambda * M %*% FacCov %*% t(M) + lambda * w %*% SpeCov %*% w)[1,1]
}

ObjFunGradient=function(w){
  M=w %*% AssetExpo-TargetExpo
  2*(1-lambda) * M %*% t(AssetExpo) + 2*lambda *M %*% FacCov %*% t(AssetExpo) +2*lambda * w %*% SpeCov
}
fullyInvested=function(w){
  sum(w)-1
}
fullyInvestedGradient=function(w){
  rep(1,length(w))
}

w0=rep(1/14,14)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=fullyInvested,eval_jac_g_eq=fullyInvestedGradient,
           lb=rep(-1,14),ub=rep(1,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-4,'maxeval'=100000000))
w_fac=opt$solution

#View(w_fac)
w_fac %*% AssetExpo


#case 2 equal exposure to each risk factor
Equal=function(w){
  return(rbind(sum(w)-1,w%*%(AssetExpo[,1]-AssetExpo[,2]),w%*%(AssetExpo[,1]-AssetExpo[,3]),
               w%*%(AssetExpo[,1]-AssetExpo[,4]),w%*%(AssetExpo[,1]-AssetExpo[,5]),w%*%(AssetExpo[,1]-AssetExpo[,6])))
}
EqualGradient=function(w){
  return(rbind(rep(1,length(w)),AssetExpo[,1]-AssetExpo[,2],AssetExpo[,1]-AssetExpo[,3],
               AssetExpo[,1]-AssetExpo[,4],AssetExpo[,1]-AssetExpo[,5],AssetExpo[,1]-AssetExpo[,6]))
}

w0=rep(1/14,14)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=Equal,eval_jac_g_eq=EqualGradient,
           lb=rep(-10,14),ub=rep(10,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-4,'maxeval'=100000000))
w_fac_riskpa=opt$solution
View(w_fac_riskpa)

w_fac_riskpa %*% AssetExpo
w_MVO  %*% AssetExpo
w_riskpa %*% AssetExpo[-14,]
w_fac %*% AssetExpo
riskexp_fac=data.frame(w_fac %*% AssetExpo)
riskexp_fac_riskpa=data.frame(w_fac_riskpa %*% AssetExpo)
riskexp_MVO=data.frame(w_MVO  %*% AssetExpo)
riskexp_riskpa=data.frame(w_riskpa  %*% AssetExpo[-14,])
riskexp=cbind(riskexp_fac,riskexp_fac_riskpa,riskexp_MVO,riskexp_riskpa)
#write.csv(riskexp,file = "C:/Users/32388/Desktop/indiproj/riskexp.csv")














