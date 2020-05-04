
library("quantmod")
library("copula")
library("fGarch")
library("rugarch")
library("LambertW")
library("forecast")

###Load data from Yahoo
getSymbols.yahoo('SPXL',env=globalenv(), from = "2009-11-25", to = "2019-11-25",periodicity = "daily")
getSymbols.yahoo('WFC',env=globalenv(), from = "2009-11-25", to = "2019-11-25",periodicity = "daily")
summary(SPXL)
summary(WFC)
dim(SPXL)
dim(WFC)

###Exploratory data analysis
plot(SPXL$SPXL.Adjusted)
plot(WFC$WFC.Adjusted)
#Stationarity test
acf(SPXL$SPXL.Adjusted , lag.max = 300, main="ACF of S&P500")
acf(WFC$WFC.Adjusted, lag.max = 300, main="ACF of Well Fargo")
Box.test(SPXL$SPXL.Adjusted,lag=600,type="Ljung-Box",fitdf=2)
Box.test(WFC$WFC.Adjusted,lag=600,type="Ljung-Box",fitdf=2)
#Normarlity test
qqnorm(SPXL$SPXL.Adjusted, main="S&P500 - QQplot")
qqline(SPXL$SPXL.Adjusted)
qqnorm(WFC$WFC.Adjusted, main="Well Fargo - QQplot")
qqline(WFC$WFC.Adjusted)
ks.test.t(WFC$WFC.Adjusted)
ks.test.t(SPXL$SPXL.Adjusted)
jarque.bera.test(SPXL$SPXL.Adjusted) 
jarque.bera.test(WFC$WFC.Adjusted)

#Skewness and kurtosis
skewness(SPXL$SPXL.Adjusted)
skewness(WFC$WFC.Adjusted)
kurtosis(SPXL$SPXL.Adjusted)
kurtosis(WFC$WFC.Adjusted)

###Log return
sp=diff(log(SPXL$SPXL.Adjusted))[2:length(SPXL$SPXL.Adjusted)]
stock=diff(log(WFC$WFC.Adjusted))[2:length(WFC$WFC.Adjusted)]
plot(sp)
plot(stock)
#Stationarity test
acf(sp,main="ACF of S&P500 log return")
acf(stock, main="ACF of Well Fargo log return")
Box.test(sp,lag=600,type="Ljung-Box")
Box.test(stock,lag=600,type="Ljung-Box")
#Normality test
qqnorm(sp, main="S&P500 log return - QQplot")
qqline(sp)
qqnorm(stock, main="Well Fargo log return - QQplot")
qqline(stock)
jarque.bera.test(sp)
jarque.bera.test(stock)

library("tseries")
#Test Unit root
adf.test(sp)
pp.test(sp)
kpss.test(sp)
adf.test(stock)
pp.test(stock)
kpss.test(stock)

#Covariance and Correlation
cov(sp,stock)
cor(sp,stock)

library("fGarch")
###AR(1)-GARCH(1,1) s&p500
fit1=garchFit(formula=~arma(1,0)+garch(1,1),data=sp,cond.dist="std")
res1=residuals(fit1)#residual of arma 
res_sd1=residuals(fit1,standardize=TRUE)#Residual of arma garch part
acf(res1,main="ACF of S&P500 ARMA residual")
acf(res1^2,main="ACF of S&P500 ARMA residual square")
Box.test(res1^2,lag=10,type="Ljung-Box",fitdf=2)# not wn
acf(res_sd1,main="ACF of S&P500 ARMA/GARCH residual")
acf(res_sd1^2,main="ACF of S&P500 ARMA/GARCH residual square")
Box.test(res_sd1^2,lag=10,type="Ljung-Box",fitdf=2)# wn
qqnorm(res_sd1, main = "Normal Q-Q Plot S&p500 ARMA/GARCH residual square") ; qqline(res_sd1)# not normal distribution

###AR(1)-GARCH(1,1) Well fargo stock
fit2=garchFit(formula=~arma(1,0)+garch(1,1),data=stock,cond.dist="std")
res2=residuals(fit2)#residual of arma 
res_sd2=residuals(fit2,standardize=TRUE)#Residual of arma garch part
acf(res2,main="ACF of Well Fargo ARMA residual")
acf(res2^2,main="ACF of Well Fargo ARMA residual square")
Box.test(res2^2,lag=10,type="Ljung-Box",fitdf=2)# not wn
acf(res_sd2,main="ACF of Well Fargo ARMA/GARCH residual")
acf(res_sd2^2,main="ACF of Well FarO ARMA/GARCH residual square")
Box.test(res_sd2^2,lag=10,type="Ljung-Box",fitdf=2)# wn
qqnorm(res_sd2, main = "Normal Q-Q Plot Well fargo ARMA/GARCH residual square") ; qqline(res_sd2)# not normal distribution



library(sn)
###fit parametric distribution families to both epsilons
#explore res_sd1 and res_sd2
skewness(res_sd1)
skewness(res_sd2)
kurtosis(res_sd1)
kurtosis(res_sd2)
hist(res_sd1)
hist(res_sd2)
ks.test.t(res_sd1)#Test student t
ks.test.t(res_sd2)

#skewed t
fit_sp_skew=sstdFit(res_sd1)
#aic
-2*(fit_sp_skew$minimum)+2*length(fit_sp_skew$estimate)
#bic
-2*(fit_sp_skew$minimum)+log(length(sp))*length(fit_sp_skew$estimate)

fit_stock_skew=sstdFit(res_sd2)
#aic
-2*(fit_stock_skew$minimum)+2*length(fit_stock_skew$estimate)
#bic
-2*(fit_stock_skew$minimum)+log(length(stock))*length(fit_stock_skew$estimate)


#student t
fit_sp=fitdistr(res_sd1,'t')
#aic
2*(fit_sp$loglik)+2*length(fit_sp$estimate)
#bic
2*(fit_sp$loglik)+log(length(sp))*length(fit_sp$estimate)
u_sp=pst(sp,df= fit_sp$estimate[3]); #transform to a variable with uniform distribution

fit_stock=fitdistr(res_sd2,'t')
#aic
2*(fit_stock$loglik)+2*length(fit_stock$estimate)
#bic
2*(fit_stock$loglik)+log(length(stock))*length(fit_stock$estimate)

u_stock=pst(stock,df= fit_stock$estimate[3]); #transform to a variable with uniform distribution


U.hat=cbind(u_sp,u_stock)
fhatU=kde(x=U.hat,H=Hscv(x=U.hat));#nonparametric density estimation 
plot(fhatU,cont=seq(10,80,10)); #contour plots
library(rgl)
plot3d(u_sp,u_stock,pch=50,col='navyblue')

###Fit parametric copula family

#KENDAL TAU
tau=as.numeric(cor.test(u_sp,u_stock,method="kendall")$estimate);
omega=sin(tau*pi/2); omega#estimator for rho
#Spearman
spear=as.numeric(cor.test(u_sp,u_stock,method="spearman")$estimate);
omega2=sin(spear*pi/2); omega2#estimator for rho
#Pearson
pearson=as.numeric(cor.test(u_sp,u_stock,method="pearson")$estimate);
omega3=sin(pearson*pi/2); omega3#estimator for rho

set.seed(1)
#fit t copula with kendal tau at starting point
Ct=fitCopula(copula=tCopula(dim=2),data=U.hat,method="ml",start=c(omega,0.1));
Ct@estimate;
loglikCopula(param=Ct@estimate,u=U.hat,copula=tCopula(dim=2));#compute loglikelihood function
-2*.Last.value+2*length(Ct@estimate);#compute 
coef(Ct)
AIC(Ct)
BIC(Ct)
#fit t copula with spearman
Ct=fitCopula(copula=tCopula(dim=2),data=U.hat,method="ml",start=c(omega2,0.1));
Ct@estimate;
loglikCopula(param=Ct@estimate,u=U.hat,copula=tCopula(dim=2));
-2*.Last.value+2*length(Ct@estimate);
AIC(Ct)
BIC(Ct)
#fit t copula with pearson
Ct=fitCopula(copula=tCopula(dim=2),data=U.hat,method="ml",start=c(omega3,0.1));
Ct@estimate;
loglikCopula(param=Ct@estimate,u=U.hat,copula=tCopula(dim=2));
-2*.Last.value+2*length(Ct@estimate);
AIC(Ct)
BIC(Ct)
#Visualize copula
persp(tCopula(dim=2,Ct@estimate[1],df=Ct@estimate[2]),dCopula)

#fit Gaussian copula with kendal tau 
Cgauss=fitCopula(copula=normalCopula(dim=2),data=U.hat,method="ml",start=c(omega));
Cgauss@estimate;
loglikCopula(param=Cgauss@estimate,u=U.hat,copula=normalCopula(dim=2));
-2*.Last.value+2*length(Cgauss@estimate); 
AIC(Cgauss)
BIC(Cgauss)
#fit Gaussian copula with spearman
Cgauss=fitCopula(copula=normalCopula(dim=2),data=U.hat,method="ml",start=c(omega2));
Cgauss@estimate;
loglikCopula(param=Cgauss@estimate,u=U.hat,copula=normalCopula(dim=2));
-2*.Last.value+2*length(Cgauss@estimate); 
AIC(Cgauss)
BIC(Cgauss)
#fit Gaussian copula with pearson
Cgauss=fitCopula(copula=normalCopula(dim=2),data=U.hat,method="ml",start=c(omega3));
Cgauss@estimate;
loglikCopula(param=Cgauss@estimate,u=U.hat,copula=normalCopula(dim=2));
-2*.Last.value+2*length(Cgauss@estimate);
AIC(Cgauss)
BIC(Cgauss)

#fit frank copula 
Cfr=fitCopula(copula=frankCopula(1,dim=2),data=U.hat,method="ml");
Cfr@estimate;
loglikCopula(param=Cfr@estimate,u=U.hat,copula=frankCopula(dim=2));
-2*.Last.value+2*length(Cfr@estimate);
AIC(Cfr)
BIC(Cfr)

#fit Gumbel copula 
ParGum<-1/(1-omega)
ParGum
Cgumbel = fitCopula(copula=gumbelCopula(ParGum,dim=2),data=U.hat,method="ml")
Cgumbel@estimate;
loglikCopula(param=Cgumbel@estimate,u=U.hat,copula=gumbelCopula(dim=2));
-2*.Last.value+2*length(Cgumbel@estimate); 
AIC(Cgumbel)
BIC(Cgumbel)
Cgumbel@estimate
#Visualize copula
persp(gumbelCopula(dim=2,Cgumbel@estimate[1]),dCopula)


#fit clayton copula
ParClay<-(2*omega)/(1-omega)
ParClay
Cclayton = fitCopula(copula = claytonCopula(ParClay, dim=2),data = U.hat, method = "ml")
Cclayton@estimate;
loglikCopula(param=Cclayton@estimate,u=U.hat,copula=claytonCopula(dim=2));
-2*.Last.value+2*length(Cclayton@estimate); 
AIC(Cclayton)
BIC(Cclayton)

### Simulation from copula
library(copula)
cop_simulation = gumbelCopula(Cgumbel@estimate[1] ,dim = 2)
set.seed(99)
rand_t_cop = rCopula(n = 10000, copula = cop_simulation)
cor(rand_t_cop)
cor.test(rand_t_cop[,1],rand_t_cop[,2])
coef(Cgumbel)

#Simulate standardized errors from random sample generated by copula
e_sp_sim = (qt(rand_t_cop[,1], df = fit_sp$estimate[3])*fit_sp$estimate[2])+fit_sp$estimate[1]
e_stock_sim = (qt(rand_t_cop[,2], df = fit_stock$estimate[3])*fit_stock$estimate[3])+fit_stock$estimate[1]

sigma_next = sqrt(fit1@fit$par[3]+fit1@fit$par[4]*(res_sd1[length(res_sd1)]^2)+fit1@fit$par[5]*fit1@sigma.t[length(res_sd1)])
sigma_bar_next =sqrt(fit2@fit$par[3]+fit2@fit$par[4]*(res_sd2[length(res_sd1)]^2)+fit2@fit$par[5]*fit2@sigma.t[length(res_sd2)])

X_next = as.numeric(sp[length(res_sd1)])
X_bar_next = as.numeric(stock[length(res_sd1)])

portfolio_sim = (fit1@fit$par[1]+fit1@fit$par[2]*(X_next-fit1@fit$par[1])+sigma_next*e_sp_sim/100)+(fit2@fit$par[1]+fit2@fit$par[2]*(X_bar_next-fit2@fit$par[1])+sigma_bar_next*e_stock_sim/100)
portfolio_sim
summary(portfolio_sim)
VAR= quantile(portfolio_sim, 0.01) ; VAR

tail = portfolio_sim[(portfolio_sim<=VAR)]
expected_shortfall= mean(tail) ;expected_shortfall

#Bootstrapping

boot = rep(0,1000)
for (i in 1:length(boot)){
  for_boot_sp_e = qstd(sample(rand_t_cop[,1], size = 1000, replace =TRUE),est_sp/100)
  boot_sp =(fit1@fit$par[1]+fit1@fit$par[2]*Xn+sigma_next*for_boot_sp_e)
  for_boot_stock_e = qstd(sample(rand_t_cop[,2], size = 1000, replace =TRUE),est_stock/100)
  boot_stock =(fit2@fit$par[1]+fit2@fit$par[2]*Xn_bar+sigma_bar_next*for_boot_stock_e)
  boot_port = boot_sp+ boot_stock
  boot[i] = -quantile(boot_port, alpha)
  
}
lower=mean(boot)-sd(boot) 
upper=mean(boot)+sd(boot)


