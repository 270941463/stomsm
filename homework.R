##problem 1
N=500
x<-rep(0,N)
c=4#choose  sup(f(x)/g(x))<c)
n=0;m=0
while (m<=500)
{
y<-runif(1)
u<-runif(1)
if (u<=dbeta(y,2,9)/c)
  {
  x[m]<- y
  m<-m+1
}
n=n+1
} 
acc=m/n#acceptance rate
hist(x,breaks = 30)#histgram
plot(density(x))#denstiy curve
curve(dbeta(x,2,9),col=2,add=T)#the line of the beta(2,9)

##problem 2
#produce de(0,1)xde(0,1)
install.packages("fMultivar")
library(fMultivar)
N=1000
x<-matrix(0,ncol=2,nrow=N)
c=2#choose  #sup f(x)/g(x)~0.25/0.183
n=0;m=0
while (m<=1000){
rde<-rep(0,2)
for (i in 1:2){
u<-runif(1)
rde[i]<-ifelse(u>=1/2,-log(2*(1-u)),log(2*u))} 
ud<-runif(1)
if (ud<=dnorm2d(rde[1],rde[2],rho = 0.5)/
    (1/4*exp(-abs(rde[1])-abs(rde[2])))  )
{
  x[m,]<- rde
  m<-m+1
}
n=n+1
} 
acc=m/n#acceptance rate
cor(x)#correlation matrix o fsample
###problem 3
x<-NULL
x[1]<-tan(pi/2*runif(1,min=-1,max=1))#choose x[1] from -Inf to Inf
#burn in 
N1<-3000
for (i in 1:N1){
u<-runif(1)
y<-ifelse(u>=1/2,x[i]-log(2*(1-u)),x[i]+log(2*u))
alpha<-min(dnorm(y)/dnorm(x),1)
t<-runif(1)
x[i+1]<-ifelse(t>alpha,x[i],y)
}
while(mean(x[i:(i-N1/3)])-mean(x[(i-N1/3):(i-2*N1/3)])>0.0001){
  u<-runif(1)
  y<-ifelse(u>=1/2,x[i]-log(2*(1-u)),x[i]+log(2*u))
  alpha<-min(dnorm(y)/dnorm(x),1)
  t<-runif(1)
  x[i+1]<-ifelse(t>alpha,x[i],y)
  i<-i+1
}
##sample 1000
N2=2000
for (j in i:(i+N2)){
  u<-runif(1)
  y<-ifelse(u>=1/2,x[j]-log(2*(1-u)),x[j]+log(2*u))
  alpha<-min(dnorm(y)/dnorm(x),1)
  t<-runif(1)
  x[j+1]<-ifelse(t>alpha,x[j],y)
}
sample<-x[j:(j-999)]
summary(sample)
var(sample)#sample variance
plot(density(sample))
curve(dnorm(x),add=T,col=2)#standard norm curve
qs<-quantile(sample,c(0.95,0.975))
##problem 4
#(1)
sigma1 <- NULL
sigma11 <- NULL
result1 <- NULL
result11 <- NULL
N1 <- NULL
for (k in 1:6) {
  N1[k] = 10 ^ k
  x <- rnorm(N1[k])
  #M<-length(x[x<=1])
  M <- as.numeric(x <= 1)
  sigma1[k] <- sd(M) / (sqrt(N1[k]))
  result1[k] <- mean(M)
  #result1[k]<-M/N1[k]
  #sigma1[k]<-sqrt(M*(1-M/N1[k])/(N1[k]^2))
}
cbind(N1,result1,sigma1)
#(2)
sigma2 <- NULL
result2 <- NULL
N2 <- NULL
for (k in 1:6) {
  N2[k] = 10 ^ k
  x <- rnorm(N2[k])
  # newdata <- subset(x, x >1 | x < -1)
  M <- as.numeric(x > 1 | x < -1)
  result2[k] <- 1 - mean(M / 2)
  sigma2[k] <- sd(M / 2) / sqrt(N2[k])
}
cbind(N2,result2,sigma2)
#3
sigma3 <- NULL
result3 <- NULL
N3 <- NULL
for (k in 1:6) {
  N3[k] = 10 ^ k
  x <- rnorm(N3[k])
  #newdata <- subset(x, x >0 & x < 1)
  M <- as.numeric(x > 0 & x < 1)
  result3[k] <- 1 / 2 + mean(M)
  sigma3[k] <- sd(M) / sqrt(N3[k])
}
cbind(N3,result3,sigma3)
#4
sigma4 <- NULL
result4 <- NULL
N4 <- NULL
for (k in 1:6) {
  N4[k] = 10 ^ k
  x <- runif(N4[k])
  M <- dnorm(x)
  result4[k] <- mean(M) + 0.5
  sigma4[k] <- sd(M) / sqrt(N4[k])
}
cbind(N4,result4,sigma4)

#5 important sampling
sigma5 <- NULL
result5 <- NULL
N5 <- NULL
A <- exp(1) / (exp(1) - 1)
for (k in 1:6) {
  N5[k] = 10 ^ k
  x <- runif(N5[k])
  y <- -log(1 - x / A)
  M <- dnorm(y) / (A * exp(-y))
  result5[k] <- mean(M) + 0.5
  sigma5[k] <- sd(M) / sqrt(N5[k])
}
cbind(N5,result5,sigma5)
#6 antithetic variates
sigma6 <- NULL
result6 <- NULL
N6 <- NULL
for (k in 1:6) {
  N6[k] = 10 ^ k
  x <- runif(N6[k] / 2)
  y <- 1 - x#antithetic
  Mx <- dnorm(x)
  My <- dnorm(y)
  result6[k] <- mean(Mx + My) / 2 + 0.5
  sigma6[k] <- sqrt(var(Mx + My) / (2 * N6[k]))
}
cbind(N6,result6,sigma6)
#7 antithetic variates+control variates
sigma7 <- NULL
result7 <- NULL
N7 <- NULL
lamda <- NULL
r <- function(x) {
  -x ^ 2 + x
}
u <- 1 / 6
for (k in 1:6) {
  N7[k] = 10 ^ k
  x <- runif(N7[k] / 2)
  y <- 1 - x#antithetic
  Mx <- dnorm(x)
  rx <- r(x)
  lamda[k] <- sum((Mx - mean(Mx)) * (rx - u)) / sum((rx - u) ^ 2)
  wx <- Mx - lamda[k] * (rx - u)
  #varwx<-var(Mx)-2*lamda[k]*cov(Mx,rx)+lamda[k]^2*var(rx)
  My <- dnorm(y)
  ry <- r(y)
  #lamda<-sum((My-mean(My))*(ry-u))/sum((ry-u)^2)
  wy <- My - lamda[k] * (ry - u)
  #varwy<-var(My)-2*lamda*cov(My,rx)+lamda^2*var(ry)
  result7[k] <- mean(wx + wy) / 2 + 0.5
  sigma7[k] <- sqrt(var(wx + wy) / (2 * N7[k]))
}
cbind(N7,result7,sigma7)
#Reserved 5 significant digits
round(cbind(
  N1,result7,sigma7,result6,sigma6,result5,sigma5,result4,sigma4
),5)
round(cbind(N1,result3,sigma3,result2,sigma2,result1,sigma1),5)

##problem 5
N <- NULL
MCresult <- NULL
MCsigma <- NULL
RBresult <- NULL
RBsigma <- NULL
for (k in 1:7) {
  N[k] <- 10 ^ k
  ##MC-extimate
  MClnx <- rnorm(N[k])
  #lny<--2+1.5*MClnx+rnorm(N[k])
  MCy <- exp(-2 + 1.5 * MClnx + rnorm(N[k]))
  #RBy<-exp(rnorm(N[k],mean=-2+1.5*lnx,sd=1))
  MCresult[k] <- mean(MCy)
  MCsigma[k] <- sd(MCy) / sqrt(N[k])
  #RB-extimate
  RBlnx <- rnorm(N[k])
  RBy <- exp(-2 + 1.5 * RBlnx)*mean(exp(rnorm(N[k])))
  RBresult[k] <- mean(RBy)
  RBsigma[k] <- sd(RBy) / sqrt(N[k])
}
cbind(N,MCresult,MCsigma,RBresult,RBsigma)
##problem 6
#1)estimate beta and the variance
setwd("G:/2015-2016xueyuan/课程/现代统计方法")
data <- read.table("Ex6.txt")
da <- as.data.frame(t(as.matrix(data[,-1])))
colnames(da) <- c("x","y")
#n<-dim(da)[1]
fit <- lm(y ~ x,data = da)
#beta
beta_hat<-fit$coefficients
##variance of beta
f<-summary(fit)
beta_hat_cov<-f$sigma^2*f$cov.unscaled
#caculate by formula
#X<-cbind(1,da[,1])
#sig<-sum((da[,2]-X%*%as.matrix(beta))^2)/(n-2)
#varb<-solve(t(X)%*%X)*sig

#2) theta and variance
theta_hat<-beta_hat[2]/beta_hat[1]
#delta method
delta_theta<-matrix(c(-beta_hat[2]/beta_hat[1]^2,1/beta_hat[1]),2,1)
theta_hat_var<-t(delta_theta)%*%beta_hat_cov%*%delta_theta

#3) bootstrap imitate theta distribution 
##
B<-10000
theta_bt<-NULL
t_bt<-NULL
for (i in 1:B){
  index_x <- sample(1:13, size = 13, replace = T)
  da_x <- da[index_x,1]
  index_y <- sample(1:13, size = 13, replace = T)
  da_y <- da[index_y,2]
  fit_bt <- summary(lm(da_y~da_x))
  beta_bt<-fit_bt$coefficients[,1]
  #theta_bt_hat
  theta_bt[i]<-beta_bt[2]/beta_bt[1]
  #theta_bt_var
  beta_bt_cov<-fit_bt$sigma^2*fit_bt$cov.unscaled
  delta_bt_theta<-matrix(c(-beta_bt[2]/beta_bt[1]^2,1/beta_bt[1]),2,1)
  theta_bt_var<-t(delta_bt_theta)%*%beta_bt_cov%*%delta_bt_theta
  t_bt[i]<-(theta_bt[i] - theta_hat)/sqrt(theta_bt_var)
}
qplot(theta_bt, geom="histogram",col=3)
##confidence interval
con_l<-quantile(theta_bt, 0.025)
con_u<-quantile(theta_bt, 0.975)
cat("[",con_l, ",",con_u,"]")

## t interval
t_l<-theta_hat+sqrt(theta_hat_var)*quantile(t_bt, 0.025)
t_u<-theta_hat+sqrt(theta_hat_var)*quantile(t_bt, 0.975)
cat("[",t_l, ",",t_u,"]")

quantile(theta_hat + t_bt*sqrt(theta_hat_var), 0.025)
quantile(theta_hat + t_bt*sqrt(theta_hat_var), 0.975)
##bias and bootstrap corrections
bias_bt <- mean(theta_bt) - theta_hat
theta_bt_correction <- 2 * theta_hat -bias_bt
