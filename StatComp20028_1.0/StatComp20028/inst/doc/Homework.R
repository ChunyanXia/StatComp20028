## ----fig.width=8,fig.height=4-------------------------------------------------
n <- 1000
u <- runif(n)
x <- 2/((1-u)^(1/2)) # F(x) = 1-4/(x^2), x>=2
hist(x, prob = TRUE, main = expression(f(x)==8*x^{-3}),xlim = range(0,50)  )
y <- seq(2, 50, .01)
lines(y, 8*y^(-3)) 

## ----fig.width=8,fig.height=4-------------------------------------------------
n <- 1000
fe <- c()
x <- runif(n, min = -1, max = 1)
y <- runif(n, min = -1, max = 1)
z <- runif(n, min = -1, max = 1)
for (i in 1 :1000){
  if(abs(z[i])>=abs(y[i]) && abs(z[i])>=abs(x[i])){fe[i] <- y[i]}else{fe[i] <- z[i]}}
hist(fe, prob = TRUE,main = (expression(f(x)==3/4*(1-x^2))))
u <- seq(-1,1, .01)
lines(u,3/4*(1-u^2))

## ----fig.width=8,fig.height=4-------------------------------------------------
n <- 1000  
u <- runif(n)
x <- 2/((1-u)^(1/4))-2 # F(y) = 1-(2/(2+y))^r, y>=0
hist(x, prob = TRUE, main = expression(f(y)==64/((2+y)^5)))
y <- seq(0, 100, .01)
lines(y, 64/((2+y)^5)) 

## -----------------------------------------------------------------------------
m <- 1e5 ; x <- runif(m, min=0,max=pi/3)  #generate the random number
esti <- mean(pi/3*sin(x))   #calculate the estimate 
print(c(esti,1-cos(pi/3)))  # give the comparision of the estimate and the true value

## -----------------------------------------------------------------------------
MC.Anti <- function( Q = 10000, anti = TRUE) {
u <- runif(Q/2)
if (!anti) v <- runif(Q/2) else
v <- 1 - u
u <- c(u, v)
g <- mean(exp(u))
g
}  #define a function for generating random number
m <- 1000
SMC <- AN <- numeric(m)
for (i in 1:m) {
SMC[i] <- MC.Anti( Q = 1000, anti = FALSE)
AN[i] <- MC.Anti( Q = 1000)
}   #SMC from simple Monte Carlo method and AN from antithetic variable approach 

print(sd(SMC))
print(sd(AN))
print((var(SMC) - var(AN))/var(SMC))


## -----------------------------------------------------------------------------
size <-10000 
the_est <- Se_the_est <- numeric(2) 
gx_function <- function(x) { (x^2)*exp(-x^2/2)/sqrt(2*pi)*(x>1) }

x1 <-rexp(size,1)   #f1=exp(-x)
f1gx <- gx_function(x1) / exp(-x1)
the_est[1] <- mean(f1gx)
Se_the_est[1] <- sd(f1gx)

x2 <- rlnorm(size,1) #f2=exp(-(log(x))^2/2)/(x*sqrt(2*pi))
f2 <- function(x) { exp(-(log(x))^2/2)/(x*sqrt(2*pi))}
f2gx <- gx_function(x2) / f2(x2)
the_est[2] <- mean(f2gx)
Se_the_est[2] <- sd(f2gx)

print(c(the_est, Se_the_est))

## -----------------------------------------------------------------------------
M <- 10000; 
k <- 5  # number of intervals
N <- 50    # size of stratified sampling
T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x) {exp(-x-log(1+x^2))*(x>0)*(x<1)}
for (i in 1:N) {
  u <- runif(M)
  u1 <- -log(1-u*(1-exp(-1)))
  fg <- g(u1)/(exp(-u1)/(1-exp(-1)))
  est[i, 1] <- mean(fg)
  for(j in 1:k) {
    u <- runif(M)
    x <- -log((exp(-(j-1)/5)-u*(exp(-(j-1)/5)-exp(-j/5))))
    f <-function(x) exp(-x)/((exp(-(j-1)/5)-exp(-j/5)))
    T2[j]<-mean(g(x)/(f(x)))
  }
    T3 <- sum(T2)
    est[i, 2] <- mean(T3)
}
apply(est,2,mean)
apply(est,2,sd)


## -----------------------------------------------------------------------------
set.seed(12345)
size <- 50  # size of estimate
alp <- 0.05   # significance level
fx <- function(x) {
  Lcl =mean(x)-sd(x)*qt(1-alp/2,size-1)/sqrt(size) 
  Ucl = mean(x)+sd(x)*qt(1-alp/2,size-1)/sqrt(size)  
  return(c(Lcl,Ucl))
}     # give the expression of confidence interval
g <- replicate(1000, expr = {x <- log(rlnorm(size))
fx(x)})   # calculate confidence intervals
mean((0<g[2,]) & (0>g[1,]))  # check confidence intervals

## -----------------------------------------------------------------------------
size <- 50  # size of estimate
alp <- 0.05   # significance level
fx <- function(x) {
  Lcl =mean(x)-sd(x)*qt(1-alp/2,size-1)/sqrt(size) 
  Ucl = mean(x)+sd(x)*qt(1-alp/2,size-1)/sqrt(size)  
  return(c(Lcl,Ucl))
}     # give the expression of confidence interval
g <- replicate(1000, expr = {x <- rchisq(size,2)
fx(x)})   # calculate confidence intervals
mean((2<g[2,]) & (2>g[1,]))  # check confidence intervals

## -----------------------------------------------------------------------------
size <-10000 
the_est <- Se_the_est <- numeric(2) 
gx_function <- function(x) { (x^2)*exp(-x^2/2)/sqrt(2*pi)*(x>1) }

x1 <-rexp(size,1)   #f1=exp(-x)
f1gx <- gx_function(x1) / exp(-x1)
the_est[1] <- mean(f1gx)
Se_the_est[1] <- sd(f1gx)

x2 <- rlnorm(size,1) #f2=exp(-(log(x))^2/2)/(x*sqrt(2*pi))
f2 <- function(x) { exp(-(log(x))^2/2)/(x*sqrt(2*pi))}
f2gx <- gx_function(x2) / f2(x2)
the_est[2] <- mean(f2gx)
Se_the_est[2] <- sd(f2gx)

print(c(the_est, Se_the_est))

## -----------------------------------------------------------------------------
M <- 10000; 
k <- 5  # number of intervals
N <- 50    # size of stratified sampling
T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x) {exp(-x-log(1+x^2))*(x>0)*(x<1)}
for (i in 1:N) {
  u <- runif(M)
  u1 <- -log(1-u*(1-exp(-1)))
  fg <- g(u1)/(exp(-u1)/(1-exp(-1)))
  est[i, 1] <- mean(fg)
  for(j in 1:k) {
    u <- runif(M)
    x <- -log((exp(-(j-1)/5)-u*(exp(-(j-1)/5)-exp(-j/5))))
    f <-function(x) exp(-x)/((exp(-(j-1)/5)-exp(-j/5)))
    T2[j]<-mean(g(x)/(f(x)))
  }
    T3 <- sum(T2)
    est[i, 2] <- mean(T3)
}
apply(est,2,mean)
apply(est,2,sd)


## -----------------------------------------------------------------------------
set.seed(12345)
size <- 50  # size of estimate
alp <- 0.05   # significance level
fx <- function(x) {
  Lcl =mean(x)-sd(x)*qt(1-alp/2,size-1)/sqrt(size) 
  Ucl = mean(x)+sd(x)*qt(1-alp/2,size-1)/sqrt(size)  
  return(c(Lcl,Ucl))
}     # give the expression of confidence interval
g <- replicate(1000, expr = {x <- log(rlnorm(size))
fx(x)})   # calculate confidence intervals
mean((0<g[2,]) & (0>g[1,]))  # check confidence intervals

## -----------------------------------------------------------------------------
size <- 50  # size of estimate
alp <- 0.05   # significance level
fx <- function(x) {
  Lcl =mean(x)-sd(x)*qt(1-alp/2,size-1)/sqrt(size) 
  Ucl = mean(x)+sd(x)*qt(1-alp/2,size-1)/sqrt(size)  
  return(c(Lcl,Ucl))
}     # give the expression of confidence interval
g <- replicate(1000, expr = {x <- rchisq(size,2)
fx(x)})   # calculate confidence intervals
mean((2<g[2,]) & (2>g[1,]))  # check confidence intervals

## -----------------------------------------------------------------------------
library(boot); data("law", package = "bootstrap");
n <- 15
b.cor <- function(x,i) cor(x[i,1],x[i,2])
x <- matrix(c(law$LSAT, law$GPA),15,2)

theta.hat <- b.cor(x,1:n)
theta.jack <- numeric(n)

for(i in 1:n){
theta.jack[i] <- b.cor(x,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),3)

## -----------------------------------------------------------------------------
library(boot); data(aircondit, package = "boot");set.seed(1234)
boot.mean <- function(x,i) mean(x[i])
n <-12; m <- 1e2 
x <- aircondit$hours; ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,1,2)
for (i in 1: n){
  de <- boot(data=x,statistic=boot.mean, R = 999)
  ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
}
ci


## -----------------------------------------------------------------------------
library(boot);data("scor", package = "bootstrap")
n <- 88
# define correlation and eigen function to simplify calculation
b.eig <- function(x,i) eigen(cor(x))$values
eig.hat <- b.eig(scor)
theta.hat <- max(eig.hat)/sum(eig.hat)

# bootstrap method
eig.jack <- matrix(NA, n,5)
theta.jack <- numeric(n)
for(i in 1:n){
  score <- scor[(1:n)[-i], ]
  eig.jack[i, ] <- b.eig(score)
  theta.jack[i] <- max(eig.jack[i, ])/sum(eig.jack[i, ])
}

bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),5)


## ----message= FALSE-----------------------------------------------------------
library(DAAG);  attach(ironslag)  #data from DAAG ironslag
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)

# fit models on leave-two-out samples of n-fold cross validation
for (k in 1:(n-1)) {
  for (i in (k+1):n) {
    y <- magnetic[-c(k,i)]
    x <- chemical[-c(k,i)]
    
  # the linear model  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1]+J1$coef[2]*chemical[c(k,i)]
  e1[k,i] <- mean((magnetic[c(k,i)]-yhat1)^2)
  
  # the quadratic model
  J2 <- lm(y~x+I(x^2))
  yhat2 <- J2$coef[1]+J2$coef[2]*chemical[c(k,i)]+
    J2$coef[3]*chemical[c(k,i)]^2
  e2[k,i] <- mean((magnetic[c(k,i)]-yhat2)^2)
  
  # the exponential model
  J3 <- lm(log(y)~x)
  logyhat3 <- J3$coef[1]+J3$coef[2]*chemical[c(k,i)]
  yhat3 <- exp(logyhat3)
  e3[k,i] <- mean((magnetic[c(k,i)]-yhat3)^2)
  
  # log-log model
  J4 <- lm(log(y)~log(x))
  logyhat4 <- J4$coef[1]+J4$coef[2]*log(chemical[c(k,i)])
  yhat4 <- exp(logyhat4)
  e4[k,i] <- mean((magnetic[c(k,i)]-yhat4)^2)
  }
}
  # compute the MSE of models
c(2*sum(e1)/(n*(n-1)),2*sum(e2)/(n*(n-1)),2*sum(e3)/(n*(n-1)),2*sum(e4)/(n*(n-1)))

## -----------------------------------------------------------------------------
 count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
 }
n1 <- 20;   n2 <- 30;  set.seed(1234)
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000
K <- 1:50

tests <- replicate(m, expr = {
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
z <- c(x,y)
kx <- sample(K, size = n1, replace = FALSE)
ky <- sample(K, size = n2, replace = FALSE)

x1 <- z[kx];y1 <- z[-ky]
x1 <- x1 - mean(x1) 
y1 <- y1 - mean(y1)
count5test(x1, y1)
} )
alphahat <- mean(tests)
alphahat 


## -----------------------------------------------------------------------------
set.seed(123)
library(RANN); library(energy);library(Ball)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}

#  1.For unequal variances and equal expectations
m <- 1e2; k<-3; p<-2; mu <- 0.3
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)

for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1),ncol=p);
  y <- cbind(rnorm(n2,0,1.2),rnorm(n2,0,1.2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
  p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow1 <- colMeans(p.values<alpha)
pow1


#  2.For unequal variances and unequal expectations

m <- 1e2; k<-3; p<-2; mu <- 0.3
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)

for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1),ncol=p);
  y <- cbind(rnorm(n2,0.1,1.2),rnorm(n2,0.1,1.2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
  p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow2 <- colMeans(p.values<alpha)
pow2


#  3. Non-normal distributions(t distribution with 1 df (heavy-tailed distribution) and bimodel distribution (mixture of two normal distributions)) 

m <- 1e2; k<-3; p<-2; mu <- 0.3
n1 <- n2 <- 10; R<-999; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)

for(i in 1:m){
 x <- matrix(rt(n1*p,1),ncol = p);
  y <- cbind(rnorm(n2),rnorm(n2,mean = 1))
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
  p.values[i,3] <-bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow3 <- colMeans(p.values<alpha)
pow3



#  4. For unbalanced samples (say, 1 case versus 10 controls) 
m <- 100; k<-3; p<-2; mu <- 0.3; set.seed(345)
n1 <- 30; n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,3.5),ncol=p);
y <- matrix(rnorm(n2*p,0,2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*1)$p.value
}
alpha <- 0.1;
pow4 <- colMeans(p.values<alpha)
pow4


## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(abs(x[i-1])-abs(y)))
x[i] <- y else {
x[i] <- x[i-1]
k <- k + 1
} }
return(list(x=x, k=k))
}



N <- 2000
sigma <- c(.05, .25, 1, 5)
x0 <- 20
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

print(c(rw1$k, rw2$k, rw3$k, rw4$k))

## -----------------------------------------------------------------------------
# function of Gelman.Rubin
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}


laplace.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for Normal(0,1)
#with Normal(X[t], sigma) proposal distribution
#and starting value X1
x <- numeric(N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
  y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(abs(x[i-1])-abs(y))) x[i] <- y else
x[i] <- x[i-1]
}
return(x)
}



GR_method_analysis <- function(sigma) {
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 1000 #burn-in length

#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- laplace.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
} 

x0 <- c(-10, -5, 5, 10)
sigma_.05 <- GR_method_analysis(0.05)
sigma_.25 <- GR_method_analysis(0.25)
sigma_1 <- GR_method_analysis(1)
sigma_5 <- GR_method_analysis(5)


## -----------------------------------------------------------------------------
ck <- function(k,a){
  sqrt(a^2*k/(k+1-a^2))
}
equ <- function(k,a){
  pt(ck(k-1,a),df = k-1) - pt(ck(k,a),df = k)
}
root.curve <- sapply(c(4:25,100,500,1000),function(k){uniroot(equ,interval = c(0.1,sqrt(k)-0.1),k=k)$root})

root.curve

## -----------------------------------------------------------------------------
set.seed(1234)
n_a. <- 444
n_b. <- 132
n_ab <- 63
n_oo <- 361
p0 <- runif(1,0,1)
q0 <- runif(1,0,1-p0)

likelihood_e <- function(prob,p0,q0){
  r0 <- 1-p0-q0 
  p <- prob[1]; q <- prob[2] ; r <- 1-p-q
  - n_a. * (2*log(p)*(p0^2/(p0^2+2*p0*r0)) + log(2*p*r)*(2*p0*r0/(p0^2+2*p0*r0))) -
    n_b. * (2*log(q)*(q0^2/(q0^2+2*q0*r0)) + log(2*q*r)*(2*q0*r0/(q0^2+2*q0*r0))) -
    n_ab * log(2*p*q) - 2*n_oo * log(r) 
}

iter <- 0;E1 <- 0;E2 <- 1; N <- 10
likelihood_trace <- p_maxtrace <- q_maxtrace <- numeric(N)
while(iter < N && abs(E1-E2)> 1e-5){
  output <- optim(par = c(0.1,0.1),likelihood_e,p0 = p0,q0 = q0)
  E1 <- E2;E2 <- output$value
  p0 <- output$par[1]
  q0 <- output$par[2]
  iter <- iter + 1
  p_maxtrace[iter+1] <- output$par[1]
  q_maxtrace[iter+1] <- output$par[2]
  likelihood_trace[iter+1] <- output$value
}
estimate <- data.frame(p0,q0,iter)
colnames(estimate) <- c("p","q","iteration times")
knitr::kable(estimate)
maxtrace <- data.frame(p_maxtrace,q_maxtrace,-likelihood_trace)
colnames(estimate) <- c("p_maxtrace","q_maxtrace","ikelihood_trace")
knitr::kable(maxtrace)


## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp +wt,
  mpg ~ I(1 / disp) + wt
)

# lapply (2 versions)
la1 <- lapply(formulas, lm, data = mtcars) 
la2 <- lapply(formulas, function(x) lm(formula = x, data = mtcars)) 
# for loop
lf1 <- vector("list", length(formulas)) 
for (i in seq_along(formulas)){ 
 lf1[[i]] <- lm(formulas[[i]], data = mtcars) 
}


## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10,10),rpois(7,10)),
  simplify = FALSE
)

round(sapply(1:100,function(i){trials[[i]]$p.value}),3)

# without anonymous function
round(sapply(trials,"[[","p.value"),3)


## -----------------------------------------------------------------------------
example_list <- list(iris, mtcars, cars) 
lapply(example_list, function(x) vapply(x, mean, numeric(1)))
map_vapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){ 
 out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X) 
 if(simplify == TRUE){return(simplify2array(out))} 
 out 
} 
map_vapply(example_list, mean, numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  
#  #include <cmath>
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  //[[Rcpp::export]]
#  double f(double x) {
#    return exp(-abs(x));
#  }
#  
#  //[[Rcpp::export]]
#  NumericVector rwMetropolis (double sigma, double x0, int N) {
#    NumericVector x(N);
#    x[0] = x0;
#    NumericVector u = runif(N);
#    for (int i = 1; i < N;i++ ) {
#      NumericVector y = rnorm(1, x[i-1], sigma);
#      if (u[i] <= (f(y[0]) / f(x[i-1]))){
#        x[i] = y[0];
#      }
#      else {
#        x[i] = x[i-1];
#      }
#    }
#    return(x);
#  }
#  

## ---- eval=FALSE--------------------------------------------------------------
#      library(Rcpp)
#      library(microbenchmark)
#      # R
#      lap_f = function(x) exp(-abs(x))
#  
#      rw.Metropolis = function(sigma, x0, N){
#      x = numeric(N)
#      x[1] = x0
#      u = runif(N)
#      k = 0
#      for (i in 2:N) {
#      y = rnorm(1, x[i-1], sigma)
#      if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y
#      else {
#      x[i] = x[i-1]
#      k = k+1
#       }
#      }
#       return(list(x = x, k = k))
#      }
#  
#      dir_cpp = 'F:/Rcpp/StatComp20028/src/'
#      sourceCpp(paste0(dir_cpp,"rwMetropolis.cpp"))
#      x0 = 25
#      N = 2000
#      sigma = 2
#      (time = microbenchmark(rwR=rw.Metropolis(sigma,x0,N),rwC=rwMetropolis(sigma,x0,N)))

## ----eval=FALSE---------------------------------------------------------------
#  
#  set.seed(12345)
#  rwR = rw.Metropolis(sigma,x0,N)$x[-(1:500)]
#  rwC = rwMetropolis(sigma,x0,N)[-(1:500)]
#  qqplot(rwR,rwC)
#  abline(a=0,b=1,col='black')

