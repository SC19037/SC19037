## ---- results = 'asis'--------------------------------------------------------
library(DAAG)
library(xtable)
data(allbacks)
print(xtable(allbacks[1:15, ]),type='html')

## -----------------------------------------------------------------------------
attach(allbacks)
plot(volume,weight)
detach(allbacks)

## -----------------------------------------------------------------------------
n <- 100000
u <- runif(n)
m <- c(2,4,6,8)
#par(mfrow=c(2,2))
for(sigma in m){
x <- sqrt(-2*(sigma^2)*log(1-u))
#density histogram of sample
hist(x, prob = TRUE, main = paste("sigma = ",sigma))
y <- seq(0, 60, .001)
#density curve f(x)
lines(y, y/(sigma^2)*exp(-(y^2)/(2*sigma^2))) 
}

## -----------------------------------------------------------------------------
n <- 1e6
l <- c(0.1,0.3,0.5,0.6,0.75,0.9)
#par(mfrow=c(2,3))
 for(p in l){
  x <- sample(c(0,1), size = n, replace = TRUE, prob = c(p, 1-p))
  y <- vector(length=n)
  for(i in 1:n){
    if(x[i] == 0) 
    y[i] <- rnorm(1, mean = 0, sd = 1)
    else y[i] <- rnorm(1, mean = 3, sd = 1)
  }
 hist(y, prob = TRUE, main = paste("p = ", p))
}

## -----------------------------------------------------------------------------
library(MASS)
n <- 10
d <-5
m <- matrix(1:d^2, d, d)
sigma <- cov(m)
X <- mvrnorm(d, mu=1:d, sigma)
Y <- t(X)
Z <- Y%*%X
Z

## -----------------------------------------------------------------------------
library(scales)
m <- 10000
x <- runif(m, min = 0, max = pi/3)
theta.hat <- mean(sin(x))*(pi/3)
print(theta.hat)
print(-cos(pi/3)+1)

## -----------------------------------------------------------------------------
library(scales)
m <- 10000
U <- runif(m)
#simple MC
T1 <- exp(-U)/(1+U^2)
#MC with antithetic variables
U1 <- runif(m/2)
U2 <- 1-U1
T2 <- (exp(-U1)/(1+U1^2)+exp(-U2)/(1+U2^2))/2
#print the variance of T1 and T2
mean(T1)
mean(T2)
var(T1)
var(T2)
(var(T1)-var(T2))/var(T1)

## -----------------------------------------------------------------------------
M <- 1e5 
k <- 5 
N <- 50
T2 <- numeric(k)
est <- matrix(0, N, 2)
g <- function(x){
  (1-exp(-1))/(1+x^2)
}
for (i in 1:N) {
#importance sampling estimate  
  u<-runif(M)
  x<--log(1-(1-exp(-1))*u)
  est[i, 1] <- mean(g(x))
  for(j in 1:k){
#stratified importance sampling estimate    
    F1<-(1-exp((1-j)/5))/(1-exp(-1))
    F2<-(1-exp(-j/5))/(1-exp(-1))
    u<- runif((F2-F1)*M)
    x<- -log(1-(1-exp(-1))*((F2-F1)*u+F1))
    T2[j]<-mean((F2-F1)*g(x))
  }
  est[i, 2] <- sum(T2)
}
#print the outcoming
apply(est,2,mean)
apply(est,2,var)
apply(est,2,sd)

## -----------------------------------------------------------------------------
#x is from a normal distribution 
library(scales)
n <- 20
alpha <- 0.05
UCL1 <- replicate(1000, expr={
  x <- rnorm(n, mean=0, sd=2)
  (n-1)*var(x) / qchisq(alpha, df=n-1)
})
sum(UCL1 > 4)
mean(UCL1 > 4)

## -----------------------------------------------------------------------------
#x is from a chi-square data.
n <- 20 
alpha <- 0.05 
UCL2 <- replicate(1000, expr = { 
  x <- rchisq(n, df = 2) 
  (n-1) * var(x) / qchisq(alpha, df = n-1) 
})
sum(UCL2 > 4)
mean(UCL2 > 4) 

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05
U1 <- replicate(1000, expr={
  x <- rchisq(n, df = 2)
  mean(x) - sqrt(var(x)/n)* qt(alpha/2, df = n-1)
})
U2 <- replicate(1000, expr={
  x <- rchisq(n, df = 2)
  mean(x) + sqrt(var(x)/n)* qt(alpha/2, df = n-1)
})
mean(U1 > 2 & U2 < 2)

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05
U3 <- replicate(1000, expr={
  x <- rnorm(n, mean = 0, sd = 2 )
  mean(x) - sqrt(var(x)/n)* qt(alpha/2, df = n-1)
})
U4 <- replicate(1000, expr={
  x <- rnorm(n, mean = 0, sd = 2 )
  mean(x) + sqrt(var(x)/n)* qt(alpha/2, df = n-1)
})
mean(U3 > 0 & U4 < 0)

## -----------------------------------------------------------------------------
sk.chisq = 4/sqrt(10)
for(i in 1:m){ 
  U <- rchisq(n, df = 5)
  de <- boot(data = U, statistic = sk, B) 
  ci <- boot.ci(de, type = c("norm", "basic", "perc")) 
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5] 
  ci.perc[i,] <- ci$percent[4:5]
  } 
cat('norm =', mean(ci.norm[,1] <= sk.chisq & ci.norm[,2] >= sk.chisq), 'basic =', mean(ci.basic[,1] <= sk.chisq & ci.basic[,2] >= sk.chisq), 'perc =', mean(ci.perc[,1] <= sk.chisq & ci.perc[,2] >= sk.chisq))

## -----------------------------------------------------------------------------
# To obtain the parameter theta
library(bootstrap)
library(boot) #for boot function
n <- nrow(scor)
theta.hat <- function(x){
  y <- cov(x)
  lambda <- eigen(y)$values
  theta <- lambda[1] / sum(lambda)
}

## -----------------------------------------------------------------------------
set.seed(0)
#set up the bootstrap
B <- 2000
n <- nrow(scor)
theta.boot <- numeric(B)

#bootstrap
for(b in 1 : B){
  i <- sample(1:n, size = n, replace = TRUE)
  x <- scor[i,]
  theta.boot[b] <- theta.hat(x)
}

# Estimate bias and standard error of theta.hat
bias.hat_boot <- mean(theta.boot) - theta.hat(scor)
se.hat_boot <- sd(theta.boot)
# print the result
bias.hat_boot
se.hat_boot

## -----------------------------------------------------------------------------
theta.jack <- numeric(n)
for(i in 1 : n){
  x <- scor[-i,]
  theta.jack[i] <- theta.hat(x)
}

# Estimate bias and standard error of theta.hat
bias.hat_jack <- (n - 1) * (mean(theta.jack) - theta.hat(scor))
se.hat_jack <- sqrt((n - 1) * mean((theta.jack - mean(theta.jack))^2))
# print the result
bias.hat_jack
se.hat_jack


## -----------------------------------------------------------------------------
library(DAAG)
magnetic <- ironslag$magnetic
chemical <- ironslag$chemical
a <- seq(10, 40, .1) #sequence for plotting fits
#par(mfrow = c(2,2))

L1 <- lm(magnetic ~ chemical) 
plot(chemical, magnetic, main = "Linear", pch = 16) 
yhat1 <- L1$coef[1] + L1$coef[2] * a 
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2)) 
plot(chemical, magnetic, main = "Quadratic", pch = 16) 
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2 
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical) 
plot(chemical, magnetic, main = "Exponential", pch = 16) 
logyhat3 <- L3$coef[1] + L3$coef[2] * a 
yhat3 <- exp(logyhat3) 
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)) 
plot(chemical, magnetic, main = "cubic", pch = 16) 
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
#Estimate the quantiles of the skewness under normality by a Monte Carlo experiment
n <- 1000
m <- 10000
sk <- numeric(n)
for (i in 1:m){
  x <- rnorm(n)
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2)
  sk[i] <- m3 / m2^1.5  
}
y1 <- unname(quantile(sk,c(0.025, 0.05, 0.95, 0.975)))
y1

## -----------------------------------------------------------------------------
#Compute the standard error of the estimates from (2.14) using the normal approximation for the density (with exact variance formula)
n <- 1000
q <- c(0.025, 0.05, 0.95, 0.975)
f <- qnorm(q, mean = 0, sd = sqrt(6*(n-2)/((n+1)*(n+3))))
sd <- sqrt(q*(1-q)/(n*f^2))
sd

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag 
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation 
# fit models on leave-one-out samples 
for (k in 1:n) { 
  y <- magnetic[-k] 
  x <- chemical[-k]
  
  J1 <- lm(y ~ x) 
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2)) 
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x) 
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] 
  yhat3 <- exp(logyhat3) 
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3)) 
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
} 
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)) 

## -----------------------------------------------------------------------------
#estimated quantiles with the quantiles of the large sample approximation N(0,6/n)
n <- 1000
y2 <- qnorm(c(0.025, 0.05, 0.95, 0.975), mean = 0, sd = sqrt(6/n))
y2

## -----------------------------------------------------------------------------
#compare the quantiles
data.frame(quantiles = c(0.025, 0.05, 0.95, 0.975), q_MC_est = y1, q_sample_approximation = y2, sd)

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag 
y <- magnetic
x <- chemical

J1 <- lm(y ~ x) 
summary(J1)$adj.r.squared# to get the adjust R^2

J2 <- lm(y ~ x + I(x^2)) 
summary(J2)$adj.r.squared

J3 <- lm(log(y) ~ x) 
summary(J3)$adj.r.squared

J4 <- lm(y ~ x + I(x^2) + I(x^3)) 
summary(J4)$adj.r.squared

## -----------------------------------------------------------------------------
L2

## -----------------------------------------------------------------------------
# define the function of "count five test"
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # get the statistic of "count five test"
  return(as.integer(max(c(outx, outy))))
}

## -----------------------------------------------------------------------------
# get two samples with the same variance while sample sizes are not equal
set.seed(12345)
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
#centered by sample mean
x <- x - mean(x) 
y <- y - mean(y)

## -----------------------------------------------------------------------------
t0 <- count5test(x,y)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:(n1+n2)
reps <- numeric(R) #storage for replicates
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n1, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- count5test(x1, y1)
}
# print the p-value
mean(reps > t0)

## -----------------------------------------------------------------------------
library(MASS)
library(Ball)
library(boot)

dCov <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d)
    M <- mean(d)
    a <- sweep(d, 1, m)
    b <- sweep(a, 2, m)
    return(b + M)
  }
  A <- Akl(x)
  B <- Akl(y)
  dCov <- sqrt(mean(A * B))
}

ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y
  p <- dims[1]
  q1 <- dims[2] + 1
  d <- p + dims[2]
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y
  return(nrow(z) * dCov(x, y)^2)
}

dc.test <- function(z, dims){
  boot.obj <- boot(data = z, statistic = ndCov2, R = R, sim = "permutation",  dims = c(2,2))
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic=ts[1], p.value=p.value)
}

## -----------------------------------------------------------------------------
set.seed(12345)
m <- 100
u <- c(0,0)
sigma <- matrix(c(1,0,0,1), 2, 2)
R <- 1e2
alpha <- 0.1
n <- seq(10,100,10)
power.dc <- power.ball <- numeric(length(n))
p.dc <- p.ball <- numeric(m)

## -----------------------------------------------------------------------------
alpha <- 1
n <- 100 
m <- 2500 
epsilon <-seq(0, 1, length.out = 30) 
N <- length(epsilon) 
power <- numeric(N) 

#computes the sample skewness coeff. 
sk <- function(x) { 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
  }

#critical value for the skewness test 
cv <- qnorm(0.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

for (j in 1:N) { 
  #for each epsilon 
  e <- epsilon[j] 
  sktests <- numeric(m)
  for (i in 1:m) {
    #for each replicate   
    Y <- sample(0:1, size = n, replace = TRUE,prob = c(1-e, e)) 
    x <- rbeta(n, alpha, alpha) * (Y == 0) + rbeta(n, alpha, 100*alpha) * (Y == 1)
    sktests[i] <- (abs(sk(x)) >= cv)
  }
  power[j] <- mean(sktests) 
} 

#plot power vs epsilon
plot(epsilon, power, type = "b", xlab = bquote(epsilon), ylim = c(0,1)) 
#abline(h = .1, lty = 3)
se <- sqrt(power * (1-power) / m)
#add standard errors
lines(epsilon, power+se, lty = 3)
lines(epsilon, power-se, lty = 3)

## -----------------------------------------------------------------------------
alpha <- 3
n <- 100 
m <- 2500 
epsilon <-seq(0, 1 ,length.out = 30) 
N <- length(epsilon) 
power <- numeric(N) 

#computes the sample skewness coeff. 
sk <- function(x) { 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
}

#critical value for the skewness test 
cv <- qnorm(0.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

for (j in 1:N) { 
  #for each epsilon 
  e <- epsilon[j] 
  sktests <- numeric(m)
  for (i in 1:m) {
    #for each replicate   
    Y <- sample(0:1, size = n, replace = TRUE, prob = c(1-e, e)) 
    x <- rbeta(n, alpha, alpha) * (Y == 0) + rt(n, 100*alpha) * (Y == 1)
    sktests[i] <- (abs(sk(x)) >= cv)
  }
  power[j] <- mean(sktests) 
} 

#plot power vs epsilon
plot(epsilon, power, type = "b", xlab = bquote(epsilon), ylim = c(0,1)) 
#abline(h = .1, lty = 3)
se <- sqrt(power * (1-power) / m)
#add standard errors
lines(epsilon, power+se, lty = 3)
lines(epsilon, power-se, lty = 3)

## -----------------------------------------------------------------------------
n <- 20 
alpha <- .05
mu0 <- 1 #null hypothesis
m <- 10000 #number of replicates 
p <- numeric(m) #storage for p-values 
for (j in 1:m) { 
  x <- rchisq(n, mu0) 
  ttest <- t.test(x, mu = mu0) 
  p[j] <- ttest$p.value 
  }
p.hat <- mean(p < alpha) 
p.hat

## -----------------------------------------------------------------------------
n <- 20 
alpha <- 0.05 
mu0 <- 1 #null hypothesis
m <- 10000 #number of replicates 
p <- numeric(m) #storage for p-values 
for (j in 1:m) { 
  x <- runif(n, min = 0, max = 2) 
  ttest <- t.test(x, mu = mu0) 
  p[j] <- ttest$p.value 
  }
p.hat <- mean(p < alpha) 
p.hat

## -----------------------------------------------------------------------------
n <- 20 
alpha <- .05 
mu0 <- 1  #null hypothesis
m <- 10000 #number of replicates 
p <- numeric(m) #storage for p-values 
for (j in 1:m) { 
  x <- rexp(n, rate = 1) 
  ttest <- t.test(x, mu = mu0) 
  p[j] <- ttest$p.value 
  }
p.hat <- mean(p < alpha) 
p.hat

## -----------------------------------------------------------------------------
library(bootstrap)
#par(mfrow = c(2, 2))
plot(scor$mec, scor$vec)
plot(scor$alg, scor$ana)
plot(scor$alg, scor$sta)
plot(scor$ana, scor$sta)

#get the sample correlation matrix
cor(scor)

## ----results='asis'-----------------------------------------------------------
#Obtain bootstrap estimates of the standard errors
#set up the bootstrap 
B <- 200 #number of replicates
n <- nrow(scor) #sample size
R12 <- R34 <- R35 <- R45 <- numeric(B) #storage for replicates

#bootstrap estimate of standard error of R 
for (b in 1:B) { 
  #randomly select the indices 
  i <- sample(1:n, size = n, replace = TRUE) 
  mec <- scor$mec[i] #i is a vector of indices 
  vec <- scor$vec[i]
  alg <- scor$alg[i]
  ana <- scor$ana[i]
  sta <- scor$sta[i]
  R12[b] <- cor(mec, vec)
  R34[b] <- cor(alg, ana)
  R35[b] <- cor(alg, sta)
  R45[b] <- cor(ana, sta)
}
# to compute standard error by using the bootstrap estimate.
se.R12 <- sd(R12)
se.R34 <- sd(R34)
se.R35 <- sd(R35)
se.R45 <- sd(R45)
# draw a table to present the result.
p <- matrix(c(se.R12, se.R34, se.R35, se.R45), 1)
colnames(p) <- c("Rho_12","Rho_34","Rho_35","Rho_45")
rownames(p) <- "Standard Errors"
knitr::kable(p)

## -----------------------------------------------------------------------------
library(boot)
n <- 10
m <- 1e3
B <- 2000
ci.norm <- ci.basic <- ci.perc <- ci.bca <- matrix(0, m, 2) 

#computes the sample skewness coeff.
sk <- function(x,i) { 
  x<-x[i]
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
}

## -----------------------------------------------------------------------------
sk.norm = 0
for(i in 1:m){ 
  U <- rnorm(n)
  de <- boot(data = U, statistic = sk, B) 
  ci <- boot.ci(de, type = c("norm", "basic", "perc")) 
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5] 
  ci.perc[i,] <- ci$percent[4:5]
  } 
cat('norm =', mean(ci.norm[,1] <= sk.norm & ci.norm[,2] >= sk.norm), 'basic =', mean(ci.basic[,1] <= sk.norm & ci.basic[,2] >= sk.norm), 'perc =', mean(ci.perc[,1] <= sk.norm & ci.perc[,2] >= sk.norm))

## -----------------------------------------------------------------------------
sk.chisq = 4/sqrt(10)
for(i in 1:m){ 
  U <- rchisq(n, df = 5)
  de <- boot(data = U, statistic = sk, B) 
  ci <- boot.ci(de, type = c("norm", "basic", "perc")) 
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5] 
  ci.perc[i,] <- ci$percent[4:5]
  } 
cat('norm =', mean(ci.norm[,1] <= sk.chisq & ci.norm[,2] >= sk.chisq), 'basic =', mean(ci.basic[,1] <= sk.chisq & ci.basic[,2] >= sk.chisq), 'perc =', mean(ci.perc[,1] <= sk.chisq & ci.perc[,2] >= sk.chisq))

## -----------------------------------------------------------------------------
# To obtain the parameter theta
library(bootstrap)
library(boot) #for boot function
n <- nrow(scor)
theta.hat <- function(x){
  y <- cov(x)
  lambda <- eigen(y)$values
  theta <- lambda[1] / sum(lambda)
}

## -----------------------------------------------------------------------------
set.seed(0)
#set up the bootstrap
B <- 2000
n <- nrow(scor)
theta.boot <- numeric(B)

#bootstrap
for(b in 1 : B){
  i <- sample(1:n, size = n, replace = TRUE)
  x <- scor[i,]
  theta.boot[b] <- theta.hat(x)
}

# Estimate bias and standard error of theta.hat
bias.hat_boot <- mean(theta.boot) - theta.hat(scor)
se.hat_boot <- sd(theta.boot)
# print the result
bias.hat_boot
se.hat_boot

## -----------------------------------------------------------------------------
theta.jack <- numeric(n)
for(i in 1 : n){
  x <- scor[-i,]
  theta.jack[i] <- theta.hat(x)
}

# Estimate bias and standard error of theta.hat
bias.hat_jack <- (n - 1) * (mean(theta.jack) - theta.hat(scor))
se.hat_jack <- sqrt((n - 1) * mean((theta.jack - mean(theta.jack))^2))
# print the result
bias.hat_jack
se.hat_jack


## -----------------------------------------------------------------------------
library(DAAG)
magnetic <- ironslag$magnetic
chemical <- ironslag$chemical
a <- seq(10, 40, .1) #sequence for plotting fits
#par(mfrow = c(2,2))

L1 <- lm(magnetic ~ chemical) 
plot(chemical, magnetic, main = "Linear", pch = 16) 
yhat1 <- L1$coef[1] + L1$coef[2] * a 
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2)) 
plot(chemical, magnetic, main = "Quadratic", pch = 16) 
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2 
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical) 
plot(chemical, magnetic, main = "Exponential", pch = 16) 
logyhat3 <- L3$coef[1] + L3$coef[2] * a 
yhat3 <- exp(logyhat3) 
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)) 
plot(chemical, magnetic, main = "cubic", pch = 16) 
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag 
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation 
# fit models on leave-one-out samples 
for (k in 1:n) { 
  y <- magnetic[-k] 
  x <- chemical[-k]
  
  J1 <- lm(y ~ x) 
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2)) 
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x) 
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] 
  yhat3 <- exp(logyhat3) 
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3)) 
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
} 
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)) 

## -----------------------------------------------------------------------------
L2

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag 
y <- magnetic
x <- chemical

J1 <- lm(y ~ x) 
summary(J1)$adj.r.squared# to get the adjust R^2

J2 <- lm(y ~ x + I(x^2)) 
summary(J2)$adj.r.squared

J3 <- lm(log(y) ~ x) 
summary(J3)$adj.r.squared

J4 <- lm(y ~ x + I(x^2) + I(x^3)) 
summary(J4)$adj.r.squared

## -----------------------------------------------------------------------------
L2

## -----------------------------------------------------------------------------
# define the function of "count five test"
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # get the statistic of "count five test"
  return(as.integer(max(c(outx, outy))))
}

## -----------------------------------------------------------------------------
# get two samples with the same variance while sample sizes are not equal
set.seed(12345)
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
#centered by sample mean
x <- x - mean(x) 
y <- y - mean(y)

## -----------------------------------------------------------------------------
t0 <- count5test(x,y)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:(n1+n2)
reps <- numeric(R) #storage for replicates
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n1, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- count5test(x1, y1)
}
# print the p-value
mean(reps > t0)

## -----------------------------------------------------------------------------
for(j in 1:length(n)){
 for(i in 1:m){
    x<-matrix(rnorm(n[j]*2),nrow = n[j])
    e<-matrix(rnorm(n[j]*2),nrow = n[j])
    y <- x/4+e
    z <- as.data.frame(cbind(x,y))
    p.dc[i] <- dc.test(z, c(2,2))$p.value
    p.ball[i] <- bcov.test(x, y, R=R, seed = i*12345)$p.value
  }
  power.dc[j] <- mean(p.dc < alpha)
  power.ball[j] <- mean(p.ball < alpha)
}

# plot the power
plot(n, power.dc, type = "o", xlab = "n", ylab = "power", main = "Model 1", col = 1)
lines(n, power.ball, type = "o", col = 2)
legend('bottomright', legend = c('dc', 'ball'), col = 1:2, lty = 1)

## -----------------------------------------------------------------------------
library(MASS)
library(Ball)
library(boot)

dCov <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d)
    M <- mean(d)
    a <- sweep(d, 1, m)
    b <- sweep(a, 2, m)
    return(b + M)
  }
  A <- Akl(x)
  B <- Akl(y)
  dCov <- sqrt(mean(A * B))
}

ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y
  p <- dims[1]
  q1 <- dims[2] + 1
  d <- p + dims[2]
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y
  return(nrow(z) * dCov(x, y)^2)
}

dc.test <- function(z, dims){
  boot.obj <- boot(data = z, statistic = ndCov2, R = R, sim = "permutation",  dims = c(2,2))
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic=ts[1], p.value=p.value)
}

## -----------------------------------------------------------------------------
set.seed(12345)
m <- 100
u <- c(0,0)
sigma <- matrix(c(1,0,0,1), 2, 2)
R <- 1e2
alpha <- 0.1
n <- seq(10,100,10)
power.dc <- power.ball <- numeric(length(n))
p.dc <- p.ball <- numeric(m)

## -----------------------------------------------------------------------------
for(j in 1:length(n)){
  for(i in 1:m){
    x<-matrix(rnorm(n[j]*2),nrow = n[j])
    e<-matrix(rnorm(n[j]*2),nrow = n[j])
    y <- x/4*e
    z <- as.data.frame(cbind(x,y))
    p.dc[i] <- dc.test(z, c(2,2))$p.value
    p.ball[i] <- bcov.test(x, y, R=R, seed = i*12345)$p.value
  }
  power.dc[j] <- mean(p.dc < alpha)
  power.ball[j] <- mean(p.ball < alpha)
}

# plot the power
plot(n, power.dc, type = "o", xlab = "n", ylab = "power", main = "Model 1", col = 1)
lines(n, power.ball, type = "o", col = 2)
legend('bottomright', legend = c('dc', 'ball'), col = 1:2, lty = 1)

## -----------------------------------------------------------------------------
laplace <- function(x, mu = 0, lambda = 1){
  return(1/(2*lambda)*exp(-abs(x-mu)/lambda))
}

rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (laplace(y, 0, 1) / laplace(x[i-1], 0, 1)))
      x[i] <- y 
    else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

#par(mfrow = c(2,2))
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
#number of candidate points rejected
accept <- 1-c(rw1$k, rw2$k, rw3$k, rw4$k)/N
print(accept)

plot(rw1$x, type="l", main = paste("accept rate = ",accept[1]), xlab = paste("sigma = ",sigma[1]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

plot(rw2$x, type="l", main = paste("accept rate = ",accept[2]), xlab = paste("sigma = ",sigma[2]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

plot(rw3$x, type="l", main = paste("accept rate = ",accept[3]), xlab = paste("sigma = ",sigma[3]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

plot(rw4$x, type="l", main = paste("accept rate = ",accept[4]), xlab = paste("sigma = ",sigma[4]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

## -----------------------------------------------------------------------------
x <- 1000
log(exp(x)) == exp(log(x))
isTRUE(all.equal(log(exp(x)), exp(log(x))))

## -----------------------------------------------------------------------------
k <- c(4:25, 100, 500, 1000)
r1 <- rep(0,length(k))
for(i in 1:length(k)){
  ki <- k[i]
  ck <- function(a, ki){
    sqrt(a^2*ki/(ki+1-a^2))
  }

  s <- function(a, ki){
    pt(ck(a, ki), df=ki) - pt(ck(a, ki-1), df=ki-1)
  }
  
  f <- function(a){
    return(s(a,ki))
  }
  r1[i] <- uniroot(f, lower = 1e-6, upper = sqrt(ki)-1e-6)$root
}
r1

## -----------------------------------------------------------------------------
k <- c(4:25)
r2 <- rep(0,length(k))
for(i in 1:length(k)){
  ki <- k[i]
  
  h <- function(a, ki){
    ck <- function(a, ki){
       sqrt(a^2*ki/(ki+1-a^2))
    } 
    
    g <- function(u){
      (1+u^2/ki)^(-(ki+1)/2)
    }
    
    # for simplily, get the log of the function
    r<-lgamma((ki+1)/2)-log(sqrt(pi*ki))-lgamma(ki/2)+log(integrate(g, lower = 0, upper=ck(a,ki), rel.tol=.Machine$double.eps^0.25)$value)
    return(r)
  }
  
  f <- function(a){
    h(a, ki) - h(a, ki-1)
  }
  r2[i] <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-0.15)$root
}
r2

## -----------------------------------------------------------------------------
ki <- 100
r3 <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-3)$root
r3

## -----------------------------------------------------------------------------
ki <- 500
r4 <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-20)$root
r4

## -----------------------------------------------------------------------------
ki <- 1000
r5 <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-26)$root
r5

## ----echo=FALSE---------------------------------------------------------------
    dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB','Sum'),
                 Frequency=c('p^2','q^2','r^2','2pr','2qr','2pq',1),
                 Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
    knitr::kable(dat,format='latex')

## -----------------------------------------------------------------------------
# define the function to compute n!
pro <- function(n){
  if(n<0) 
    print("warning: your input is wrong!")
  else if(n==0) 
    print(1)
  else {
    s <- length(n)
    s[1] <- 1
    for(i in 2:n){
      s[i] <- i*s[i-1]
    }
    return(s[n])
  }
}

## -----------------------------------------------------------------------------
# Use EM algorithm to solve MLE of p and q 
blood <- function(P, n.obs){
  n <- sum(n.obs)
  na <- n.obs[1]
  nb <- n.obs[2]
  nab <- n.obs[3]
  no <- n.obs[4]
  noo <- no
  
  ## E step
#  cat(P, "\n")  # P is the vector of original probability c(p,q,r)
  p <- q <- r <- rep(0, m)
  lnf <- rep(0, m-1)
  p[1] <- P[1]
  q[1] <- P[2]
  r[1] <- 1 - P[1] - P[2]
  for(i in 2:m){
    # instead old value by new value
    p.old <- p[i-1]  
    q.old <- q[i-1]
    r.old <- r[i-1]
    
    ## M step
    den1 <- p.old^2 + 2 * p.old * r.old
    naa <- na * p.old^2 / den1
    nao <- 2 * na * p.old * r.old / den1
    den2 <- q.old^2 + 2 * q.old * r.old
    nbb <- nb * q.old^2 / den2
    nbo <- 2 * nb * q.old * r.old / den2
    p[i] <- (2*naa + nao + nab)/(2*n)
    q[i] <- (2*nbb + nbo + nab)/(2*n)
    r[i] <- (2*noo + nao + nbo)/(2*n)
    lnf[i-1] <- naa*log(p[i]^2)+nao*log(2*p[i]*r[i])+nbb*log(q[i]^2)+nbo*log(2*q[i]*r[i])+nab*log(p[i]*q[i])+noo*log(r[i]^2)+log(pro(n)/(pro(naa)*pro(nao)*pro(nbb)*pro(nbo)*pro(nab)*pro(noo)))
  }
  
  return(list(P = c(p[m], q[m], r[m]), lh = lnf[1:50]))  # return the first 50 values
}
n.obs <- c(28, 24, 70, 41)
m <- 1000 #repicate times
P <- c(1/3, 1/3)  # set the original value

c <- blood(P, n.obs)
c$P
#c$lh
plot(c$lh, xlab="n", ylab="likelihood", type="b", cex=.5, col=2)

## -----------------------------------------------------------------------------
formulas <- list(    
  mpg ~ disp,     
  mpg ~ I(1 / disp),    
  mpg ~ disp + wt,    
  mpg ~ I(1 / disp) + wt    
) 

# use for loops
r1 <- list()
for (i in 1:4){
  r1[[i]]  <- lm(formula = formulas[[i]], data = mtcars)
}
r1

# use lapply()
r2 <- lapply(formulas, function(formulas){
  lm(formula = formulas, data = mtcars)
})
r2

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {    
  rows <- sample(1: nrow(mtcars), rep = TRUE)     
  mtcars[rows, ]     
}) 

# use for loops
r3 <- list()
for (i in 1:10){
r3[[i]] <- lm(formula = mpg ~ disp, data = bootstraps[[i]])
}
r3

# use lapply
r4 <- lapply(bootstraps, function(x) lm(formula = mpg ~ disp, data = x))
r4

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared 
R1 <- lapply(r1, rsq)
R1
R2 <- lapply(r2, rsq)
R2
R3 <- lapply(r3, rsq)
R3
R4 <- lapply(r4, rsq)
R4

## -----------------------------------------------------------------------------
trials <- replicate(100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE) 
# use sapply
p1 <- sapply(trials, function(x) x$p.value)
p1

# use for loops
p2 <- numeric(100)
for (i in 1:100){
p2[i] <- trials[[i]]$p.value
}

## -----------------------------------------------------------------------------
library(parallel)
f <- function(x) c(sum(x), mean(x), sd(x))
x <- list(1:10, c(2,3,4), c(TRUE,FALSE,TRUE,TRUE))

# use mcsapply
mcsapply <- function(x, f){
  cl <- makeCluster(4)
  result <- parLapply(cl, x, f)
  stopCluster(cl) 
  return(unlist(result))
}
# print the system.time
system.time(mcsapply(x, f))  

# use mcvapply
mcvapply <- function(x, f){
  cl <- makeCluster(4)
  result <- vapply(x,f,c("a" = 0, "b" = 0, "c" = 0))
  stopCluster(cl) 
  return(result)
}
# print the system.time
system.time(mcvapply(x, f))

## -----------------------------------------------------------------------------
for(j in 1:length(n)){
 for(i in 1:m){
    x<-matrix(rnorm(n[j]*2),nrow = n[j])
    e<-matrix(rnorm(n[j]*2),nrow = n[j])
    y <- x/4+e
    z <- as.data.frame(cbind(x,y))
    p.dc[i] <- dc.test(z, c(2,2))$p.value
    p.ball[i] <- bcov.test(x, y, R=R, seed = i*12345)$p.value
  }
  power.dc[j] <- mean(p.dc < alpha)
  power.ball[j] <- mean(p.ball < alpha)
}

# plot the power
plot(n, power.dc, type = "o", xlab = "n", ylab = "power", main = "Model 1", col = 1)
lines(n, power.ball, type = "o", col = 2)
legend('bottomright', legend = c('dc', 'ball'), col = 1:2, lty = 1)

## -----------------------------------------------------------------------------
for(j in 1:length(n)){
  for(i in 1:m){
    x<-matrix(rnorm(n[j]*2),nrow = n[j])
    e<-matrix(rnorm(n[j]*2),nrow = n[j])
    y <- x/4*e
    z <- as.data.frame(cbind(x,y))
    p.dc[i] <- dc.test(z, c(2,2))$p.value
    p.ball[i] <- bcov.test(x, y, R=R, seed = i*12345)$p.value
  }
  power.dc[j] <- mean(p.dc < alpha)
  power.ball[j] <- mean(p.ball < alpha)
}

# plot the power
plot(n, power.dc, type = "o", xlab = "n", ylab = "power", main = "Model 1", col = 1)
lines(n, power.ball, type = "o", col = 2)
legend('bottomright', legend = c('dc', 'ball'), col = 1:2, lty = 1)

## -----------------------------------------------------------------------------
laplace <- function(x, mu = 0, lambda = 1){
  return(1/(2*lambda)*exp(-abs(x-mu)/lambda))
}

rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (laplace(y, 0, 1) / laplace(x[i-1], 0, 1)))
      x[i] <- y 
    else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

#par(mfrow = c(2,2))
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
#number of candidate points rejected
accept <- 1-c(rw1$k, rw2$k, rw3$k, rw4$k)/N
print(accept)

plot(rw1$x, type="l", main = paste("accept rate = ",accept[1]), xlab = paste("sigma = ",sigma[1]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

plot(rw2$x, type="l", main = paste("accept rate = ",accept[2]), xlab = paste("sigma = ",sigma[2]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

plot(rw3$x, type="l", main = paste("accept rate = ",accept[3]), xlab = paste("sigma = ",sigma[3]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

plot(rw4$x, type="l", main = paste("accept rate = ",accept[4]), xlab = paste("sigma = ",sigma[4]), ylab = "x")
abline(h = c(log(0.05/2), -log(0.05/2)), col = "red")

## -----------------------------------------------------------------------------
x <- 1000
log(exp(x)) == exp(log(x))
isTRUE(all.equal(log(exp(x)), exp(log(x))))

## -----------------------------------------------------------------------------
k <- c(4:25, 100, 500, 1000)
r1 <- rep(0,length(k))
for(i in 1:length(k)){
  ki <- k[i]
  ck <- function(a, ki){
    sqrt(a^2*ki/(ki+1-a^2))
  }

  s <- function(a, ki){
    pt(ck(a, ki), df=ki) - pt(ck(a, ki-1), df=ki-1)
  }
  
  f <- function(a){
    return(s(a,ki))
  }
  r1[i] <- uniroot(f, lower = 1e-6, upper = sqrt(ki)-1e-6)$root
}
r1

## -----------------------------------------------------------------------------
k <- c(4:25)
r2 <- rep(0,length(k))
for(i in 1:length(k)){
  ki <- k[i]
  
  h <- function(a, ki){
    ck <- function(a, ki){
       sqrt(a^2*ki/(ki+1-a^2))
    } 
    
    g <- function(u){
      (1+u^2/ki)^(-(ki+1)/2)
    }
    
    # for simplily, get the log of the function
    r<-lgamma((ki+1)/2)-log(sqrt(pi*ki))-lgamma(ki/2)+log(integrate(g, lower = 0, upper=ck(a,ki), rel.tol=.Machine$double.eps^0.25)$value)
    return(r)
  }
  
  f <- function(a){
    h(a, ki) - h(a, ki-1)
  }
  r2[i] <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-0.15)$root
}
r2

## -----------------------------------------------------------------------------
ki <- 100
r3 <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-3)$root
r3

## -----------------------------------------------------------------------------
ki <- 500
r4 <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-20)$root
r4

## -----------------------------------------------------------------------------
ki <- 1000
r5 <- uniroot(f, lower = 1e-7, upper = sqrt(ki)-26)$root
r5

## ----echo=FALSE---------------------------------------------------------------
    dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB','Sum'),
                 Frequency=c('p^2','q^2','r^2','2pr','2qr','2pq',1),
                 Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
    knitr::kable(dat,format='latex')

## -----------------------------------------------------------------------------
# define the function to compute n!
pro <- function(n){
  if(n<0) 
    print("warning: your input is wrong!")
  else if(n==0) 
    print(1)
  else {
    s <- length(n)
    s[1] <- 1
    for(i in 2:n){
      s[i] <- i*s[i-1]
    }
    return(s[n])
  }
}

## -----------------------------------------------------------------------------
# Use EM algorithm to solve MLE of p and q 
blood <- function(P, n.obs){
  n <- sum(n.obs)
  na <- n.obs[1]
  nb <- n.obs[2]
  nab <- n.obs[3]
  no <- n.obs[4]
  noo <- no
  
  ## E step
#  cat(P, "\n")  # P is the vector of original probability c(p,q,r)
  p <- q <- r <- rep(0, m)
  lnf <- rep(0, m-1)
  p[1] <- P[1]
  q[1] <- P[2]
  r[1] <- 1 - P[1] - P[2]
  for(i in 2:m){
    # instead old value by new value
    p.old <- p[i-1]  
    q.old <- q[i-1]
    r.old <- r[i-1]
    
    ## M step
    den1 <- p.old^2 + 2 * p.old * r.old
    naa <- na * p.old^2 / den1
    nao <- 2 * na * p.old * r.old / den1
    den2 <- q.old^2 + 2 * q.old * r.old
    nbb <- nb * q.old^2 / den2
    nbo <- 2 * nb * q.old * r.old / den2
    p[i] <- (2*naa + nao + nab)/(2*n)
    q[i] <- (2*nbb + nbo + nab)/(2*n)
    r[i] <- (2*noo + nao + nbo)/(2*n)
    lnf[i-1] <- naa*log(p[i]^2)+nao*log(2*p[i]*r[i])+nbb*log(q[i]^2)+nbo*log(2*q[i]*r[i])+nab*log(p[i]*q[i])+noo*log(r[i]^2)+log(pro(n)/(pro(naa)*pro(nao)*pro(nbb)*pro(nbo)*pro(nab)*pro(noo)))
  }
  
  return(list(P = c(p[m], q[m], r[m]), lh = lnf[1:50]))  # return the first 50 values
}
n.obs <- c(28, 24, 70, 41)
m <- 1000 #repicate times
P <- c(1/3, 1/3)  # set the original value

c <- blood(P, n.obs)
c$P
#c$lh
plot(c$lh, xlab="n", ylab="likelihood", type="b", cex=.5, col=2)

## -----------------------------------------------------------------------------
formulas <- list(    
  mpg ~ disp,     
  mpg ~ I(1 / disp),    
  mpg ~ disp + wt,    
  mpg ~ I(1 / disp) + wt    
) 

# use for loops
r1 <- list()
for (i in 1:4){
  r1[[i]]  <- lm(formula = formulas[[i]], data = mtcars)
}
r1

# use lapply()
r2 <- lapply(formulas, function(formulas){
  lm(formula = formulas, data = mtcars)
})
r2

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {    
  rows <- sample(1: nrow(mtcars), rep = TRUE)     
  mtcars[rows, ]     
}) 

# use for loops
r3 <- list()
for (i in 1:10){
r3[[i]] <- lm(formula = mpg ~ disp, data = bootstraps[[i]])
}
r3

# use lapply
r4 <- lapply(bootstraps, function(x) lm(formula = mpg ~ disp, data = x))
r4

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared 
R1 <- lapply(r1, rsq)
R1
R2 <- lapply(r2, rsq)
R2
R3 <- lapply(r3, rsq)
R3
R4 <- lapply(r4, rsq)
R4

## -----------------------------------------------------------------------------
trials <- replicate(100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE) 
# use sapply
p1 <- sapply(trials, function(x) x$p.value)
p1

# use for loops
p2 <- numeric(100)
for (i in 1:100){
p2[i] <- trials[[i]]$p.value
}

## -----------------------------------------------------------------------------
library(parallel)
f <- function(x) c(sum(x), mean(x), sd(x))
x <- list(1:10, c(2,3,4), c(TRUE,FALSE,TRUE,TRUE))

# use mcsapply
mcsapply <- function(x, f){
  cl <- makeCluster(4)
  result <- parLapply(cl, x, f)
  stopCluster(cl) 
  return(unlist(result))
}
# print the system.time
system.time(mcsapply(x, f))  

# use mcvapply
mcvapply <- function(x, f){
  cl <- makeCluster(4)
  result <- vapply(x,f,c("a" = 0, "b" = 0, "c" = 0))
  stopCluster(cl) 
  return(result)
}
# print the system.time
system.time(mcvapply(x, f))

## ----eval=TRUE----------------------------------------------------------------
# R function

laplace <- function(x, mu = 0, lambda = 1){
  return(1/(2*lambda)*exp(-abs(x-mu)/lambda))
}

rw <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (laplace(y, 0, 1) / laplace(x[i-1], 0, 1)))
      x[i] <- y 
    else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(x)
}


# Cpp function
library(microbenchmark)
library(SC19037)
# define variables
N <- 10000
sigma <- c(0.5, 2, 8, 16)
x0 <- 5
#par(mfrow = c(2,2))
for (i in 1:4){
  # qqplot:
  r <- rw(sigma[i],x0,N)
  R <- as.vector(r)
  C <- rwC(sigma[i],x0,N)
  qqplot(R,C, main = paste("sigma = ",sigma[i]))
  points(x=seq(-8,6), y=seq(-8,6), type = "l", col = "red")
}

for (i in 1:4){
  # compare two function:
  ts <- microbenchmark(R=rw(sigma[i],x0,N), Cpp=rwC(sigma[i],x0,N))
  print(summary(ts)[,c(1,3,5,6)])
}


