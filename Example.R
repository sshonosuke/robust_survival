###------------------------------------------------------###
###           Code for example of applying               ###
###      Robust Bayesian inference Survival Models       ###
###------------------------------------------------------###

## preparation 
rm(list=ls())
set.seed(123)

source("RSV-function.R")
library(ggamma)    # for generation of generalized gamma random variable 

# regression function 
fun_reg <- function(x, al, be, ga){  
  p <- ga
  d <- ga * al
  rate <- al * exp( - sum(x * be))
  a <- (1/rate)^(1/p)
  reg <- a * gamma((d + 1)/p) / gamma(d/p)
  return(reg)
}


## settings 
om1 <- 0.1    # proportion of too long survival time 
om2 <- 0.1  # proportion of too short survival time

n <- 200   # sample size 
p <- 3   # number of covariate
out <- 100    # magnitude of outlier 

# true parameters of data generating distribution
alpha <- 10
gam <- 1      # gamma distribution if gam=1

cens <- T    # censored (T) or not censored (F) 
mc <- 4000     # length of MCMC
bn <- 2000     # length of burn-in 
th <- 5        # number of thinning
thin <- seq(from=1, to=mc-bn, by=th)



## data generation 
x1 <- runif(n, 0, 2)
x2 <- runif(n, -2, 2)
X <- cbind(1, x1, x2)
Beta <- c(0.5, 2, -0.5)
tY <- rggamma(n, a=(alpha*exp(-as.vector(X%*%Beta)))^(-1/gam), b=gam, k=alpha) 
hist(tY)
mean(tY)


## true regression point
reg1 <- fun_reg(x=c(1, 0.5, -1), al=alpha, be=Beta, ga=gam)
reg2 <- fun_reg(x=apply(X, 2, mean), al=alpha, be=Beta, ga=gam)
reg3 <- fun_reg(x=c(1, 1.5, 1), al=alpha, be=Beta, ga=gam)

## outlier 
ch1 <- rbinom(n, 1, 2*om1)*ifelse(x2>0.5, 1, 0)
ch2 <- rbinom(n, 1, om2)*ifelse(ch1==0, 1, 0)
YY <- tY + out*ch1
YY[ch2==1] <- 1

## censoring 
max.C <- max(max(tY), 50) 
C <- runif(n, 10, max.C)    # censoring time
delta <- ifelse(YY<C, 1, 0)   # censoring indicator
delta[ch1==1] <- 1
mean(1-delta)   # censoring rate
Time <- YY
Time[delta==0] <- C[delta==0]    # observed survival time

left0 <- right0 <- rep(FALSE, n)
right0[delta==0] <- TRUE


## estimation 
Fit <- list()

# RGG (robust generalized gamma model)
Fit[[1]] <- GG_robust(y=Time, X=X, left=left0, right=right0, mc=mc, burn=bn,
                      est_la=T, est_al=T, est_ga=T)

# GG (non-robust generalized gamma model)
Fit[[2]] <- GG_robust(y=Time, X=X, left=left0, right=right0, mc=mc, burn=bn,
                      est_la=F, est_al=T, est_ga=T)

# RGA (robust gamma model)
Fit[[3]] <- GG_robust(y=Time, X=X, left=left0, right=right0, mc=mc, burn=bn,
                      est_la=T, est_al=T, est_ga=F)

# GA (non-robust gamma model)
Fit[[4]] <- GG_robust(y=Time, X=X, left=left0, right=right0, mc=mc, burn=bn,
                      est_la=F, est_al=T, est_ga=F)

# RWB (robust Weibul model)
Fit[[5]] <- GG_robust(y=Time, X=X, left=left0, right=right0, mc=mc, burn=bn,
                      est_la=T, est_al=F, est_ga=T)

# WB (non-robust Weibul model)
Fit[[6]] <- GG_robust(y=Time, X=X, left=left0, right=right0, mc=mc, burn=bn,
                      est_la=F, est_al=F, est_ga=T)



## Result 
meth <- c("RGG", "GG", "RGA", "GA", "RWB", "WB")
L <- length(meth)

# posterior of regression line
fn1 <- function(para){ fun_reg(x=c(1, 0.5, -1), al=para[1], be=para[2:4], ga=para[5]) }
fn2 <- function(para){ fun_reg(x=apply(X, 2, mean), al=para[1], be=para[2:4], ga=para[5]) }
fn3 <- function(para){ fun_reg(x=c(1, 1.5, 1), al=para[1], be=para[2:4], ga=para[5]) }

fn1_pos <- fn2_pos <- fn3_pos <- matrix(NA, mc-bn, L)
dimnames(fn1_pos)[[2]] <- dimnames(fn2_pos)[[2]] <- dimnames(fn3_pos)[[2]] <- meth
for(k in 1:6){
  pos <- cbind(Fit[[k]]$al, (-1)*Fit[[k]]$be, Fit[[k]]$ga)
  fn1_pos[,k] <- apply(pos, 1, fn1)
  fn2_pos[,k] <- apply(pos, 1, fn2)
  fn3_pos[,k] <- apply(pos, 1, fn3)
}

# posterior mean
est1 <- apply(fn1_pos, 2, mean)
est2 <- apply(fn2_pos, 2, mean)
est3 <- apply(fn3_pos, 2, mean)

# credible interval 
CI1 <- apply(fn1_pos, 2, quantile, prob=c(0.025, 0.975))
CI2 <- apply(fn2_pos, 2, quantile, prob=c(0.025, 0.975))
CI3 <- apply(fn3_pos, 2, quantile, prob=c(0.025, 0.975))

# evaluation 
est1   # posterior mean  
CI1    # 95% credible interval 
reg1   # true regression point 

est2   # posterior mean  
CI2    # 95% credible interval 
reg2   # true regression point 

est3   # posterior mean  
CI3    # 95% credible interval 
reg3   # true regression point 

