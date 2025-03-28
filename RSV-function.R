###------------------------------------------------###
###      Code for Robust Bayesian inference        ###
###        for Censored Survival Models            ###
###------------------------------------------------###

## packages
library(expint)   # incomplete gamma function
library(MASS)
library(truncdist)
library(progress)   # for progress bar


## likelihood functions
log_lik_GG <- function(y, X, alpha, beta, gamma, lambda){
  log( (gamma * alpha^alpha / gamma(alpha)) * (exp(X%*%beta)/lambda)^alpha * y^(alpha*gamma-1) * exp(-alpha * (exp(X%*%beta)/lambda) * y^gamma) )
}

log_surv_fun_GG <- function(y, X, alpha, beta, gamma, lambda){
  log( gammainc(alpha*(exp(X%*%beta)/lambda)*y^gamma, alpha) / gamma(alpha) )
}


## customized function for generating gamma random variables
rgamma01 <- function(n, shape, rate){
  sr <- cbind(shape, rate)
  qgamma(shape=shape, rate=rate, p=runif(n)*pgamma(shape=shape, rate=rate, q=1))
}

rgamma1Inf <- function(n, shape, rate){
  prob <- pgamma(shape=shape, rate=rate, q=1)
  qgamma(shape=shape, rate=rate, p=runif(n)*(1-prob)+prob)
}


### Main function 
## INPUT
# y: censored survival time (vector)
# X: covariate matrix
# left: vector of left censoring indicator ("T" means censored)
# right: vector of right censoring indicator ("T" means censored)
# mc: length of MCMC
# burn: length of burn-in period
# est_la: if "T", the robust model is applied 
# est_al: if "T", shape parameter of generalized gamma is estimated 
#        (if "F", the model reduces to Weibul model)
# est_ga: if "T", shape parameter of generalized gamma is estimated 
#        (if "F", the model reduces to gamma model)

## OUTPUT
# Posterior samples of the model parameters 

GG_robust <- function(y, X, left=NULL, right=NULL, mc=5000, burn=2000, 
                      est_la=T, est_al=T, est_ga=T){
  ## preparation 
  if(is.null(left)){  left <- rep(F, dim(X)[1])  }
  if(is.null(right)){  right <- rep(F, dim(X)[1])  }
  p <- dim(X)[2]
  n <- dim(X)[1]

  ## tuning parameters
  a_al <- b_al <- a_ga <- b_ga <- 1/10
  c_al <- c_ga <- 0
  c_EH <- 1
  mu <- rep(0, 3)
  psi <- rep(1/10^3, 3)
  a_s <- b_s <- 1
  ep <- ep_bet <- 10^(-8)
  M <- 10
  M_bet <- 50
  G <- 100
  Psi <- diag(psi, nrow=p)
  
  C_left <- y[left]
  C_right <- y[right]
  
  ## initial values
  la <- rep(1, n)
  etat <- rep(1, n)
  z <- rep(0, n)
  u <- v <- w <- rep(1, n)
  alt <- 1
  ga <- 1
  bet <- rep(0, p)
  s <- rep(1/2, n)
  
  ## objects for storing posterior samples
  Al <- rep(NA, mc)
  Be <- matrix(NA, mc, p)
  Ga <- rep(NA, mc)
  Alt <- rep(NA, mc)    # for transformed parameters
  Bet <- matrix(NA, mc, p)   # for transformed parameters
  Gat <- rep(NA, mc)    # for transformed parameters
  log_likelihood <- rep(NA, mc)
  ss <- rep(NA, mc)
  zz <- matrix(NA, mc, n)
  La <- matrix(NA, mc, n)
  Etat <- matrix(NA, mc, n)
  
  ac <- c(al=0, be=0, ga=0)   # acceptance ratio 
  iter_bet <- rep(NA, mc)
  
  ## MCMC
  pb <- progress_bar$new(total=mc)   # progress bar 
  if(est_al | (!est_ga)){
    # MCMC for fixed alpha
    for (item in 1:mc){
      # s 
      s <- rbeta(1, shape1=sum(z) + a_s, shape2=n - sum(z) + b_s)
      
      # eta 
      u <- rgamma(n, shape=1 + c_EH, rate=1 + log(1 + log(1 + etat^(1 - z + z * ga))))
      v <- rgamma(n, shape=1 + u, rate=1 + log(1 + etat^(1 - z + z * ga)))
      w <- rgamma(n, shape=1 + v, rate=1 + 1 / etat^(1 - z + z * ga))
      shape <- v + alt / ga^2
      rate <- w + (alt / ga^2) * (exp(c(X %*% bet)) * y)^ga
      nonout <- (z == 0)
      if (any(!nonout)){
        etat[!nonout] <- (1 / rgamma(n=sum(!nonout), shape=shape[!nonout], rate=rate[!nonout]))^(1 / ga)
      }
      rand <- rbeta(n=sum(nonout), shape1=1, shape2=c_EH)
      rand <- rand / (1 - rand)
      rand <- exp(rand) - 1
      rand <- exp(rand) - 1
      etat[nonout] <- rand    
      etat[etat > 10^8] <- 10^8
      
      # z 
      log_den_1 <- log(s) + log(ga) + (ga - 1) * log(etat) + log(c_EH) - log(1 + etat^ga) - log(1 + log(1 + etat^ga)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat^ga))) - (alt / ga^2) * log(etat^ga) - (alt / ga^2) * (exp(c(X %*% bet)) * y)^ga / etat^ga
      log_den_0 <- log(1 - s) + log(c_EH) - log(1 + etat) - log(1 + log(1 + etat)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat))) - (alt / ga^2) * (exp(c(X %*% bet)) * y)^ga
      z <- ifelse(test=(runif(n) <= 1 / (1 + exp(log_den_0 - log_den_1))), 
                  yes=1, no=0)
      if(!est_la){ z <- rep(0, n) }
      
      # la 
      la <- etat^(z*ga)
      
      # al 
      if(est_al){
        S <- c(X %*% bet) * ga + log(y^ga / la)
        S <- S - exp(S)
        S <- sum(S)
        A <- B <- 1
        for (j in 1:M) {
          alt_j <- A / B
          A <- c_al + 2 * b_al / (alt_j / ga^2) - n * (alt_j / ga^2) * (1 - (alt_j / ga^2) * trigamma(alt_j / ga^2))
          B <- (A - c_al) / alt_j + a_al / ga^2 - b_al / (alt_j / ga)^2 - (n / ga^2) * (log(alt_j / ga^2) + 1 - digamma(alt_j / ga^2)) - S / ga^2
          if (abs(alt_j / (A / B) - 1) < ep) {
            break
          }
        }
        ABj <- c(A=A, B=B, j=j)
        alt_proposal <- rgamma(1, shape=ABj[1], rate=ABj[2])
        log_MHratio_gamma <- ((ABj[1] - 1) * log(alt_proposal) - ABj[2] * alt_proposal) - ((ABj[1] - 1) * log(alt) - ABj[2] * alt)
        log_MHratio_target <- (c_al-1)*log(alt_proposal) - a_al*(alt_proposal/ga^2) - b_al/(alt_proposal/ga^2) + n*((alt_proposal/ga^2) * log(alt_proposal/ga^2) - lgamma(alt_proposal/ga^2)) + (alt_proposal/ga^2)*S
        log_MHratio_target <- log_MHratio_target - ((c_al-1)*log(alt) - a_al*(alt/ga^2) - b_al/(alt/ga^2) + n*((alt/ga^2) * log(alt/ga^2) - lgamma(alt/ga^2)) + (alt/ga^2)*S)
        if (log(runif(1)) <= log_MHratio_target - log_MHratio_gamma){
          alt <- alt_proposal
          ac["al"] <- ac["al"] + 1
        }else{
          alt <- alt
        }
      }else{
        alt <- alt
      }
      
      # ga 
      if(est_ga){
        gatgat <- seq(0, 1, by=1 / G)
        gatgat[1] <- 1 / (2 * G)
        gatgat[length(gatgat)] <- 1 - 1 / (2 * G)
        gaga <- gatgat / (1 - gatgat)
        log_den <- (-2) * log(1 - gatgat) + (p - 2) * log(gaga) + (c_al - 1) * log(alt / gaga^2) - a_al * (alt / gaga^2) - b_al / (alt / gaga^2) - (1 / 2) * colSums(psi * (outer(bet, gaga) - mu)^2) + (c_ga - 1) * log(gaga) - a_ga * gaga - b_ga / gaga
        inner <- c(X %*% bet)
        SS <- outer(inner, gaga) + outer(log(y) - z * log(etat), gaga)
        SS <- SS - exp(SS)
        SS <- colSums(SS)
        log_den <- log_den + n * log(gaga) + n * ((alt / gaga^2) * log(alt / gaga^2) - lgamma(alt / gaga^2)) + (alt / gaga^2) * SS
        etat_out <- etat[z == 1]
        log_den <- log_den + sum(z) * log(gaga) + any(z == 1) * colSums(
          outer(log(etat_out), gaga - 1) - log(1 + exp(outer(log(etat_out), gaga))) - log(1 + log(1 + exp(outer(log(etat_out), gaga)))) - (1 + c_EH) * log(1 + log(1 + log(1 + exp(outer(log(etat_out), gaga)))))
        )
        middle <- log_den[-(c(1, length(gatgat)))]
        middle_max <- pmax(middle[-length(middle)], middle[-1])
        middle_min <- pmin(middle[-length(middle)], middle[-1])
        log_area_middle <- log(1 / 2) + middle_max + log(1 + exp(middle_min - middle_max))
        log_area_middle[middle_min == middle_max] <- (log(1 / 2) + middle_max + log(1 + exp(0)))[middle_min == middle_max]
        log_area <- c(log_den[1], log_area_middle, log_den[length(gatgat)])
        prob <- exp(log_area - max(log_area))
        prob <- prob / sum(prob)
        g <- sample(x=1:G, size=1, replace=T, prob=prob)
        if(g==1){
          gat_proposal <- runif(1)/G
        }else if(g==G){
          gat_proposal <- 1 - runif(1) / G
        }else{
          a_g <- (1/G) * (g-1)
          b_g <- (1/G) * g
          log_den_f_1_t_G <- log_den[-1]
          log_c_g <- log_den_f_1_t_G[g - 1]
          log_d_g <- log_den_f_1_t_G[g]
          ratio <- 1 / (exp(log_d_g - log_c_g) - 1)
          sign_d_gminusc_g <- sign(log_d_g - log_c_g)
          difference <- b_g - a_g
          if(abs(exp(log_d_g-log_c_g) - 1) < 10^(-8)){
            gat_proposal <- runif(n=1, min=a_g, max=b_g)
          }else{
            gat_proposal <- a_g - difference * ratio + sign_d_gminusc_g * sqrt((ratio * difference - a_g)^2 + a_g * (2 * difference * ratio - a_g) + runif(1) * difference^2 * (2 * ratio + 1))
          }
        }
        ga_proposal <- gat_proposal / (1 - gat_proposal)
        S_proposal <- c(X %*% bet) * ga_proposal + ga_proposal * log(y / etat^z)
        S_proposal <- S_proposal - exp(S_proposal)
        S_proposal <- sum(S_proposal)
        log_MHratio_target_proposal <- (-2) * log(1 - gat_proposal) + (p - 2) * log(ga_proposal) + (c_al - 1) * log(alt / ga_proposal^2) - a_al * (alt / ga_proposal^2) - b_al / (alt / ga_proposal^2) - (1 / 2) * c(t(ga_proposal * bet - mu) %*% Psi %*% (ga_proposal * bet - mu)) + (c_ga - 1) * log(ga_proposal) - a_ga * ga_proposal - b_ga / ga_proposal + n * log(ga_proposal) + n * ((alt / ga_proposal^2) * log(alt / ga_proposal^2) - lgamma(alt / ga_proposal^2)) + (alt / ga_proposal^2) * S_proposal
        log_MHratio_target_proposal <- log_MHratio_target_proposal + any(z == 1) * sum(
          log(ga_proposal) + (ga_proposal - 1) * log(etat_out) - log(1 + etat_out^ga_proposal) - log(1 + log(1 + etat_out^ga_proposal)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat_out^ga_proposal)))
        )
        gat <- ga / (1 + ga)
        S <- c(X %*% bet) * ga + ga * log(y / etat^z)
        S <- S - exp(S)
        S <- sum(S)
        log_MHratio_target_old <- (-2) * log(1 - gat) + (p - 2) * log(ga) + (c_al - 1) * log(alt / ga^2) - a_al * (alt / ga^2) - b_al / (alt / ga^2) - (1 / 2) * c(t(ga * bet - mu) %*% Psi %*% (ga * bet - mu)) + (c_ga - 1) * log(ga) - a_ga * ga - b_ga / ga + n * log(ga) + n * ((alt / ga^2) * log(alt / ga^2) - lgamma(alt / ga^2)) + (alt / ga^2) * S
        log_MHratio_target_old <- log_MHratio_target_old + any(z == 1) * sum(
          log(ga) + (ga - 1) * log(etat_out) - log(1 + etat_out^ga) - log(1 + log(1 + etat_out^ga)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat_out^ga)))
        )
        if(g==1){
          log_MHratio_linear_proposal <- log_den[1]
        }else if(g == G){
          log_MHratio_linear_proposal <- log_den[length(log_den)]
        }else{
          if(abs(exp(log_d_g - log_c_g) - 1) < 10^(-8)) {
            log_MHratio_linear_proposal <- log_c_g
          }else{
            log_MHratio_linear_proposal <- log_c_g + log(1 + (gat_proposal - a_g) / (difference * ratio))
          }
        }
        g_old <- gat%/%(1/G) + 1
        if(g_old==1){
          log_MHratio_linear_old <- log_den[1]
        }else if(g_old==G){
          log_MHratio_linear_old <- log_den[length(log_den)]
        }else{
          a_g_old <- (1 / G) * (g_old - 1)
          b_g_old <- (1 / G) * g_old
          log_den_f_1_t_G <- log_den[-1]
          log_c_g_old <- log_den_f_1_t_G[g_old - 1]
          log_d_g_old <- log_den_f_1_t_G[g_old]
          ratio_old <- 1 / (exp(log_d_g_old - log_c_g_old) - 1)
          difference_old <- b_g_old - a_g_old
          if(abs(exp(log_d_g_old - log_c_g_old) - 1) < 10^(-8)){
            log_MHratio_linear_old <- log_c_g_old
          }else{
            log_MHratio_linear_old <- log_c_g_old + log(1 + (gat - a_g_old) / (difference_old * ratio_old))
          }
        }
        if(log(runif(1)) <= (log_MHratio_target_proposal - log_MHratio_target_old) - (log_MHratio_linear_proposal - log_MHratio_linear_old)) {
          gat_new <- gat_proposal
          ac["ga"] <- ac["ga"] + 1
        }else{
          gat_new <- gat
        }
        ga <- gat_new / (1 - gat_new)
      }else{
        ga <- ga
      }
      ga <- min(100, ga)
      
      # la 
      la <- etat^(z * ga)
      
      # be 
      mu0 <- rep(0, p)
      Psi0 <- diag(p)
      for(j in 1:M_bet){
        bet_j <- mu0
        Psi0 <- ga^2 * Psi + ((alt / ga^2) * ga^2) * t(X) %*% (((y^ga / la) * exp(c(X %*% (ga * bet_j)))) * X)
        Psi0 <- (Psi0 + t(Psi0)) / 2
        mu0 <- ga * c(Psi %*% (ga * bet_j - mu)) - (alt / ga^2) * colSums(X) * ga + (alt / ga^2) * ga * c(t(X) %*% ((y^ga / la) * exp(c(X %*% (ga * bet_j)))))
        mu0 <- bet_j - c(solve(Psi0) %*% mu0)
        if(sum(abs(bet_j - mu0)) < ep_bet){  break  }
      }
      iter_bet[item] <- j
      bet_proposal <- mvrnorm(n=1, mu=mu0, Sigma=solve(Psi0))
      log_MHratio_target_proposal <- ((-1) / 2) * c(t(ga * bet_proposal - mu) %*% Psi %*% (ga * bet_proposal - mu)) + (alt / ga^2) * ga * sum(colSums(X) * bet_proposal) - (alt / ga^2) * sum((y^ga / la) * exp(c(X %*% (ga * bet_proposal))))
      log_MHratio_target_old <- ((-1) / 2) * c(t(ga * bet - mu) %*% Psi %*% (ga * bet - mu)) + (alt / ga^2) * ga * sum(colSums(X) * bet) - (alt / ga^2) * sum((y^ga / la) * exp(c(X %*% (ga * bet))))
      log_MHratio_q_proposal <- ((-1) / 2) * c(t(bet_proposal - mu0) %*% Psi0 %*% (bet_proposal - mu0))
      log_MHratio_q_old <- ((-1) / 2) * c(t(bet - mu0) %*% Psi0 %*% (bet - mu0))
      log_prob <- (log_MHratio_target_proposal - log_MHratio_target_old) - (log_MHratio_q_proposal - log_MHratio_q_old)
      if(log(runif(1))<log_prob){
        bet <- bet_proposal
        ac["be"] <- ac["be"] + 1
      }else{
        bet <- bet
      }
      
      # y (imputation of censored observations)
      if(any(right)){
        rate <- C_right^ga * (alt / ga^2) * exp(c(X[right,,drop=F] %*% (bet*ga))) / la[right]
        rand_right <- rgamma1Inf(sum(right), shape=alt/ga^2, rate=rate)
        rand_right <- rand_right * C_right^ga
        rand_right <- rand_right^(1 / ga)
        rand_right[rand_right==Inf] <- 10^8
        y[right] <- rand_right
      }
      if(any(left)){
        rate <- C_left^ga * (alt / ga^2) * exp(c(X[left,,drop=F] %*% (bet*ga))) / la[left]
        rand_left <- rgamma01(sum(left), shape=alt/ga^2, rate=rate)
        rand_left <- rand_left * C_left^ga
        rand_left <- rand_left^(1 / ga)
        y[left] <- rand_left
      }
      y[y>10^8] <- 10^8    # upper bound to avoid numerical problems
      
      # store posterior samples
      Al[item] <- alt / ga^2
      Be[item, ] <- bet * ga
      Ga[item] <- ga
      log_likelihood[item] <- sum( as.numeric(right) * log_lik_GG(y=y, X=X, alpha=alt / ga^2, beta=bet * ga, gamma=ga, lambda=la) + (1-as.numeric(right)) * log_surv_fun_GG(y=y, X=X, alpha=alt / ga^2, beta=bet * ga, gamma=ga, lambda=la) )
      ss[item] <- s
      Alt[item] <- alt
      Bet[item, ] <- bet
      Gat[item] <- ga / (1 + ga)
      zz[item, ] <- z
      La[item, ] <- la
      Etat[item, ] <- etat
      
      pb$tick()
    }
  }else{
    # MCMC for estimated alpha
    alt <- NA
    al <- 1
    for (item in 1:mc) {
      # mixing proportion
      s <- rbeta(1, shape1=sum(z)+a_s, shape2=n-sum(z)+b_s)
      
      # eta 
      u <- rgamma(n, shape=1 + c_EH, rate=1 + log(1 + log(1 + etat^(1 - z + z * ga))))
      v <- rgamma(n, shape=1 + u, rate=1 + log(1 + etat^(1 - z + z * ga)))
      w <- rgamma(n, shape=1 + v, rate=1 + 1 / etat^(1 - z + z * ga))
      
      shape <- v + al
      rate <- w + al * (exp(c(X %*% bet)) * y)^ga
      nonout <- (z==0)
      if (any(!nonout)) {
        etat[!nonout] <- (1 / rgamma(n=sum(!nonout), shape=shape[!nonout], rate=rate[!nonout]))^(1 / ga)
      }
      rand <- rbeta(n=sum(nonout), shape1=1, shape2=c_EH)
      rand <- rand / (1 - rand)
      rand <- exp(rand) - 1
      rand <- exp(rand) - 1
      etat[nonout] <- rand    
      etat[etat > 10^8] <- 10^8
    
      # z 
      log_den_1 <- log(s) + log(ga) + (ga - 1) * log(etat) + log(c_EH) - log(1 + etat^ga) - log(1 + log(1 + etat^ga)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat^ga))) - al * log(etat^ga) - al * (exp(c(X %*% bet)) * y)^ga / etat^ga
      log_den_0 <- log(1 - s) + log(c_EH) - log(1 + etat) - log(1 + log(1 + etat)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat))) - al * (exp(c(X %*% bet)) * y)^ga
      z <- ifelse(test=(runif(n) <= 1 / (1 + exp(log_den_0 - log_den_1))), 
                  yes=1, no=0)
      if(!est_la){ z <- rep(0, n) }
      
      # la 
      la <- etat^(z*ga)
      
      # ga 
      if(est_ga){
        gatgat <- seq(0, 1, by=1 / G)
        gatgat[1] <- 1 / (2 * G)
        gatgat[length(gatgat)] <- 1 - 1 / (2 * G)
        gaga <- gatgat / (1 - gatgat)
        log_den <- (-2) * log(1 - gatgat) + p * log(gaga) - (1 / 2) * colSums(psi * (outer(bet, gaga) - mu)^2) + (c_ga - 1) * log(gaga) - a_ga * gaga - b_ga / gaga
        inner <- c(X %*% bet)
        SS <- outer(inner, gaga) + outer(log(y) - z * log(etat), gaga)
        SS <- SS - exp(SS)
        SS <- colSums(SS)
        log_den <- log_den + n * log(gaga) + al * SS
        etat_out <- etat[z == 1]
        log_den <- log_den + sum(z) * log(gaga) + any(z == 1) * colSums(
          outer(log(etat_out), gaga - 1) - log(1 + exp(outer(log(etat_out), gaga))) - log(1 + log(1 + exp(outer(log(etat_out), gaga)))) - (1 + c_EH) * log(1 + log(1 + log(1 + exp(outer(log(etat_out), gaga)))))
        )
        middle <- log_den[-(c(1, length(gatgat)))]
        middle_max <- pmax(middle[-length(middle)], middle[-1])
        middle_min <- pmin(middle[-length(middle)], middle[-1])
        log_area_middle <- log(1 / 2) + middle_max + log(1 + exp(middle_min - middle_max))
        log_area_middle[middle_min == middle_max] <- (log(1 / 2) + middle_max + log(1 + exp(0)))[middle_min == middle_max]
        log_area <- c(log_den[1], log_area_middle, log_den[length(gatgat)])
        prob <- exp(log_area - max(log_area))
        prob <- prob / sum(prob)
        g <- sample(x=1:G, size=1, replace=T, prob=prob)
        if(g==1){
          gat_proposal <- runif(1) / G
        }else if(g==G){
          gat_proposal <- 1 - runif(1) / G
        }else{
          a_g <- (1 / G) * (g - 1)
          b_g <- (1 / G) * g
          log_den_f_1_t_G <- log_den[-1]
          log_c_g <- log_den_f_1_t_G[g - 1]
          log_d_g <- log_den_f_1_t_G[g]
          ratio <- 1 / (exp(log_d_g - log_c_g) - 1)
          sign_d_gminusc_g <- sign(log_d_g - log_c_g)
          difference <- b_g - a_g
          if(abs(exp(log_d_g - log_c_g) - 1) < 10^(-8)){
            gat_proposal <- runif(n=1, min=a_g, max=b_g)
          }else{
            gat_proposal <- a_g - difference * ratio + sign_d_gminusc_g * sqrt((ratio * difference - a_g)^2 + a_g * (2 * difference * ratio - a_g) + runif(1) * difference^2 * (2 * ratio + 1))
          }
        }
        ga_proposal <- gat_proposal / (1 - gat_proposal)
        S_proposal <- c(X %*% bet) * ga_proposal + ga_proposal * log(y / etat^z)
        S_proposal <- S_proposal - exp(S_proposal)
        S_proposal <- sum(S_proposal)
        log_MHratio_target_proposal <- (-2) * log(1 - gat_proposal) + p * log(ga_proposal) - (1 / 2) * c(t(ga_proposal * bet - mu) %*% Psi %*% (ga_proposal * bet - mu)) + (c_ga - 1) * log(ga_proposal) - a_ga * ga_proposal - b_ga / ga_proposal + n * log(ga_proposal) + al * S_proposal
        log_MHratio_target_proposal <- log_MHratio_target_proposal + any(z == 1) * sum(
          log(ga_proposal) + (ga_proposal - 1) * log(etat_out) - log(1 + etat_out^ga_proposal) - log(1 + log(1 + etat_out^ga_proposal)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat_out^ga_proposal)))
        )
        gat <- ga / (1 + ga)
        S <- c(X %*% bet) * ga + ga * log(y / etat^z)
        S <- S - exp(S)
        S <- sum(S)
        log_MHratio_target_old <- (-2) * log(1 - gat) + p * log(ga) - (1 / 2) * c(t(ga * bet - mu) %*% Psi %*% (ga * bet - mu)) + (c_ga - 1) * log(ga) - a_ga * ga - b_ga / ga + n * log(ga) + al * S
        log_MHratio_target_old <- log_MHratio_target_old + any(z == 1) * sum(
          log(ga) + (ga - 1) * log(etat_out) - log(1 + etat_out^ga) - log(1 + log(1 + etat_out^ga)) - (1 + c_EH) * log(1 + log(1 + log(1 + etat_out^ga)))
        )
        if(g==1){
          log_MHratio_linear_proposal <- log_den[1]
        }else if(g==G){
          log_MHratio_linear_proposal <- log_den[length(log_den)]
        }else{
          if(abs(exp(log_d_g - log_c_g) - 1) < 10^(-8)) {
            log_MHratio_linear_proposal <- log_c_g
          }else{
            log_MHratio_linear_proposal <- log_c_g + log(1 + (gat_proposal - a_g) / (difference * ratio))
          }
        }
        g_old <- gat %/% (1/G) + 1
        if(g_old == 1){
          log_MHratio_linear_old <- log_den[1]
        }else if(g_old == G){
          log_MHratio_linear_old <- log_den[length(log_den)]
        }else{
          a_g_old <- (1 / G) * (g_old - 1)
          b_g_old <- (1 / G) * g_old
          log_den_f_1_t_G <- log_den[-1]
          log_c_g_old <- log_den_f_1_t_G[g_old - 1]
          log_d_g_old <- log_den_f_1_t_G[g_old]
          ratio_old <- 1 / (exp(log_d_g_old - log_c_g_old) - 1)
          difference_old <- b_g_old - a_g_old
          if(abs(exp(log_d_g_old - log_c_g_old) - 1) < 10^(-8)){
            log_MHratio_linear_old <- log_c_g_old
          }else{
            log_MHratio_linear_old <- log_c_g_old + log(1 + (gat - a_g_old) / (difference_old * ratio_old))
          }
        }
        if(log(runif(1)) <= (log_MHratio_target_proposal - log_MHratio_target_old) - (log_MHratio_linear_proposal - log_MHratio_linear_old)) {
          gat_new <- gat_proposal
          ac["ga"] <- ac["ga"] + 1
        }else{
          gat_new <- gat
        }
        ga <- gat_new / (1 - gat_new)
      }
      ga <- min(100, ga)
    
      # lambda (local parameter for absorbing outier effects)
      la <- etat^(z*ga)
      
      # beta (Muller's method)
      mu0 <- rep(0, p)
      Psi0 <- diag(p)
      for(j in 1:M_bet){
        bet_j <- mu0
        Psi0 <- ga^2 * Psi + (al * ga^2) * t(X) %*% (((y^ga / la) * exp(c(X %*% (ga * bet_j)))) * X)
        Psi0 <- (Psi0 + t(Psi0)) / 2
        mu0 <- ga * c(Psi %*% (ga * bet_j - mu)) - al * colSums(X) * ga + al * ga * c(t(X) %*% ((y^ga / la) * exp(c(X %*% (ga * bet_j)))))
        mu0 <- bet_j - c(solve(Psi0) %*% mu0)
        if (sum(abs(bet_j - mu0)) < ep_bet) {
          break
        }
      }
      iter_bet[item] <- j
      bet_proposal <- mvrnorm(n=1, mu=mu0, Sigma=solve(Psi0))
      log_MHratio_target_proposal <- ((-1) / 2) * c(t(ga * bet_proposal - mu) %*% Psi %*% (ga * bet_proposal - mu)) + al * ga * sum(colSums(X) * bet_proposal) - al * sum((y^ga / la) * exp(c(X %*% (ga * bet_proposal))))
      log_MHratio_target_old <- ((-1) / 2) * c(t(ga * bet - mu) %*% Psi %*% (ga * bet - mu)) + al * ga * sum(colSums(X) * bet) - al * sum((y^ga / la) * exp(c(X %*% (ga * bet))))
      log_MHratio_q_proposal <- ((-1) / 2) * c(t(bet_proposal - mu0) %*% Psi0 %*% (bet_proposal - mu0))
      log_MHratio_q_old <- ((-1) / 2) * c(t(bet - mu0) %*% Psi0 %*% (bet - mu0))
      if(log(runif(1)) <= (log_MHratio_target_proposal - log_MHratio_target_old) - (log_MHratio_q_proposal - log_MHratio_q_old)) {
        bet <- bet_proposal
        ac["be"] <- ac["be"] + 1
      }else{
        bet <- bet
      }
      
      # imputation of censoring data 
      if(any(right)){
        rand_right <- rgamma1Inf(sum(right), shape=al, 
                                 rate=C_right^ga * al * exp(c(X[right, , drop=F] %*% (bet * ga))) / la[right])
        rand_right <- rand_right * C_right^ga
        rand_right <- rand_right^(1 / ga)
        rand_right[rand_right == Inf] <- 10^8
        y[right] <- rand_right
      }
      
      if(any(left)){
        rand_left <- rgamma01(sum(left), shape=al, 
                              rate=C_left^ga*al*exp(c(X[left,,drop=F] %*% (bet * ga))) / la[left])
        rand_left <- rand_left * C_left^ga
        rand_left <- rand_left^(1 / ga)
        y[left] <- rand_left
      }
      # upper bound to avoid numerical instability 
      y[y>10^8] <- 10^8    
      
      # store the current values
      Al[item] <- al
      Be[item, ] <- bet * ga
      Ga[item] <- ga
      log_likelihood[item] <- sum( as.numeric(right) * log_lik_GG(y=y, X=X, alpha=al, beta=bet * ga, gamma=ga, lambda=la) + (1-as.numeric(right)) * log_surv_fun_GG(y=y, X=X, alpha=alt / ga^2, beta=bet * ga, gamma=ga, lambda=la) )
      ss[item] <- s
      Alt[item] <- alt
      Bet[item, ] <- bet
      Gat[item] <- ga / (1 + ga)
      zz[item, ] <- z
      La[item, ] <- la
      Etat[item, ] <- etat
      
      pb$tick()
    }
  }
  
  ## Summary 
  Result <- list(al=Al[-(1:burn)], be=Be[-(1:burn),,drop=F], ga=Ga[-(1:burn)], 
                alt=Alt[-(1:burn)], bet=Bet[-(1:burn),,drop=F], gat=Gat[-(1:burn)], 
                s=ss[-(1:burn)], zz=zz[-(1:burn),,drop=F],
                la=La[-(1:burn),,drop=F], etat=Etat[-(1:burn)], 
                ac=ac, iter_bet=iter_bet, ll=log_likelihood[-(1:burn)])
  return(Result)
}







