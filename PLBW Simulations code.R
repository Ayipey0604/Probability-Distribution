###############################################################
#        MONTE CARLO SIMULATION FOR THE PLBW DISTRIBUTION
###############################################################

rm(list = ls())
set.seed(123)

###############################################################
# Required Packages
###############################################################

library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)

###############################################################
# True Parameter Values
###############################################################

alpha.true <- 2
beta.true  <- 2
c.true     <- 1.5

###############################################################
# PLBW Probability Density Function
###############################################################

dPLBW <- function(x, alpha, beta, c){
  
  const <- (c * beta^2 * alpha^(1 + 1/beta)) /
    gamma(1/beta)
  
  pdf <- const *
    x^(beta*c + c - 1) *
    exp(-alpha * x^(beta*c))
  
  pdf[x <= 0] <- 0
  
  return(pdf)
  
}

###############################################################
# PLBW Distribution Function
###############################################################

pPLBW <- function(x, alpha, beta, c){
  
  F <- pgamma(alpha * x^(beta*c),
              shape = 1 + 1/beta,
              lower.tail = TRUE)
  
  F[x <= 0] <- 0
  
  return(F)
  
}

###############################################################
# Numerical Quantile Function
###############################################################

qPLBW <- function(p, alpha, beta, c){
  
  sapply(p,function(pp){
    
    uniroot(function(x){
      
      pPLBW(x,alpha,beta,c)-pp
      
    },
    interval=c(1e-8,50),
    tol=1e-10)$root
    
  })
  
}

###############################################################
# Random Number Generator
###############################################################

rPLBW <- function(n,alpha,beta,c){
  
  u <- runif(n)
  
  qPLBW(u,alpha,beta,c)
  
}

###############################################################
# Example
###############################################################

set.seed(123)

sample.data <- rPLBW(10,
                     alpha.true,
                     beta.true,
                     c.true)

sample.data

######Part 2######

###############################################################
# LOG-LIKELIHOOD FUNCTION
###############################################################

loglik_PLBW <- function(par, x){
  
  alpha <- par[1]
  beta  <- par[2]
  c     <- par[3]
  
  # Parameter constraints
  if(alpha <= 0 || beta <= 0 || c <= 0)
    return(1e10)
  
  pdf <- dPLBW(x, alpha, beta, c)
  
  # Avoid log(0)
  pdf[pdf <= 0] <- 1e-15
  
  # Negative log-likelihood
  return(-sum(log(pdf)))
  
}

###############################################################
# PARAMETER ESTIMATION FUNCTION
###############################################################

estimate_PLBW <- function(x){
  
  start <- c(
    alpha = runif(1,0.5,3),
    beta  = runif(1,0.5,3),
    c      = runif(1,0.5,3)
  )
  
  fit <- try(
    
    optim(
      par = start,
      fn = loglik_PLBW,
      x = x,
      method = "L-BFGS-B",
      lower = c(0.001,0.001,0.001),
      upper = c(20,20,20)
    ),
    
    silent = TRUE
    
  )
  
  #############################################################
  # Retry if optimization fails
  #############################################################
  
  if(inherits(fit,"try-error") || fit$convergence != 0){
    
    start2 <- c(
      alpha = runif(1,1,4),
      beta  = runif(1,1,4),
      c      = runif(1,1,4)
    )
    
    fit <- try(
      
      optim(
        par = start2,
        fn = loglik_PLBW,
        x = x,
        method = "L-BFGS-B",
        lower = c(0.001,0.001,0.001),
        upper = c(20,20,20)
      ),
      
      silent = TRUE
      
    )
    
  }
  
  #############################################################
  # Return estimates
  #############################################################
  
  if(inherits(fit,"try-error") || fit$convergence != 0){
    
    return(rep(NA,3))
    
  }
  
  return(fit$par)
  
}


######Part 3######

###############################################################
# PERFORMANCE MEASURES
###############################################################

performance <- function(estimates, true){
  
  estimates <- estimates[!is.na(estimates)]
  
  Mean <- mean(estimates)
  
  Bias <- Mean - true
  
  SD <- sd(estimates)
  
  MSE <- mean((estimates - true)^2)
  
  RMSE <- sqrt(MSE)
  
  data.frame(
    
    Mean = Mean,
    Bias = Bias,
    SD = SD,
    MSE = MSE,
    RMSE = RMSE
    
  )
  
}

###############################################################
# MONTE CARLO SIMULATION
###############################################################

MonteCarlo <- function(n,
                       R,
                       alpha,
                       beta,
                       c){
  
  estimates <- matrix(NA,R,3)
  
  colnames(estimates) <- c("alpha","beta","c")
  
  for(i in 1:R){
    
    x <- rPLBW(n,
               alpha,
               beta,
               c)
    
    estimates[i,] <- estimate_PLBW(x)
    
    if(i %% 100 == 0)
      cat("Replication:",i,"\n")
    
  }
  
  estimates <- as.data.frame(estimates)
  
  result.alpha <- performance(estimates$alpha,alpha)
  
  result.beta <- performance(estimates$beta,beta)
  
  result.c <- performance(estimates$c,c)
  
  results <- rbind(
    
    cbind(Parameter="alpha",result.alpha),
    
    cbind(Parameter="beta",result.beta),
    
    cbind(Parameter="c",result.c)
    
  )
  
  return(list(
    
    Estimates=estimates,
    
    Summary=results
    
  ))
  
}


######Part 4#####

###############################################################
# SAMPLE SIZES
###############################################################

sample.sizes <- c(
  30,
  50,
  100,
  150,
  200,
  300,
  400,
  600
)

###############################################################
# NUMBER OF REPLICATIONS
###############################################################

R <- 10000

###############################################################
# TRUE PARAMETERS
###############################################################

alpha.true <- 2

beta.true <- 2

c.true <- 1.5

###############################################################
# STORAGE
###############################################################

All.Results <- list()

###############################################################
# RUN SIMULATION
###############################################################

for(n in sample.sizes){
  
  cat("\n")
  cat("=============================\n")
  cat("Sample Size =",n,"\n")
  cat("=============================\n")
  
  out <- MonteCarlo(
    
    n=n,
    R=R,
    alpha=alpha.true,
    beta=beta.true,
    c=c.true
    
  )
  
  print(out$Summary)
  
  All.Results[[paste0("n_",n)]] <- out$Summary
  
}

#####Part 5####

###############################################################
# COMBINE RESULTS
###############################################################

Final.Results <- data.frame()

for(i in seq_along(All.Results)){
  
  temp <- All.Results[[i]]
  
  temp$n <- sample.sizes[i]
  
  Final.Results <- rbind(
    Final.Results,
    temp
  )
  
}

###############################################################
# DISPLAY RESULTS
###############################################################

print(Final.Results)

###############################################################
# SAVE RESULTS
###############################################################

write.csv(
  
  Final.Results,
  
  "PLBW_MonteCarlo_Results.csv",
  
  row.names = FALSE
  
)

cat("\n")
cat("Simulation Completed Successfully.\n")


######Part 6####

###############################################################
# BIAS PLOT
###############################################################

ggplot(Final.Results,
       aes(n,Bias,
           colour=Parameter,
           group=Parameter))+
  
  geom_line(size=1)+
  
  geom_point(size=2)+
  
  theme_bw()+
  
  labs(
    
    title="Bias of the PLBW Parameter Estimates",
    
    x="Sample Size",
    
    y="Bias"
    
  )


#####Part 7###

###############################################################
# MSE PLOT
###############################################################

ggplot(Final.Results,
       aes(n,MSE,
           colour=Parameter,
           group=Parameter))+
  
  geom_line(size=1)+
  
  geom_point(size=2)+
  
  theme_bw()+
  
  labs(
    
    title="Mean Squared Error of the PLBW Parameter Estimates",
    
    x="Sample Size",
    
    y="MSE"
    
  )


#####Part 8###

###############################################################
# RMSE PLOT
###############################################################

ggplot(Final.Results,
       aes(n,RMSE,
           colour=Parameter,
           group=Parameter))+
  
  geom_line(size=1)+
  
  geom_point(size=2)+
  
  theme_bw()+
  
  labs(
    
    title="Root Mean Squared Error of the PLBW Parameter Estimates",
    
    x="Sample Size",
    
    y="RMSE"
    
  )