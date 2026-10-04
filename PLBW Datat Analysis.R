###############################################################################
#                     PART 1
#      POWER LENGTH-BIASED WEIBULL (PLBW) DISTRIBUTION
#          DATA + PDF + CDF + QUANTILE + RANDOM GENERATOR
###############################################################################

rm(list = ls())

###############################################################################
# REQUIRED PACKAGES
###############################################################################

library(stats)
library(numDeriv)
library(goftest)

###############################################################################
# FAILURE TIME DATA (Murthy et al., 2004)
###############################################################################

data <- c(
  0.036,0.058,0.061,0.074,0.078,0.086,0.102,0.103,0.114,0.116,
  0.148,0.183,0.192,0.254,0.262,0.379,0.381,0.538,0.570,0.574,
  0.590,0.618,0.645,0.961,1.228,1.600,2.006,2.054,2.804,3.058,
  3.076,3.147,3.625,3.704,3.931,4.037,4.393,4.534,4.893,6.274,
  6.816,7.896,7.904,8.022,9.337,10.940,11.020,13.880,14.730,15.080
)

n <- length(data)

###############################################################################
# POWER LENGTH-BIASED WEIBULL PDF
#
# f(y)=
# [c β² α^(1+1/β) y^(βc+c−1) exp(-αy^(βc))]
# -------------------------------------------
#            Γ(1/β)
###############################################################################

dplbw <- function(x, alpha, beta, c){
  
  f <- numeric(length(x))
  
  ind <- which(x > 0)
  
  if(length(ind)==0) return(f)
  
  xx <- x[ind]
  
  f[ind] <-
    (c *
       beta^2 *
       alpha^(1 + 1/beta) *
       xx^(beta*c + c - 1) *
       exp(-alpha * xx^(beta*c))) /
    gamma(1/beta)
  
  f[!is.finite(f)] <- 1e-300
  
  f[f <= 0] <- 1e-300
  
  return(f)
  
}

###############################################################################
# POWER LENGTH-BIASED WEIBULL CDF
#
# F(y)=γ(1+1/β, αy^(βc))/Γ(1+1/β)
###############################################################################

pplbw <- function(x, alpha, beta, c){
  
  F <- numeric(length(x))
  
  ind <- which(x > 0)
  
  if(length(ind)==0) return(F)
  
  xx <- x[ind]
  
  F[ind] <- pgamma(
    
    alpha * xx^(beta*c),
    
    shape = 1 + 1/beta,
    
    rate = 1
    
  )
  
  F[F < 0] <- 0
  
  F[F > 1] <- 1
  
  return(F)
  
}

###############################################################################
# SURVIVAL FUNCTION
###############################################################################

splbw <- function(x, alpha, beta, c){
  
  1 - pplbw(x, alpha, beta, c)
  
}

###############################################################################
# HAZARD FUNCTION
###############################################################################

hplbw <- function(x, alpha, beta, c){
  
  f <- dplbw(x, alpha, beta, c)
  
  s <- splbw(x, alpha, beta, c)
  
  s[s <= 1e-12] <- 1e-12
  
  f/s
  
}

###############################################################################
# CUMULATIVE HAZARD FUNCTION
###############################################################################

Hplbw <- function(x, alpha, beta, c){
  
  -log(splbw(x, alpha, beta, c))
  
}

###############################################################################
# QUANTILE FUNCTION
#
# Numerical inversion using uniroot()
###############################################################################

qplbw <- function(u, alpha, beta, c){
  
  sapply(u,function(p){
    
    uniroot(
      
      function(y)
        
        pplbw(y, alpha, beta, c) - p,
      
      lower = 1e-10,
      
      upper = 100
      
    )$root
    
  })
  
}

###############################################################################
# RANDOM NUMBER GENERATOR
###############################################################################

rplbw <- function(n, alpha, beta, c){
  
  u <- runif(n)
  
  qplbw(
    
    u,
    
    alpha,
    
    beta,
    
    c
    
  )
  
}

###############################################################################
# BASIC CHECKS
###############################################################################

cat("=================================================\n")
cat("POWER LENGTH-BIASED WEIBULL DISTRIBUTION\n")
cat("=================================================\n\n")

cat("Sample Size :", n, "\n\n")

cat("First Five Density Values\n")
print(
  dplbw(
    data[1:5],
    alpha = 2,
    beta = 2,
    c = 1.5
  )
)

cat("\n")

cat("First Five CDF Values\n")
print(
  pplbw(
    data[1:5],
    alpha = 2,
    beta = 2,
    c = 1.5
  )
)

cat("\n")

cat("First Five Hazard Values\n")
print(
  hplbw(
    data[1:5],
    alpha = 2,
    beta = 2,
    c = 1.5
  )
)

cat("\nPart 1 Loaded Successfully.\n")


###############################################################################
# 4. NEGATIVE LOG-LIKELIHOOD FUNCTION
###############################################################################

negloglik_plbw <- function(par, x){
  
  alpha <- par[1]
  beta  <- par[2]
  c     <- par[3]
  
  if(any(par <= 0))
    return(1e10)
  
  dens <- dplbw(x, alpha, beta, c)
  
  dens[dens <= 0] <- 1e-300
  
  if(any(!is.finite(dens)))
    return(1e10)
  
  -sum(log(dens))
  
}

###############################################################################
# 5. MULTI-START MAXIMUM LIKELIHOOD ESTIMATION
###############################################################################

starting.values <- list(
  
  c(1,1,1),
  
  c(2,2,1.5),
  
  c(0.5,0.5,0.8),
  
  c(3,1,2),
  
  c(1,3,1),
  
  c(5,2,2),
  
  c(2,5,1),
  
  c(4,4,2)
  
)

best.fit <- NULL

best.value <- Inf

for(start in starting.values){
  
  fit <- try(
    
    optim(
      
      par=start,
      
      fn=negloglik_plbw,
      
      x=data,
      
      method="L-BFGS-B",
      
      lower=c(0.001,0.001,0.001),
      
      upper=c(20,20,20),
      
      control=list(maxit=10000),
      
      hessian=TRUE
      
    ),
    
    silent=TRUE
    
  )
  
  if(!inherits(fit,"try-error")){
    
    if(fit$value < best.value){
      
      best.value <- fit$value
      
      best.fit <- fit
      
    }
    
  }
  
}

fit <- best.fit

###############################################################################
# ESTIMATED PARAMETERS
###############################################################################

alpha_hat <- fit$par[1]

beta_hat <- fit$par[2]

c_hat <- fit$par[3]

cat("Estimated Parameters\n")

print(fit$par)


###############################################################################
# 6. STANDARD ERRORS
###############################################################################

vcov.mat <- solve(fit$hessian)

SE <- sqrt(diag(vcov.mat))

z.value <- fit$par/SE

p.value <- 2*(1-pnorm(abs(z.value)))

Lower95 <- fit$par-1.96*SE

Upper95 <- fit$par+1.96*SE

results <- data.frame(
  
  Parameter=c("alpha","beta","c"),
  
  Estimate=round(fit$par,6),
  
  Std.Error=round(SE,6),
  
  z.value=round(z.value,4),
  
  p.value=round(p.value,5),
  
  Lower95=round(Lower95,6),
  
  Upper95=round(Upper95,6)
  
)

print(results)

###############################################################################
# 7. MODEL SELECTION CRITERIA
###############################################################################

logLik <- -fit$value

k <- 3

n <- length(data)

AIC <- -2*logLik+2*k

BIC <- -2*logLik+k*log(n)

AICc <- AIC+(2*k*(k+1))/(n-k-1)

cat("\n")

cat("Log-Likelihood :",round(logLik,5),"\n")

cat("AIC :",round(AIC,5),"\n")

cat("BIC :",round(BIC,5),"\n")

cat("AICc :",round(AICc,5),"\n")


###############################################################################
# 8. HISTOGRAM WITH FITTED PLBW DENSITY
###############################################################################

xgrid <- seq(
  
  min(data),
  
  max(data),
  
  length.out=1000
  
)

pdf.fit <- dplbw(
  
  xgrid,
  
  alpha_hat,
  
  beta_hat,
  
  c_hat
  
)

hist(
  
  data,
  
  probability=TRUE,
  
  breaks=10,
  
  col="grey85",
  
  border="white",
  
  main="Histogram with Fitted PLBW Density",
  
  xlab="Failure Time"
  
)

lines(
  
  xgrid,
  
  pdf.fit,
  
  col="red",
  
  lwd=3
  
)

rug(data)

legend(
  
  "topright",
  
  legend="PLBW",
  
  col="red",
  
  lwd=3,
  
  bty="n"
  
)


###############################################################################
# 9. EMPIRICAL VS FITTED CDF
###############################################################################

plot(
  
  ecdf(data),
  
  verticals=TRUE,
  
  do.points=FALSE,
  
  main="Empirical vs Fitted PLBW CDF",
  
  xlab="Failure Time",
  
  ylab="CDF"
  
)

lines(
  
  xgrid,
  
  pplbw(
    
    xgrid,
    
    alpha_hat,
    
    beta_hat,
    
    c_hat
    
  ),
  
  col="blue",
  
  lwd=3
  
)

legend(
  
  "bottomright",
  
  legend=c("Empirical","PLBW"),
  
  col=c("black","blue"),
  
  lwd=2,
  
  bty="n"
  
)


###############################################################################
# 10. Q-Q PLOT
###############################################################################

u <- ppoints(n)

theoretical.quantiles <- qplbw(
  u,
  alpha_hat,
  beta_hat,
  c_hat
)

observed.quantiles <- sort(data)

plot(
  theoretical.quantiles,
  observed.quantiles,
  pch = 19,
  col = "blue",
  cex = 1.1,
  xlab = "Theoretical Quantiles",
  ylab = "Observed Quantiles",
  main = "Q-Q Plot for the PLBW Distribution"
)

abline(0,1,
       col="red",
       lwd=2,
       lty=2)

grid()



###############################################################################
# 11. P-P PLOT
###############################################################################

empirical.prob <- (1:n - 0.5)/n

fitted.prob <- pplbw(
  sort(data),
  alpha_hat,
  beta_hat,
  c_hat
)

plot(
  fitted.prob,
  empirical.prob,
  pch=19,
  col="darkgreen",
  cex=1.1,
  xlab="Theoretical Probabilities",
  ylab="Empirical Probabilities",
  main="P-P Plot for the PLBW Distribution"
)

abline(0,1,
       col="red",
       lwd=2,
       lty=2)

grid()


library(goftest)


###############################################################################
# 12. KOLMOGOROV-SMIRNOV TEST
###############################################################################

KS <- ks.test(
  data,
  function(x)
    pplbw(
      x,
      alpha_hat,
      beta_hat,
      c_hat
    )
)

print(KS)


###############################################################################
# 13. ANDERSON-DARLING TEST
###############################################################################

AD <- ad.test(
  data,
  null=function(x)
    pplbw(
      x,
      alpha_hat,
      beta_hat,
      c_hat
    )
)

print(AD)


###############################################################################
# 14. CRAMER-VON MISES TEST
###############################################################################

CVM <- cvm.test(
  data,
  null=function(x)
    pplbw(
      x,
      alpha_hat,
      beta_hat,
      c_hat
    )
)

print(CVM)



###############################################################################
# 15. PARAMETER SUMMARY
###############################################################################

summary.table <- data.frame(
  
  Parameter=c("alpha","beta","c"),
  
  Estimate=round(c(alpha_hat,beta_hat,c_hat),6),
  
  Std.Error=round(SE,6),
  
  z.value=round(z.value,4),
  
  p.value=round(p.value,6),
  
  Lower95=round(Lower95,6),
  
  Upper95=round(Upper95,6)
  
)

print(summary.table)



###############################################################################
# 16. FINAL OUTPUT
###############################################################################

cat("\n")
cat("===========================================================\n")
cat("        POWER LENGTH-BIASED WEIBULL DISTRIBUTION\n")
cat("===========================================================\n\n")

print(summary.table)

cat("\n")

cat("Sample Size       :",n,"\n")
cat("Log-Likelihood    :",round(logLik,6),"\n")
cat("AIC               :",round(AIC,6),"\n")
cat("BIC               :",round(BIC,6),"\n")
cat("AICc              :",round(AICc,6),"\n")

cat("\n")
cat("=============== GOODNESS OF FIT ===============\n")

cat("KS Statistic      :",round(KS$statistic,6),"\n")
cat("KS p-value        :",KS$p.value,"\n\n")

cat("AD Statistic      :",round(AD$statistic,6),"\n")
cat("AD p-value        :",AD$p.value,"\n\n")

cat("CVM Statistic     :",round(CVM$statistic,6),"\n")
cat("CVM p-value       :",CVM$p.value,"\n")

cat("===========================================================\n")




###############################################################################
# 17. SAVE RESULTS
###############################################################################

write.csv(
  summary.table,
  "PLBW_Parameter_Estimates.csv",
  row.names=FALSE
)

fit.indices <- data.frame(
  
  Statistic=c(
    "LogLik",
    "AIC",
    "BIC",
    "AICc",
    "KS",
    "KS_pvalue",
    "AD",
    "AD_pvalue",
    "CVM",
    "CVM_pvalue"
  ),
  
  Value=c(
    logLik,
    AIC,
    BIC,
    AICc,
    KS$statistic,
    KS$p.value,
    AD$statistic,
    AD$p.value,
    CVM$statistic,
    CVM$p.value
  )
  
)

write.csv(
  fit.indices,
  "PLBW_Model_Fit.csv",
  row.names=FALSE
)

cat("\nResults successfully saved.\n")



describe(data)
summary(data)
