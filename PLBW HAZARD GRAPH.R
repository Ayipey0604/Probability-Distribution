##############################################################################
# Final Hazard Function Plot for the PLBW Distribution
# Professional Version
##############################################################################

rm(list = ls())

# ---------------------------------------------------------------------------
# Hazard Function of the PLBW Distribution
# ---------------------------------------------------------------------------
h_plbw <- function(x, alpha, beta, c){
  
  pdf <- (c * beta^2 * alpha^(1 + 1/beta) /
            gamma(1/beta)) *
    x^(beta*c + c - 1) *
    exp(-alpha * x^(beta*c))
  
  cdf <- pgamma(alpha * x^(beta*c),
                shape = 1 + 1/beta,
                lower.tail = TRUE)
  
  surv <- 1 - cdf
  
  surv[surv < 1e-12] <- 1e-12
  
  h <- pdf / surv
  
  return(h)
}

# ---------------------------------------------------------------------------
# Parameter Configurations (Three Parameters)
# ---------------------------------------------------------------------------
configs <- list(
  list(alpha=0.5, beta=0.8, c=0.8),
  list(alpha=1.0, beta=1.0, c=1.0),
  list(alpha=1.5, beta=1.2, c=1.1),
  list(alpha=2.0, beta=1.5, c=1.2),
  list(alpha=2.5, beta=1.8, c=1.4),
  list(alpha=3.0, beta=2.0, c=1.5)
)

# ---------------------------------------------------------------------------
# x-grid
# ---------------------------------------------------------------------------
x <- seq(0.001, 4, length.out = 1500)

# ---------------------------------------------------------------------------
# Compute Hazard Curves
# ---------------------------------------------------------------------------
hazards <- lapply(configs, function(cfg){
  
  h_plbw(x,
         alpha = cfg$alpha,
         beta  = cfg$beta,
         c     = cfg$c)
  
})

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
cols <- c("#1b9e77",
          "#d95f02",
          "#7570b3",
          "#e7298a",
          "#66a61e",
          "#e6ab02")

# ---------------------------------------------------------------------------
# Plot Limits
# ---------------------------------------------------------------------------
ymax <- max(unlist(hazards), na.rm = TRUE)

ylim <- c(0, ymax*1.10)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
par(mar=c(5.2,5.2,4.5,2.5))

plot(x,
     hazards[[1]],
     type="l",
     lwd=2.6,
     col=cols[1],
     ylim=ylim,
     xlim=range(x),
     xlab="x",
     ylab=expression(h(x)),
     main="Hazard Function of the PLBW Distribution",
     las=1)

grid(col="gray90", lty="solid")

for(i in 2:length(hazards)){
  lines(x,
        hazards[[i]],
        col=cols[i],
        lwd=2.6)
}

# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
legend("topleft",
       legend=c(
         expression(alpha==0.5~","~beta==0.8~","~c==0.8),
         expression(alpha==1.0~","~beta==1.0~","~c==1.0),
         expression(alpha==1.5~","~beta==1.2~","~c==1.1),
         expression(alpha==2.0~","~beta==1.5~","~c==1.2),
         expression(alpha==2.5~","~beta==1.8~","~c==1.4),
         expression(alpha==3.0~","~beta==2.0~","~c==1.5)
       ),
       col=cols,
       lwd=2.6,
       cex=0.80,
       bty="o",
       bg=adjustcolor("white",0.92),
       box.col="gray50")

# ---------------------------------------------------------------------------
# Subtitle
# ---------------------------------------------------------------------------
mtext("Power Length-Biased Weibull (PLBW) Distribution",
      side=3,
      line=0.4,
      cex=0.90,
      col="darkblue")

cat("\nHazard Function plot generated successfully.\n")