# ============================================================================
# CDF PLOTS FOR THE PLBW DISTRIBUTION
# ============================================================================

rm(list = ls())

# ---------------------------------------------------------------------------
# CDF OF THE PLBW DISTRIBUTION
# ---------------------------------------------------------------------------

F_plbw <- function(y, alpha, beta, c)
{
  pgamma(alpha * y^(beta * c),
         shape = 1 + 1/beta,
         scale = 1)
}

# ---------------------------------------------------------------------------
# PARAMETER CONFIGURATIONS
# ---------------------------------------------------------------------------

configs <- list(
  
  list(alpha = 0.5, beta = 0.8, c = 0.7),
  
  list(alpha = 1.0, beta = 1.0, c = 1.0),
  
  list(alpha = 1.5, beta = 1.2, c = 1.0),
  
  list(alpha = 2.0, beta = 1.0, c = 1.5),
  
  list(alpha = 2.5, beta = 1.5, c = 1.2),
  
  list(alpha = 3.0, beta = 2.0, c = 1.5)
  
)

# ---------------------------------------------------------------------------
# GRID OF y VALUES
# ---------------------------------------------------------------------------

y <- seq(0.001, 5, length.out = 3000)

# ---------------------------------------------------------------------------
# COMPUTE CDF CURVES
# ---------------------------------------------------------------------------

curve_list <- lapply(configs, function(cfg)
{
  F_plbw(y,
         alpha = cfg$alpha,
         beta  = cfg$beta,
         c     = cfg$c)
})

# ---------------------------------------------------------------------------
# COLORS
# ---------------------------------------------------------------------------

cols <- c(
  "#1b9e77",
  "#d95f02",
  "#7570b3",
  "#e7298a",
  "#66a61e",
  "#e6ab02"
)

# ---------------------------------------------------------------------------
# LEGEND EXPRESSIONS
# ---------------------------------------------------------------------------

labels_expr <- as.expression(
  lapply(seq_along(configs), function(i)
  {
    cfg <- configs[[i]]
    
    bquote(
      C[.(i)] * ":" ~
        alpha == .(cfg$alpha) * "," ~
        beta == .(cfg$beta) * "," ~
        c == .(cfg$c)
    )
  })
)

# ---------------------------------------------------------------------------
# PLOT
# ---------------------------------------------------------------------------

par(mar = c(5,5,4,2)+0.1)

plot(y,
     curve_list[[1]],
     type = "n",
     xlab = "y",
     ylab = expression(F(y)),
     main = "CDF of the PLBW Distribution",
     ylim = c(-0.02,1.02),
     xlim = range(y),
     las = 1,
     cex.lab = 1.2,
     cex.axis = 1.1,
     cex.main = 1.3)

grid(col = "gray90", lty = "dotted")

abline(h = c(0,1),
       lty = "dotted",
       col = "gray80")

for(i in seq_along(curve_list))
{
  lines(y,
        curve_list[[i]],
        col = cols[i],
        lwd = 3)
}

legend("bottomright",
       legend = labels_expr,
       col = cols,
       lwd = 3,
       bty = "o",
       bg = "white",
       cex = 0.82)

box()

mtext("Power Length-Biased Weibull (PLBW) Distribution",
      side = 3,
      line = 0.5,
      col = "darkblue",
      cex = 0.95)

cat("\n✔ PLBW CDF plot generated successfully.\n")
