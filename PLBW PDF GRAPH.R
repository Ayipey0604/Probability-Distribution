# ============================================================================
# PDF PLOTS FOR THE PLBW DISTRIBUTION
# ============================================================================

rm(list = ls())

# ---------------------------------------------------------------------------
# PDF OF THE PLBW DISTRIBUTION
# ---------------------------------------------------------------------------

f_plbw <- function(y, alpha, beta, c)
{
  const <- (c * beta^2 * alpha^(1 + 1/beta)) /
    gamma(1/beta)
  
  const *
    y^(beta*c + c - 1) *
    exp(-alpha * y^(beta*c))
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
# GRID OF Y VALUES
# ---------------------------------------------------------------------------

y <- seq(0.001, 5, length.out = 3000)

# ---------------------------------------------------------------------------
# COMPUTE PDF CURVES
# ---------------------------------------------------------------------------

curve_list <- lapply(configs, function(cfg)
{
  f_plbw(y,
         alpha = cfg$alpha,
         beta  = cfg$beta,
         c     = cfg$c)
})

# ---------------------------------------------------------------------------
# DETERMINE AXIS LIMITS
# ---------------------------------------------------------------------------

all_vals <- unlist(curve_list)

ymax <- max(all_vals[is.finite(all_vals)]) * 1.15

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
     ylim = c(0, ymax),
     xlim = range(y),
     xlab = "y",
     ylab = expression(f(y)),
     main = "PDF of the PLBW Distribution",
     las = 1,
     cex.lab = 1.2,
     cex.axis = 1.1,
     cex.main = 1.3)

grid(col = "gray90")

for(i in seq_along(curve_list))
{
  lines(y,
        curve_list[[i]],
        col = cols[i],
        lwd = 3)
}

legend("topright",
       legend = labels_expr,
       col = cols,
       lwd = 3,
       bty = "n",
       cex = 0.85)

box()

mtext("Power Length-Biased Weibull (PLBW) Distribution",
      side = 3,
      line = 0.5,
      col = "darkblue",
      cex = 0.95)

cat("\n✔ PLBW PDF plot generated successfully.\n")
