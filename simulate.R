# Generate deterministic solutions and add noise
simulate_curve <- function(A_level, p) {
  p$A_0 <- A_level; times <- seq(0, p$tmax, by = p$dt)
  out <- deSolve::ode(y = p$N0, times = times, func = rhs, parms = p, method = "lsoda")
  data.frame(time = out[, 1], density = out[, 2], concentration = A_level, replicate = 0)
}

make_synthetic <- function(df, p) {
  df$density <- pmax(df$density + rnorm(nrow(df), 0, p$noise_sd), 0); df
}
