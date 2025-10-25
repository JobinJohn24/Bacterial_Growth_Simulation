# ODE: dN/dt = r·N·(1 - N/K) - k_max · [A(t)^h / (EC50^h + A(t)^h)] · N
dose_A <- function(t, p) {
  if (p$dose_schedule == "constant") p$A_0 
  else if ((t %% p$pulse_interval) < p$pulse_duration) p$A_0 else 0
}
rhs <- function(t, y, p) {
  A_t <- dose_A(t, p); kill <- p$kmax * (A_t^p$h) / (p$EC50^p$h + A_t^p$h)
  list(c(p$r * y * (1 - y / p$K) - kill * y))
}
