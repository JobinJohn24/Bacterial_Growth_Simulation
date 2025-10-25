# Fits model parameters via least squares
fit_params <- function(df, p_init, fixed_params) {
  obj <- function(theta) {
    p <- p_init; names(theta) <- names(p_init)[!(names(p_init) %in% fixed_params)]
    for (nm in names(theta)) p[[nm]] <- theta[nm]
    pred <- NULL
    for (A_lv in unique(df$concentration)) {
      p$A_0 <- A_lv; times <- sort(unique(df$time[df$concentration == A_lv]))
      out <- deSolve::ode(y = p$N0, times = times, func = rhs, parms = p, method = "lsoda")
      pred <- rbind(pred, data.frame(time = out[, 1], density = out[, 2], concentration = A_lv))
    }
    merged <- merge(df[, c("time", "density", "concentration")], pred, 
                    by = c("time", "concentration"), suffixes = c(".obs", ".pred"))
    sum((merged$density.obs - merged$density.pred)^2, na.rm = TRUE)
  }
  free_nms <- names(p_init)[!(names(p_init) %in% fixed_params)]
  theta0 <- p_init[free_nms]; result <- optim(theta0, obj, method = "Nelder-Mead")
  p_init[free_nms] <- result$par; p_init
}
