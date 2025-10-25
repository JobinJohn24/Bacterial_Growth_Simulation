source("model.R"); source("simulate.R"); source("fit.R"); source("visualize.R"); source("io.R")
p <- read_params("sim_params.yaml"); set.seed(p$seed)
dir.create("plots", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)

df_user <- maybe_load_user_data()
if (is.null(df_user)) {
  df <- NULL; for (A_lv in p$A_levels) df <- rbind(df, make_synthetic(simulate_curve(A_lv, p), p))
} else { df <- df_user }

p_init <- list(N0 = p$N0, r = p$r, K = p$K, kmax = p$kmax, EC50 = p$EC50, h = p$h,
  tmax = p$tmax, dt = p$dt, seed = p$seed, A_0 = 1, dose_schedule = p$dose_schedule,
  pulse_interval = p$pulse_interval, pulse_duration = p$pulse_duration, noise_sd = p$noise_sd)

fixed <- setdiff(names(p_init), p$fit_free)
p_fit <- tryCatch(fit_params(df, p_init, fixed), error = function(e) { 
  cat("Fitting warning:", e$message, "\n"); p_init })

df_pred <- NULL
for (A_lv in p$A_levels) {
  p_tmp <- p_fit; p_tmp$A_0 <- A_lv; times <- seq(0, p$tmax, by = p$dt)
  out <- deSolve::ode(y = p$N0, times = times, func = rhs, parms = p_tmp, method = "lsoda")
  df_pred <- rbind(df_pred, data.frame(time = out[, 1], density = out[, 2], 
                                        concentration = A_lv, replicate = 0))
}

metrics <- compute_metrics(df_pred, p_fit)
plot_timecourses(df_pred, p$dose_schedule); plot_dose_response(metrics)
write_metrics(metrics, "results/metrics.csv")

best_auc_idx <- which.min(metrics$AUC[metrics$concentration > 0])
best_conc <- metrics$concentration[metrics$concentration > 0][best_auc_idx]
cat("\n=== RESULTS ===\nBest concentration:", best_conc, "μg/mL (AUC:", 
    sprintf("%.2e", metrics$AUC[metrics$concentration == best_conc]), ")\n")
cat("EC50:", sprintf("%.2f", p_fit$EC50), "μg/mL; Hill:", sprintf("%.2f", p_fit$h), "\n")
cat("Results saved to results/metrics.csv\n")
