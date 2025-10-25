#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
# Bacterial Growth Under Antibiotics - STANDALONE VERSION
# Complete pipeline in one file (105 lines of core code)
# Run: Rscript bacterial_growth_standalone.R
# ════════════════════════════════════════════════════════════════════════════

# ── MODEL LAYER ──────────────────────────────────────────────────────────────
# ODE: dN/dt = r·N·(1 - N/K) - k_max · [A(t)^h / (EC50^h + A(t)^h)] · N
dose_A <- function(t, p) {
  if (p$dose_schedule == "constant") p$A_0 
  else if ((t %% p$pulse_interval) < p$pulse_duration) p$A_0 else 0
}
rhs <- function(t, y, p) {
  A_t <- dose_A(t, p); kill <- p$kmax * (A_t^p$h) / (p$EC50^p$h + A_t^p$h)
  list(c(p$r * y * (1 - y / p$K) - kill * y))
}

# ── SIMULATION LAYER ─────────────────────────────────────────────────────────
simulate_curve <- function(A_level, p) {
  p$A_0 <- A_level; times <- seq(0, p$tmax, by = p$dt)
  out <- deSolve::ode(y = p$N0, times = times, func = rhs, parms = p, method = "lsoda")
  data.frame(time = out[, 1], density = out[, 2], concentration = A_level, replicate = 0)
}
make_synthetic <- function(df, p) {
  df$density <- pmax(df$density + rnorm(nrow(df), 0, p$noise_sd), 0); df
}

# ── FITTING LAYER ────────────────────────────────────────────────────────────
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

# ── VISUALIZATION LAYER ──────────────────────────────────────────────────────
plot_timecourses <- function(df, schedule) {
  g <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = density, color = factor(concentration))) +
    ggplot2::geom_line(size = 0.7) + ggplot2::facet_wrap(~ concentration) +
    ggplot2::labs(title = paste("Timecourses -", schedule), x = "Time (hr)", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(sprintf("plots/timecourses_%s.png", schedule), g, width = 10, height = 6)
}
plot_dose_response <- function(metrics) {
  ec50 <- metrics$estimated_EC50[1]
  g <- ggplot2::ggplot(metrics, ggplot2::aes(x = concentration, y = AUC)) +
    ggplot2::geom_point(size = 2) + ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = ec50, linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = ec50, y = max(metrics$AUC) * 0.95, 
                      label = sprintf("EC50=%.2f", ec50), hjust = -0.1) +
    ggplot2::labs(title = "Dose-Response Curve", x = "Antibiotic Conc", y = "AUC") +
    ggplot2::theme_minimal()
  ggplot2::ggsave("plots/dose_response.png", g, width = 7, height = 5)
}

# ── I/O LAYER ────────────────────────────────────────────────────────────────
read_params <- function(path = "sim_params.yaml") yaml::read_yaml(path)
write_metrics <- function(metrics, path) write.csv(metrics, path, row.names = FALSE)
maybe_load_user_data <- function(path = "data/growth.csv") 
  if (file.exists(path)) read.csv(path) else NULL

compute_metrics <- function(df, p_fit) {
  metrics <- data.frame()
  for (A_lv in unique(df$concentration)) {
    sub <- df[df$concentration == A_lv, ]; ord <- order(sub$time)
    sub <- sub[ord, ]
    auc <- sum(diff(sub$time) * (sub$density[-nrow(sub)] + sub$density[-1]) / 2, na.rm = TRUE)
    idx_50 <- which(sub$density <= 0.5 * p_fit$N0)[1]
    t_50 <- if (is.na(idx_50)) NA else sub$time[idx_50]
    metrics <- rbind(metrics, data.frame(concentration = A_lv, AUC = auc, 
      max_density = max(sub$density, na.rm = TRUE), time_to_50pct_reduction = t_50,
      estimated_r = p_fit$r, estimated_K = p_fit$K, estimated_kmax = p_fit$kmax,
      estimated_EC50 = p_fit$EC50, estimated_h = p_fit$h))
  }
  metrics
}

# ════════════════════════════════════════════════════════════════════════════
# ORCHESTRATION - MAIN PIPELINE
# ════════════════════════════════════════════════════════════════════════════

# Check for required packages
required_pkgs <- c("yaml", "deSolve", "ggplot2")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("Installing required packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, quietly = TRUE)
}
library(yaml); library(deSolve); library(ggplot2)

# Load configuration
if (!file.exists("sim_params.yaml")) {
  cat("Creating default sim_params.yaml...\n")
  cat("seed: 42\ntmax: 48\ndt: 0.1\nN0: 1.0e6\nr: 0.9\nK: 1.0e9\nkmax: 2.0\n",
      "EC50: 2.0\nh: 1.5\nA_levels: [0, 0.5, 1, 2, 5, 10]\n",
      "dose_schedule: 'constant'\npulse_interval: 6\npulse_duration: 1\n",
      "noise_sd: 5.0e6\nfit_free: ['r', 'K', 'EC50']\n", 
      file = "sim_params.yaml")
}

p <- read_params("sim_params.yaml")
set.seed(p$seed)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# Load or generate data
df_user <- maybe_load_user_data()
if (is.null(df_user)) {
  df <- NULL; for (A_lv in p$A_levels) df <- rbind(df, make_synthetic(simulate_curve(A_lv, p), p))
} else { df <- df_user }

# Prepare parameters
p_init <- list(N0 = p$N0, r = p$r, K = p$K, kmax = p$kmax, EC50 = p$EC50, h = p$h,
  tmax = p$tmax, dt = p$dt, seed = p$seed, A_0 = 1, dose_schedule = p$dose_schedule,
  pulse_interval = p$pulse_interval, pulse_duration = p$pulse_duration, noise_sd = p$noise_sd)

# Fit parameters
fixed <- setdiff(names(p_init), p$fit_free)
p_fit <- tryCatch(fit_params(df, p_init, fixed), error = function(e) { 
  cat("Fitting warning:", e$message, "\n"); p_init })

# Generate predictions
df_pred <- NULL
for (A_lv in p$A_levels) {
  p_tmp <- p_fit; p_tmp$A_0 <- A_lv; times <- seq(0, p$tmax, by = p$dt)
  out <- deSolve::ode(y = p$N0, times = times, func = rhs, parms = p_tmp, method = "lsoda")
  df_pred <- rbind(df_pred, data.frame(time = out[, 1], density = out[, 2], 
                                        concentration = A_lv, replicate = 0))
}

# Compute metrics
metrics <- compute_metrics(df_pred, p_fit)

# Generate plots
plot_timecourses(df_pred, p$dose_schedule)
plot_dose_response(metrics)

# Write results
write_metrics(metrics, "results/metrics.csv")

# Print summary
best_auc_idx <- which.min(metrics$AUC[metrics$concentration > 0])
best_conc <- metrics$concentration[metrics$concentration > 0][best_auc_idx]
cat("\n════════════════════════════════════════════════════════════════\n")
cat("Best antibiotic concentration:", best_conc, "μg/mL (AUC:", 
    sprintf("%.2e", metrics$AUC[metrics$concentration == best_conc]), ")\n")
cat("Estimated EC50:", sprintf("%.2f", p_fit$EC50), "μg/mL; Hill:", 
    sprintf("%.2f", p_fit$h), "\n")
cat("Results saved to results/metrics.csv\n")
cat("Plots saved to plots/\n")
cat("════════════════════════════════════════════════════════════════\n\n")
