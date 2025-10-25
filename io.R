# I/O and metrics computation
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
