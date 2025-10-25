# Create time-course and dose-response plots
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
