## Analysis: Supplementary Figure 4 bootstrap validation of AUC and cut-off
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "ggplot2", "patchwork", "pROC", "tidyr"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
event_col <- "t2el"
marker <- "ostial_area_mm2"
require_columns(data, c(event_col, marker), "Supplementary Figure 4")
d <- data.frame(event = as_binary(data[[event_col]]), marker = as.numeric(data[[marker]])) |>
  tidyr::drop_na()
if (length(unique(d$event)) != 2) stop("Both outcome classes are required.")

original_roc <- pROC::roc(d$event, d$marker, levels = c(0, 1), direction = "<", quiet = TRUE)
original_auc <- as.numeric(pROC::auc(original_roc))
original_auc_ci <- as.numeric(pROC::ci.auc(original_roc, method = "delong"))
original_best <- pROC::coords(
  original_roc, "best", best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
)

bootstrap <- data.frame(
  replicate = seq_len(cfg$bootstrap_repetitions), AUC = NA_real_, threshold = NA_real_,
  sensitivity = NA_real_, specificity = NA_real_
)
for (b in seq_len(cfg$bootstrap_repetitions)) {
  index <- sample(seq_len(nrow(d)), nrow(d), replace = TRUE)
  sample_data <- d[index, , drop = FALSE]
  if (length(unique(sample_data$event)) != 2) next
  roc <- tryCatch(
    pROC::roc(sample_data$event, sample_data$marker, levels = c(0, 1), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc)) next
  best <- pROC::coords(
    roc, "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
  )
  bootstrap$AUC[b] <- as.numeric(pROC::auc(roc))
  bootstrap$threshold[b] <- as.numeric(best$threshold[1])
  bootstrap$sensitivity[b] <- as.numeric(best$sensitivity[1])
  bootstrap$specificity[b] <- as.numeric(best$specificity[1])
}
valid <- bootstrap[stats::complete.cases(bootstrap[, c("AUC", "threshold")]), , drop = FALSE]
if (nrow(valid) < 0.9 * cfg$bootstrap_repetitions) warning("More than 10% of bootstrap resamples were invalid.")

summary <- data.frame(
  Metric = c("AUC", "Youden cut-off, mm²", "Sensitivity", "Specificity"),
  Original = c(original_auc, as.numeric(original_best$threshold[1]),
               as.numeric(original_best$sensitivity[1]), as.numeric(original_best$specificity[1])),
  Bootstrap_median = c(median(valid$AUC), median(valid$threshold),
                       median(valid$sensitivity), median(valid$specificity)),
  Bootstrap_CI_lower = c(stats::quantile(valid$AUC, 0.025), stats::quantile(valid$threshold, 0.025),
                         stats::quantile(valid$sensitivity, 0.025), stats::quantile(valid$specificity, 0.025)),
  Bootstrap_CI_upper = c(stats::quantile(valid$AUC, 0.975), stats::quantile(valid$threshold, 0.975),
                         stats::quantile(valid$sensitivity, 0.975), stats::quantile(valid$specificity, 0.975)),
  Valid_repetitions = nrow(valid)
)

auc_plot <- ggplot2::ggplot(valid, ggplot2::aes(AUC)) +
  ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 35,
                          fill = "#DCEAF7", color = "grey55") +
  ggplot2::geom_density(color = "#2C7BB6", linewidth = 1) +
  ggplot2::geom_vline(xintercept = original_auc, color = "#C44E52", linetype = "dashed") +
  ggplot2::labs(x = "AUC", y = "Density") + theme_rsna(10)
threshold_plot <- ggplot2::ggplot(valid, ggplot2::aes(threshold)) +
  ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 35,
                          fill = "#EADBED", color = "grey55") +
  ggplot2::geom_density(color = "#9B59B6", linewidth = 1) +
  ggplot2::geom_vline(xintercept = as.numeric(original_best$threshold[1]), color = "#C44E52", linetype = "dashed") +
  ggplot2::labs(x = "Youden-derived cut-off (mm²)", y = "Density") + theme_rsna(10)
combined <- auc_plot + threshold_plot + patchwork::plot_annotation(tag_levels = "A")

dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
save_publication_plot(combined, "supp_figure4_bootstrap_validation", io$output_dir, 10, 4.8)
utils::write.csv(summary, file.path(io$output_dir, "supp_figure4_bootstrap_summary.csv"), row.names = FALSE)
utils::write.csv(valid, file.path(io$output_dir, "supp_figure4_bootstrap_replicates.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
