## Analysis: Supplementary Figure 1 Bland-Altman agreement plots
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "ggplot2", "irr", "patchwork"))

io <- resolve_io(cfg$operator_file)
data <- read_analysis_data(io$input_file)
results <- agreement_comparisons(data)

make_panel <- function(result, title) {
  points <- data.frame(average = result$average, difference = result$difference)
  label <- sprintf(
    "ICC = %.2f (95%% CI: %.2f–%.2f)\nBias = %.2f\n95%% LoA: %.2f to %.2f",
    result$icc, result$icc_lower, result$icc_upper,
    result$bias, result$loa_lower, result$loa_upper
  )
  ggplot2::ggplot(points, ggplot2::aes(average, difference)) +
    ggplot2::geom_point(shape = 21, fill = "#78A7CC", color = "black", size = 2) +
    ggplot2::geom_hline(yintercept = result$bias, color = "#A84255", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(result$loa_lower, result$loa_upper), linetype = "dashed", linewidth = 0.7) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = label, hjust = 1.05, vjust = 1.2, size = 3) +
    ggplot2::labs(title = title, x = "Mean of paired measurements (mm²)", y = "Difference (mm²)") +
    theme_rsna(10)
}

panel_titles <- c(
  "Observer W: W1 vs W2", "Observer Z: Z1 vs Z2",
  "Inter-rater at time 1: W1 vs Z1", "Inter-rater at time 2: W2 vs Z2"
)
plots <- Map(make_panel, results, panel_titles)
combined <- patchwork::wrap_plots(plots, ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
save_publication_plot(combined, "supp_figure1_bland_altman", io$output_dir, 10, 8)
statistics <- agreement_results_frame(results)
utils::write.csv(statistics, file.path(io$output_dir, "supp_figure1_agreement_statistics.csv"), row.names = FALSE)
print(statistics, row.names = FALSE)
