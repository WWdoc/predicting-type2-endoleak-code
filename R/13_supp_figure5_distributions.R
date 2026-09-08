## Analysis: Supplementary Figure 5 threshold-defined anatomical distributions
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "ggplot2", "patchwork", "tidyr"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
variables <- c("ostial_area_mm2", "ima_diameter_mm", "patent_la_number", "la_ge2_number")
require_columns(data, variables, "Supplementary Figure 5")

bar_panel <- function(variable, threshold, x_label, y_label) {
  d <- data.frame(value = as.numeric(data[[variable]])) |>
    tidyr::drop_na() |>
    dplyr::mutate(group = factor(ifelse(value <= threshold, "Lower", "Higher"), levels = c("Lower", "Higher")))
  summaries <- d |>
    dplyr::group_by(group) |>
    dplyr::summarise(mean = mean(value), sd = stats::sd(value), .groups = "drop")
  p <- stats::t.test(value ~ group, data = d)$p.value
  labels <- c(Lower = paste0("≤ ", threshold), Higher = paste0("> ", threshold))
  ggplot2::ggplot(summaries, ggplot2::aes(group, mean, fill = group)) +
    ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_fill_manual(values = c(Lower = "#91ADD0", Higher = "#286C98")) +
    ggplot2::annotate("text", x = 1.5, y = max(summaries$mean + summaries$sd) * 1.05,
                      label = paste0("P ", ifelse(p < 0.001, "< .001", paste0("= ", format_p_rsna(p)))), size = 3) +
    ggplot2::labs(x = x_label, y = y_label) + theme_rsna(9) +
    ggplot2::theme(legend.position = "none")
}

violin_panel <- function(variable, threshold, x_label, y_label) {
  d <- data.frame(value = as.numeric(data[[variable]])) |>
    tidyr::drop_na() |>
    dplyr::mutate(group = factor(ifelse(value <= threshold, "Lower", "Higher"), levels = c("Lower", "Higher")))
  p <- stats::wilcox.test(value ~ group, data = d, exact = FALSE)$p.value
  labels <- c(Lower = paste0("≤ ", threshold), Higher = paste0("> ", threshold))
  ggplot2::ggplot(d, ggplot2::aes(group, value, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.45, color = "black") +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.08, size = 0.8, alpha = 0.4) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_fill_manual(values = c(Lower = "#91ADD0", Higher = "#286C98")) +
    ggplot2::annotate("text", x = 1.5, y = max(d$value) * 1.05,
                      label = paste0("P ", ifelse(p < 0.001, "< .001", paste0("= ", format_p_rsna(p)))), size = 3) +
    ggplot2::labs(x = x_label, y = y_label) + theme_rsna(9) +
    ggplot2::theme(legend.position = "none")
}

plots <- list(
  bar_panel("ostial_area_mm2", cfg$thresholds[["ostial_area_mm2"]],
            "Branch-vessel ostial-area category", "Branch-vessel ostial area (mm²)"),
  bar_panel("ima_diameter_mm", cfg$thresholds[["ima_diameter_mm"]],
            "IMA diameter category", "IMA diameter (mm)"),
  violin_panel("patent_la_number", cfg$thresholds[["patent_la_number"]],
               "Number of patent lumbar arteries category", "Number of patent lumbar arteries"),
  violin_panel("la_ge2_number", cfg$thresholds[["la_ge2_number"]],
               "Number of lumbar arteries ≥2 mm category", "Number of lumbar arteries ≥2 mm")
)
combined <- patchwork::wrap_plots(plots, ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
save_publication_plot(combined, "supp_figure5_anatomical_distributions", io$output_dir, 8, 8)
