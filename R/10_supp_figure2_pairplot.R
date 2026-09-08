## Analysis: Supplementary Figure 2 Spearman correlation pair plot
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "GGally", "ggplot2", "scales"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
variables <- c(
  ostial_area_mm2 = "Ostial Area", ima_diameter_mm = "IMA Diameter",
  patent_la_number = "LA Number", la_ge2_number = "No. LA ≥2 mm",
  infrarenal_angle_deg = "Infrarenal Angle", thrombus_area_mm2 = "Thrombus Area",
  thrombus_fraction = "Thrombus Proportion", aaa_max_diameter_mm = "AAA Max Diameter",
  aaa_volume_ml = "AAA Volume"
)
require_columns(data, names(variables), "Supplementary Figure 2")
plot_data <- data[, names(variables), drop = FALSE]
names(plot_data) <- unname(variables)

upper_panel <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  keep <- stats::complete.cases(x, y)
  if (sum(keep) < 3 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    rho <- p <- NA_real_
  } else {
    test <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    rho <- unname(test$estimate)
    p <- test$p.value
  }
  fill <- if (is.na(rho)) "white" else scales::col_numeric(
    palette = c("#BFD7E8", "white", "#B985C5"), domain = c(-1, 1)
  )(rho)
  labels <- data.frame(
    x = 1, y = c(1.12, 0.88),
    label = c(ifelse(is.na(rho), "NA", sprintf("%.2f", rho)), paste0("P", ifelse(is.na(p), " = NA", paste0(" ", format_p_rsna(p)))))
  )
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = data.frame(x = 1, y = 1), ggplot2::aes(x, y), fill = fill, color = "grey40") +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label), family = "Arial", size = 3.4) +
    ggplot2::xlim(0.5, 1.5) + ggplot2::ylim(0.5, 1.5) + ggplot2::theme_void()
}

lower_panel <- function(data, mapping, ...) {
  ggplot2::ggplot(data, mapping) +
    ggplot2::geom_point(color = "#9B59B6", alpha = 0.35, size = 0.8) +
    ggplot2::theme_bw(base_family = "Arial", base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
}

diagonal_panel <- function(data, mapping, ...) {
  ggplot2::ggplot(data, mapping) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 18,
                            fill = "#DCE8F5", color = "grey55", linewidth = 0.25) +
    ggplot2::geom_density(color = "#2C7BB6", linewidth = 0.7, na.rm = TRUE) +
    ggplot2::theme_bw(base_family = "Arial", base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
}

pair_plot <- GGally::ggpairs(
  plot_data, upper = list(continuous = upper_panel), lower = list(continuous = lower_panel),
  diag = list(continuous = diagonal_panel), progress = FALSE
) + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 8))
save_publication_plot(pair_plot, "supp_figure2_spearman_pairplot", io$output_dir, 12, 12)

correlations <- list()
for (i in seq_len(ncol(plot_data) - 1L)) {
  for (j in (i + 1L):ncol(plot_data)) {
    x <- plot_data[[i]]
    y <- plot_data[[j]]
    keep <- stats::complete.cases(x, y)
    test <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    correlations[[length(correlations) + 1L]] <- data.frame(
      Variable_1 = names(plot_data)[i], Variable_2 = names(plot_data)[j], N = sum(keep),
      Spearman_rho = unname(test$estimate), CI_lower = NA_real_, CI_upper = NA_real_, P_value = test$p.value
    )
  }
}
correlation_table <- dplyr::bind_rows(correlations)
bootstrap_spearman_ci <- function(variable_1, variable_2, repetitions = cfg$bootstrap_repetitions) {
  x <- plot_data[[variable_1]]
  y <- plot_data[[variable_2]]
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  estimates <- replicate(repetitions, {
    index <- sample(seq_along(x), length(x), replace = TRUE)
    suppressWarnings(stats::cor(x[index], y[index], method = "spearman"))
  })
  as.numeric(stats::quantile(estimates, c(0.025, 0.975), na.rm = TRUE))
}
bootstrap_cis <- Map(
  bootstrap_spearman_ci, correlation_table$Variable_1, correlation_table$Variable_2
)
correlation_table$CI_lower <- vapply(bootstrap_cis, `[[`, numeric(1), 1)
correlation_table$CI_upper <- vapply(bootstrap_cis, `[[`, numeric(1), 2)
correlation_table$P_value_Holm <- stats::p.adjust(correlation_table$P_value, method = "holm")
utils::write.csv(correlation_table, file.path(io$output_dir, "supp_figure2_spearman_values.csv"), row.names = FALSE)
print(correlation_table, row.names = FALSE)
