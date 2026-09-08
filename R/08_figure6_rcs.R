## Analysis: Figure 6 restricted cubic spline analyses
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("broom", "dplyr", "ggplot2", "patchwork", "splines", "survival"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
time_col <- "t2el_followup_years"
event_col <- "t2el"
exposures <- c(
  ostial_area_mm2 = "Total branch-vessel ostial area (mm²)",
  ima_diameter_mm = "IMA diameter (mm)",
  patent_la_number = "Number of patent lumbar arteries",
  la_ge2_number = "Number of lumbar arteries ≥2 mm"
)
adjustment <- cfg$comprehensive_covariates
require_columns(data, unique(c(time_col, event_col, names(exposures), adjustment)), "Figure 6 RCS analyses")
data[[event_col]] <- as_binary(data[[event_col]])

formula_from_rhs <- function(rhs) {
  stats::as.formula(paste0(
    "survival::Surv(", backtick_names(time_col), ", ", backtick_names(event_col), ") ~ ", rhs
  ))
}

reference_value <- function(x) {
  if (is.factor(x)) return(factor(names(which.max(table(x))), levels = levels(x)))
  if (is.character(x)) return(names(which.max(table(x))))
  stats::median(as.numeric(x), na.rm = TRUE)
}

draw_rcs <- function(exposure, label) {
  columns <- unique(c(time_col, event_col, exposure, adjustment))
  d <- data[stats::complete.cases(data[, columns, drop = FALSE]), columns, drop = FALSE]
  d <- d[d[[time_col]] > 0, , drop = FALSE]
  if (length(unique(d[[exposure]])) < 5) stop("At least five unique values are required for ", exposure)

  adjustment_rhs <- if (length(adjustment)) paste(backtick_names(adjustment), collapse = " + ") else "1"
  spline_term <- paste0("splines::ns(", backtick_names(exposure), ", df = 3)")
  spline_rhs <- paste(c(spline_term, if (adjustment_rhs != "1") adjustment_rhs), collapse = " + ")
  linear_rhs <- paste(c(backtick_names(exposure), if (adjustment_rhs != "1") adjustment_rhs), collapse = " + ")
  fit_spline <- survival::coxph(formula_from_rhs(spline_rhs), data = d, x = TRUE)
  fit_linear <- survival::coxph(formula_from_rhs(linear_rhs), data = d, x = TRUE)
  fit_null <- survival::coxph(formula_from_rhs(adjustment_rhs), data = d, x = TRUE)
  p_overall <- stats::anova(fit_null, fit_spline, test = "LRT")[["Pr(>|Chi|)"]][2]
  p_nonlinear <- stats::anova(fit_linear, fit_spline, test = "LRT")[["Pr(>|Chi|)"]][2]
  linear_row <- broom::tidy(fit_linear, exponentiate = TRUE, conf.int = TRUE)
  linear_row <- linear_row[linear_row$term == exposure, , drop = FALSE]

  base_row <- d[1, , drop = FALSE]
  for (column in adjustment) base_row[[column]] <- reference_value(d[[column]])
  threshold <- unname(cfg$thresholds[exposure])
  base_row[[exposure]] <- threshold
  grid_values <- seq(
    stats::quantile(d[[exposure]], 0.02), stats::quantile(d[[exposure]], 0.98), length.out = 200
  )
  grid <- base_row[rep(1, length(grid_values)), , drop = FALSE]
  grid[[exposure]] <- grid_values
  terms_no_response <- stats::delete.response(stats::terms(fit_spline))
  x_grid <- stats::model.matrix(terms_no_response, grid)
  x_reference <- stats::model.matrix(terms_no_response, base_row)
  coefficients <- stats::coef(fit_spline)
  x_grid <- x_grid[, names(coefficients), drop = FALSE]
  x_reference <- x_reference[, names(coefficients), drop = FALSE]
  contrast <- sweep(x_grid, 2, x_reference[1, ], FUN = "-")
  variance <- stats::vcov(fit_spline)
  log_hr <- as.numeric(contrast %*% coefficients)
  se <- sqrt(rowSums((contrast %*% variance) * contrast))
  predictions <- data.frame(
    exposure = exposure, value = grid_values,
    HR = exp(log_hr), lower = exp(log_hr - 1.96 * se), upper = exp(log_hr + 1.96 * se)
  )

  annotation <- sprintf(
    "Linear HR %.2f (95%% CI %.2f–%.2f)\nP overall %s\nP nonlinear %s",
    linear_row$estimate, linear_row$conf.low, linear_row$conf.high,
    ifelse(p_overall < 0.001, "< .001", paste0("= ", format_p_rsna(p_overall))),
    ifelse(p_nonlinear < 0.001, "< .001", paste0("= ", format_p_rsna(p_nonlinear)))
  )
  plot <- ggplot2::ggplot(predictions, ggplot2::aes(value, HR)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#AFC3DF", alpha = 0.65) +
    ggplot2::geom_line(color = "#2878A8", linewidth = 1.1) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed") +
    ggplot2::annotate("text", x = Inf, y = Inf, label = annotation, hjust = 1.05, vjust = 1.2, size = 3) +
    ggplot2::labs(x = label, y = "Hazard ratio for T2EL") +
    theme_rsna(10)

  list(
    plot = plot, predictions = predictions,
    statistics = data.frame(
      Exposure = exposure, N = nrow(d), Events = sum(d[[event_col]]),
      Linear_HR = linear_row$estimate, Linear_CI_lower = linear_row$conf.low,
      Linear_CI_upper = linear_row$conf.high, P_overall = p_overall,
      P_nonlinear = p_nonlinear, Reference = threshold
    )
  )
}

analyses <- Map(draw_rcs, names(exposures), unname(exposures))
names(analyses) <- names(exposures)
combined <- patchwork::wrap_plots(lapply(analyses, `[[`, "plot"), ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A")
dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
save_publication_plot(combined, "figure6_rcs", io$output_dir, 10, 8)
statistics <- dplyr::bind_rows(lapply(analyses, `[[`, "statistics"))
predictions <- dplyr::bind_rows(lapply(analyses, `[[`, "predictions"))
utils::write.csv(statistics, file.path(io$output_dir, "figure6_rcs_statistics.csv"), row.names = FALSE)
utils::write.csv(predictions, file.path(io$output_dir, "figure6_rcs_predictions.csv"), row.names = FALSE)
print(statistics, row.names = FALSE)
