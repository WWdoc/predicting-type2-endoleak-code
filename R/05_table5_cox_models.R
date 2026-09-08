## Analysis: Table 5 sequential Cox proportional-hazards models
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("broom", "dplyr", "flextable", "survival"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
time_col <- "t2el_followup_years"
event_col <- "t2el"
exposure <- "ostial_area_mm2"
required <- unique(c(time_col, event_col, exposure, unlist(cfg$cox_models)))
require_columns(data, required, "Table 5 Cox models")
data[[event_col]] <- as_binary(data[[event_col]])

fit_one <- function(model_name, covariates) {
  columns <- unique(c(time_col, event_col, exposure, covariates))
  analysis_data <- data[stats::complete.cases(data[, columns, drop = FALSE]), columns, drop = FALSE]
  analysis_data <- analysis_data[analysis_data[[time_col]] > 0, , drop = FALSE]
  events <- sum(analysis_data[[event_col]] == 1)
  n_parameters <- length(c(exposure, covariates))
  if (events == 0) stop("No T2EL events are available for ", model_name)
  if (events / n_parameters < 10) {
    warning(model_name, " has fewer than 10 events per candidate parameter; interpret Wald CIs cautiously.")
  }

  formula <- make_survival_formula(time_col, event_col, c(exposure, covariates))
  fit <- survival::coxph(formula, data = analysis_data, x = TRUE, y = TRUE)
  tidy <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  exposure_row <- tidy[tidy$term == exposure, , drop = FALSE]
  if (nrow(exposure_row) != 1) stop("Could not uniquely identify the ostial-area coefficient.")
  ph <- survival::cox.zph(fit, transform = "km")
  ph_table <- as.data.frame(ph$table)
  exposure_ph_row <- grep(exposure, rownames(ph_table), fixed = TRUE)
  ph_exposure <- if (length(exposure_ph_row)) ph_table[exposure_ph_row[1], "p"] else NA_real_
  ph_global <- ph_table["GLOBAL", "p"]
  concordance <- summary(fit)$concordance

  list(
    fit = fit, analysis_data = analysis_data, ph = ph,
    result = data.frame(
      Model = model_name,
      Adjustment = unname(cfg$model_adjustment_labels[model_name]),
      N = nrow(analysis_data), Events = events,
      HR = exposure_row$estimate, CI_lower = exposure_row$conf.low,
      CI_upper = exposure_row$conf.high, P_value = exposure_row$p.value,
      C_index = unname(concordance[1]), C_index_SE = unname(concordance[2]),
      Ostial_area_PH_P = ph_exposure, Global_PH_P = ph_global,
      stringsAsFactors = FALSE
    )
  )
}

fits <- Map(fit_one, names(cfg$cox_models), cfg$cox_models)
names(fits) <- names(cfg$cox_models)
results <- dplyr::bind_rows(lapply(fits, `[[`, "result")) |>
  dplyr::mutate(
    `HR per 1-mm² increase (95% CI)` = sprintf("%.2f (%.2f, %.2f)", HR, CI_lower, CI_upper),
    `P value` = format_p_rsna(P_value),
    `C-index (95% CI)` = sprintf(
      "%.3f (%.3f, %.3f)", C_index,
      pmax(0, C_index - 1.96 * C_index_SE), pmin(1, C_index + 1.96 * C_index_SE)
    )
  )

dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, file.path(io$output_dir, "table5_cox_models_detailed.csv"), row.names = FALSE)
publication <- results[, c("Model", "Adjustment", "N", "Events", "HR per 1-mm² increase (95% CI)", "P value")]
utils::write.csv(publication, file.path(io$output_dir, "table5_cox_models.csv"), row.names = FALSE)
ft <- flextable::flextable(publication) |>
  flextable::theme_booktabs() |>
  flextable::autofit()
flextable::save_as_docx(`Table 5` = ft, path = file.path(io$output_dir, "table5_cox_models.docx"))

reverse_km <- survival::survfit(
  survival::Surv(data[[time_col]], 1 - data[[event_col]]) ~ 1,
  data = data[stats::complete.cases(data[, c(time_col, event_col)]), , drop = FALSE]
)
followup_summary <- summary(reverse_km)$table
writeLines(
  c(
    "Table 5 reproducibility summary",
    paste0("Median follow-up by reverse Kaplan-Meier: ", format(followup_summary[["median"]], digits = 4), " years"),
    "Cox proportional-hazards assumptions were evaluated with Schoenfeld residuals.",
    "Model-specific fitted N and event counts are included in the detailed CSV."
  ),
  file.path(io$output_dir, "table5_methods_notes.txt")
)

print(publication, row.names = FALSE)
