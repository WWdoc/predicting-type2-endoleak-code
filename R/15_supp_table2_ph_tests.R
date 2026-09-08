## Analysis: Supplementary Table 2 proportional-hazards assumption tests
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "flextable", "survival"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
time_col <- "t2el_followup_years"
event_col <- "t2el"
exposure <- "ostial_area_mm2"
require_columns(data, unique(c(time_col, event_col, exposure, unlist(cfg$cox_models))), "Supplementary Table 2")
data[[event_col]] <- as_binary(data[[event_col]])

rows <- Map(function(model_name, covariates) {
  columns <- unique(c(time_col, event_col, exposure, covariates))
  d <- data[stats::complete.cases(data[, columns, drop = FALSE]) & data[[time_col]] > 0, columns, drop = FALSE]
  fit <- survival::coxph(make_survival_formula(time_col, event_col, c(exposure, covariates)), data = d, x = TRUE)
  ph <- as.data.frame(survival::cox.zph(fit, transform = "km")$table)
  exposure_row <- grep(exposure, rownames(ph), fixed = TRUE)
  data.frame(
    `Cox model` = model_name, N = nrow(d), Events = sum(d[[event_col]]),
    `Total branch-vessel ostial area, P` = if (length(exposure_row)) ph[exposure_row[1], "p"] else NA_real_,
    `Global test, P` = ph["GLOBAL", "p"], check.names = FALSE
  )
}, names(cfg$cox_models), cfg$cox_models)
results <- dplyr::bind_rows(rows)
publication <- results |>
  dplyr::mutate(
    `Total branch-vessel ostial area, P` = format_p_rsna(`Total branch-vessel ostial area, P`),
    `Global test, P` = format_p_rsna(`Global test, P`)
  )
dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, file.path(io$output_dir, "supp_table2_ph_tests_detailed.csv"), row.names = FALSE)
utils::write.csv(publication, file.path(io$output_dir, "supp_table2_ph_tests.csv"), row.names = FALSE)
ft <- flextable::flextable(publication) |> flextable::theme_booktabs() |> flextable::autofit()
flextable::save_as_docx(`Supplementary Table 2` = ft,
                        path = file.path(io$output_dir, "supp_table2_ph_tests.docx"))
print(publication, row.names = FALSE)
