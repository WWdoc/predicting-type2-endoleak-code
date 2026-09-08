## Analysis: Supplementary Table 1 observer reproducibility
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("flextable", "irr"))

io <- resolve_io(cfg$operator_file)
data <- read_analysis_data(io$input_file)
results <- agreement_results_frame(agreement_comparisons(data))

publication <- results[, c(
  "Comparison", "N", "Measurement_1_mean_SD", "Measurement_2_mean_SD",
  "Bias_95CI", "LoA", "ICC_95CI", "P_value_bias_display", "P_value_bias_Holm_display"
)]
names(publication) <- c(
  "Comparison", "N", "Measurement 1, mean ± SD", "Measurement 2, mean ± SD",
  "Bias (95% CI)", "95% limits of agreement", "ICC (95% CI)",
  "P value for bias", "Holm-adjusted P value"
)
dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, file.path(io$output_dir, "supp_table1_agreement_detailed.csv"), row.names = FALSE)
utils::write.csv(publication, file.path(io$output_dir, "supp_table1_agreement.csv"), row.names = FALSE)
ft <- flextable::flextable(publication) |>
  flextable::theme_booktabs() |>
  flextable::autofit()
flextable::save_as_docx(`Supplementary Table 1` = ft,
                        path = file.path(io$output_dir, "supp_table1_agreement.docx"))
print(publication, row.names = FALSE)
