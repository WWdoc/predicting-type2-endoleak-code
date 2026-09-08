## Analysis: run all public manuscript workflows
## Date: 2026-09-08
## Random seed: 42

set.seed(42)

scripts <- sprintf("R/%02d_%s.R", 1:16, c(
  "table1_baseline", "table2_procedural", "table3_outcomes", "table4_anatomy",
  "table5_cox_models", "figures3_4_survival", "figure5_roc", "figure6_rcs",
  "supp_figure1_agreement", "supp_figure2_pairplot", "supp_figure3_cox_lasso",
  "supp_figure4_bootstrap", "supp_figure5_distributions", "supp_table1_agreement",
  "supp_table2_ph_tests", "supp_table3_diagnostics"
))

for (script in scripts) {
  message("\n===== Running ", script, " =====")
  status <- system2("Rscript", script)
  if (!identical(status, 0L)) stop("Analysis failed: ", script)
}

message("All analyses completed.")
