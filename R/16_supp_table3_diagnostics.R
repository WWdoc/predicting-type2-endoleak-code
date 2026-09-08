## Analysis: Supplementary Table 3 diagnostic performance of 17 methods
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "flextable", "pROC", "tidyr"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
required <- c(
  "t2el", "ostial_area_mm2", "ima_diameter_mm", "patent_la_number",
  "la_ge2_number", "any_la_ge2", "thrombus_area_mm2", "thrombus_fraction"
)
require_columns(data, required, "Supplementary Table 3")

d <- data.frame(
  outcome = as_binary(data$t2el),
  ostial_area = as.numeric(data$ostial_area_mm2),
  ima = as.numeric(data$ima_diameter_mm),
  patent_la = as.numeric(data$patent_la_number),
  la_ge2_count = as.numeric(data$la_ge2_number),
  any_la_ge2 = as_binary(data$any_la_ge2),
  thrombus_area = as.numeric(data$thrombus_area_mm2),
  thrombus_fraction = as.numeric(data$thrombus_fraction)
)

methods <- data.frame(
  Method = sprintf("Method %d", 1:17),
  Item = c(
    "Total branch-vessel ostial area (continuous)",
    "Guideline risk: IMA ≥3 mm OR >3 patent LAs OR any LA ≥2 mm",
    "IMA diameter (continuous)", "IMA diameter ≥3 mm", "IMA diameter ≥2.5 mm",
    "Any lumbar artery ≥2 mm", "Number of patent lumbar arteries (continuous)",
    "Number of lumbar arteries ≥2 mm (continuous)", ">3 patent lumbar arteries",
    ">3 patent LAs OR any LA ≥2 mm", ">3 patent LAs AND any LA ≥2 mm",
    "IMA ≥3 mm OR >3 patent LAs", "IMA ≥3 mm AND >3 patent LAs",
    "IMA ≥3 mm OR any LA ≥2 mm", "IMA ≥3 mm AND any LA ≥2 mm",
    "Thrombus cross-sectional area (continuous)", "Thrombus proportion (continuous)"
  ),
  stringsAsFactors = FALSE
)

scores <- data.frame(
  Method_1 = d$ostial_area,
  Method_2 = as.integer(d$ima >= 3 | d$patent_la > 3 | d$any_la_ge2 == 1),
  Method_3 = d$ima,
  Method_4 = as.integer(d$ima >= 3),
  Method_5 = as.integer(d$ima >= 2.5),
  Method_6 = d$any_la_ge2,
  Method_7 = d$patent_la,
  Method_8 = d$la_ge2_count,
  Method_9 = as.integer(d$patent_la > 3),
  Method_10 = as.integer(d$patent_la > 3 | d$any_la_ge2 == 1),
  Method_11 = as.integer(d$patent_la > 3 & d$any_la_ge2 == 1),
  Method_12 = as.integer(d$ima >= 3 | d$patent_la > 3),
  Method_13 = as.integer(d$ima >= 3 & d$patent_la > 3),
  Method_14 = as.integer(d$ima >= 3 | d$any_la_ge2 == 1),
  Method_15 = as.integer(d$ima >= 3 & d$any_la_ge2 == 1),
  Method_16 = d$thrombus_area,
  Method_17 = d$thrombus_fraction,
  check.names = FALSE
)

analyse_method <- function(index) {
  score <- scores[[index]]
  keep <- stats::complete.cases(d$outcome, score)
  outcome <- d$outcome[keep]
  score <- score[keep]
  roc <- pROC::roc(outcome, score, levels = c(0, 1), direction = "auto", quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  best <- pROC::coords(
    roc, "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
  )[1, , drop = FALSE]
  threshold <- as.numeric(best$threshold)
  prediction <- if (roc$direction == "<") as.integer(score >= threshold) else as.integer(score <= threshold)
  tp <- sum(prediction == 1 & outcome == 1)
  fn <- sum(prediction == 0 & outcome == 1)
  tn <- sum(prediction == 0 & outcome == 0)
  fp <- sum(prediction == 1 & outcome == 0)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  sensitivity_ci <- wilson_ci(tp, tp + fn)
  specificity_ci <- wilson_ci(tn, tn + fp)
  list(
    roc = roc,
    result = data.frame(
      Method = methods$Method[index], Item = methods$Item[index], N = length(outcome),
      Events = sum(outcome), AUC = as.numeric(pROC::auc(roc)),
      AUC_CI_lower = auc_ci[1], AUC_CI_upper = auc_ci[3], Direction = roc$direction,
      Threshold = threshold, TP = tp, FN = fn, TN = tn, FP = fp,
      Sensitivity = sensitivity, Sensitivity_CI_lower = sensitivity_ci[1],
      Sensitivity_CI_upper = sensitivity_ci[2], Specificity = specificity,
      Specificity_CI_lower = specificity_ci[1], Specificity_CI_upper = specificity_ci[2]
    )
  )
}

analyses <- lapply(seq_len(nrow(methods)), analyse_method)
detailed <- dplyr::bind_rows(lapply(analyses, `[[`, "result"))

reference_p <- rep(NA_real_, nrow(methods))
for (index in 2:nrow(methods)) {
  score_reference <- scores[[1]]
  score_comparator <- scores[[index]]
  keep <- stats::complete.cases(d$outcome, score_reference, score_comparator)
  roc_reference <- pROC::roc(d$outcome[keep], score_reference[keep], levels = c(0, 1), direction = "auto", quiet = TRUE)
  roc_comparator <- pROC::roc(d$outcome[keep], score_comparator[keep], levels = c(0, 1), direction = "auto", quiet = TRUE)
  reference_p[index] <- pROC::roc.test(
    roc_reference, roc_comparator, method = "delong", paired = TRUE
  )$p.value
}
detailed$P_value_vs_ostial_area <- reference_p
detailed$P_value_vs_ostial_area_Holm <- c(NA_real_, stats::p.adjust(reference_p[-1], method = "holm"))

publication <- detailed |>
  dplyr::transmute(
    Method, Item,
    `AUC (95% CI)` = sprintf("%.2f (%.2f, %.2f)", AUC, AUC_CI_lower, AUC_CI_upper),
    `Sensitivity (95% CI)` = sprintf(
      "%.1f%% (%d/%d; %.1f%%, %.1f%%)", 100 * Sensitivity, TP, TP + FN,
      100 * Sensitivity_CI_lower, 100 * Sensitivity_CI_upper
    ),
    `Specificity (95% CI)` = sprintf(
      "%.1f%% (%d/%d; %.1f%%, %.1f%%)", 100 * Specificity, TN, TN + FP,
      100 * Specificity_CI_lower, 100 * Specificity_CI_upper
    ),
    `P value vs total ostial area (Holm-adjusted)` = ifelse(
      is.na(P_value_vs_ostial_area_Holm), "Reference", format_p_rsna(P_value_vs_ostial_area_Holm)
    )
  )

dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(detailed, file.path(io$output_dir, "supp_table3_diagnostic_metrics_detailed.csv"), row.names = FALSE)
utils::write.csv(publication, file.path(io$output_dir, "supp_table3_diagnostic_metrics.csv"), row.names = FALSE)
ft <- flextable::flextable(publication) |>
  flextable::theme_booktabs() |>
  flextable::autofit()
flextable::save_as_docx(`Supplementary Table 3` = ft,
                        path = file.path(io$output_dir, "supp_table3_diagnostic_metrics.docx"))
print(publication, row.names = FALSE)
