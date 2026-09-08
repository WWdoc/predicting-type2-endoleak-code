## Analysis: Figure 5 ROC curves for four anatomical markers
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "ggplot2", "pROC"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
markers <- c(
  ostial_area_mm2 = "Total branch-vessel ostial area",
  ima_diameter_mm = "IMA diameter",
  patent_la_number = "Number of patent lumbar arteries",
  la_ge2_number = "Number of lumbar arteries ≥2 mm"
)
require_columns(data, c("t2el", names(markers)), "Figure 5 ROC analysis")
data$t2el <- as_binary(data$t2el)

roc_objects <- list()
curve_rows <- list()
result_rows <- list()
for (marker in names(markers)) {
  d <- data[stats::complete.cases(data[, c("t2el", marker)]), c("t2el", marker), drop = FALSE]
  if (length(unique(d$t2el)) != 2) stop("Both outcome classes are required for ", marker)
  roc <- pROC::roc(d$t2el, d[[marker]], levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_objects[[marker]] <- roc
  ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  best <- pROC::coords(
    roc, x = "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
  )
  coords <- data.frame(
    false_positive_rate = 1 - roc$specificities,
    sensitivity = roc$sensitivities,
    marker = unname(markers[marker])
  )
  curve_rows[[marker]] <- coords
  result_rows[[marker]] <- data.frame(
    Marker = unname(markers[marker]), N = nrow(d), Events = sum(d$t2el),
    AUC = as.numeric(pROC::auc(roc)), CI_lower = ci[1], CI_upper = ci[3],
    Youden_threshold = as.numeric(best$threshold),
    Sensitivity = as.numeric(best$sensitivity), Specificity = as.numeric(best$specificity)
  )
}
results <- dplyr::bind_rows(result_rows)
curves <- dplyr::bind_rows(curve_rows)

raw_p <- vapply(setdiff(names(roc_objects), "ostial_area_mm2"), function(marker) {
  keep <- stats::complete.cases(data[, c("t2el", "ostial_area_mm2", marker)])
  reference <- pROC::roc(
    data$t2el[keep], data$ostial_area_mm2[keep],
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  comparator <- pROC::roc(
    data$t2el[keep], data[[marker]][keep],
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  pROC::roc.test(reference, comparator, paired = TRUE, method = "delong")$p.value
}, numeric(1))
adjusted_p <- stats::p.adjust(raw_p, method = "holm")
comparisons <- data.frame(
  Comparator = unname(markers[names(raw_p)]),
  Reference = unname(markers["ostial_area_mm2"]),
  P_value_DeLong = raw_p,
  P_value_Holm = adjusted_p,
  stringsAsFactors = FALSE
)

legend_labels <- setNames(
  sprintf("%s: AUC %.2f (95%% CI %.2f–%.2f)", results$Marker, results$AUC, results$CI_lower, results$CI_upper),
  results$Marker
)
figure5 <- ggplot2::ggplot(curves, ggplot2::aes(false_positive_rate, sensitivity, color = marker)) +
  ggplot2::geom_step(linewidth = 1) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::scale_color_manual(
    values = c("#4C72B0", "#DD8452", "#55A868", "#C44E52"), labels = legend_labels
  ) +
  ggplot2::labs(x = "1 − Specificity", y = "Sensitivity", color = NULL) +
  theme_rsna(10) +
  ggplot2::theme(legend.position = c(0.62, 0.24), legend.text = ggplot2::element_text(size = 7))

dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)
save_publication_plot(figure5, "figure5_roc_curves", io$output_dir, 6.5, 6.2)
utils::write.csv(results, file.path(io$output_dir, "figure5_roc_statistics.csv"), row.names = FALSE)
utils::write.csv(comparisons, file.path(io$output_dir, "figure5_paired_delong_comparisons.csv"), row.names = FALSE)
print(results, row.names = FALSE)
print(comparisons, row.names = FALSE)
