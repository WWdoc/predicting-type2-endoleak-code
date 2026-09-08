## Analysis: Table 4 aneurysm anatomy measurements
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "gtsummary", "flextable", "knitr", "rlang"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)

mean_vars <- c("neck_short_axis_mm", "neck_long_axis_mm")
median_vars <- c(
  "neck_length_mm", "infrarenal_angle_deg", "aaa_max_diameter_mm", "aaa_volume_ml",
  "ostial_area_mm2", "ima_diameter_mm", "patent_la_number", "la_ge2_number",
  "thrombus_area_mm2", "thrombus_fraction"
)
binary_vars <- c("ima_ge3", "ima_ge2_5", "any_la_ge2")
required_base <- c("t2el", mean_vars, median_vars, "any_la_ge2")
require_columns(data, required_base, "Table 4")
data$ima_ge3 <- as.integer(data$ima_diameter_mm >= 3)
data$ima_ge2_5 <- as.integer(data$ima_diameter_mm >= 2.5)
data <- prepare_binary_factors(data, binary_vars)
data <- prepare_t2el_group(data)

labels <- c(
  neck_length_mm = "Neck length, mm", infrarenal_angle_deg = "Infrarenal angle, degrees",
  neck_short_axis_mm = "Neck short-axis diameter, mm",
  neck_long_axis_mm = "Neck long-axis diameter, mm",
  aaa_max_diameter_mm = "AAA maximum diameter, mm", aaa_volume_ml = "AAA volume, mL",
  ostial_area_mm2 = "Total branch-vessel ostial area, mm²",
  ima_diameter_mm = "IMA diameter, mm", ima_ge3 = "IMA diameter ≥3 mm",
  ima_ge2_5 = "IMA diameter ≥2.5 mm", patent_la_number = "Number of patent lumbar arteries",
  any_la_ge2 = "At least one lumbar artery ≥2 mm",
  la_ge2_number = "Number of lumbar arteries ≥2 mm",
  thrombus_area_mm2 = "Thrombus cross-sectional area, mm²",
  thrombus_fraction = "Thrombus proportion"
)

table4 <- data |>
  dplyr::select(dplyr::all_of(c("t2el", mean_vars, median_vars, binary_vars))) |>
  gtsummary::tbl_summary(
    by = "t2el", label = label_formulas(labels),
    type = formula_list(binary_vars, "dichotomous"),
    value = formula_list(binary_vars, "Yes"),
    statistic = c(
      formula_list(mean_vars, "{mean} ± {sd}"),
      formula_list(median_vars, "{median} ({p25}, {p75})"),
      formula_list(binary_vars, "{n} ({p}%)")
    ),
    digits = list(gtsummary::all_continuous() ~ 1, gtsummary::all_categorical() ~ c(0, 1)),
    missing = "ifany", missing_text = "Missing"
  ) |>
  gtsummary::add_overall(last = FALSE) |>
  gtsummary::add_p(
    test = c(
      formula_list(mean_vars, "t.test"), formula_list(median_vars, "wilcox.test"),
      formula_list(binary_vars, "chisq.test.no.correct")
    ),
    pvalue_fun = format_p_rsna
  ) |>
  gtsummary::bold_labels() |>
  gtsummary::modify_caption("**Table 4. Aneurysm Anatomy Measurements**")

export_gtsummary(table4, "table4_anatomy", io$output_dir, "Table 4")
print(table4)
