## Analysis: Table 2 procedural and hospitalization details
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "gtsummary", "flextable", "knitr", "rlang"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)

median_vars <- c(
  "main_body_diameter_mm", "operation_time_min", "icu_stay_hours",
  "ventilation_time_hours", "hospital_stay_days"
)
categorical_vars <- c("asa_class", "anesthesia", "stent_graft", "femoral_approach", "iia_preservation")
binary_vars <- c("iia_embolization", "technical_success")
required <- c("t2el", categorical_vars, binary_vars, median_vars)
require_columns(data, required, "Table 2")

data <- prepare_binary_factors(data, binary_vars)
data <- prepare_t2el_group(data)
labels <- c(
  asa_class = "ASA classification", anesthesia = "Anesthesia",
  stent_graft = "Type of stent-graft",
  main_body_diameter_mm = "Main body proximal diameter, mm",
  femoral_approach = "Femoral artery approach",
  iia_embolization = "Embolization of the internal iliac arteries",
  iia_preservation = "Internal iliac artery preservation",
  technical_success = "Technical success", operation_time_min = "Operation time, min",
  icu_stay_hours = "ICU stay, hours", ventilation_time_hours = "Mechanical ventilation, hours",
  hospital_stay_days = "Hospital stay, days"
)

table2 <- data |>
  dplyr::select(dplyr::all_of(required)) |>
  gtsummary::tbl_summary(
    by = "t2el", label = label_formulas(labels),
    type = c(formula_list(binary_vars, "dichotomous")),
    value = formula_list(binary_vars, "Yes"),
    statistic = c(
      formula_list(median_vars, "{median} ({p25}, {p75})"),
      formula_list(c(categorical_vars, binary_vars), "{n} ({p}%)")
    ),
    digits = list(gtsummary::all_continuous() ~ 1, gtsummary::all_categorical() ~ c(0, 1)),
    missing = "ifany", missing_text = "Missing"
  ) |>
  gtsummary::add_overall(last = FALSE) |>
  gtsummary::add_p(
    test = c(
      formula_list(median_vars, "wilcox.test"),
      formula_list(categorical_vars, "chisq.test.no.correct"),
      formula_list(binary_vars, "chisq.test.no.correct")
    ),
    pvalue_fun = format_p_rsna
  ) |>
  gtsummary::bold_labels() |>
  gtsummary::modify_caption("**Table 2. Procedural and Hospitalization Details**")

export_gtsummary(table2, "table2_procedural", io$output_dir, "Table 2")
print(table2)
