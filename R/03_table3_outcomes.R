## Analysis: Table 3 clinical outcomes
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "gtsummary", "flextable", "knitr", "rlang"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)

time_var <- "time_to_detected_t2el_years"
binary_vars <- c(
  "aneurysm_rupture", "mortality_30d", "all_cause_death", "reintervention_30d",
  "any_reintervention", "atrial_fibrillation_30d", "acute_kidney_injury_30d",
  "acute_limb_ischemia_30d", "access_site_complication_30d", "stent_graft_stenosis"
)
required <- c("t2el", time_var, binary_vars)
require_columns(data, required, "Table 3")
data <- prepare_binary_factors(data, binary_vars)
data <- prepare_t2el_group(data)

labels <- c(
  time_to_detected_t2el_years = "Time from EVAR to first CTA-detected T2EL, years",
  aneurysm_rupture = "Aneurysm rupture", mortality_30d = "30-day mortality",
  all_cause_death = "All-cause mortality", reintervention_30d = "30-day reintervention",
  any_reintervention = "Any reintervention",
  atrial_fibrillation_30d = "Atrial fibrillation (30-day)",
  acute_kidney_injury_30d = "Acute kidney injury (30-day)",
  acute_limb_ischemia_30d = "Acute limb ischemia (30-day)",
  access_site_complication_30d = "Access-site complication (30-day)",
  stent_graft_stenosis = "Stent-graft stenosis"
)

table3 <- data |>
  dplyr::select(dplyr::all_of(required)) |>
  gtsummary::tbl_summary(
    by = "t2el", label = label_formulas(labels),
    type = formula_list(binary_vars, "dichotomous"),
    value = formula_list(binary_vars, "Yes"),
    statistic = c(
      formula_list(time_var, "{median} ({p25}, {p75})"),
      formula_list(binary_vars, "{n} ({p}%)")
    ),
    digits = list(gtsummary::all_continuous() ~ 2, gtsummary::all_categorical() ~ c(0, 1)),
    missing = "ifany", missing_text = "Missing"
  ) |>
  gtsummary::add_overall(last = FALSE) |>
  gtsummary::add_p(
    include = -dplyr::all_of(time_var),
    test = formula_list(binary_vars, "chisq.test.no.correct"),
    pvalue_fun = format_p_rsna
  ) |>
  gtsummary::bold_labels() |>
  gtsummary::modify_caption("**Table 3. Clinical Outcomes by T2EL Status**")

export_gtsummary(table3, "table3_clinical_outcomes", io$output_dir, "Table 3")
print(table3)
