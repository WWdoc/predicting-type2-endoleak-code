## Analysis: Table 1 baseline demographics and clinical characteristics
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "gtsummary", "flextable", "knitr", "rlang"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)

mean_vars <- c("age_years", "bmi")
median_vars <- c("preoperative_ef", "followup_years")
binary_vars <- c(
  "male", "smoking", "hypertension", "diabetes", "antidiabetic_drugs",
  "hypercholesterolemia", "lipid_lowering_drugs", "cerebrovascular_disease",
  "coronary_artery_disease", "renal_cysts", "concomitant_iliac_aneurysm",
  "preoperative_aortic_repair", "symptom_or_rupture", "anticoagulants",
  "antiplatelet_agents"
)
required <- c("t2el", mean_vars, median_vars, binary_vars)
require_columns(data, required, "Table 1")

data <- prepare_binary_factors(data, binary_vars)
data <- prepare_t2el_group(data)

labels <- c(
  age_years = "Age, years", male = "Male sex", bmi = "Body mass index, kg/m²",
  smoking = "Smoking", hypertension = "Hypertension", diabetes = "Diabetes",
  antidiabetic_drugs = "Antidiabetic drugs",
  hypercholesterolemia = "Hypercholesterolemia",
  lipid_lowering_drugs = "Lipid-lowering drugs",
  cerebrovascular_disease = "Cerebrovascular disease",
  coronary_artery_disease = "Coronary artery disease", renal_cysts = "Renal cysts",
  concomitant_iliac_aneurysm = "Concomitant iliac aneurysm",
  preoperative_aortic_repair = "Preoperative aortic repair",
  symptom_or_rupture = "Symptom/rupture", anticoagulants = "Anticoagulants",
  antiplatelet_agents = "Antiplatelet agents",
  preoperative_ef = "Preoperative ejection fraction",
  followup_years = "Follow-up, years"
)

table1 <- data |>
  dplyr::select(dplyr::all_of(required)) |>
  gtsummary::tbl_summary(
    by = "t2el",
    label = label_formulas(labels),
    type = c(formula_list(binary_vars, "dichotomous")),
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
      formula_list(mean_vars, "t.test"),
      formula_list(median_vars, "wilcox.test"),
      formula_list(binary_vars, "chisq.test.no.correct")
    ),
    pvalue_fun = format_p_rsna
  ) |>
  gtsummary::bold_labels() |>
  gtsummary::modify_caption("**Table 1. Baseline Demographics and Clinical Characteristics**")

export_gtsummary(table1, "table1_baseline", io$output_dir, "Table 1")
print(table1)
