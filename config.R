## Public analysis configuration. No patient-level data or private paths.

cfg <- list(
  seed = 42,
  input_file = "data/analysis_data.csv",
  operator_file = "data/operator_variability.csv",
  output_dir = "outputs",
  bootstrap_repetitions = 1000L,
  cv_folds = 10L,
  risk_times = c(0, 3, 6, 9, 12, 15),
  thresholds = c(
    ostial_area_mm2 = 15.55,
    ima_diameter_mm = 2.4,
    patent_la_number = 3,
    la_ge2_number = 2
  )
)

cfg$baseline_covariates <- c(
  "age_years", "male", "bmi", "smoking", "hypertension", "diabetes"
)

cfg$nonbranch_anatomical_covariates <- c(
  "concomitant_iliac_aneurysm", "neck_length_mm", "neck_short_axis_mm",
  "neck_long_axis_mm", "infrarenal_angle_deg", "aaa_max_diameter_mm",
  "thrombus_area_mm2", "thrombus_fraction", "aaa_volume_ml"
)

cfg$comprehensive_covariates <- unique(c(
  cfg$baseline_covariates,
  "antidiabetic_drugs", "hypercholesterolemia", "lipid_lowering_drugs",
  "cerebrovascular_disease", "coronary_artery_disease", "renal_cysts",
  "concomitant_iliac_aneurysm", "preoperative_aortic_repair",
  "symptom_or_rupture", "anticoagulants", "antiplatelet_agents",
  "preoperative_ef", cfg$nonbranch_anatomical_covariates
))

cfg$cox_models <- list(
  `Model 1` = character(0),
  `Model 2` = cfg$baseline_covariates,
  `Model 3` = cfg$nonbranch_anatomical_covariates,
  `Model 4` = unique(c(cfg$baseline_covariates, cfg$nonbranch_anatomical_covariates)),
  `Model 5` = cfg$comprehensive_covariates
)

cfg$model_adjustment_labels <- c(
  `Model 1` = "Unadjusted",
  `Model 2` = "Clinically relevant baseline covariates",
  `Model 3` = "Non-branch anatomical covariates",
  `Model 4` = "Baseline plus non-branch anatomical covariates",
  `Model 5` = "Comprehensive preoperative covariates"
)
