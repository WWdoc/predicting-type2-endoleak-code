# Standardized analysis data contract

The scripts expect one de-identified patient per row. Column names below are public analysis aliases, not a disclosure of any original database schema.

## Core outcomes and grouping

| Column | Type | Definition |
|---|---|---|
| `t2el` | 0/1 | Post-EVAR type II endoleak indicator |
| `t2el_followup_years` | numeric | Time from EVAR to T2EL or censoring |
| `time_to_detected_t2el_years` | numeric | Time to first CTA-detected T2EL; leave missing for patients without T2EL |
| `followup_years` | numeric | Time to death or last clinical follow-up |
| `all_cause_death` | 0/1 | All-cause mortality event |
| `reintervention_followup_years` | numeric | Time to any reintervention or censoring |
| `any_reintervention` | 0/1 | Any reintervention event |
| `t2el_reintervention_followup_years` | numeric | Time to T2EL-related reintervention or censoring |
| `t2el_related_reintervention` | 0/1 | T2EL-related reintervention event |
| `sac_enlargement` | 0/1 | Aneurysm sac volume increase of at least 10% |

## Baseline and procedural variables

Baseline aliases: `age_years`, `male`, `bmi`, `smoking`, `hypertension`, `diabetes`, `antidiabetic_drugs`, `hypercholesterolemia`, `lipid_lowering_drugs`, `cerebrovascular_disease`, `coronary_artery_disease`, `renal_cysts`, `concomitant_iliac_aneurysm`, `preoperative_aortic_repair`, `symptom_or_rupture`, `anticoagulants`, `antiplatelet_agents`, and `preoperative_ef`.

Procedural aliases: `asa_class`, `anesthesia`, `stent_graft`, `main_body_diameter_mm`, `femoral_approach`, `iia_embolization`, `iia_preservation`, `technical_success`, `operation_time_min`, `icu_stay_hours`, `ventilation_time_hours`, and `hospital_stay_days`.

Clinical outcome aliases: `aneurysm_rupture`, `mortality_30d`, `reintervention_30d`, `atrial_fibrillation_30d`, `acute_kidney_injury_30d`, `acute_limb_ischemia_30d`, `access_site_complication_30d`, and `stent_graft_stenosis`.

## Anatomical variables

| Column | Unit/encoding |
|---|---|
| `neck_length_mm` | mm |
| `infrarenal_angle_deg` | degrees |
| `neck_short_axis_mm` | mm |
| `neck_long_axis_mm` | mm |
| `aaa_max_diameter_mm` | mm |
| `aaa_volume_ml` | mL |
| `ostial_area_mm2` | mm2, total patent branch-vessel ostial area |
| `ima_diameter_mm` | mm |
| `patent_la_number` | count |
| `la_ge2_number` | count of lumbar arteries with diameter at least 2 mm |
| `any_la_ge2` | 0/1 |
| `thrombus_area_mm2` | mm2 |
| `thrombus_fraction` | proportion from 0 to 1 |

## Observer-variability input

`operator_variability.csv` contains only four numeric columns: `W1`, `W2`, `Z1`, and `Z2`, representing two repeated measurements by each of two observers. No identifier is required by the analysis.

Missing values should be coded as blank/`NA`. Binary aliases must use 0 for absence and 1 for presence. Time variables must be positive and use years.
