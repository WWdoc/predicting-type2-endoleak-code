# Expected analysis outputs

When `run_all.R` is executed with de-identified local inputs, the `outputs/` directory receives:

- Tables 1–5 as CSV; Tables 1–5 as editable Word where applicable.
- Figures 3–6 as PDF and 300-dpi PNG, with numerical statistics in adjacent CSV files.
- Supplementary Figures 1–5 as PDF and 300-dpi PNG.
- Supplementary Tables 1–3 as CSV and editable Word.
- Detailed model, bootstrap, correlation, ROC, and proportional-hazards results as CSV.

Patient-level input files are never copied into `outputs/`. Some detailed derived files, such as bootstrap replicates and spline prediction grids, are ignored by Git because the entire `outputs/` directory is excluded by `.gitignore`.
