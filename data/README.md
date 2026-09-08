# Local data folder

Do not commit patient-level data to the public repository.

Expected local files:

- `analysis_data.csv`: one row per patient, de-identified.
- `operator_variability.csv`: one row per measured patient with columns `W1`, `W2`, `Z1`, and `Z2`.

See `../DATA_DICTIONARY.md` for the standardized analysis columns. This folder is ignored by Git except for this file.
