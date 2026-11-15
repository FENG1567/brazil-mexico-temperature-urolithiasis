# Processed analysis inputs

`data/processed_inputs/inventory.json` lists the six preserved input files, original byte counts and SHA-256 values, row counts, column names and storage types. These are geographic daily panels or stratified daily counts, rather than individual hospital records. Their source windows can exceed the fitted study windows: the archived model code applies the relevant eligibility and lag filters. Input rows are not admission counts.

| Input | Observational level and purpose |
| --- | --- |
| `panel_composite_humidity.parquet` | Brazil health-region/calendar-day exposure and aggregate outcome panel. Preserved as eight binary parts. |
| `mexico_composite_state_day_panel.parquet` | Mexico state/calendar-day exposure, matching strata and aggregate outcome panel. |
| `clinical_code_day_aggregate.parquet` | Brazil counts and sums grouped by region, date and administrative diagnostic/admission/procedure categories. `IDENT` is an authorization-type category, not an individual identifier. |
| `sp_linked_daily_outcomes.csv` | Brazil region/day counts of linked procedure-record categories. No admission/procedure identifiers are included. |
| `mexico_boundary_daily_counts.csv` | Mexico state/day aggregate admission and record-defined procedure counts used in the admission-window correction. |
| `mex_state_day_2015_2023.csv` | Mexico state/day ERA5 temperature series, including calendar fields. |

The clinical aggregate construction is preserved in `analysis/work/upgrade_20261112/audit_inputs.py`: `beddays` sums SIH-SUS `DIAS_PERM` and `nominal_cost` sums `VAL_TOT`, grouped by region, admission date, diagnosis and administrative categories. The service models then select `IDENT=1`, strict `N20`, and `CAR_INT=02`.

`hrid` and `state` identify geographic units; `date` is a calendar date. Diagnosis fields contain administrative ICD categories. Admission and procedure outcomes are event counts; `beddays`/`*_bd` are summed hospital days. The scripts define filtering, zero filling, analysis windows and denominator rules. Temperature variables are in degrees Celsius; fire-PM2.5 concentration is in micrograms per cubic metre. Other fields retain their source coding and units; inspect the archived script and analysis specification rather than inferring their scale from storage type.

Reassemble the large Brazil panel without changing any bytes:

```sh
python src/data/reassemble.py --verify-only
python src/data/reassemble.py
```

The helper validates every part and the reconstructed file against the original checksum. It writes the full panel only to ignored `artifacts/processed_inputs/`. Each part is at most 16 MiB and below GitHub's single-file limit. The parts do not require Git LFS.

See [Model runner](model_runner.md) for the optional isolated revision-stage entrypoints. These are distinct from the default figure/table rebuild. They have only been prepared and checked here; statistical models have not been refitted in this packaging task. No complete original source-to-panel pipeline is asserted. Source-specific access, attribution, and reuse terms continue to apply to these provider-derived inputs; the MIT code license does not replace them.
