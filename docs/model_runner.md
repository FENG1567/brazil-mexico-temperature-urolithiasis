# Optional model replay

`src/analysis/run_stage.py` prepares the archived service or nonlinear model in
an isolated directory under `artifacts/model_workspace/`. It copies the source
scripts while rewriting only their input/output path constants and the retired
`work/pylibs` path. It does not change estimation code or install packages.
The Brazil-panel default preserves the historical v1.0.0 behavior. For v1.1.0,
pass the corrected panel path with --brazil-panel; the chosen path is recorded in
the preparation manifest.

Preparation is the default. It copies the scripts, checks the expected inputs,
and writes `preparation_manifest.json`; it does not import or fit a model.
The first two commands below use the historical default. The next two prepare
the same stages with the corrected v1.1.0 panel.

```bash
python src/analysis/run_stage.py service --prepare-only
python src/analysis/run_stage.py nonlinear --prepare-only
```

~~~bash
python src/analysis/run_stage.py service --prepare-only --brazil-panel data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet
python src/analysis/run_stage.py nonlinear --prepare-only --brazil-panel data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet
~~~

The corrected-panel commands above only prepare workspaces. After reviewing
the manifest and making the required inputs available, opt in to historical-input
execution with the default commands below, or to corrected-panel execution with
the explicitly selected v1.1.0 commands:

```bash
python src/analysis/run_stage.py service --run
python src/analysis/run_stage.py nonlinear --run
```

~~~bash
python src/analysis/run_stage.py service --run --brazil-panel data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet
python src/analysis/run_stage.py nonlinear --run --brazil-panel data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet
~~~

Each stage runs from its own workspace. Copied scripts and the stage manifest
are under `artifacts/model_workspace/<stage>/`; model outputs are written under
that workspace's `outputs/` directory. The archived sources remain under
`analysis/work/`.

## Inputs and source keys

The preparation manifest reports each source key, expected path, and whether
the file is present. The Brazil panel rows below show the historical default;
select the corrected file with --brazil-panel for v1.1.0. Service inputs are:

| Source key | Repository path |
| --- | --- |
| `panel_composite_humidity` | Historical default: `artifacts/processed_inputs/panel_composite_humidity.parquet`; v1.1.0: `data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet` |
| `clinical_code_day_aggregate` | `data/processed_inputs/clinical_code_day_aggregate.parquet` |
| `ipca_monthly` | `data/aggregate_results/sensitivity_results/ipca_monthly.csv` |
| `sp_linked_daily_outcomes` | `data/processed_inputs/sp_linked_daily_outcomes.csv` |
| `holiday_calendar_dates` | `data/aggregate_results/sensitivity_results/holiday_calendar_dates.csv` |
| `acute_repair_model_results` | `data/aggregate_results/sensitivity_results/acute_repair_model_results.csv` |

Nonlinear inputs are:

| Source key | Repository path |
| --- | --- |
| `panel_composite_humidity` | Historical default: `artifacts/processed_inputs/panel_composite_humidity.parquet`; v1.1.0: `data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet` |
| `clinical_code_day_aggregate` | `data/processed_inputs/clinical_code_day_aggregate.parquet` |
| `mexico_composite_state_day_panel` | `data/processed_inputs/mexico_composite_state_day_panel.parquet` |
| `mex_state_day_2015_2023` | `data/processed_inputs/mex_state_day_2015_2023.csv` |
| `mexico_boundary_daily_counts` | `data/processed_inputs/mexico_boundary_daily_counts.csv` |
| `holiday_calendar_dates` | `data/aggregate_results/sensitivity_results/holiday_calendar_dates.csv` |
| `n20_composite_continuous_lag` | `data/aggregate_results/alternative_model_results/n20_composite_continuous_lag.csv` |
| `n20_composite_groups_corrected` | `data/aggregate_results/alternative_model_results/n20_composite_groups_corrected.csv` |

The nonlinear runner also copies
`analysis/specifications/prespecified_specification.json` and
`analysis/work/upgrade_20261112/mexico_strict_phenotype.py`. The archived M3
script calls the helper's Brazil loader. The helper's Mexico clinical and NASA
temperature constants are redirected to `data/processed_inputs/`, but are
optional for the current M3 entry point; the helper's Mexico-specific functions
need those files if called separately.

The service runner copies `analysis/specifications/analysis_specification.json`
to the service workspace because the archived output manifest records that
specification path.

For the archived `input_manifest()` check, preparation also materializes
`data/processed_inputs/clinical_code_day_aggregate.parquet` at
`artifacts/model_workspace/nonlinear/work/upgrade_20261112/clinical_code_day_aggregate.parquet`.
It uses a hard link when the filesystem allows it and falls back to a byte-for-
byte copy otherwise. The manifest records the method and path.

## Panel reassembly and preflight

The Brazil composite panel is distributed in
`data/processed_inputs/panel_composite_humidity.parquet.parts/` with its input
inventory at `data/processed_inputs/inventory.json`. Reassemble it before `--run`:

```bash
python src/data/reassemble.py
```

The helper materializes the historical panel at
`artifacts/processed_inputs/panel_composite_humidity.parquet`. The model runner
does not invoke reassembly implicitly. If you select the corrected v1.1.0 panel,
use its packaged file directly; legacy panel reassembly is not required. If the selected panel or another required input
is absent, `--run` stops before executing a model and lists the missing paths.

The nonlinear historical group table is an aggregate reference input; the
runner requires the packaged `n20_composite_groups_corrected.csv` source key and
does not rebuild that table from lower-level records. No model fit is claimed
by preparation or preflight alone. Execution also requires the scientific
Python environment used by the archived scripts, including DuckDB, NumPy,
pandas, pyfixest, and SciPy.
