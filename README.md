# Short-Term Temperature and Urolithiasis Hospitalizations in Brazil and Mexico: Exploratory Heat–Smoke Comparisons With Cholelithiasis and Cholecystitis in Brazil


Code release: **v1.1.0**  
Corrected-results identity: **TEMP-20261115-C1**  
Historical release: **v1.0.0** remains available unchanged.  
Authors: Shuo Feng, Chao Shi, Tianming Liu, Haoran Ding, Fenghe Hang, Zhenghao Huang, and Xianzhi Meng.

This v1.1.0 release packages the corrected Brazil temperature panel, fixed-grid weights, and corrected-results identity TEMP-20261115-C1. It rebuilds five figures and three main tables from the packaged aggregate results; this build does not retrieve raw data, reconstruct weather inputs, or fit models. See [Temperature-input correction](docs/temperature_input_correction.md), [Reproduction](docs/reproduction.md), and [Scope and limitations](docs/limitations.md) for the supported scope and prerequisites.

## Rebuild the released artifacts

Use Python 3.12.3, the version used for the tested plotting environment. In an activated environment, install the plotting and table-export dependencies and run:

```sh
python -m pip install -r requirements-figures.txt
python reproduce.py --verify-only
python reproduce.py
```

`--verify-only` checks the manifest-declared input files by byte count and SHA-256. The default command verifies those files, then rebuilds figures and tables without refitting models. Generated files go under `artifacts/` and are ignored by Git. See [Reproduction](docs/reproduction.md) for outputs and details.

## Contents

- `data/aggregate_results/`: supplementary aggregate outputs and the released `main_tables/table1.csv` through `table3.csv`.
- `data/figure_inputs/`: the frozen and revision aggregate inputs mapped to the five figures.
- `manifest.csv`: repository-relative payload paths, SHA-256 values, byte counts, and roles.
- `src/`: figure and main-table builders called by `reproduce.py`.
- `data/processed_inputs/`: geographic-day inputs, including the retained historical Brazil panel and the compact corrected Brazil panel; see [Processed inputs](docs/processed-inputs.md) and [Temperature-input correction](docs/temperature_input_correction.md).
- `analysis/work/upgrade_20261112/` and `analysis/work/revision_20261113/`: archival model and analysis scripts; these are not the supported aggregate rebuild path.
- `reference/figures/`: optional reference PDFs, if included. They are not read by the builders.

The optional [model runner](docs/model_runner.md) prepares isolated service and nonlinear stages. For v1.1.0, select the corrected panel explicitly with the --brazil-panel option described in [Temperature-input correction](docs/temperature_input_correction.md). Preparation does not fit models; the full upstream raw-weather reconstruction remains outside this package.

## Data, code, and reuse

The source families include Brazil SIH/DATASUS and CNES, Mexico DGIS/SAEH SSA, ERA5-Land and ERA5 single-level reanalysis for the corrected Brazil temperature panel, NASA POWER, IPCA, and the GFED-driven fire PM2.5 concentration product described by Hu et al. (2025; DOI: [10.5194/essd-17-3741-2025](https://doi.org/10.5194/essd-17-3741-2025)). See [Data sources](docs/data-sources.md) and [Temperature-input correction](docs/temperature_input_correction.md) for retrieval limits and correction scope.

Author-written code and associated code documentation are released under the [MIT License](LICENSE). Third-party data, provider-derived inputs, and software dependencies retain their applicable terms; see [LICENSE_STATUS.md](LICENSE_STATUS.md).

## Citation

Use [CITATION.cff](CITATION.cff) for the citation metadata. This prepared corrected code release is v1.1.0 (results identity TEMP-20261115-C1). The historical v1.0.0 release remains unchanged. No DOI has been assigned.

## Display revision

The v1.1.0 figures and main tables are rebuilt from the corrected aggregate inputs packaged in this release. The build step only renders those results; model execution is separate and requires the input and software prerequisites documented in the repository.

## Code availability

The v1.1.0 release package, including the corrected Brazil geographic-day panel and fixed-grid weight table, is prepared for https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis (results identity TEMP-20261115-C1). The repository instructions rebuild five figures and three main tables from packaged results and describe the separate, explicitly selected model-rerun path. Author-written code is licensed under the MIT License; third-party data retain their providers’ terms. The previous v1.0.0 release remains available as the historical package.
