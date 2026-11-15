# Short-Term Temperature and Urolithiasis Hospitalizations in Brazil and Mexico: Exploratory Heat–Smoke Comparisons With Cholelithiasis and Cholecystitis in Brazil


Code release: **v1.0.0**  
Results/display identity: **SRR-20261114-C3-display-3**  
Authors: Shuo Feng, Chao Shi, Tianming Liu, Haoran Ding, Fenghe Hang, Zhenghao Huang, and Xianzhi Meng.

This repository preserves the released aggregate results and rebuilds five figures and three main tables from those packaged results. The build does not retrieve raw data, reconstruct the analytic panels, or fit statistical models. See [Reproduction](docs/reproduction.md) and [Scope and limitations](docs/limitations.md) before interpreting it as a reproducibility package.

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
- `data/processed_inputs/`: six geographic-day or stratified daily input files, including a losslessly split Brazil panel; see [Processed inputs](docs/processed-inputs.md).
- `analysis/work/upgrade_20261112/` and `analysis/work/revision_20261113/`: archival model and analysis scripts; these are not the supported aggregate rebuild path.
- `reference/figures/`: optional reference PDFs, if included. They are not read by the builders.

The release package contains the 121 supplementary aggregate files, one additional frozen historical reference table, 19 figure-input copies, and three main-table CSVs within the aggregate collection. The optional [model runner](docs/model_runner.md) prepares isolated service and nonlinear stages from existing analysis panels. Its input checks and path substitutions have been checked without refitting models; they do not establish a complete raw-data pipeline.

## Data, code, and reuse

The source families include Brazil SIH/DATASUS and CNES, Mexico DGIS/SAEH SSA, ERA5-Land, NASA POWER, IPCA, and the GFED-driven fire PM2.5 concentration product described by Hu et al. (2025; DOI: [10.5194/essd-17-3741-2025](https://doi.org/10.5194/essd-17-3741-2025)). See [Data sources](docs/data-sources.md) for the distinction between named source locations and verified downloads.

Author-written code and associated code documentation are released under the [MIT License](LICENSE). Third-party data, provider-derived inputs, and software dependencies retain their applicable terms; see [LICENSE_STATUS.md](LICENSE_STATUS.md).

## Citation

Use the citation metadata in [CITATION.cff](CITATION.cff). The repository is [https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis](https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis); the code release is `v1.0.0`. No DOI has been assigned.

## Display revision

The display revision retains the numerical results and updates the five figures and three main tables. Complete pre-revision main tables remain in `data/aggregate_results/full_model_tables/`; five themed secondary-table CSVs are in `display_tables/`. No model was refitted.

## Code availability

Code and supporting aggregate inputs are available at https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis (version v1.0.0). The repository includes analysis scripts, processed geographic-day inputs, and instructions for rebuilding the five main figures and three main tables from released aggregate results. Author-written code is licensed under the MIT License; third-party data retain their providers’ terms. The repository documentation describes the scope and prerequisites for model reruns.
