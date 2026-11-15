# Brazil temperature-input correction

Release: v1.1.0  
Corrected-results identity: TEMP-20261115-C1  
Historical release: v1.0.0 remains available unchanged.

## What changed

The Brazil geographic-day panel in v1.1.0 adds the ERA5 single-level weighted contribution to the ERA5-Land contribution using the supplied fixed grid weights. The daily temperature mean and maximum therefore represent the recorded full fixed-weight sum across the two products. The calculation does not renormalize weights or replace a missing contribution with zero.

The panel has 2,629,800 rows for 450 health regions from 2008-01-01 through 2023-12-31. Temperature thresholds are recalculated by health region from corrected daily mean temperature over the original 2008-2020 reference window. Heat, cold, and heat-smoke composite indicators are recalculated under the existing strict threshold rules. Outcomes, dates, and unrelated covariates remain unchanged.

The fixed-weight file contains 66,495 grid records covering the 450 health regions. Within each region, the supplied weights across ERA5-Land and ERA5 single-level grid support sum to one. This documents the supplied support and assembly; it does not independently reconstruct the population raster or prove the population-support lineage.

Relative humidity is included for descriptive use only. The shared hourly/grid validity mask has not been independently verified, and the original missing RH values in 2023 are preserved. This release does not authorize RH-adjusted inferential models.

## Released inputs and checksums

| File | Rows / bytes | SHA-256 |
| --- | ---: | --- |
| data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet | 2,629,800 rows; 83,382,645 bytes | b837328e149714d631257aa43eaffbb95506510b6b9275a35a6fd892d977d597 |
| src/data/temperature_repair/final_production_exposure_weights.csv | 6,935,134 bytes | 2c5311c83064aa9ea4845f941163ad508d7a91305807247290549b26eeda51ae |

The historical Brazil panel and its split parts remain in the package for v1.0.0 compatibility. The corrected panel is a separate file; it does not replace the historical input.

Check the released payload hashes from the repository root:

~~~sh
sha256sum data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet src/data/temperature_repair/final_production_exposure_weights.csv
~~~

The packaged root input gate and input-validation report record passing checks for source identity, fixed-weight accounting, Celsius units, source dates, source-specific eligible mass, component-to-total reconciliation, unique region-day keys, absence of missing-as-zero assembly, and unchanged outcomes and unrelated fields. These checks cover the supplied fixed health-region/grid support and recorded contribution assembly. They do not independently rebuild hourly weather data or the population-support raster. Passing the input gate is not a model-fit check; the historical old-input replication evidence does not establish refitted estimates for every corrected-input model.

## Rebuild the released figures and tables

The package-level rebuild uses the released aggregate result files and does not fit models:

~~~sh
python -m pip install -r requirements-figures.txt
python reproduce.py --verify-only
python reproduce.py
~~~

The first command checks manifest-listed payloads. The second recreates the five figures and three main tables under artifacts/. These commands do not reconstruct the corrected weather panel.

## Prepare or run corrected-panel model stages

The service and nonlinear stage runner keeps the historical v1.0.0 input as its default. Select the v1.1.0 panel explicitly:

~~~sh
python src/analysis/run_stage.py service --prepare-only --brazil-panel data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet
python src/analysis/run_stage.py nonlinear --prepare-only --brazil-panel data/processed_inputs/panel_composite_humidity_exposure_repaired.parquet
~~~

Preparation copies the archived model source into artifacts/model_workspace/ and records the selected panel path in its preparation manifest; it does not fit models. After checking each manifest and ensuring that the required inputs and archived-script software dependencies are available, opt in to model execution with the same command and --run instead of --prepare-only. The override only routes the Brazil panel used by these service and nonlinear stage runners. Other archived analysis scripts are not globally redirected by this option.

## Rebuild the corrected panel from source contributions

The repair runner now accepts paths instead of relying on the original server directory. Reassemble the historical Brazil panel first:

~~~sh
python src/data/reassemble.py
~~~

With the exact original contribution files available locally, run the fail-closed repair from the repository root. The assembly manifest and fixed-weight CSV default to their released paths; the old panel defaults to the reassembled file:

~~~sh
python src/data/temperature_repair/repair_exposure.py --primary-weather /private/path/primary_D_signed_rain_weather.parquet --mixed-contributions /private/path/mixed_region_day_contributions.parquet --old-panel artifacts/processed_inputs/panel_composite_humidity.parquet --output artifacts/temperature_repair/rebuilt_panel.parquet --report artifacts/temperature_repair/rebuild_validation.json
~~~

Use --validate-only to run the input and assembly gates without writing the corrected panel. The runner checks the expected SHA-256 values and byte sizes for all source inputs. It refuses to overwrite an existing panel output.

The primary and mixed contribution parquet files and the 948 MB upstream source-weather files are private and are not distributed in this repository. As a result, a public checkout can use the released corrected panel for supported analysis reruns and verify its checksum, but cannot reconstruct that panel from raw weather inputs using only the repository contents. No public weather downloader or raw-data pipeline is claimed.

## Reproducibility and licensing scope

The package includes author-written code under the MIT License, the released corrected geographic-day panel, the fixed-weight table, and the frozen results used by the figure/table builder. Provider data and provider-derived inputs retain their applicable terms. The corrected exposure is a mean over the supplied fixed health-region/grid support, not an independently reconstructed population-weighted exposure. See LICENSE_STATUS.md and docs/data-sources.md.
