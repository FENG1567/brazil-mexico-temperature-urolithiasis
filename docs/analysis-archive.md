# Archived analysis code and environment

The folders `analysis/work/upgrade_20261112/` and `analysis/work/revision_20261113/` preserve selected scripts used in the upgrade and revision stages. They are provided as historical research records. The supported public rebuild is `python reproduce.py`, which uses already-produced aggregate outputs and performs no model fitting.

The recorded original model environment was Python 3.12.14 with NumPy 2.5.4, pandas 3.1.0, DuckDB 2.0.0, SciPy 1.18.1, pyfixest 0.60.0, pyreaddbc 2.0.4, dbfread 2.0.7, and holidays 0.108. These pins are listed in `requirements-analysis.txt` as a version reference, not as a complete environment lock.

The archive includes the clinical aggregate construction in `analysis/work/upgrade_20261112/audit_inputs.py` and the IPCA retrieval and index construction in `analysis/work/upgrade_20261112/retrieve_prices_inspect_sp.py`. The clinical query sums SIH-SUS `DIAS_PERM` and `VAL_TOT` as `beddays` and `nominal_cost`. The IPCA code compounds monthly percentage changes and divides the December 2022 index by the index for each payment month. The recorded IPCA retrieval metadata are in `data/aggregate_results/sensitivity_results/ipca_download_provenance.json`. These scripts were recovered as source records; they were not executed during manuscript editing.

The archived scripts retain original machine-specific paths and require external source or intermediate files. The original upstream pipeline for the frozen Brazil/Mexico panels and route6 results was not recovered. Consequently, the archives are not claimed to run on a clean checkout, and they do not provide a verified raw-data reproduction.
