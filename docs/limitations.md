# Reproducibility scope and limitations

## Aggregate rebuild only

The supported rebuild starts from packaged aggregate outputs. It regenerates five figures and three main tables but does not refit any model, recalculate estimates, fetch raw administrative records, or recreate the analysis panels. A successful `reproduce.py` run therefore verifies the packaged payload and aggregate-to-artifact step only.

The original upstream pipeline that produced the frozen Brazil and Mexico panels and original route6 results was not recovered. The exact source-generation scripts for those results are not included as a complete, runnable pipeline. Full raw-data reproduction is unverified.

## Archived analytical scripts

The `analysis/work/upgrade_20261112/` and `analysis/work/revision_20261113/` trees are archival snapshots. Some scripts retain machine-specific paths and depend on inputs outside this package. `requirements-analysis.txt` records versions from the original model environment; it is not a complete lockfile, does not configure those paths, and does not make the archived scripts portable or demonstrably runnable.

## Unresolved source limitations

Packaging does not resolve incomplete exposure-source provenance and period-boundary documentation, or uncertainty in Mexico's historical procedure-code coverage and semantics. These affect interpretation of source and period coverage and remain unresolved by the aggregate rebuild.

## Optional reference and processed files

Reference PDFs under `reference/figures/`, if included, are contextual files and are not builder inputs. Any separately included processed analysis inputs should be treated as archival unless a builder explicitly reads them; their presence does not establish a raw-data reproduction route or public redistribution rights.
