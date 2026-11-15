# Data sources and retrieval status

This page identifies source systems and products associated with the analyses. It distinguishes named sources from verified downloads: `reproduce.py` does not retrieve data, and the package manifest verifies local payload bytes rather than validating how an upstream file was obtained.

| Source family | Source location or product identified | Package and verification status |
| --- | --- | --- |
| Brazil hospital records | SIH/DATASUS | The public package is an aggregate-result release. The raw-record pipeline and a full raw-data recreation are not verified here. |
| Mexico hospital records | DGIS/SAEH SSA | The public package is an aggregate-result release. The raw-record pipeline and a full raw-data recreation are not verified here. |
| Weather and environmental inputs | ERA5-Land and NASA POWER | These source families are identified in the analysis record. No downloader or source-query replay is part of the package build. |
| Fire-related PM2.5 | GFED-driven fire PM2.5 concentration product described by Hu et al. (2025), DOI [10.5194/essd-17-3741-2025](https://doi.org/10.5194/essd-17-3741-2025) | The DOI identifies the cited product description. The package does not replay product retrieval or establish redistribution rights. |
| Health-system and price information | CNES and IPCA | The IPCA retrieval and index-construction script is archived at `analysis/work/upgrade_20261112/retrieve_prices_inspect_sp.py`. It uses Banco Central do Brasil SGS series 433 for monthly IPCA changes; retrieval metadata and the applied monthly factors are packaged. The default rebuild does not rerun these downloads or the CNES source retrieval. |

## What a raw-data reconstruction would need to record

For any future source-to-result reconstruction, document the exact portal or catalog entry, query parameters, date and geographic coverage, source release or code vintage, retrieval date, original filenames, and checksums. Then document the cleaning, linkage, geography, and exposure processing that creates analytic inputs. These steps are guidance for transparent provenance; they are not a claim that this package has completed or validated those downloads.

The included manifest records path, SHA-256, byte count, and role for local package payloads. Figure provenance additionally maps plotted values to packaged CSV/JSON records. Neither record substitutes for a complete upstream raw-file acquisition log.

Data-provider terms may differ by source and file. See [LICENSE_STATUS.md](../LICENSE_STATUS.md); inclusion in this package does not itself establish redistribution permission.
