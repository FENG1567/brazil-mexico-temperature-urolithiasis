# FROZEN RESULTS: TEMP-20261115-C1

Corrected Brazilian temperature input; code version v1.1.0. Historical v1.0.0 outputs remain archived. No Mexico refits were performed.

## Submitted tables

### data/aggregate_results/main_tables/table1.csv

SHA-256: 88d748b650c6291563ab44009726bac013f5182a783d7bd82140279640682494

Row 1: Analysis stream | Period and diagnostic scope | Source frame and admissions | Model rows / geographic units | Exposure assignment
Row 2: Brazil primary | 2008–2022; N20–N23 | 1,132,163 admissions; 450 health regions; 2,465,550 region-days | 1,481,520 model rows | Daily Brazilian temperature from ERA5-Land and ERA5 single-level contributions; fire-source PM₂.₅
Row 3: Brazil harmonized broad series | 2018–2022; N20–N23 | 431,374 unique AIH authorizations | 46,387 model rows; 27 states | Population-weighted health-region temperature aggregated to state
Row 4: Brazil strict-diagnosis state series | 2018–2022; strict N20 | 321,201 admissions | 45,436 model rows; 27 states | Population-weighted health-region temperature aggregated to state
Row 5: Brazil urgent clinical / procedure extensions | 2018–2022; urgent strict N20 | 226,402 admissions; 84,764 with any drainage record; 71,833 with any stone-extraction/lithotripsy record | 450 health regions; outcome-specific model rows | Categorical heat/smoke contrasts and continuous lag 0–7 models
Row 6: Mexico baseline state series | 2018–2022 admission window; SAEH/SSA | 40,415 N20–N23 and 32,417 strict-N20 admissions | 43,483 broad and 39,359 strict-N20 model rows; 32 states | ERA5-Land state-area weighting; NASA POWER state-capital sensitivity
Row 7: Mexico source-expanded state series | 2018–2022 admission window; discharge sources through 2023 | 40,448 N20–N23 and 32,446 strict-N20 admissions | 43,503 broad and 39,390 strict-N20 model rows; 32 states | ERA5-Land state-area weighting
Row 8: Mexico municipality sensitivity | 2018–2022; N20–N23; exposure-complete municipality frame | Raw 40,907 admissions; complete exposure 39,809 | 133,775 model rows; 1,830 municipalities | Same-day NASA POWER temperature at state-capital coordinates and municipality fire PM₂.₅

### data/aggregate_results/main_tables/table2.csv

SHA-256: b29ed764e2351bfe369221d08ae28c82eccd792b4a8ae22d9641c4e74f2331bd

Row 1: Population / outcome | Exposure contrast, unit, lag, and reference | RR (95% CI) | Model rows
Row 2: Brazil, N20–N23, national region-day series | Temperature, per 1 °C; lag 0–7 cumulative; region/date clusters | 1.0144 (1.0108–1.0180) | 1,481,520
Row 3: Brazil, N20–N23, national region-day series | Fire-source PM₂.₅, per 10 µg/m³; lag 0–7 cumulative; region/date clusters | 1.0555 (1.0318–1.0797) | 1,481,520
Row 4: Brazil, N20–N23, harmonized state series | Temperature, per 1 °C; lag 0–7 cumulative | 1.0143 (1.0096–1.0189) | 46,387
Row 5: Mexico, N20–N23, harmonized state series | ERA5-Land temperature, per 1 °C; lag 0–7 cumulative | 1.0289 (1.0174–1.0406) | 43,483
Row 6: Brazil, strict N20, calendar-adjusted | Temperature, per 1 °C; lag 0–7 cumulative; national public dates and federal elections; two-way state/date clusters | 1.0115 (1.0052–1.0179) | 45,436
Row 7: Mexico, strict N20, source-expanded | ERA5-Land temperature, per 1 °C; lag 0–7 cumulative; calendar controls; two-way state/date clusters | 1.0202 (1.0073–1.0333) | 39,390
Row 8: Brazil, N20–N23, compound days | Same-day compound vs neither; temperature > P90 and fire PM₂.₅ > P75; region/date clusters | 1.0683 (1.0515–1.0853) | 1,484,094
Row 9: Mexico, N20–N23, compound days | Same-day compound vs neither; temperature > P90 and fire PM₂.₅ > P75 | 0.996 (0.936–1.060) | 43,671
Row 10: Brazil, N20–N23, multiplicative interaction | Same-day heat × smoke; ratio beyond multiplicative main associations; region/date clusters | 1.0308 (1.0100–1.0521) | 1,484,094
Row 11: Brazil, K85 acute pancreatitis | Compound day vs neither day; temperature > P90 and fire PM₂.₅ > P75 | 0.9610 (0.9285–0.9947) | 570,696
Row 12: Brazil, K81 cholecystitis | Compound day vs neither day; temperature > P90 and fire PM₂.₅ > P75 | 0.9987 (0.9840–1.0137) | 1,561,954
Row 13: Brazil, K80 gallstones | Compound day vs neither day; temperature > P90 and fire PM₂.₅ > P75 | 1.0082 (0.9983–1.0181) | 1,735,546

### data/aggregate_results/main_tables/table3.csv

SHA-256: 3a3a28e32d275894d8ec03c9e856ecebbd3858a70a6b43727eee92dcde98a581

Row 1: Outcome | Contrast / estimand | 2018–2022 ratio (95% CI) | Excluding 2020–2021 ratio (95% CI)
Row 2: Urgent/elective RR ratio, N20–N23 | Heat-only vs neither; calendar-adjusted | 1.0546 (1.0097–1.1015) | 1.0655 (1.0067–1.1278)
Row 3: Urgent/elective RR ratio, strict N20 | Heat-only vs neither; calendar-adjusted | 1.0409 (0.9924–1.0918) | 1.0486 (0.9846–1.1168)
Row 4: Urgent admissions, strict N20 | Heat-only vs neither; calendar-adjusted | 1.0509 (1.0259–1.0765) | 1.0272 (0.9933–1.0623)
Row 5: Admission-attributed hospital-days, strict N20 | Heat-only vs neither; calendar-adjusted | 1.0486 (1.0043–1.0948) | 1.0428 (0.9874–1.1012)
Row 6: Total reimbursement, strict N20 (Dec 2022 BRL) | Heat-only vs neither; calendar-adjusted | 1.0371 (1.0035–1.0719) | 1.0179 (0.9720–1.0661)
Row 7: Reimbursement per admission, strict N20 | Heat-only vs neither; geometric mean ratio | 0.9964 (0.9729–1.0205) | —
Row 8: Drainage/stent records, strict N20 | Temperature per 1 °C; cumulative lag 0–7 | 1.0119 (1.0046–1.0193) | 1.0122 (1.0033–1.0212)
Row 9: Stone-extraction/lithotripsy records, strict N20† | Temperature per 1 °C; cumulative lag 0–7 | 1.0062 (0.9991–1.0134) | —

### data/aggregate_results/display_tables/tableS1a.csv

SHA-256: 3222739d5a29484269b2b76d744febc66f2547c2261c76aee36d6dd4444e129c

Row 1: Analysis stream | Diagnostic scope and period | Source events / exposure frame | Complete-case and fitted sample
Row 2: Brazil primary series | N20–N23; 2008–2022 | 1,132,163 admissions; 450 health regions; 2,465,550 region-days | 1,130,926 events after lag 0–7 completion; 1,481,520 fitted region-day rows; 227,283 zero-event strata removed 980,880 rows
Row 3: Brazil harmonized broad state series | N20–N23; 2018–2022 | 431,374 unique AIH authorizations | 46,387 rows; 27 states; original state-clustered comparison
Row 4: Brazil strict-N20 state series | N20; 2018–2022 | 321,201 admissions | 320,389 events after lag completion; 45,436 fitted rows; 856 zero-event strata removed 3,677 rows; 27 states
Row 5: Brazil urgent-N20 service models | Urgent N20; 2018–2022 | 226,402 admissions; 805,746 admission-attributed hospital-days; BRL 194,858,997 in December 2022 reimbursement | Admissions and reimbursement: 366,584 rows; hospital-days: 356,593 rows; 450 regions and 1,826 dates
Row 6: Brazil service pandemic-exclusion sensitivity | Urgent N20; 2018–2019 and 2022 outcome dates | 136,094 admissions; 497,269 admission-attributed hospital-days; BRL 117,395,551 reimbursement | Admissions and reimbursement: 220,501 rows; hospital-days: 215,201 rows; 450 regions and 1,095 dates
Row 7: Brazil drainage-record lag model | Urgent N20 admissions with any linked drainage record; 2018–2022 | 84,764 full-window records; 84,572 events after lag 0–7 completion | 180,241 fitted rows; 443 retained regions and 1,819 dates
Row 8: Mexico baseline state series | N20–N23 / N20; 2018–2022 admission window | 40,415 / 32,417 valid-state admissions | 43,483 / 39,359 rows; 32 states in the original state-series models
Row 9: Mexico source-expanded state series | N20–N23 / N20; 2018–2022 admission window, discharge sources through 2023 | 40,448 / 32,446 admissions | 43,503 / 39,390 rows; 32 states
Row 10: Mexico municipality sensitivity | Broad N20–N23; 2018–2022; separate exposure-complete branch | 40,907 raw records; 492 invalid/out-of-range state codes; 40,415 valid-state records; 39,809 complete-exposure records across 1,830 municipalities | 606 valid-state records lacked municipality-date exposure; 133,775 model rows
Row 11: Mexico strict-N20 spline sensitivity | Strict N20; source-expanded state frame | 32,446 admissions; 32,352 events after lag completion | 39,390 fitted rows; 4,372 zero-event strata removed 18,818 rows; 32 states and 1,819 dates

### data/aggregate_results/display_tables/tableS1b.csv

SHA-256: dc2c18db1d27c01e07a1024bae25b36297d3f7ae2bcf7203c7ac5b3c2fd61a06

Row 1: Country / period | Dimension | Category | Admissions | % of cohort
Row 2: Brazil, 2008–2022 | sex_code | 1 (male) | 567,503 | 50.126
Row 3: Brazil, 2008–2022 | sex_code | 3 (female) | 564,660 | 49.874
Row 4: Brazil, 2008–2022 | age_group | 0–14 | 27,381 | 2.418
Row 5: Brazil, 2008–2022 | age_group | 15–44 | 598,260 | 52.842
Row 6: Brazil, 2008–2022 | age_group | 45–64 | 388,975 | 34.357
Row 7: Brazil, 2008–2022 | age_group | 65+ | 117,547 | 10.383
Row 8: Mexico, 2018–2022 baseline | sex_code | 1 (male) | 17,131 | 42.388
Row 9: Mexico, 2018–2022 baseline | sex_code | 2 (female) | 23,281 | 57.605
Row 10: Mexico, 2018–2022 baseline | sex_code | 9 (unspecified) | 3 | 0.007
Row 11: Mexico, 2018–2022 baseline | age_category | 0–14 years | 1,208 | 2.989
Row 12: Mexico, 2018–2022 baseline | age_category | 15–44 years | 18,069 | 44.709
Row 13: Mexico, 2018–2022 baseline | age_category | 45–64 years | 15,505 | 38.364
Row 14: Mexico, 2018–2022 baseline | age_category | 65+ years | 3,943 | 9.756
Row 15: Mexico, 2018–2022 baseline | age_category | 2019-source age, 15–29 years (shown separately) | 1,689 | 4.179
Row 16: Mexico, 2018–2022 baseline | age_category | unspecified age (CVEEDAD=9) | 1 | 0.002

### data/aggregate_results/display_tables/tableS1c.csv

SHA-256: 407d136bd39549f9083fec5c8c7b2f7635fe7c972f09604a1fd856c11262609b

Row 1: Year | Brazil N20–N23 admissions | Mexico valid-state N20–N23 admissions
Row 2: 2008 | 63,648 | Outside study window
Row 3: 2009 | 66,786 | Outside study window
Row 4: 2010 | 68,912 | Outside study window
Row 5: 2011 | 66,357 | Outside study window
Row 6: 2012 | 67,653 | Outside study window
Row 7: 2013 | 67,369 | Outside study window
Row 8: 2014 | 71,956 | Outside study window
Row 9: 2015 | 73,628 | Outside study window
Row 10: 2016 | 75,422 | Outside study window
Row 11: 2017 | 79,058 | Outside study window
Row 12: 2018 | 83,057 | 9,403
Row 13: 2019 | 89,490 | 10,517
Row 14: 2020 | 78,175 | 4,608
Row 15: 2021 | 82,258 | 6,672
Row 16: 2022 | 98,394 | 9,215

### data/aggregate_results/display_tables/tableS2a.csv

SHA-256: cbe08e0f94c59681629b78d81c121fc4a1bb2b9bfe78a8db4cdf7589af8e7484

Row 1: Analysis | Definition | Interpretation
Row 2: Continuous temperature | Daily mean temperature; per 1 °C; distributed lag 0–7 unless otherwise labelled | Cumulative RR is the exponentiated sum of separate lag coefficients
Row 3: Fire-related PM₂.₅ | Daily concentration; per 10 µg/m³; lag 0–7 in specified models | Not included in the source-expanded Mexican broad/strict state calendar models
Row 4: Categorical heat and smoke | Heat > area-specific P90; smoke > area-specific P75; reference group neither; Brazil thresholds use 2008–2020 and Mexico 2018–2020 | Heat-only, smoke-only, and compound contrasts are same-day categorical comparisons, not per-unit effects
Row 5: National calendar control | 49 Brazilian dates and 38 Mexican dates in 2018–2022: recurring public dates plus national federal-election dates | Does not capture every local holiday, optional day, local election, or hospital schedule
Row 6: Same-day temperature quintiles | Pooled area-day cut points: Brazil 2008–2023; Mexico 2018–2022; Q1 reference | Quintiles are not cut locally. Q5-versus-Q1 RRs were 1.094 (1.075–1.114) in Brazil and 1.142 (1.047–1.247) in Mexico; Mexico Q4 was 1.144 (1.056–1.238), so the pattern is not monotonic
Row 7: Mexico state same-day joint sensitivity | NASA POWER state-capital temperature and fire-related PM₂.₅, jointly modeled on the same day; 43,671 rows | Temperature RR 1.0165 (1.0091–1.0240) per 1 °C; fire-PM₂.₅ RR 0.874 (0.814–0.939) per 10 µg/m³. The below-one fire estimate is not interpreted as protective
Row 8: Mexico municipality-stratified sensitivity | Broad N20–N23; same-day NASA POWER temperature at state-capital coordinates and municipality fire PM₂.₅ | This is not a lag 0–7 model; municipality refers to the matching stratum and the temperature remains a state-capital point assignment

### data/aggregate_results/display_tables/tableS2b.csv

SHA-256: 83afe5057e5fab38e46e8f82a706a235fe7eca8c61bda84159ca11a3a664755c

Row 1: Frame | Exposure | Mean (SD) | Median (P25–P75) | Missing / area-days
Row 2: Brazil, 450 health regions | tmean (°C) | 23.36 (3.92) | 24.08 (21.17–26.13) | 0 / 2,465,550
Row 3: Brazil, 450 health regions | fire_pm25 (µg/m³) | 1.89 (3.65) | 0.72 (0.26–2.01) | 0 / 2,465,550
Row 4: Brazil, 450 health regions | rh (%) | 72.81 (13.80) | 75.30 (64.78–83.34) | 0 / 2,465,550
Row 5: Mexico, valid-state municipalities | tmean (°C) | 19.79 (5.15) | 19.33 (16.22–23.14) | 0 / 3,390,882
Row 6: Mexico, valid-state municipalities | fire_pm25 (µg/m³) | 1.63 (2.74) | 0.62 (0.18–1.89) | 49,302 / 3,390,882
Row 7: Mexico, valid-state municipalities | rh (%) | 65.47 (17.23) | 68.94 (54.11–79.33) | 0 / 3,390,882

### data/aggregate_results/display_tables/tableS3.csv

SHA-256: 742127f428c2ebad30aef76ff863937b81515a93b51e27243ed97ab2e374578e

Row 1: Country and endpoint | Population and record source | Code definition | Full-window records in the final analysis
Row 2: Brazil: drainage/stent-related | Urgent strict N20; main-operation and/or SP detail record | 0409010103, 0409010162, 0409010170, 0409010278, 0409010286, 0409010294 | 84,764
Row 3: Brazil: stone extraction/lithotripsy-related | Urgent strict N20; main-operation and/or SP detail record | 0309030102, 0309030110, 0309030129, 0309030137, 0409010146, 0409010154, 0409010189, 0409010227, 0409010235, 0409010316, 0409010391, 0409010561, 0409010596 | 71,833
Row 4: Mexico: nephrostomy | Strict N20; PROMED | 5502, 5503, 5504 | 1,760
Row 5: Mexico: ureteral catheterization | Strict N20; PROMED | 598X | 3,313
Row 6: Mexico: urinary stone fragmentation/lithotripsy | Strict N20; PROMED | 5504, 5995, 9851 | 2,107
Row 7: Mexico: any recorded surgery | Strict N20; any TIPO=Q row | Surgery flag, not a code-specific endpoint | 10,498
Row 8: Mexico: no recorded TIPO=Q surgery | Strict N20; no TIPO=Q row | Complement of the recorded surgery flag | 21,948

### data/aggregate_results/display_tables/tableS4.csv

SHA-256: 344b1f984c127e0bd913c8737ff6df1fe6e48ae94d44808c354162bfca8dd8b8

Row 1: Analysis stream | Outcome and exposure specification | Matching strata | Cluster-robust inference
Row 2: Brazil regional primary | Daily N20–N23 count; linear temperature and fire-PM₂.₅ lags 0–7; no humidity or holiday term in this model | Health region × year-month × weekday | Primary linear model (Table 2): health-region/date two-way CRV1; matched spline: health-region/date two-way. Region-only clustering is an alternative inference specification
Row 3: Brazil strict-N20 / Mexico strict-N20 | Daily N20 count; matched linear and three-df natural-cubic temperature lags 0–7 plus national calendar design | State × year-month × weekday; national-calendar indicator in the mean model | Two-way state and date, including intersection; t df=min(retained cluster counts)-1
Row 4: Brazil / Mexico common-period state-level state comparisons | Daily N20–N23 count; temperature lags 0–7; exposure product and calendar as labelled | State × year-month × weekday | State for the common-period state-level estimates
Row 5: Mexico source-expanded broad and strict models | N20–N23 or N20 count; ERA5-Land temperature lags 0–7 plus national calendar indicators | State × year-month × weekday | Two-way state and date; small-cluster t reference
Row 6: Brazil service models | Urgent N20 admissions, admission-attributed hospital-days, and 2022-BRL total reimbursement; heat-only, smoke-only, compound, holiday | Health region × year-month × weekday | Two-way health region/date; 450 retained region clusters and 1,826 date clusters; df=449
Row 7: Brazil pandemic-exclusion service sensitivity | Same outcomes and model as service model; outcome dates limited to 2018–2019 and 2022 | Health region × year-month × weekday | Two-way health region/date; 450 retained region clusters and 1,095 date clusters; df=449
Row 8: Brazil urgency | Joint urgent/elective counts for N20–N23 or N20; heat-only urgent-to-elective RR ratio plus national holiday | Health region × year-month × weekday × outcome | Two-way health region/date; joint covariance used for ratio interval
Row 9: Brazil drainage-record models | Urgent strict-N20 any-drainage outcome; temperature and fire-PM₂.₅ lags 0–7; calendar-adjusted main model and pandemic sensitivity | Health region × year-month × weekday | Two-way health region/date; model-specific retained cluster counts and t reference
Row 10: Brazil procedure-record categorical models | Urgent strict-N20 record counts; temperature/fire categories and model specifications as labelled | Health region × year-month × weekday | Health-region clustering; q families as labelled
Row 11: Mexico specific procedures | Strict-N20 PROMED code-group counts; temperature and fire lags 0–7 plus 38-date calendar | State × year-month × weekday | Two-way state/date with model-specific t reference; current intervals use the explicit PSD sensitivity
Row 12: Mexico TIPO=Q / no-TIPO=Q | Strict-N20 groups; same-day temperature and fire terms from original state-day panel | State × year-month × weekday | State; no formal interaction test
Row 13: Mexico municipality sensitivity | Broad N20–N23; same-day NASA POWER state-capital temperature and municipality fire PM₂.₅ | Municipality × year-month × weekday | State
Row 14: Brazil reconstructed hospital-days | Resident-associated daily hospital-days; same-day temperature and fire PM₂.₅ | Health region × year-month × weekday | Two-way health region/date

### data/aggregate_results/display_tables/tableS5A.csv

SHA-256: 383c3fda3db1245348712382c2ea2b8e71c9ad5cf95bdd22323a684fc5a3b81e

Row 1: Outcome or comparison | Contrast / model | Estimate (95% CI) | p / q / sample size / notes
Row 2: Brazil urgent strict-N20 admissions | Heat-only vs neither, holiday + two-way model | RR 1.0509 (1.0259–1.0765) | 226,402 input events; 366,584 fitted rows; 450 regions and 1,826 dates
Row 3: Brazil strict-N20 admission-attributed hospital-days | Heat-only vs neither, holiday + two-way model | RR 1.0486 (1.0043–1.0948) | 805,746 input days; 356,593 fitted rows; p=0.0312; 356,593 fitted rows; 805,746 outcome total
Row 4: Brazil strict-N20 total reimbursement | Heat-only vs neither, holiday + two-way model | RR 1.0371 (1.0035–1.0719) | 194,858,997 BRL in December 2022 prices; p=0.0303; 366,584 fitted rows; 194,858,997 outcome total
Row 5: Brazil urgent strict-N20 admissions | Heat-only vs neither, excluding 2020–2021 | RR 1.0272 (0.9933–1.0623) | 136,094 input events; interval includes one; p=0.1168; 220,501 fitted rows; 136,094 total
Row 6: Brazil admission-attributed hospital-days | Heat-only vs neither, excluding 2020–2021 | RR 1.0428 (0.9874–1.1012) | 497,269 input days; interval includes one; p=0.1322; 215,201 fitted rows; 497,269 total
Row 7: Brazil strict-N20 total reimbursement | Heat-only vs neither, excluding 2020–2021 | RR 1.0179 (0.9720–1.0661) | 117,395,551 BRL in December 2022 prices; interval includes one; p=0.4500; 220,501 fitted rows; 117,395,551 total
Row 8: Brazil N20–N23 urgent-to-elective ratio | Heat-only urgent RR divided by elective RR, holiday + two-way | 1.0546 (1.0097–1.1015) | p=0.0167; 431,374 joint input events; p=0.0167; 658,844 fitted rows; 431,374 outcome total
Row 9: Brazil N20–N23 urgent-to-elective ratio | Same contrast, excluding 2020–2021 | 1.0655 (1.0067–1.1278) | p=0.0286
Row 10: Brazil strict-N20 urgent-to-elective ratio | Heat-only urgent RR divided by elective RR, holiday + two-way | 1.0409 (0.9924–1.0918) | p=0.0997; interval includes one; p=0.0997; 543,323 fitted rows; 321,201 outcome total
Row 11: Brazil strict-N20 urgent-to-elective ratio | Same contrast, excluding 2020–2021 | 1.0486 (0.9846–1.1168) | p=0.1395; interval includes one
Row 12: Brazil urgent strict-N20 admissions | Heat-only reclassified as neither; holiday + two-way model | +881 (parameter interval, 450–1,296) | Difference is 0.39% of the observed five-year input total; conditions on stratum totals
Row 13: Brazil admission-attributed hospital-days | Same heat-only reassignment | +3,004 (266–5,616) | Difference is 0.37% of the observed five-year input total
Row 14: Brazil strict-N20 reimbursement | Same heat-only reassignment; December 2022 BRL | +579,066 (55,561–1,094,968) | Difference is 0.30% of the observed five-year input total
Row 15: Brazil urgent strict-N20 admission-attributed hospital-days | Heat-only vs neither, health-region-clustered model | RR 1.0503 (1.0075–1.0948) | Health-region-clustered model; calendar-adjusted two-way results are listed separately
Row 16: Brazil strict-N20 reimbursement per admission | Conditional geometric mean ratio; December 2022 BRL | 0.9964 (0.9729–1.0205) | Positive-cost area-days; not patient-weighted arithmetic mean; p=0.7698; 100,805 fitted rows; 147,285 positive-cost days

### data/aggregate_results/display_tables/tableS5B.csv

SHA-256: 82b6e26474ba75508fae20b411258a74287296e722ad209c7e4f1bcc0bd1d18c

Row 1: Outcome or comparison | Contrast / model | Estimate (95% CI) | p / q / sample size / notes
Row 2: Brazil drainage-record admissions | Cumulative temperature lag 0–7 per 1 °C, holiday + two-way | RR 1.0119 (1.0046–1.0193) | 84,572 events; 180,241 fitted rows; 443 regions and 1,819 dates; p=0.0015; 180,241 fitted rows; 84,572 outcome total
Row 3: Brazil drainage-record admissions | Cumulative fire-PM₂.₅ lag 0–7 per 10 µg/m³, holiday + two-way | RR 1.0601 (0.9939–1.1307) | Interval includes one
Row 4: Brazil drainage-record admissions | Cumulative temperature lag 0–7 per 1 °C, excluding 2020–2021 | RR 1.0122 (1.0033–1.0212) | 49,201 events; 105,097 fitted rows; 437 regions and 1,088 dates; p=0.0072; 105,097 fitted rows; 49,201 events
Row 5: Brazil drainage-record admissions | Cumulative fire-PM₂.₅ lag 0–7 per 10 µg/m³, excluding 2020–2021 | RR 1.0534 (0.9686–1.1457) | Interval includes one
Row 6: Brazil drainage/stent record | Temperature, cumulative lag 0–7, health-region-clustered model | RR 1.0123 (1.0052–1.0194) | Health-region-clustered model; calendar-adjusted two-way results are listed separately; p=0.0007; q=0.0026 within 4 continuous contrasts; 180,241 rows; 84,572 admissions
Row 7: Brazil stone-extraction/lithotripsy record | Temperature, cumulative lag 0–7, health-region-clustered model | RR 1.0069 (0.9998–1.0141) | Health-region-clustered estimate; outcome records may overlap with drainage group; p=0.0579; q=0.0579 within 4 continuous contrasts; 158,873 rows; 71,684 admissions
Row 8: Mexico nephrostomy record | Temperature, cumulative lag 0–7 | RR 1.0134 (0.9382–1.0948) | 1,760 admissions; q=0.7490 (explicit PSD sensitivity) in six temperature/fire contrasts; p=0.7265; q=0.7490 in 6 contrasts; 5,325 rows; 1,760 events
Row 9: Mexico ureteral-catheterization record | Temperature, cumulative lag 0–7 | RR 0.9897 (0.9558–1.0249) | 3,313 admissions; q=0.7490 (explicit PSD sensitivity); p=0.5513; q=0.7490 in 6 contrasts; 9,061 rows; 3,313 events
Row 10: Mexico fragmentation/lithotripsy record | Temperature, cumulative lag 0–7 | RR 1.0085 (0.9559–1.0640) | 2,107 admissions; q=0.7490 (explicit PSD sensitivity); p=0.7490; q=0.7490 in 6 contrasts; 6,084 rows; 2,107 events
Row 11: Mexico TIPO=Q surgery | Same-day temperature per 1 °C | RR 1.0113 (0.9986–1.0242) | No formal comparison with no-TIPO=Q group; p=0.0820; 22,380 rows; 10,498 events
Row 12: Mexico no recorded TIPO=Q surgery | Same-day temperature per 1 °C | RR 1.0180 (1.0067–1.0294) | No formal comparison with TIPO=Q group; p=0.0017; 34,267 rows; 21,948 events
Row 13: Stone-extraction/lithotripsy record, urgent strict N20 [health-region-clustered model] | Temperature per 1 °C; cumulative lag 0–7 with national calendar control | RR 1.0062 (0.9991–1.0134) | Outcome events unchanged; p=0.0875

### data/aggregate_results/display_tables/tableS5C.csv

SHA-256: 0417511308ff54a2e4251161c966fe3cbcab4756a1265252b939c68f99ee6ea2

Row 1: Outcome or comparison | Contrast / model | Estimate (95% CI) | p / q / sample size / notes
Row 2: Mexico state compound-day group | Compound vs neither, same-day P90 heat and P75 fire thresholds | RR 0.996 (0.936–1.060) | Heat-only RR 0.983 (0.894–1.081); smoke-only RR 0.945 (0.897–0.996); interaction RR 1.072 (0.974–1.179)
Row 3: Brazil compound-day N20–N23 | Heat-only vs neither, two-way model | RR 1.0206 (1.0044–1.0370) | 450 health regions and 5,479 retained dates; exploratory categorical contrast
Row 4: Brazil compound-day N20–N23 | Smoke-only vs neither, two-way model | RR 1.0155 (1.0071–1.0239) | Same fitted model
Row 5: Brazil compound-day N20–N23 | Compound vs neither, two-way model | RR 1.0683 (1.0515–1.0853) | Same fitted model; no cross-country contrast estimated
Row 6: Brazil compound-day N20–N23 | Multiplicative heat-by-smoke interaction, two-way model | RR 1.0308 (1.0100–1.0521) | p=0.0037; contrast is compound − heat-only − smoke-only on the log scale
Row 7: Mexico state compound-day N20–N23 | Compound vs neither | RR 0.996 (0.936–1.060) | No formal cross-country contrast; municipality three-group model not estimable
Row 8: Brazil K85 acute pancreatitis | Heat-only vs neither | RR 0.9710 (0.9513–0.9911) | Health-region-clustered categorical comparison; exploratory
Row 9: Brazil K85 acute pancreatitis | Smoke-only vs neither | RR 0.9837 (0.9679–0.9998) | Health-region-clustered categorical comparison; exploratory
Row 10: Brazil K85 acute pancreatitis | Compound vs neither | RR 0.9610 (0.9285–0.9947) | Nominal p=0.0235; three-outcome BH q=0.0704; exploratory
Row 11: Brazil K81 cholecystitis | Heat-only vs neither | RR 0.9661 (0.9553–0.9771) | Health-region-clustered categorical comparison; exploratory
Row 12: Brazil K81 cholecystitis | Smoke-only vs neither | RR 1.0094 (1.0021–1.0167) | Health-region-clustered categorical comparison; exploratory
Row 13: Brazil K81 cholecystitis | Compound vs neither | RR 0.9987 (0.9840–1.0137) | Nominal p=0.8640; three-outcome BH q=0.8640; exploratory
Row 14: Brazil K80 gallstones | Heat-only vs neither | RR 0.9710 (0.9610–0.9811) | Health-region-clustered categorical comparison; exploratory
Row 15: Brazil K80 gallstones | Smoke-only vs neither | RR 1.0042 (0.9991–1.0094) | Health-region-clustered categorical comparison; exploratory
Row 16: Brazil K80 gallstones | Compound vs neither | RR 1.0082 (0.9983–1.0181) | Nominal p=0.1037; three-outcome BH q=0.1555; exploratory

### data/aggregate_results/display_tables/tableS5D.csv

SHA-256: b682a4fba956b0a770b5bf3866a46f8da709408cf0d74cb063bb1e26435502b3

Row 1: Outcome or comparison | Contrast / model | Estimate (95% CI) | p / q / sample size / notes
Row 2: Brazil ICU admission or death, urgent strict N20 | Heat-only vs neither | RR 0.971 (0.842–1.120) | 4,967 events; imprecise; not evidence of equivalence; p=0.6869; 20,463 model rows; 4,967 events; limited precision
Row 3: Brazil occupancy, inclusive dates | Resident-associated hospital-days; same-day temperature | 1,031,207 days; RR 1.0031 (1.0008–1.0053) per °C | Fire-PM₂.₅ RR 1.0205 (0.9977–1.0438) per 10 µg/m³; p=0.0075; p=0.0789; 619,735 model rows
Row 4: Brazil occupancy, discharge day excluded | Resident-associated hospital-days; same-day temperature | 804,968 days; RR 1.0034 (1.0011–1.0057) per °C | Fire-PM₂.₅ RR 1.0219 (0.9972–1.0472); p=0.0039; p=0.0823; 576,627 model rows
Row 5: Brazil bed-supply interaction | Heat by baseline tertile; temperature by tertile | p=0.8121; p=0.9730 | No evidence of modification under evaluated ecological measures
Row 6: Brazil monthly bed-supply interaction | Heat RR ratio per 100 beds/100,000 residents | 0.9846 (0.9541–1.0160) | p=0.3323
Row 7: Brazil cataract comparator | Continuous temperature per 1 °C | RR 0.9991 (0.9944–1.0038) | Negative-control comparison, not evidence against all confounding; p=0.697; 335,928 model rows
Row 8: Brazil fracture comparator | Continuous temperature per 1 °C | RR 1.0048 (1.0040–1.0056) | Heat-sensitive comparator; a heat–fall pathway is plausible; p<0.0001; 2,356,578 model rows

### data/aggregate_results/display_tables/tableS5E.csv

SHA-256: f6239e5d22493817302b25daa4c32582b212cdecdc879200c78889bb1e7f7c0a

Row 1: Outcome or comparison | Contrast / model | Estimate (95% CI) | p / q / sample size / notes
Row 2: Mexico state same-day sensitivity | NASA POWER temperature per 1 °C, jointly modeled with fire-PM₂.₅ | RR 1.0165 (1.0091–1.0240) | 43,671 model rows
Row 3: Mexico state same-day sensitivity | Fire-PM₂.₅ per 10 µg/m³, jointly modeled with temperature | RR 0.874 (0.814–0.939) | Below-one association is not interpreted as protective
Row 4: Mexico municipality sensitivity | Same-day state-capital NASA POWER temperature / municipality fire PM₂.₅ | RRs 1.0165 (1.0099–1.0231) / 0.881 (0.813–0.955) | Broad N20–N23; 133,775 model rows; not lag 0–7
Row 5: Mexico, N20–N23, state harmonized sensitivity | NASA POWER temperature, per 1 °C; lag 0–7 cumulative | 1.0264 (1.0148–1.0380) | Model rows: 43,483
Row 6: Brazil, strict N20, 2018–2022 baseline frame | Temperature, per 1 °C; lag 0–7 cumulative; two-way state/date clusters | 1.0125 (1.0059–1.0191) | Model rows: 45,436
Row 7: Mexico, strict N20, 2018–2022 baseline frame | ERA5-Land temperature, per 1 °C; lag 0–7 cumulative; two-way state/date clusters | 1.0265 (1.0127–1.0404) | Model rows: 39,359
Row 8: Mexico, broad N20–N23, source-expanded | ERA5-Land temperature, per 1 °C; lag 0–7 cumulative; calendar controls; two-way state/date clusters | 1.0231 (1.0098–1.0366) | Model rows: 43,503
Row 9: Brazil, strict N20, 2018–2022 | ERA5-Land land-normalized temperature, per 1 °C; cumulative lag 0–7; calendar controls; state/date clusters | 1.0115 (1.0052–1.0179) | Model rows: 45,436
Row 10: Brazil, regional N20–N23, 2008–2022 | ERA5-Land land-normalized temperature, per 1 °C; cumulative lag 0–7; joint fire terms; region/date clusters | 1.0143 (1.0107–1.0179) | Model rows: 1,481,520

### data/aggregate_results/display_tables/tableS6a.csv

SHA-256: d5d721e370dc26e44fbefe56d1a9ba4b19a29e46f22cce0de6ae61de08d7fd26

Row 1: Frame | Matched sample flow | Knot set in °C (P5, P33.33, P66.67, P95) | Raw covariance test | Explicit PSD sensitivity test
Row 2: Brazil strict N20, 2018–2022 | 49,302 source rows; 49,113 complete lag rows; 320,389 events; 856 zero-event strata removed 3,677 rows; 45,436 fitted rows; 27 state clusters and 1,819 date clusters; t df=26 | 17.91, 23.96, 26.21, 28.46 | F(16,26)=13.689; p=8.12×10⁻⁹ | F(16,26)=7.143; p=6.32×10⁻⁶
Row 3: Mexico strict N20, 2018–2022 admission dates, expanded discharge source | 58,432 source rows; 58,208 complete lag rows; 32,352 events; 4,372 zero-event strata removed 18,818 rows; 39,390 fitted rows; 32 states and 1,819 dates; t df=31 | 12.32, 18.63, 23.07, 28.41 | Not interpreted: restriction covariance was non-positive-semidefinite (minimum eigenvalue −1.37×10⁻⁷) | F(16,31)=0.649; p=0.819
Row 4: Brazil regional primary, 2008–2022 | 2,465,550 source rows; 2,462,400 complete lag rows; 1,130,926 events; 227,283 zero-event strata removed 980,880 rows; 1,481,520 fitted rows; 450 region clusters and 5,472 date clusters; t df=449 | 15.99, 22.28, 25.48, 28.57 | F(16,449)=4.292; p=7.71×10⁻⁸ | F(16,449)=4.292; p=7.71×10⁻⁸

### data/aggregate_results/display_tables/tableS6b.csv

SHA-256: be53cac017e4f1534e0d668031ce879bcb74ce01de98bbbd51d84e55a89a4522

Row 1: Frame | Contrast | Temperature reference (°C) | Raw RR (95% CI) | PSD sensitivity RR (95% CI)
Row 2: Brazil strict N20 | P10 vs P50 | 19.96 vs 25.35 | 0.930 (0.887–0.976) | 0.930 (0.886–0.977)
Row 3: Brazil strict N20 | P90 vs P50 | 27.75 vs 25.35 | 1.078 (1.033–1.125) | 1.078 (1.031–1.127)
Row 4: Brazil strict N20 | P95 vs P50 | 28.46 vs 25.35 | 1.107 (1.044–1.174) | 1.107 (1.042–1.177)
Row 5: Mexico strict N20 | P10 vs P50 | 13.94 vs 20.90 | 0.859 (0.777–0.949) | 0.859 (0.776–0.950)
Row 6: Mexico strict N20 | P90 vs P50 | 27.20 vs 20.90 | 1.075 (0.968–1.193) | 1.075 (0.962–1.201)
Row 7: Mexico strict N20 | P95 vs P50 | 28.41 vs 20.90 | 1.107 (0.976–1.255) | 1.107 (0.967–1.266)
Row 8: Brazil regional primary | P10 vs P50 | 17.99 vs 24.08 | 0.936 (0.915–0.957) | 0.936 (0.915–0.957)
Row 9: Brazil regional primary | P90 vs P50 | 27.66 vs 24.08 | 1.114 (1.085–1.143) | 1.114 (1.085–1.143)
Row 10: Brazil regional primary | P95 vs P50 | 28.57 vs 24.08 | 1.155 (1.114–1.196) | 1.155 (1.114–1.196)

## Figure source bindings

Every plotted estimate, confidence interval, curve, and displayed count is selected in data/figure_input_map.json and recorded on rebuild in artifacts/provenance/figure_sources.json. The saved-fit verification record is docs/figure_fit_binding_verification.json.
