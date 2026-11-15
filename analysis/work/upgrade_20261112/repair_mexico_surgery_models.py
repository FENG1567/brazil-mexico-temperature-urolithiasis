"""Correct the annual source linkage in the two existing surgery subgroups.

Use corrected aggregate counts. Preserve the original same-day formula,
state clustering and 1.96-normal reporting convention for comparability.
"""
import os
for variable in ['OMP_NUM_THREADS', 'NUMBA_NUM_THREADS', 'OPENBLAS_NUM_THREADS']:
    os.environ[variable] = '8'
from pathlib import Path
import hashlib
import json
import sys
sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import norm

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / 'outputs/upgrade_20261112'
DAILY = OUT / 'mexico_specific_daily_counts.csv'
SOURCE = Path('F:/analysis/climate_biliary_20261102/work/mexico_composite_state_day_panel.parquet')
PRIOR = Path('F:/analysis/climate_biliary_20261102/outputs/route6_Brazil_Mexico_combined/tables/mexico_surgery_fixed_results.csv')


def sha(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


daily = pd.read_csv(DAILY, dtype={'state': str})
daily.state = daily.state.str.zfill(2)
daily.date = pd.to_datetime(daily.date).dt.normalize()
if not {'admissions', 'any_surgery'}.issubset(daily.columns):
    raise ValueError('Corrected surgical-record counts missing from source aggregate')
if int(daily.admissions.sum()) != 32417 or not daily.any_surgery.between(0, daily.admissions).all():
    raise ValueError('Corrected surgery endpoint counts inconsistent')
daily['n20_surg'] = daily.any_surgery.astype(int)
daily['n20_nosurg'] = daily.admissions.astype(int) - daily.n20_surg
panel = duckdb.connect().execute('SELECT * FROM read_parquet(?)', [str(SOURCE)]).df()
panel.date = pd.to_datetime(panel.date).dt.normalize()
panel.state = panel.state.astype(str).str.zfill(2)
panel = panel[panel.date.dt.year.between(2018, 2022)].copy()
if len(panel) != 58432 or panel.duplicated(['state','date']).any():
    raise ValueError('Original state panel calendar changed')
panel = panel.merge(daily[['state','date','admissions','n20_surg','n20_nosurg']],
                    on=['state','date'], how='left', validate='one_to_one')
for column in ['admissions','n20_surg','n20_nosurg']:
    panel[column] = panel[column].fillna(0).astype(int)
if not (panel.n20_surg + panel.n20_nosurg).eq(panel.admissions).all():
    raise ValueError('Surgery subgroups do not partition daily admissions')
panel['fire10'] = panel.fire_pm25 / 10
if panel[['tmean','fire10','stratum']].isna().any().any():
    raise ValueError('Original model covariate coverage changed')
old = pd.read_csv(PRIOR)
rows = []
models = {}
for outcome in ['n20_surg', 'n20_nosurg']:
    formula = f'{outcome} ~ tmean + fire10 | stratum'
    print('FIT_SURGERY_REPAIR', outcome, int(panel[outcome].sum()), flush=True)
    fitted = pf.fepois(formula, data=panel, vcov={'CRV1':'state'})
    names = list(fitted.coef().index)
    beta = fitted.coef().to_numpy()
    covariance = np.asarray(fitted._vcov)
    if not np.isfinite(covariance).all():
        raise ValueError('Surgery covariance invalid')
    models[outcome] = {'terms': names, 'coef': beta.tolist(), 'vcov': covariance.tolist(),
                       'model_rows': int(fitted._N), 'formula': formula, 'cluster': 'state',
                       'events_full_window': int(panel[outcome].sum())}
    for index, term in enumerate(names):
        coefficient = float(beta[index])
        se = float(np.sqrt(covariance[index,index]))
        prior = old[(old.outcome == outcome) & (old.term == term)]
        if len(prior) != 1:
            raise ValueError('Original frozen comparator row missing')
        prior = prior.iloc[0]
        rows.append({'outcome':outcome, 'term':term, 'RR':float(np.exp(coefficient)),
                     'CI_low':float(np.exp(coefficient - 1.96*se)),
                     'CI_high':float(np.exp(coefficient + 1.96*se)),
                     'p':float(2*norm.sf(abs(coefficient/se))), 'beta':coefficient, 'se':se,
                     'n':int(fitted._N), 'events_full_window':int(panel[outcome].sum()),
                     'frozen_prior_RR':float(prior.RR), 'frozen_prior_CI_low':float(prior.CI_low),
                     'frozen_prior_CI_high':float(prior.CI_high),
                     'exposure_definition':'same-day temperature per C; same-day fire PM2.5 per10',
                     'surgery_definition':'any recorded TIPO=Q during admission, source-file-year linkage',
                     'inference':'state CRV1; original 1.96 normal interval',
                     'interpretation':'group-specific association; between-group specificity requires a joint contrast'})
result = pd.DataFrame(rows)
result.to_csv(OUT/'mexico_surgery_linkage_repaired_results.csv', index=False, encoding='utf-8-sig')
(OUT/'mexico_surgery_linkage_repaired_covariance.json').write_text(
    json.dumps({'models':models, 'daily_sha256':sha(DAILY), 'original_model_panel_sha256':sha(SOURCE),
                'frozen_comparator_sha256':sha(PRIOR), 'linkage_repair':'ID plus source file year'},
               ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(result.to_string(index=False), flush=True)
print('MEXICO_SURGERY_LINKAGE_REPAIR_COMPLETE', flush=True)
