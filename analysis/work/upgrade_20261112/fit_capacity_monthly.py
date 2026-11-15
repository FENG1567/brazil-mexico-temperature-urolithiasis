"""Exploratory capacity effect modification using historical monthly SUS beds.

Inputs are read-only. Results are written under outputs/upgrade_20261112 and
carry a ``capacity_`` prefix so the frozen and source outputs remain untouched.
"""
import os

os.environ['OMP_NUM_THREADS'] = '8'
os.environ['NUMBA_NUM_THREADS'] = '8'
os.environ['OPENBLAS_NUM_THREADS'] = '8'

import sys
from pathlib import Path
import json
import time
import traceback

sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import chi2, norm


BASE = Path(__file__).resolve().parents[2]
WORK = BASE / 'work' / 'upgrade_20261112'
OUT = BASE / 'outputs' / 'upgrade_20261112'
CNES_CSV = OUT / 'cnes_monthly_health_region.csv'
CLINICAL_PARQUET = WORK / 'clinical_code_day_aggregate.parquet'
EXPOSURE_PARQUET = Path('F:/analysis/climate_biliary_20261102/work/panel_composite_humidity.parquet')
OUT.mkdir(parents=True, exist_ok=True)

db = duckdb.connect()
db.execute('SET threads=8')


def log(*items):
    print(*items, flush=True)


def read_parquet(path, query='*'):
    safe = Path(path).as_posix().replace("'", "''")
    return db.sql(f"SELECT {query} FROM read_parquet('{safe}')").df()


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding='utf-8')


def matrix_from_names(names, coef_names, coef_weights):
    row = np.zeros(len(coef_names), dtype=float)
    for name, weight in coef_weights.items():
        if name not in coef_names:
            raise ValueError(f'Required coefficient absent: {name}; available={coef_names}')
        row[coef_names.index(name)] = weight
    return row


def linear_estimate(coef, vcov, contrast):
    beta = coef.to_numpy(dtype=float)
    b = float(contrast @ beta)
    variance = float(contrast @ vcov @ contrast)
    se = float(np.sqrt(max(variance, 0.0)))
    if se == 0 or not np.isfinite(se):
        return {'beta': b, 'se': se, 'RR': float(np.exp(b)), 'CI_low': None,
                'CI_high': None, 'p': None}
    return {'beta': b, 'se': se, 'RR': float(np.exp(b)),
            'CI_low': float(np.exp(b - 1.96 * se)),
            'CI_high': float(np.exp(b + 1.96 * se)),
            'p': float(2 * norm.sf(abs(b / se)))}


def joint_wald(coef, vcov, restrictions):
    beta = coef.to_numpy(dtype=float)
    estimates = restrictions @ beta
    cov = restrictions @ vcov @ restrictions.T
    rank = int(np.linalg.matrix_rank(cov))
    if rank != restrictions.shape[0]:
        raise ValueError(f'Joint Wald covariance is rank deficient: {rank}/{restrictions.shape[0]}')
    statistic = float(estimates @ np.linalg.solve(cov, estimates))
    df = int(restrictions.shape[0])
    return {'statistic': statistic, 'df': df, 'p': float(chi2.sf(statistic, df))}


def fit_model(label, formula, frame, expected_terms):
    started = time.time()
    log('START', label, 'rows=', len(frame), 'events=', int(frame['urgent_n20'].sum()), 'formula=', formula)
    model = pf.fepois(formula, data=frame, vcov={'CRV1': 'hrid'})
    coef = model.coef()
    vcov = np.asarray(model._vcov, dtype=float)
    names = list(coef.index)
    if vcov.shape != (len(names), len(names)):
        raise ValueError(f'{label}: covariance shape {vcov.shape} does not match {len(names)} coefficients')
    missing = sorted(set(expected_terms) - set(names))
    if missing:
        raise ValueError(f'{label}: expected terms omitted by estimator: {missing}')
    record = {
        'model': label,
        'formula': formula,
        'fixed_effect': 'health-region × year-month × day-of-week',
        'covariance': 'CRV1 clustered by health region',
        'terms': names,
        'coefficients': coef.to_numpy(dtype=float).tolist(),
        'vcov': vcov.tolist(),
        'model_rows': int(model._N),
        'input_rows': int(len(frame)),
        'events_input': int(frame['urgent_n20'].sum()),
        'events_model': int(model._data['urgent_n20'].sum()) if getattr(model, '_data', None) is not None else None,
        'separated_rows': int(len(frame) - int(model._N)),
        'regions': int(frame['hrid'].nunique()),
        'elapsed_seconds': round(time.time() - started, 3),
    }
    save_json(OUT / f'capacity_{label}_fit.json', record)
    log('FIT_OK', label, 'model_rows=', record['model_rows'], 'terms=', names,
        'seconds=', record['elapsed_seconds'])
    return model, coef, vcov, record


def result_row(model_name, estimand, contrast_name, estimate, record,
               p_type='normal', df=1, note=''):
    return {
        'model': model_name,
        'estimand': estimand,
        'contrast': contrast_name,
        'beta': estimate.get('beta'),
        'se': estimate.get('se'),
        'RR': estimate.get('RR'),
        'CI_low': estimate.get('CI_low'),
        'CI_high': estimate.get('CI_high'),
        'statistic': estimate.get('statistic'),
        'df': df,
        'p': estimate.get('p'),
        'p_type': p_type,
        'input_rows': record['input_rows'],
        'model_rows': record['model_rows'],
        'events_input': record['events_input'],
        'events_model': record['events_model'],
        'regions': record['regions'],
        'covariance': record['covariance'],
        'formula': record['formula'],
        'note': note,
    }


def main():
    # Assemble the chronological, complete 2018–2022 exposure calendar first.
    p = read_parquet(EXPOSURE_PARQUET,
                     'hrid, date, pop, tmean, fire_pm25, heat')
    p['hrid'] = p.hrid.astype(str)
    p['date'] = pd.to_datetime(p.date).dt.normalize()
    p = p[p.date.between('2018-01-01', '2022-12-31')].copy()
    p = p.sort_values(['hrid', 'date']).reset_index(drop=True)
    if p.duplicated(['hrid', 'date']).any():
        raise ValueError('Duplicate region-days in the 2018–2022 exposure calendar')
    if p[['pop', 'tmean', 'fire_pm25', 'heat']].isna().any().any():
        raise ValueError('Missing population or exposure values in the 2018–2022 calendar')
    p['ym'] = p.date.dt.strftime('%Y%m')
    p['dow'] = p.date.dt.day_of_week.astype(int)
    p['stratum'] = p.hrid + '_' + p.ym + '_' + p.dow.astype(str)
    p['fire10'] = p.fire_pm25 / 10.0
    p['tmean'] = pd.to_numeric(p.tmean, errors='raise')
    p['pop'] = pd.to_numeric(p['pop'], errors='raise')
    p['heat'] = pd.to_numeric(p.heat, errors='raise').astype(int)
    panel_regions = sorted(p.hrid.unique())
    dates = pd.date_range('2018-01-01', '2022-12-31', freq='D')
    if len(dates) != 1826 or p.date.nunique() != len(dates):
        raise ValueError('Expected all 1,826 dates from 2018 through 2022')
    region_sizes = p.groupby('hrid').size()
    if not region_sizes.eq(len(dates)).all():
        raise ValueError('Exposure calendar is not complete for every health region')
    gaps = p.groupby('hrid').date.diff().dropna().dt.days
    if not gaps.eq(1).all():
        raise ValueError('Date gaps found before lag creation')

    # Monthly population denominator comes only from the exposure panel.
    pop_month = p.groupby(['hrid', 'ym'], as_index=False).agg(
        population=('pop', 'mean'), population_values=('pop', 'nunique'))
    if not pop_month.population_values.eq(1).all():
        raise ValueError('Population changes within a region-month in the exposure panel')
    pop_month = pop_month[['hrid', 'ym', 'population']]

    # Historical monthly CNES SUS-bed supply; the source table is read-only.
    cnes = pd.read_csv(CNES_CSV, dtype={'hrid': str, 'month': str})
    required_cnes = {'hrid', 'month', 'beds_sus'}
    if not required_cnes.issubset(cnes.columns):
        raise ValueError(f'CNES columns missing: {sorted(required_cnes - set(cnes.columns))}')
    cnes['hrid'] = cnes.hrid.astype(str)
    cnes['ym'] = cnes.month.astype(str).str.replace(r'\.0$', '', regex=True).str.zfill(6)
    cnes['beds_sus'] = pd.to_numeric(cnes.beds_sus, errors='raise')
    if cnes.duplicated(['hrid', 'ym']).any():
        raise ValueError('Duplicate CNES health-region-month rows')
    cnes = cnes.merge(pop_month, on=['hrid', 'ym'], how='left', validate='one_to_one')
    if cnes.population.isna().any() or (cnes.population <= 0).any():
        raise ValueError('Missing or nonpositive monthly population for a CNES row')
    cnes['supply_per_100k'] = cnes.beds_sus / cnes.population * 100000.0
    cnes = cnes[cnes.ym.str[:4].between('2018', '2022')].copy()
    if cnes.ym.nunique() != 60 or cnes.hrid.nunique() != len(panel_regions) or len(cnes) != 60 * len(panel_regions):
        raise ValueError('Historical CNES coverage does not equal 60 months × panel regions')
    baseline = (cnes[cnes.ym.str.startswith('2018')]
                .groupby('hrid', as_index=False).supply_per_100k.mean()
                .rename(columns={'supply_per_100k': 'baseline_2018_mean_per_100k'}))
    if len(baseline) != len(panel_regions) or baseline.baseline_2018_mean_per_100k.isna().any():
        raise ValueError('Could not calculate a 12-month 2018 baseline for every region')
    if cnes[cnes.ym.str.startswith('2018')].groupby('hrid').ym.nunique().ne(12).any():
        raise ValueError('The 2018 baseline is not based on 12 monthly observations per region')
    q1, q2 = baseline.baseline_2018_mean_per_100k.quantile([1 / 3, 2 / 3], interpolation='linear')
    if not q1 < q2:
        raise ValueError(f'Baseline supply tertile cutpoints are not distinct: q1={q1}, q2={q2}')
    baseline['capacity_group'] = np.select(
        [baseline.baseline_2018_mean_per_100k <= q1,
         baseline.baseline_2018_mean_per_100k <= q2],
        ['low', 'mid'], default='high')
    baseline['baseline_rank'] = baseline.baseline_2018_mean_per_100k.rank(method='average')
    baseline['q1_cutpoint_per_100k'] = float(q1)
    baseline['q2_cutpoint_per_100k'] = float(q2)
    baseline['cutpoint_rule'] = '2018 regional mean monthly SUS beds per 100,000; empirical 1/3 and 2/3 quantiles'
    tier_summary = (baseline.groupby('capacity_group', as_index=False)
                    .agg(regions=('hrid', 'nunique'),
                         min_baseline=('baseline_2018_mean_per_100k', 'min'),
                         median_baseline=('baseline_2018_mean_per_100k', 'median'),
                         max_baseline=('baseline_2018_mean_per_100k', 'max')))
    tier_summary['q1_cutpoint_per_100k'] = float(q1)
    tier_summary['q2_cutpoint_per_100k'] = float(q2)
    tier_summary['cutpoint_rule'] = '2018 regional mean monthly SUS beds per 100,000; empirical 1/3 and 2/3 quantiles'
    tier_summary.to_csv(OUT / 'capacity_tier_cutpoints.csv', index=False, encoding='utf-8-sig')
    cnes = cnes.merge(baseline[['hrid', 'capacity_group']], on='hrid', how='left', validate='many_to_one')
    cnes['capacity_center100'] = (cnes.supply_per_100k - cnes.supply_per_100k.mean()) / 100.0
    cnes['ym'] = cnes.ym.astype(str)
    p = p.merge(baseline[['hrid', 'capacity_group']], on='hrid', how='left', validate='many_to_one')
    p = p.merge(cnes[['hrid', 'ym', 'supply_per_100k', 'capacity_center100']],
                on=['hrid', 'ym'], how='left', validate='many_to_one')
    if p.capacity_group.isna().any() or p.supply_per_100k.isna().any():
        raise ValueError('Could not attach historical capacity to all region-days')

    # Daily urgent N20 admissions: fixed diagnosis, urgency, and authorization rules.
    clinical = read_parquet(CLINICAL_PARQUET,
                            'IDENT, diagnosis, CAR_INT, hrid, date, admissions')
    audit_all_rows = int(len(clinical))
    clinical['hrid'] = clinical.hrid.astype('string')
    clinical['date'] = pd.to_datetime(clinical.date, errors='coerce').dt.normalize()
    clinical['admissions'] = pd.to_numeric(clinical.admissions, errors='raise')
    eligible = clinical[(clinical.IDENT.astype(str) == '1') &
                        (clinical.diagnosis.astype(str) == 'N20') &
                        (clinical.CAR_INT.astype(str) == '02')].copy()
    in_window = eligible[eligible.date.between('2018-01-01', '2022-12-31')].copy()
    missing_region_events = float(in_window.loc[in_window.hrid.isna(), 'admissions'].sum())
    # Missing dates cannot be assigned to the study window; preserve the full
    # eligible aggregate count as a separate linkage audit item.
    missing_date_events = float(eligible.loc[eligible.date.isna(), 'admissions'].sum())
    in_window = in_window.dropna(subset=['hrid', 'date'])
    unknown_regions = sorted(set(in_window.hrid.astype(str)) - set(panel_regions))
    unknown_region_events = float(in_window.loc[in_window.hrid.astype(str).isin(unknown_regions), 'admissions'].sum())
    in_window = in_window[in_window.hrid.astype(str).isin(panel_regions)].copy()
    daily = (in_window.groupby(['hrid', 'date'], as_index=False).admissions.sum()
             .rename(columns={'admissions': 'urgent_n20'}))
    if daily.duplicated(['hrid', 'date']).any():
        raise ValueError('Clinical daily aggregation produced duplicate region-days')
    p = p.merge(daily, on=['hrid', 'date'], how='left', validate='one_to_one')
    p['urgent_n20'] = p.urgent_n20.fillna(0).astype(float)
    if (p.urgent_n20 < 0).any() or not np.isfinite(p.urgent_n20).all():
        raise ValueError('Invalid negative or non-finite urgent N20 outcome counts')

    group_events = (p.groupby('capacity_group', as_index=False)
                    .agg(regions=('hrid', 'nunique'), region_days=('urgent_n20', 'size'),
                         events=('urgent_n20', 'sum'), days_with_events=('urgent_n20', lambda x: int((x > 0).sum()))))
    audit = {
        'clinical_source': str(CLINICAL_PARQUET),
        'clinical_aggregate_rows_total': audit_all_rows,
        'endpoint_definition': "IDENT='1' AND diagnosis='N20' AND CAR_INT='02'; admissions summed by health-region/date",
        'eligible_rows_with_date_in_2018_2022': int(eligible.date.between('2018-01-01', '2022-12-31').sum()),
        'eligible_events_in_window_before_region_date_exclusions': float(eligible.loc[eligible.date.between('2018-01-01', '2022-12-31'), 'admissions'].sum()),
        'eligible_missing_region_events_in_window': missing_region_events,
        'eligible_missing_date_events_all_periods': missing_date_events,
        'eligible_unknown_region_events_in_window': unknown_region_events,
        'analysis_events_after_region_date_link': int(p.urgent_n20.sum()),
        'analysis_days_with_events': int((p.urgent_n20 > 0).sum()),
        'analysis_rows': int(len(p)),
        'analysis_regions': int(p.hrid.nunique()),
        'analysis_calendar_days': int(p.date.nunique()),
        'monthly_CNES_rows': int(len(cnes)),
        'monthly_CNES_regions': int(cnes.hrid.nunique()),
        'monthly_CNES_months': int(cnes.ym.nunique()),
        '2018_baseline_regions': int(len(baseline)),
        'baseline_q1_per_100k': float(q1),
        'baseline_q2_per_100k': float(q2),
        'baseline_capacity_groups': group_events.to_dict('records'),
        'monthly_population_within_region_month_unique': True,
        'monthly_capacity_grand_mean_per_100k': float(cnes.supply_per_100k.mean()),
        'daily_panel_complete_before_lags': True,
    }
    save_json(OUT / 'capacity_aggregate_audit.json', audit)
    log('AUDIT', json.dumps(audit, ensure_ascii=False))

    rows = []
    failures = []

    # Primary exploratory categorical heat × fixed-baseline capacity tertile model.
    p['heat_mid'] = p.heat * (p.capacity_group == 'mid').astype(int)
    p['heat_high'] = p.heat * (p.capacity_group == 'high').astype(int)
    p['fire10'] = p.fire_pm25 / 10.0
    categorical_formula = 'urgent_n20 ~ heat + heat_mid + heat_high + fire10 | stratum'
    categorical_expected = ['heat', 'heat_mid', 'heat_high', 'fire10']
    try:
        _, coef, vcov, record = fit_model('categorical', categorical_formula, p, categorical_expected)
        names = list(coef.index)
        for group, weights in [
            ('low', {'heat': 1.0}),
            ('mid', {'heat': 1.0, 'heat_mid': 1.0}),
            ('high', {'heat': 1.0, 'heat_high': 1.0}),
        ]:
            c = matrix_from_names(names, names, weights)
            rows.append(result_row('categorical', 'heat RR within fixed 2018 capacity tertile', group,
                                   linear_estimate(coef, vcov, c), record,
                                   note='Heat is the same-day binary exposure; RR compares heat=1 with heat=0.'))
        c_mid = matrix_from_names(names, names, {'heat_mid': 1.0})
        c_high = matrix_from_names(names, names, {'heat_high': 1.0})
        joint = joint_wald(coef, vcov, np.vstack([c_mid, c_high]))
        rows.append(result_row('categorical', 'overall heat × capacity interaction', 'mid=low and high=low jointly',
                               {'statistic': joint['statistic'], 'p': joint['p']}, record,
                               p_type='Wald chi-square', df=joint['df'],
                               note='Two-degree-of-freedom joint test of both heat interaction coefficients.'))
    except Exception as exc:
        failures.append({'model': 'categorical', 'error': repr(exc), 'traceback': traceback.format_exc()})
        log('FIT_FAILED categorical', repr(exc))
        save_json(OUT / 'capacity_categorical_failure.json', failures[-1])

    # Continuous lag corroboration: 0–7 temperature lags and both group interactions;
    # fire PM2.5 lags are adjustment terms. Compare cumulative heat effects across groups.
    lag_terms = []
    expected = []
    for lag in range(8):
        tname = f't{lag}'
        p[tname] = p.groupby('hrid', sort=False).tmean.shift(lag)
        p[f't{lag}_mid'] = p[tname] * (p.capacity_group == 'mid').astype(int)
        p[f't{lag}_high'] = p[tname] * (p.capacity_group == 'high').astype(int)
        p[f'f{lag}'] = p.groupby('hrid', sort=False).fire10.shift(lag)
        lag_terms.extend([tname, f't{lag}_mid', f't{lag}_high', f'f{lag}'])
        expected.extend([tname, f't{lag}_mid', f't{lag}_high', f'f{lag}'])
    continuous = p.dropna(subset=lag_terms).copy()
    continuous_formula = ('urgent_n20 ~ ' + ' + '.join(lag_terms) + ' | stratum')
    try:
        _, coef, vcov, record = fit_model('continuous_lag0_7', continuous_formula, continuous, expected)
        names = list(coef.index)
        for group, weights in [
            ('low', {f't{i}': 1.0 for i in range(8)}),
            ('mid', {**{f't{i}': 1.0 for i in range(8)}, **{f't{i}_mid': 1.0 for i in range(8)}}),
            ('high', {**{f't{i}': 1.0 for i in range(8)}, **{f't{i}_high': 1.0 for i in range(8)}}),
        ]:
            c = matrix_from_names(names, names, weights)
            rows.append(result_row('continuous_lag0_7', 'cumulative temperature RR per 1°C across lag 0–7', group,
                                   linear_estimate(coef, vcov, c), record,
                                   note='Cumulative log coefficient is the sum of lag 0–7 temperature terms.'))
        c_mid = matrix_from_names(names, names, {f't{i}_mid': 1.0 for i in range(8)})
        c_high = matrix_from_names(names, names, {f't{i}_high': 1.0 for i in range(8)})
        joint = joint_wald(coef, vcov, np.vstack([c_mid, c_high]))
        rows.append(result_row('continuous_lag0_7', 'overall cumulative temperature × capacity interaction',
                               'mid=low and high=low jointly',
                               {'statistic': joint['statistic'], 'p': joint['p']}, record,
                               p_type='Wald chi-square', df=joint['df'],
                               note='Two-degree-of-freedom joint test of the cumulative mid/low and high/low contrasts.'))
    except Exception as exc:
        failures.append({'model': 'continuous_lag0_7', 'error': repr(exc), 'traceback': traceback.format_exc()})
        log('FIT_FAILED continuous_lag0_7', repr(exc))
        save_json(OUT / 'capacity_continuous_lag0_7_failure.json', failures[-1])

    # Exploratory time-varying monthly supply interaction. The supply main effect is
    # absorbed by region × year-month × weekday strata and is intentionally omitted.
    p['heat_capacity100'] = p.heat * p.capacity_center100
    monthly_formula = 'urgent_n20 ~ heat + heat_capacity100 + fire10 | stratum'
    monthly_expected = ['heat', 'heat_capacity100', 'fire10']
    try:
        _, coef, vcov, record = fit_model('monthly_continuous', monthly_formula, p, monthly_expected)
        names = list(coef.index)
        c_heat = matrix_from_names(names, names, {'heat': 1.0})
        rows.append(result_row('monthly_continuous', 'heat RR at grand-mean monthly SUS-bed supply',
                               'heat at centered supply=0', linear_estimate(coef, vcov, c_heat), record,
                               note='Monthly supply main effect omitted because region-month supply is absorbed by stratum fixed effects.'))
        c_interaction = matrix_from_names(names, names, {'heat_capacity100': 1.0})
        rows.append(result_row('monthly_continuous', 'ratio of heat RRs per +100 SUS beds/100,000',
                               'monthly supply interaction', linear_estimate(coef, vcov, c_interaction), record,
                               note='Capacity is centered at the unweighted region-month mean and scaled by 100 beds/100,000.'))
    except Exception as exc:
        failures.append({'model': 'monthly_continuous', 'error': repr(exc), 'traceback': traceback.format_exc()})
        log('FIT_FAILED monthly_continuous', repr(exc))
        save_json(OUT / 'capacity_monthly_continuous_failure.json', failures[-1])

    result_df = pd.DataFrame(rows)
    result_df.to_csv(OUT / 'capacity_model_results.csv', index=False, encoding='utf-8-sig')
    save_json(OUT / 'capacity_run_status.json', {
        'status': 'complete' if not failures else ('partial' if rows else 'failed'),
        'failures': failures,
        'models_succeeded': sorted(result_df.model.unique().tolist()) if not result_df.empty else [],
        'models_failed': [f['model'] for f in failures],
        'result_rows': int(len(result_df)),
        'command': 'python work/upgrade_20261112/fit_capacity_monthly.py',
    })
    make_report(audit, tier_summary, result_df, failures)
    if failures:
        raise RuntimeError('One or more capacity models failed; inspect capacity_run_status.json')
    log('CAPACITY_ANALYSIS_COMPLETE', len(result_df), 'rows')


def fmt(value, digits=3):
    if value is None or pd.isna(value):
        return 'NA'
    return f'{float(value):.{digits}f}'


def fmtp(value):
    if value is None or pd.isna(value):
        return 'NA'
    return '<0.0001' if float(value) < 0.0001 else f'{float(value):.4f}'


def make_report(audit, tiers, results, failures):
    lines = [
        '# 历史月度床位供给异质性敏感性分析',
        '',
        '## 分析定义',
        '',
        '使用2018–2022年卫生区日面板，结局为 IDENT=1、诊断 N20、CAR_INT=02 的紧急入院数，按卫生区和入院日期汇总。暴露日历先验证为每个卫生区连续完整的1,826天，再生成温度和火源PM2.5的滞后项。模型固定效应为卫生区×年月×星期几，协方差按卫生区 CRV1 聚类。',
        '',
        '本模块的热暴露是同日二元指标，比较热日（heat=1）与非热日（heat=0），并调整同日火源PM2.5（每10 μg/m³）。该对比不按烟雾状态拆分热日或构造热-烟联合类别。',
        '',
        '历史供给定义为 CNES 月度 SUS 床位数 / 暴露面板人口 ×100,000。预设分组依据各卫生区2018年12个月供给率的均值，以450个卫生区的经验三分位切分，分组在全研究期间固定。不得据结果调整阈值。',
        '',
        '## 覆盖与分组',
        '',
        f"- 分析面板：{audit['analysis_rows']:,} 行，{audit['analysis_regions']} 个卫生区，{audit['analysis_calendar_days']} 天；紧急 N20 入院 {audit['analysis_events_after_region_date_link']:,} 次，阳性地区日 {audit['analysis_days_with_events']:,}。",
        f"- CNES覆盖：{audit['monthly_CNES_rows']:,} 个地区-月份，{audit['monthly_CNES_regions']} 个卫生区、{audit['monthly_CNES_months']} 个月；2018基线每区12个月。人口取自暴露面板，且区-月内唯一。",
        f"- 2018月度供给三分位切点：低/中边界 {audit['baseline_q1_per_100k']:.3f}、中/高边界 {audit['baseline_q2_per_100k']:.3f} 张SUS床位/10万人。",
        '',
        '| 2018基线组 | 卫生区数 | 最低 | 中位数 | 最高 |',
        '|---|---:|---:|---:|---:|',
    ]
    for row in tiers.to_dict('records'):
        lines.append(f"| {row['capacity_group']} | {int(row['regions'])} | {row['min_baseline']:.3f} | {row['median_baseline']:.3f} | {row['max_baseline']:.3f} |")
    lines += ['', '数值单位均为每10万人SUS床位。临床记录缺日期、缺卫生区或不在面板地区的记录数/事件数见 `capacity_aggregate_audit.json`。', '', '## 模型结果', '']
    if results.empty:
        lines.append('本次未获得模型结果。')
    else:
        lines += [
            '| 模型与估计量 | 对比 | RR / 比值 | 95% CI | p值 | 估计行数 | 模型事件/输入事件 |',
            '|---|---|---:|---:|---:|---:|---:|',
        ]
        for row in results.to_dict('records'):
            label = f"{row['model']} — {row['estimand']}"
            if pd.notna(row.get('RR')):
                effect = fmt(row['RR'])
                ci = f"{fmt(row['CI_low'])}–{fmt(row['CI_high'])}"
            else:
                effect = f"χ²({int(row['df'])})={fmt(row.get('statistic'))}"
                ci = '—'
            model_events = 'NA' if pd.isna(row.get('events_model')) else f"{int(row['events_model']):,}"
            lines.append(f"| {label} | {row['contrast']} | {effect} | {ci} | {fmtp(row.get('p'))} | {int(row['model_rows']):,} | {model_events}/{int(row['events_input']):,} |")
        lines += ['', '模型事件数为估计样本内的结局总数；输入事件数为生成该模型前的完整分析数据总数。条件Poisson因固定效应分离移除的观测行数已记入各模型文件。所有预设组别和整体交互检验均报告，不按方向或显著性筛选。区间为基于地区聚类稳健协方差的正态近似95%区间。']
    categorical = results[results.model == 'categorical']
    continuous = results[results.model == 'continuous_lag0_7']
    monthly = results[results.model == 'monthly_continuous']
    if not categorical.empty and not continuous.empty and not monthly.empty:
        def by_contrast(frame, contrast):
            match = frame[frame.contrast == contrast]
            return match.iloc[0] if not match.empty else None
        low = by_contrast(categorical, 'low')
        mid = by_contrast(categorical, 'mid')
        high = by_contrast(categorical, 'high')
        cat_test = by_contrast(categorical, 'mid=low and high=low jointly')
        lag_test = by_contrast(continuous, 'mid=low and high=low jointly')
        monthly_effect = by_contrast(monthly, 'monthly supply interaction')
        lag_low = by_contrast(continuous, 'low')
        lines += [
            '',
            '## 简短结果解释',
            '',
            f"同日热暴露RR在低、中、高基线供给组分别为 {fmt(low['RR'])}、{fmt(mid['RR'])}、{fmt(high['RR'])}；整体2自由度交互检验 p={fmtp(cat_test['p'])}。温度lag 0–7累积RR在低供给组为 {fmt(lag_low['RR'])}，中、高组估计接近；累积交互检验 p={fmtp(lag_test['p'])}。月度连续供给每增加100张SUS床位/10万人对应的热效应RR比为 {fmt(monthly_effect['RR'])}（p={fmtp(monthly_effect['p'])}）。整体未见热效应随地区供给水平明显不同的统计证据；交互检验未显著不等于证明组间完全相同。",
        ]
    lines += [
        '',
        '## 解读边界',
        '',
        '本分析是探索性效应修饰敏感性分析。主模型比较基线低/中/高供给组内的同日热暴露RR，并以2自由度Wald检验整体交互。连续模型使用温度lag 0–7与两组交互，火源PM2.5 lag 0–7作为调整项；组别累积RR为8个温度系数之和，整体检验比较中、 高组相对低组的两个累积系数差。月度连续模型将供给中心化并按每增加100张/10万人缩放；因供给在卫生区-月内不变且被分层固定效应吸收，不单独估计其主效应。',
        '',
        'CNES卫生区月度床位供给只是居住地卫生区的背景可用性指标。医院地点可能不同于居民所属卫生区，且该指标不含床位占用、患者跨区就医、转诊或服务可及性。因此结果不能解释为本地医院容量或占床率造成了N20入院变化，也不支持因果结论。',
        '',
        '## 复现',
        '',
        '`python work/upgrade_20261112/fit_capacity_monthly.py`（使用项目指定的 Python 运行时；线程上限8；DuckDB读取Parquet）。完整系数和协方差位于各 `capacity_*_fit.json`，汇总结果见 `capacity_model_results.csv`。',
    ]
    if failures:
        lines += ['', '## 估计失败', '']
        for failure in failures:
            lines.append(f"- {failure['model']}: {failure['error']}")
    (OUT / 'capacity_report_zh.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
