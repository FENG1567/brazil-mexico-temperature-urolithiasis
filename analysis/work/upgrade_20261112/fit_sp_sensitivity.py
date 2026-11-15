"""SP-linked daily outcome sensitivity models; run only after full SP linkage is ready.

This script consumes the completed aggregate outcome link. It does not scan raw
SP files. Outputs use the ``sp_model_`` prefix and contain no person-level data.
"""
import os

os.environ['OMP_NUM_THREADS'] = '8'
os.environ['NUMBA_NUM_THREADS'] = '8'
os.environ['OPENBLAS_NUM_THREADS'] = '8'

import sys
from pathlib import Path
import hashlib
import json
import time
import traceback

sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')

import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import norm


BASE = Path(__file__).resolve().parents[2]
OUT = BASE / 'outputs' / 'upgrade_20261112'
SP_LINK = OUT / 'sp_linked_daily_outcomes.csv'
EXPOSURE_PARQUET = Path('F:/analysis/climate_biliary_20261102/work/panel_composite_humidity.parquet')
OUT.mkdir(parents=True, exist_ok=True)

DB = duckdb.connect()
DB.execute('SET threads=8')
OUTCOMES = ['any_drainage', 'any_removal', 'rd_drainage', 'rd_removal']
MODELED_OUTCOMES = ['any_drainage', 'any_removal']
MIN_ANY_EVENTS = {'any_drainage': 47918, 'any_removal': 35976}
REFERENCE_RD_URGENT_N20 = 226402
STUDY_DATES = pd.date_range('2018-01-01', '2022-12-31', freq='D')


def log(*items):
    print(*items, flush=True)


def read_parquet(path, columns):
    safe = Path(path).as_posix().replace("'", "''")
    return DB.sql(f"SELECT {columns} FROM read_parquet('{safe}')").df()


def read_linked_csv(path):
    safe = Path(path).as_posix().replace("'", "''")
    return DB.sql(
        "SELECT hrid, date, any_drainage, any_removal, rd_drainage, rd_removal "
        f"FROM read_csv_auto('{safe}', header=true, all_varchar=true)"
    ).df()


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str), encoding='utf-8')


def linear_estimate(coef, vcov, weights):
    names = list(coef.index)
    c = np.zeros(len(names), dtype=float)
    for term, weight in weights.items():
        if term not in names:
            raise ValueError(f'Required coefficient absent: {term}; available={names}')
        c[names.index(term)] = weight
    b = float(c @ coef.to_numpy(dtype=float))
    variance = float(c @ vcov @ c)
    se = float(np.sqrt(max(variance, 0.0)))
    if not np.isfinite(se) or se == 0:
        return {'beta': b, 'se': se, 'RR': float(np.exp(b)), 'CI_low': None,
                'CI_high': None, 'p': None}
    return {'beta': b, 'se': se, 'RR': float(np.exp(b)),
            'CI_low': float(np.exp(b - 1.96 * se)),
            'CI_high': float(np.exp(b + 1.96 * se)),
            'p': float(2 * norm.sf(abs(b / se)))}


def fit_model(outcome, spec, formula, frame, expected_terms):
    label = f'{outcome}_{spec}'
    started = time.time()
    log('START', label, 'input_rows=', len(frame), 'input_events=', int(frame[outcome].sum()), 'formula=', formula)
    model = pf.fepois(formula, data=frame, vcov={'CRV1': 'hrid'})
    coef = model.coef()
    vcov = np.asarray(model._vcov, dtype=float)
    terms = list(coef.index)
    if vcov.shape != (len(terms), len(terms)):
        raise ValueError(f'{label}: covariance shape {vcov.shape} does not match {len(terms)} coefficients')
    missing = sorted(set(expected_terms) - set(terms))
    if missing:
        raise ValueError(f'{label}: expected terms omitted by estimator: {missing}')
    model_events = int(model._data[outcome].sum()) if getattr(model, '_data', None) is not None else None
    record = {
        'model': label,
        'outcome': outcome,
        'formula': formula,
        'fixed_effect': 'health-region × year-month × day-of-week',
        'covariance': 'CRV1 clustered by health region',
        'terms': terms,
        'coefficients': coef.to_numpy(dtype=float).tolist(),
        'vcov': vcov.tolist(),
        'input_rows': int(len(frame)),
        'model_rows': int(model._N),
        'separated_rows': int(len(frame) - int(model._N)),
        'input_events': int(frame[outcome].sum()),
        'model_events': model_events,
        'regions_input': int(frame.hrid.nunique()),
        'elapsed_seconds': round(time.time() - started, 3),
    }
    save_json(OUT / f'sp_model_{label}_fit.json', record)
    log('FIT_OK', label, 'model_rows=', record['model_rows'], 'model_events=', model_events,
        'terms=', terms, 'seconds=', record['elapsed_seconds'])
    return coef, vcov, record


def result_row(record, estimand, contrast, estimate, p_type='normal', df=1, note=''):
    return {
        'model': record['model'],
        'outcome': record['outcome'],
        'estimand': estimand,
        'contrast': contrast,
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
        'input_events': record['input_events'],
        'model_events': record['model_events'],
        'regions': record['regions_input'],
        'covariance': record['covariance'],
        'formula': record['formula'],
        'note': note,
    }


def prepare_panel():
    # Exposure panel supplies the complete calendar and all exposure variables.
    p = read_parquet(EXPOSURE_PARQUET, 'hrid, date, tmean, fire_pm25, heat, smoke')
    p['hrid'] = p.hrid.astype(str)
    p['date'] = pd.to_datetime(p.date, errors='raise').dt.normalize()
    p = p[p.date.between('2018-01-01', '2022-12-31')].copy()
    p = p.sort_values(['hrid', 'date']).reset_index(drop=True)
    if p.duplicated(['hrid', 'date']).any():
        raise ValueError('Duplicate region-days in exposure calendar')
    if p[['tmean', 'fire_pm25', 'heat', 'smoke']].isna().any().any():
        raise ValueError('Missing exposure or heat/smoke flags in 2018–2022 calendar')
    if not np.isfinite(p[['tmean', 'fire_pm25', 'heat', 'smoke']].to_numpy(dtype=float)).all():
        raise ValueError('Non-finite exposure or heat/smoke flag in 2018–2022 calendar')
    if p.date.nunique() != 1826 or len(p) != 450 * 1826 or p.hrid.nunique() != 450:
        raise ValueError(f'Expected complete 450-region × 1,826-day exposure panel, got {len(p)} rows')
    if not p.groupby('hrid').size().eq(1826).all():
        raise ValueError('Exposure calendar is incomplete for one or more regions')
    if not p.groupby('hrid').date.diff().dropna().dt.days.eq(1).all():
        raise ValueError('Calendar has gaps before lag creation')
    p['ym'] = p.date.dt.strftime('%Y%m')
    p['dow'] = p.date.dt.day_of_week.astype(int)
    p['stratum'] = p.hrid + '_' + p.ym + '_' + p.dow.astype(str)
    p['fire10'] = pd.to_numeric(p.fire_pm25, errors='raise') / 10.0
    p['tmean'] = pd.to_numeric(p.tmean, errors='raise')
    p['heat'] = pd.to_numeric(p.heat, errors='raise').astype(int)
    p['smoke'] = pd.to_numeric(p.smoke, errors='raise').astype(int)
    if not set(p.heat.unique()).issubset({0, 1}) or not set(p.smoke.unique()).issubset({0, 1}):
        raise ValueError('Heat and smoke indicators must each be binary')
    p['h_only'] = ((p.heat == 1) & (p.smoke == 0)).astype(int)
    p['s_only'] = ((p.heat == 0) & (p.smoke == 1)).astype(int)
    p['both'] = ((p.heat == 1) & (p.smoke == 1)).astype(int)

    # The root-created, completed daily link is the only SP input. No raw SP scan.
    if not SP_LINK.exists():
        raise FileNotFoundError(f'Completed SP link is not present: {SP_LINK}; wait for full scan and link completion')
    linked = read_linked_csv(SP_LINK)
    linked['hrid'] = linked.hrid.astype('string').str.strip()
    date_text = linked.date.astype('string').str.strip().str.replace(r'\.0$', '', regex=True)
    linked['date'] = pd.to_datetime(date_text, format='%Y%m%d', errors='coerce')
    if linked[['hrid', 'date']].isna().any().any():
        raise ValueError('SP linked outcomes contain missing or invalid region/date keys')
    if linked.duplicated(['hrid', 'date']).any():
        raise ValueError('SP linked daily outcomes contain duplicate region-days')
    if not linked.date.between('2018-01-01', '2022-12-31').all():
        raise ValueError('SP linked outcomes contain dates outside 2018–2022')
    outside = linked.merge(p[['hrid', 'date']], on=['hrid', 'date'], how='left', indicator=True)
    if (outside._merge != 'both').any():
        raise ValueError('SP linked outcome region-days include keys outside the exposure panel')
    for outcome in OUTCOMES:
        linked[outcome] = pd.to_numeric(linked[outcome], errors='coerce')
        if linked[outcome].isna().any():
            raise ValueError(f'SP linked outcome contains missing/non-numeric values: {outcome}')
        if (linked[outcome] < 0).any():
            raise ValueError(f'SP linked outcome contains negative counts: {outcome}')
        if not np.allclose(linked[outcome], np.round(linked[outcome])):
            raise ValueError(f'SP linked outcome must be integer-valued: {outcome}')

    source_rows = int(len(linked))
    source_event_totals = {name: int(linked[name].sum()) for name in OUTCOMES}
    p = p.merge(linked[['hrid', 'date'] + OUTCOMES], on=['hrid', 'date'],
                how='left', validate='one_to_one')
    zero_fill_rows = int(p[OUTCOMES].isna().all(axis=1).sum())
    if p[OUTCOMES].isna().any().any():
        # A source region-day must have all four measures together or none.
        partial = p[OUTCOMES].isna().any(axis=1) & ~p[OUTCOMES].isna().all(axis=1)
        if partial.any():
            raise ValueError('Partially missing SP outcome row after daily link')
        p[OUTCOMES] = p[OUTCOMES].fillna(0)
    for outcome in OUTCOMES:
        p[outcome] = pd.to_numeric(p[outcome], errors='raise').astype(float)
    total_events = {name: int(p[name].sum()) for name in OUTCOMES}
    if total_events != source_event_totals:
        raise ValueError('Event totals changed when zero region-days were restored')
    for name, minimum in MIN_ANY_EVENTS.items():
        if total_events[name] < minimum:
            raise ValueError(f'{name} total {total_events[name]} is below required RD reference {minimum}')
    if len(p) != 821700 or p[['hrid', 'date']].duplicated().any():
        raise ValueError('Linked analysis panel must contain exactly 821,700 unique region-days')

    audit = {
        'analysis_type': 'additional observational sensitivity specified after the RD procedure results',
        'SP_input': str(SP_LINK),
        'SP_input_sha256': hashlib.sha256(SP_LINK.read_bytes()).hexdigest(),
        'SP_input_bytes': int(SP_LINK.stat().st_size),
        'SP_linked_rows': source_rows,
        'zero_region_days_restored': zero_fill_rows,
        'analysis_rows': int(len(p)),
        'analysis_regions': int(p.hrid.nunique()),
        'analysis_calendar_days': int(p.date.nunique()),
        'reference_urgent_N20_RD_unique_admissions': REFERENCE_RD_URGENT_N20,
        'SP_and_RD_endpoint_overlap': True,
        'endpoint_definition': 'any group is RD principal procedure OR any SP procedure claim; daily group outcomes can overlap',
        'event_totals_before_link': source_event_totals,
        'event_totals_analysis': total_events,
        'zero_outcome_region_days': int(((p.any_drainage == 0) & (p.any_removal == 0)).sum()),
        'positive_region_days': {name: int((p[name] > 0).sum()) for name in OUTCOMES},
        'all_zero_days_preserved': True,
        'continuous_calendar_checked_before_lags': True,
        'fixed_effect': 'health-region × year-month × day-of-week',
        'covariance': 'CRV1 clustered by health region',
    }
    save_json(OUT / 'sp_model_audit.json', audit)
    log('INPUT_AUDIT', json.dumps(audit, ensure_ascii=False))
    return p, audit


def run_models(p):
    results = []
    failures = []
    for lag in range(8):
        p[f't{lag}'] = p.groupby('hrid', sort=False).tmean.shift(lag)
        p[f'f{lag}'] = p.groupby('hrid', sort=False).fire10.shift(lag)
    lag_terms = [f't{i}' for i in range(8)] + [f'f{i}' for i in range(8)]
    continuous = p.dropna(subset=lag_terms).copy()
    if len(continuous) != 818550:
        raise ValueError(f'Expected 818,550 continuous-lag rows after 7-day warmup; got {len(continuous)}')
    continuous_formula_terms = lag_terms
    categorical_terms = ['h_only', 's_only', 'both']

    for outcome in MODELED_OUTCOMES:
        # Continuous temperature/fire lag model; cumulative contrasts use full covariance.
        label = f'{outcome}_continuous_lag0_7'
        formula = outcome + ' ~ ' + ' + '.join(continuous_formula_terms) + ' | stratum'
        try:
            coef, vcov, record = fit_model(outcome, 'continuous_lag0_7', formula,
                                           continuous, continuous_formula_terms)
            for estimand, weights in [
                ('cumulative temperature RR per 1°C across lag 0–7', {f't{i}': 1.0 for i in range(8)}),
                ('cumulative fire PM2.5 RR per 10 μg/m³ across lag 0–7', {f'f{i}': 1.0 for i in range(8)}),
            ]:
                estimate = linear_estimate(coef, vcov, weights)
                results.append(result_row(record, estimand, 'lag 0–7 cumulative sum', estimate,
                                          note='Point estimate and CI use the full fitted coefficient covariance.'))
        except Exception as exc:
            failure = {'model': label, 'error': repr(exc), 'traceback': traceback.format_exc()}
            failures.append(failure)
            save_json(OUT / f'sp_model_{label}_failure.json', failure)
            log('FIT_FAILED', label, repr(exc))

        # Three indicator contrasts: heat-only, smoke-only, and both, each vs neither.
        label = f'{outcome}_categorical'
        formula = outcome + ' ~ ' + ' + '.join(categorical_terms) + ' | stratum'
        try:
            coef, vcov, record = fit_model(outcome, 'categorical', formula, p, categorical_terms)
            for term, estimand in [
                ('h_only', 'heat-only vs neither'),
                ('s_only', 'smoke-only vs neither'),
                ('both', 'heat-and-smoke vs neither'),
            ]:
                estimate = linear_estimate(coef, vcov, {term: 1.0})
                results.append(result_row(record, estimand, term, estimate,
                                          note='Three-indicator model; each category is compared with neither.'))
        except Exception as exc:
            failure = {'model': label, 'error': repr(exc), 'traceback': traceback.format_exc()}
            failures.append(failure)
            save_json(OUT / f'sp_model_{label}_failure.json', failure)
            log('FIT_FAILED', label, repr(exc))

    result_df = pd.DataFrame(results)
    result_df.to_csv(OUT / 'sp_model_results.csv', index=False, encoding='utf-8-sig')
    status = {
        'status': 'complete' if not failures else ('partial' if results else 'failed'),
        'models_succeeded': sorted(result_df.model.unique().tolist()) if not result_df.empty else [],
        'models_failed': [item['model'] for item in failures],
        'failure_count': len(failures),
        'result_rows': int(len(result_df)),
        'instruction': 'Run only after root confirms the full 1,944-file SP scan/link is complete.',
    }
    save_json(OUT / 'sp_model_run_status.json', status)
    return result_df, failures


def fmt(value, digits=3):
    if value is None or pd.isna(value):
        return 'NA'
    return f'{float(value):.{digits}f}'


def fmtp(value):
    if value is None or pd.isna(value):
        return 'NA'
    return '<0.0001' if float(value) < 0.0001 else f'{float(value):.4f}'


def make_report(audit, results, failures):
    lines = [
        '# SP关联操作结局的追加敏感性分析',
        '',
        '## 设计与边界',
        '',
        '本模块是临床RD操作结果之后追加的观察性敏感性分析，使用已完成链接的地区-日汇总，不扫描原始SP文件。结局分别为any_drainage与any_removal：均包括RD主操作或任一SP操作记录，两个组可重叠；操作只在住院期间记录，不能据此确定操作发生日期。参考RD主操作紧急N20唯一住院数为226,402。',
        '',
        '同一临床/暴露地区-日面板包括所有0结局日。模型固定效应为卫生区×年月×星期几，协方差按卫生区CRV1聚类。连续模型纳入温度和火源PM2.5各自lag 0–7，累积估计以系数和及完整协方差计算。分类模型将heat-only、smoke-only、heat-and-smoke分别与neither比较。',
        '',
        'any类结果高于RD主操作事件数的覆盖审计、零值补齐和事件数见`sp_model_audit.json`。any_drainage与any_removal可能重叠；不据两类模型的单独p值或显著性差异主张操作类别特异性。',
        '',
        '## 覆盖与估计',
        '',
        f"- 完整面板：{audit['analysis_rows']:,}行，{audit['analysis_regions']}个卫生区，{audit['analysis_calendar_days']}天；SP链接源记录{audit['SP_linked_rows']:,}行，补回全零地区日{audit['zero_region_days_restored']:,}行。",
        f"- 结局总数：any_drainage {audit['event_totals_analysis']['any_drainage']:,}；any_removal {audit['event_totals_analysis']['any_removal']:,}；RD drainage {audit['event_totals_analysis']['rd_drainage']:,}；RD removal {audit['event_totals_analysis']['rd_removal']:,}。",
        '',
        '| 结局 | 模型 | 对比 | RR | 95% CI | p值 | 估计行数 | 模型事件/输入事件 |',
        '|---|---|---|---:|---:|---:|---:|---:|',
    ]
    for row in results.to_dict('records'):
        rr = fmt(row['RR'])
        ci = f"{fmt(row['CI_low'])}–{fmt(row['CI_high'])}"
        me = 'NA' if pd.isna(row.get('model_events')) else f"{int(row['model_events']):,}"
        lines.append(f"| {row['outcome']} | {row['model']} | {row['estimand']} | {rr} | {ci} | {fmtp(row['p'])} | {int(row['model_rows']):,} | {me}/{int(row['input_events']):,} |")
    lines += [
        '',
        '所有模型对比均列出，不按方向或显著性筛选。区间为地区聚类稳健协方差的正态近似95%区间。结果不比较any组间的显著性差异，也不作特异性或因果解释。',
        '',
        '复现命令：`python work/upgrade_20261112/fit_sp_sensitivity.py`；仅在完整SP来源扫描与链接完成、且收到执行确认后运行。完整系数和协方差保存在各`sp_model_*_fit.json`。',
    ]
    if failures:
        lines += ['', '## 估计失败', '']
        lines.extend(f"- {item['model']}: {item['error']}" for item in failures)
    (OUT / 'sp_model_report_zh.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    p, audit = prepare_panel()
    result_df, failures = run_models(p)
    make_report(audit, result_df, failures)
    if failures:
        raise RuntimeError('One or more SP sensitivity models failed; inspect sp_model_run_status.json')
    log('SP_SENSITIVITY_COMPLETE', len(result_df), 'result rows')


if __name__ == '__main__':
    main()
