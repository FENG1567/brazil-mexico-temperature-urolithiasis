"""Explicit ISO parsing of cached Mexican admission calendar dates."""
from pathlib import Path
import hashlib
import json
import sys
sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / 'outputs/upgrade_20261112'
ROOT = Path('F:/analysis/climate_biliary_20261102/work')
checks = []
for name in ['mexico_n20_id_proced.parquet', 'mexico_n20_proced.parquet']:
    path = ROOT / name
    frame = duckdb.connect().execute('SELECT INGRE,date FROM read_parquet(?)', [str(path)]).df()
    raw = frame.INGRE.astype('string').str.strip()
    iso = raw.str.fullmatch(r'\d{4}-\d{2}-\d{2}(?:[ T].*)?', na=False)
    explicit = pd.to_datetime(raw.str[:10], format='%Y-%m-%d', errors='coerce')
    cached = pd.to_datetime(frame.date, errors='coerce').dt.normalize()
    differences = int((explicit != cached).sum())
    with path.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    counts = explicit.dt.year.value_counts(dropna=False).sort_index()
    checks.append({'file': str(path), 'sha256': digest, 'records': len(frame),
                   'ISO_format_records': int(iso.sum()), 'raw_date_parse_failures': int(explicit.isna().sum()),
                   'cached_date_missing': int(cached.isna().sum()), 'explicit_vs_cached_date_mismatches': differences,
                   'admission_year_counts': {str(int(year)): int(count) for year,count in counts.items() if pd.notna(year)},
                   'passed': bool(iso.all() and explicit.notna().all() and cached.notna().all() and differences == 0)})
report = {'passed': all(check['passed'] for check in checks),
          'date_semantics': 'Explicit YYYY-MM-DD admission calendar date versus normalized cached date',
          'checks': checks, 'individual_records_emitted': False}
(OUT / 'mexico_admission_date_audit.json').write_text(
    json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
print(json.dumps(report), flush=True)
if not report['passed']:
    raise ValueError('Admission date audit failed')
