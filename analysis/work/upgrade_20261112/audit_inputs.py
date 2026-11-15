"""Read-only input audit. Emit aggregate metadata and public procedure dictionary only."""
import sys
from pathlib import Path
import json
import zipfile
import hashlib
import re

PROJECT = Path('F:/analysis/climate_biliary_20261102')
sys.path.insert(0, str(PROJECT/'work/pylibs'))
import duckdb
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
OUT = BASE/'outputs/upgrade_20261112'
WORK = BASE/'work/upgrade_20261112'
OUT.mkdir(parents=True, exist_ok=True)
RAW = Path('D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/private/health_partitions/rd')
SIG = Path('D:/数据库集合/DATASUS数据库/02_reference/SIGTAP')
db = duckdb.connect()
db.execute('SET threads=8')
db.execute("SET memory_limit='12GB'")
db.execute("SET temp_directory='"+(WORK/'ducktmp').as_posix()+"'")
rows, manifest = [], []
for year in range(2018, 2024):
    archives = sorted((SIG/str(year)).glob('*.zip'))
    if len(archives) != 12:
        raise ValueError(f'{year} SIGTAP expected 12 archives, found {len(archives)}')
    for path in archives:
        month = re.search(r'(20\d{4})', path.name).group(1)
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith('tb_procedimento.txt'))
            for line in z.read(name).decode('latin1').splitlines():
                rows.append({'month': month, 'code': line[:10], 'name': line[10:260].strip()})
        manifest.append({'path': str(path), 'bytes': path.stat().st_size,
                         'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
mapping = pd.DataFrame(rows).drop_duplicates()
mapping.to_csv(OUT/'sigtap_monthly_dictionary.csv', index=False, encoding='utf-8-sig')
pd.DataFrame(manifest).to_csv(OUT/'sigtap_source_manifest.csv', index=False, encoding='utf-8-sig')
subset = mapping[mapping.name.str.contains('LITOT|NEFROST|DUPLO J|URETER|CALCULOSE|CALCULO ', case=False, regex=True)]
subset.groupby(['code', 'name']).agg(months=('month','nunique')).reset_index().to_csv(
    OUT/'sigtap_candidate_codes.csv', index=False, encoding='utf-8-sig')
print('SIGTAP months', mapping.month.nunique(), 'candidate codes', subset.code.nunique(), flush=True)
print(subset[['code','name']].drop_duplicates().to_string(index=False), flush=True)

source = (RAW/'*/*/*.parquet').as_posix()
date_audit = db.sql(f"""SELECT count(*) records,
 count(*) FILTER(WHERE try_strptime(replace(DT_INTER,'-',''),'%Y%m%d') IS NULL) invalid_admission_date
 FROM read_parquet('{source}') WHERE substr(DIAG_PRINC,1,3) IN ('N20','N21','N22','N23')
 AND ANO_CMPT BETWEEN '2018' AND '2023'""").df()
date_audit.to_csv(OUT/'source_date_audit_before_filter.csv',index=False,encoding='utf-8-sig')
if date_audit.invalid_admission_date.iloc[0] != 0:
    raise ValueError('Admission dates invalid before admission-year eligibility filtering')
db.execute(f"""CREATE TABLE stone AS SELECT N_AIH, IDENT, ANO_CMPT, MES_CMPT,
    DT_INTER, DT_SAIDA, DIAG_PRINC, MUNIC_RES, PROC_REA, CAR_INT,
    MORTE, UTI_MES_TO, UTI_INT_TO, DIAS_PERM, VAL_TOT, __SOURCE_SHA256,
    __ROW_NUMBER FROM read_parquet('{source}')
    WHERE substr(DIAG_PRINC,1,3) IN ('N20','N21','N22','N23')
    AND substr(replace(DT_INTER,'-',''),1,4) BETWEEN '2018' AND '2022'""")
audit = {}
queries = {
 'diagnosis_admission_type': "SELECT substr(DIAG_PRINC,1,3) diagnosis,CAR_INT,IDENT,count(*) n FROM stone GROUP BY ALL ORDER BY ALL",
 'duplicate_authorizations': "SELECT count(*) records,count(DISTINCT N_AIH) authorizations,count(*) FILTER(WHERE N_AIH IS NULL OR N_AIH='') missing_authorization FROM stone",
 'authorization_multiplicity': "SELECT n,count(*) authorizations FROM (SELECT N_AIH,count(*) n FROM stone GROUP BY N_AIH) GROUP BY n ORDER BY n",
 'duplicate_source_rows': "SELECT count(*) records,count(DISTINCT (__SOURCE_SHA256,__ROW_NUMBER)) source_rows FROM stone",
 'dates': "SELECT min(DT_INTER),max(DT_INTER),min(ANO_CMPT),max(ANO_CMPT) FROM stone",
 'top_urgent_n20_codes': "SELECT PROC_REA,count(*) n FROM stone WHERE substr(DIAG_PRINC,1,3)='N20' AND CAR_INT='02' GROUP BY 1 ORDER BY n DESC LIMIT 40",
 'duplicate_variation': "SELECT count(*) repeated_authorizations,sum(records) records,sum(date_variants>1) multiple_admission_dates,sum(proc_variants>1) multiple_procedures,sum(ident_variants>1) multiple_ident FROM (SELECT N_AIH,count(*) records,count(DISTINCT DT_INTER) date_variants,count(DISTINCT PROC_REA) proc_variants,count(DISTINCT IDENT) ident_variants FROM stone GROUP BY 1 HAVING count(*)>1)",
}
for label, query in queries.items():
    frame = db.sql(query).df()
    if label == 'duplicate_variation':
        # SUM over an empty set denotes zero repeated records in this audit.
        frame = frame.fillna(0)
    frame.to_csv(OUT/(label+'.csv'), index=False, encoding='utf-8-sig')
    audit[label] = frame.to_dict('records')
    print(label, audit[label], flush=True)
auth = audit['duplicate_authorizations'][0]
if auth['records'] != auth['authorizations'] or auth['missing_authorization']:
    raise ValueError('Authorization uniqueness gate failed')
origin = audit['duplicate_source_rows'][0]
if origin['records'] != origin['source_rows']:
    raise ValueError('Duplicate source rows gate failed')
(OUT/'input_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2, default=str, allow_nan=False), encoding='utf-8')

# Save deidentified daily code-level aggregates only, keeping all diagnosis/type/IDENT strata.
crosswalk = Path('D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/reference_step03/municipality_sih6_ibge7_health_region_2023_v2023_1.csv')
db.execute(f"""CREATE TABLE xw AS SELECT CAST(municipality_sih6 AS VARCHAR) municipality,
 CAST(health_region_2023 AS VARCHAR) hrid FROM read_csv_auto('{crosswalk.as_posix()}')""")
print('crosswalk', db.sql('SELECT count(*),count(DISTINCT municipality) FROM xw').fetchall(), flush=True)
crosswalk_counts = db.sql('SELECT count(*),count(DISTINCT municipality) FROM xw').fetchall()[0]
if crosswalk_counts[0] != crosswalk_counts[1]:
    raise ValueError('Crosswalk municipality uniqueness gate failed')
db.sql("""SELECT s.IDENT,substr(s.DIAG_PRINC,1,3) diagnosis,s.CAR_INT,s.PROC_REA,
 lpad(s.ANO_CMPT,4,'0')||lpad(s.MES_CMPT,2,'0') AS competence,x.hrid,
 try_strptime(replace(s.DT_INTER,'-',''),'%Y%m%d')::DATE date,
 count(*) admissions,sum(s.DIAS_PERM) beddays,sum(s.VAL_TOT) nominal_cost
 FROM stone s LEFT JOIN xw x ON s.MUNIC_RES=x.municipality GROUP BY ALL""").write_parquet(
 str(WORK/'clinical_code_day_aggregate.parquet'))
print('aggregate_saved', flush=True)
