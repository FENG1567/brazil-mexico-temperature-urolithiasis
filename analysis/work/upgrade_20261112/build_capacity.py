"""Validate historical monthly CNES LT and aggregate SUS beds to health region."""
import sys
from pathlib import Path
import hashlib
import json
import shutil
from collections import defaultdict
sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')
from pyreaddbc import dbc2dbf
from dbfread import DBF
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
OUT = BASE/'outputs/upgrade_20261112'
WORK = BASE/'work/upgrade_20261112'
RAW = Path('D:/数据库集合/DATASUS数据库/01_raw/CNES/LT')
XW = Path('D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/reference_step03/municipality_sih6_ibge7_health_region_2023_v2023_1.csv')
x = pd.read_csv(XW, dtype=str)
assert not x.municipality_sih6.duplicated().any()
cross = dict(zip(x.municipality_sih6, x.health_region_2023))
totals = defaultdict(lambda: [0.0,0.0,0])
man = []
for year in range(2018,2023):
    files = sorted((RAW/str(year)).rglob('LT*.dbc'))
    if len(files) != 324:
        raise ValueError(f'{year}: expected324 files, found{len(files)}')
    for i, f in enumerate(files):
        target = WORK/'cnes_tmp.dbf'
        entry = {'path': str(f), 'bytes': f.stat().st_size,
                 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()}
        try:
            file_totals = defaultdict(lambda: [0.0,0.0,0])
            ascii_source = WORK/'cnes_tmp.dbc'
            shutil.copyfile(f,ascii_source)
            dbc2dbf(str(ascii_source),str(target))
            table = DBF(str(target), encoding='latin1')
            if not {'CODUFMUN','COMPETEN','QT_EXIST','QT_SUS'}.issubset(table.field_names):
                raise ValueError('Required LT fields missing')
            if not man:
                print('CNES fields',table.field_names, flush=True)
            n, unmapped, months, duplicate = 0,0,set(),0
            seen = set()
            for r in table:
                n += 1
                month = str(r['COMPETEN']).strip()
                months.add(month)
                key = tuple(str(r.get(c,'')) for c in table.field_names)
                if key in seen:
                    duplicate += 1
                    continue
                seen.add(key)
                municipality = str(r['CODUFMUN']).strip().zfill(6)
                hrid = cross.get(municipality)
                if hrid is None:
                    unmapped += 1
                    continue
                b = float(r['QT_EXIST'] or 0)
                s = float(r['QT_SUS'] or 0)
                if b<0 or s<0 or s>b:
                    raise ValueError('Invalid SUS/total bed counts')
                total = file_totals[(hrid,month)]
                total[0] += b
                total[1] += s
                total[2] += 1
            expected = '20'+f.stem[-4:]
            if months != {expected}:
                raise ValueError(f'Competence mismatch {months}, expected {expected}')
            for key, values in file_totals.items():
                for position, value in enumerate(values):
                    totals[key][position] += value
            entry.update(status='ok',records=n,unmapped_rows=unmapped,duplicate_rows=duplicate,month=expected)
        except Exception as error:
            entry.update(status='failed',error=str(error))
        finally:
            (WORK/'cnes_tmp.dbc').unlink(missing_ok=True)
            target.unlink(missing_ok=True)
        man.append(entry)
        if (i+1)%81 == 0:
            print(year,i+1,'of',len(files),'failed',sum(r['status']=='failed' for r in man),flush=True)
    pd.DataFrame(man).to_csv(OUT/'cnes_source_manifest.csv', index=False,encoding='utf-8-sig')
if any(entry['status'] != 'ok' for entry in man):
    raise ValueError('CNES source failures recorded. Final capacity panel publication stopped.')
frame = pd.DataFrame([{'hrid':h,'month':m,'beds':v[0],'beds_sus':v[1],'bed_type_records':v[2]} for (h,m),v in totals.items()])
frame.to_csv(OUT/'cnes_monthly_health_region.csv',index=False,encoding='utf-8-sig')
summary = {'files':len(man),'failed':sum(r['status']=='failed' for r in man),
           'months':frame.month.nunique(),'regions':frame.hrid.nunique(),
           'unmapped_rows':sum(r.get('unmapped_rows',0) for r in man),
           'duplicate_rows':sum(r.get('duplicate_rows',0) for r in man),
           'rows':len(frame)}
(OUT/'cnes_audit.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(summary, flush=True)
