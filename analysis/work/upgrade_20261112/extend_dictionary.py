import sys
from pathlib import Path
import zipfile
import re
import hashlib
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import pandas as pd
out=Path(__file__).resolve().parents[2]/'outputs/upgrade_20261112'
dictionary=pd.read_csv(out/'sigtap_monthly_dictionary.csv',dtype=str)
manifest=pd.read_csv(out/'sigtap_source_manifest.csv',dtype=str)
rows=[]
sources=[]
for f in sorted(Path('D:/数据库集合/DATASUS数据库/02_reference/SIGTAP/2023').glob('*.zip')):
    month=re.search(r'(20\d{4})',f.name).group(1)
    with zipfile.ZipFile(f) as z:
        name=next(n for n in z.namelist() if n.endswith('tb_procedimento.txt'))
        for line in z.read(name).decode('latin1').splitlines():
            rows.append({'month':month,'code':line[:10],'name':line[10:260].strip()})
    sources.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
dictionary=pd.concat([dictionary,pd.DataFrame(rows)],ignore_index=True).drop_duplicates()
dictionary.to_csv(out/'sigtap_monthly_dictionary.csv',index=False,encoding='utf-8-sig')
pd.concat([manifest,pd.DataFrame(sources)],ignore_index=True).drop_duplicates('path').to_csv(
 out/'sigtap_source_manifest.csv',index=False,encoding='utf-8-sig')
print('dictionary_extended_for_delayed_claims',dictionary.month.nunique(),flush=True)
