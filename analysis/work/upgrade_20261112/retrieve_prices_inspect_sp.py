import sys
from pathlib import Path
import urllib.request
import json
import hashlib
from datetime import datetime,timezone
import shutil
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import pandas as pd
from pyreaddbc import dbc2dbf
from dbfread import DBF
base=Path(__file__).resolve().parents[2]
out=base/'outputs/upgrade_20261112'
work=base/'work/upgrade_20261112'
url='https://api.bcb.gov.br/dados/serie/bcdata.sgs.433/dados?formato=json&dataInicial=01/01/2018&dataFinal=31/12/2023'
try:
    with urllib.request.urlopen(url,timeout=40) as response:
        body=response.read()
    data=json.loads(body)
    frame=pd.DataFrame(data)
    frame['month']=pd.to_datetime(frame.data,format='%d/%m/%Y').dt.strftime('%Y%m')
    frame['ipca_pct']=pd.to_numeric(frame.valor)
    frame=frame.sort_values('month')
    assert len(frame)==72 and frame.month.nunique()==72
    assert frame.ipca_pct.between(-5,10).all()
    frame['index']=(1+frame.ipca_pct/100).cumprod()
    ref=frame.loc[frame.month=='202212','index'].iloc[0]
    frame['factor_to_2022_12']=ref/frame['index']
    frame.to_csv(out/'ipca_monthly.csv',index=False,encoding='utf-8-sig')
    (out/'ipca_api_response.json').write_bytes(body)
    (out/'ipca_download_provenance.json').write_text(json.dumps({'url':url,
        'definition_source':'https://www3.bcb.gov.br/sgspub/consultarvalores/consultarValoresSeries.do?method=consultarSeries&series=433',
        'series':'433 IPCA monthly percent','access_utc':datetime.now(timezone.utc).isoformat(),
        'rows':len(frame),'expected_rows':72,'sha256':hashlib.sha256(body).hexdigest(),
        'reference_prices':'December2022','authentication':'anonymous'},indent=2),encoding='utf-8')
    print('IPCA_DOWNLOADED',len(frame),flush=True)
except Exception as error:
    (out/'ipca_download_failure.txt').write_text(repr(error),encoding='utf-8')
    print('IPCA_DOWNLOAD_FAILED',repr(error),flush=True)

root=Path('D:/数据库集合/DATASUS数据库/01_raw/SIHSUS/SP')
files=[]
for year in range(2018,2024):
    files.extend(sorted((root/str(year)).rglob('*.dbc')))
meta={'files_2018_2023':len(files),'bytes':sum(f.stat().st_size for f in files)}
if files:
    f=files[0]
    src=work/'sp_sample.dbc'
    target=work/'sp_sample.dbf'
    shutil.copyfile(f,src)
    dbc2dbf(str(src),str(target))
    dbf=DBF(str(target),encoding='latin1')
    meta.update(sample_path=str(f),sample_rows=len(dbf),fields=dbf.field_names)
    src.unlink(missing_ok=True)
    target.unlink(missing_ok=True)
(out/'sp_feasibility.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
print('SP_METADATA',meta,flush=True)
