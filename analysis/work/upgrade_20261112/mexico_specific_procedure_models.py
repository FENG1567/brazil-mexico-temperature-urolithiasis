"""Three exploratory, catalogue-defined Mexico procedure-record outcomes."""
import os
for name in ['OMP_NUM_THREADS','NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS']:
    os.environ[name]='8'
import sys
from pathlib import Path
import hashlib
import json
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
import pandas as pd
import numpy as np
import pyfixest as pf
from scipy.stats import t
import mexico_procedure_feasibility as feasibility
import mexico_strict_phenotype as country_inputs

BASE=Path(__file__).resolve().parents[2]
OUT=BASE/'outputs/upgrade_20261112'
GROUPS=['nephrostomy','ureteral_catheterization','urinary_stone_fragmentation_lithotripsy']
db=duckdb.connect()
db.execute('SET threads=8')
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as file:
        for block in iter(lambda:file.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()

OUT.mkdir(parents=True,exist_ok=True)
feasibility.preserve_pre_repair_outputs()
cohort,ids_by_source_year,audit=feasibility.load_cohort()
coverage,source_file_years,matches=feasibility.scan_procedure_tables(ids_by_source_year)
feasibility_result=feasibility.assemble_result(
    cohort,audit,coverage,source_file_years,matches)
feasibility.write_feasibility_outputs(feasibility_result)
flags=matches['any_recorded_matching_procedure']
clinical=cohort.copy()
clinical['date']=pd.to_datetime(clinical.date).dt.normalize()
clinical['ID']=clinical.ID.astype(str).str.strip()
if len(clinical)!=feasibility.EXPECTED_STRICT_COHORT_TOTAL \
        or clinical.duplicated(['ID','source_file_year']).any():
    raise ValueError('Specific-procedure clinical admission key/count mismatch')
clinical['state']=clinical.state_number.astype(int).astype(str).str.zfill(2)
keys=list(zip(clinical.ID,clinical.source_file_year.astype(int)))
for group in GROUPS:
    clinical[group]=[int(key in flags[group]) for key in keys]
clinical['any_surgery']=[int(key in matches['any_surgery']) for key in keys]
cohort_keys=set(keys)
for group in GROUPS:
    if not flags[group].issubset(cohort_keys):
        raise ValueError('Procedure flag reconciliation failed '+group)
if not matches['any_surgery'].issubset(cohort_keys):
    raise ValueError('TIPO=Q flag reconciliation failed')
daily=clinical.groupby(['state','date']).agg(
    admissions=('ID','size'),nephrostomy=('nephrostomy','sum'),
    ureteral_catheterization=('ureteral_catheterization','sum'),
    urinary_stone_fragmentation_lithotripsy=('urinary_stone_fragmentation_lithotripsy','sum'),
    any_surgery=('any_surgery','sum')).reset_index()
if int(daily.admissions.sum())!=len(clinical) \
        or int(daily.any_surgery.sum())!=len(matches['any_surgery']):
    raise ValueError('Daily admission or TIPO=Q aggregation did not reconcile')
for group in GROUPS:
    if int(daily[group].sum())!=len(flags[group]):
        raise ValueError('Daily procedure-group aggregation did not reconcile '+group)
daily.to_csv(OUT/'mexico_specific_daily_counts.csv',index=False,encoding='utf-8-sig')

panel,era5,_=country_inputs.load_mexico()
p=panel.merge(era5,on=['state','date'],validate='one_to_one')
fire=db.execute('SELECT state,date,fire_pm25 FROM read_parquet(?)',
                [str(country_inputs.MEXICO_PANEL_PATH)]).df()
fire.state=fire.state.astype(str).str.zfill(2)
fire.date=pd.to_datetime(fire.date).dt.normalize()
p=p.merge(fire,on=['state','date'],validate='one_to_one').merge(daily,on=['state','date'],how='left',validate='one_to_one')
for column in GROUPS+['admissions','any_surgery']:
    p[column]=p[column].fillna(0).astype(int)
if not p.admissions.eq(p.strict_N20).all() or p.fire_pm25.isna().any() \
        or int(p.any_surgery.sum())!=len(matches['any_surgery']):
    raise ValueError('Daily source reconciliation or fire coverage failed')
p=p.sort_values(['state','date']).reset_index(drop=True)
calendar=pd.read_csv(OUT/'holiday_calendar_dates.csv')
days=set(pd.to_datetime(calendar.loc[calendar.country=='MX','date']))
if len(days)!=38:
    raise ValueError('Expected final38 Mexico calendar control dates')
p['holiday']=p.date.isin(days).astype(int)
p['date_id']=p.date.dt.strftime('%Y%m%d')
p['fire10']=p.fire_pm25/10
terms=[]
for lag in range(8):
    for prefix,column in [('t','temperature'),('f','fire10')]:
        term=f'{prefix}{lag}'
        p[term]=p.groupby('state')[column].shift(lag)
        terms.append(term)
p=p.dropna(subset=terms).copy()
results=[]
models={}
for outcome in GROUPS:
    informative=p.groupby('stratum')[outcome].transform('sum')>0
    frame=p.loc[informative].copy()
    G_state=int(frame.state.nunique());G_date=int(frame.date_id.nunique())
    df=min(G_state,G_date)-1
    formula=outcome+' ~ '+' + '.join(terms+['holiday'])+' | stratum'
    print('START',outcome,'records',int(clinical[outcome].sum()),
          'model_rows',len(frame),'states',G_state,flush=True)
    model=pf.fepois(formula,data=frame,vcov={'CRV1':'state + date_id'})
    names=list(model.coef().index)
    beta=model.coef().to_numpy();V=np.asarray(model._vcov)
    if model._N!=len(frame) or not np.isfinite(V).all():
        raise ValueError('Additional separation or invalid covariance requires review')
    models[outcome]={'terms':names,'coef':beta.tolist(),'vcov':V.tolist(),'formula':formula,
                     'model_rows':int(model._N),'state_clusters':G_state,'date_clusters':G_date,'t_df':df}
    for prefix,label in [('t','temperature_per_C_lag0_7'),('f','fire_per10_lag0_7')]:
        selected={f'{prefix}{lag}' for lag in range(8)}
        if not selected.issubset(names):
            raise ValueError('Distributed lag term absent')
        z=np.array([float(term in selected) for term in names])
        b=float(z@beta);variance=float(z@V@z)
        if variance<=0:
            raise ValueError('Nonpositive contrast variance')
        se=float(np.sqrt(variance));critical=float(t.ppf(.975,df))
        row={'outcome':outcome,'contrast':label,'RR':float(np.exp(b)),
             'CI_low':float(np.exp(b-critical*se)),'CI_high':float(np.exp(b+critical*se)),
             'p':float(2*t.sf(abs(b/se),df)),'beta':b,'se':se,'t_df':df,
             'events_full_window':int(clinical[outcome].sum()),'events_after_lag_window':int(p[outcome].sum()),
             'input_rows_after_lag':len(p),'model_rows':int(model._N),
             'state_clusters':G_state,'date_clusters':G_date,'formula':formula,
             'cluster':'state + calendar date','inference':'t(min(G_state,G_date)-1) approximation',
             'status':'additional exploratory procedure-record outcome; all strict N20 admissions'}
        results.append(row)
        print('RESULT',row,flush=True)
    pd.DataFrame(results).to_csv(OUT/'mexico_specific_model_results.csv',index=False,encoding='utf-8-sig')
    (OUT/'mexico_specific_full_covariance.json').write_text(json.dumps(models,indent=2,allow_nan=False),encoding='utf-8')

frame=pd.DataFrame(results)
order=np.argsort(frame.p.to_numpy())
adjusted=np.minimum.accumulate((frame.p.to_numpy()[order]*len(frame)/np.arange(1,len(frame)+1))[::-1])[::-1]
q=np.empty(len(frame));q[order]=np.minimum(adjusted,1)
frame['BH_q_six_continuous_contrasts']=q
frame.to_csv(OUT/'mexico_specific_model_results.csv',index=False,encoding='utf-8-sig')
catalog=OUT/'source_mexico_cie9mc_catalog_20260310.csv'
book=pd.read_csv(catalog,dtype=str)
codes=set().union(*(feasibility.CODE_GROUPS[group]['codes'] for group in GROUPS))
chosen=book[book.CATALOG_KEY.isin(codes)]
if set(chosen.CATALOG_KEY)!=codes:
    raise ValueError('Official catalogue code coverage incomplete')
chosen.to_csv(OUT/'mexico_specific_codebook.csv',index=False,encoding='utf-8-sig')
provenance={
 'catalogue_resource':'https://www.datos.gob.mx/es/dataset/catalogo_cie_9_mc/resource/536cbef9-80bc-4d77-b5fe-ecadb9011f57/view/923cb431-8e5d-4782-ab1e-ae0d37852fd6',
 'catalogue_file':catalog.name,'catalogue_sha256':sha(catalog),'catalogue_rows':len(book),
 'catalogue_version':'Government resource metadata10March2026; recorded codes occur in2018–2022 tables',
 'field_guide':'https://riisp.insp.mx/nada/files_nesstar/SAEH_SSA/Referencias/NOTAS_TECNICAS/GII-B002-SAEH-02-02.pdf',
 'procedure_inputs':[{'file':str(path),'sha256':sha(path)} for path,_ in feasibility.PROCEDURE_FILES.values()],
 'clinical_file':str(feasibility.COHORT_PATH),'clinical_sha256':sha(feasibility.COHORT_PATH),
 'coverage':coverage,'cohort':audit,'linkage_repair_audit':feasibility_result['linkage_repair_audit'],
 'source_partition_reconciliation':audit['source_partition_reconciliation'],
 'groups':{g:{'label':feasibility.CODE_GROUPS[g]['label'],
 'codes':sorted(feasibility.CODE_GROUPS[g]['codes'])} for g in GROUPS},
 'date_definition':'Exposure indexed to admission date; source_file_year selects procedure partition; no verified operation performance date',
 'interpretation':'Presence of recorded procedure codes during an admission; procedure groups can overlap',
 'no_linked_nonblank_PROMED_record':coverage['cohort_admissions_with_no_nonblank_PROMED'],
 'any_surgery_definition':'Strict N20 admission with at least one linked procedure row whose TIPO equals Q, linked by ID + source_file_year',
 'missing_record_interpretation':'Absence of a linked target record does not establish absence of clinical treatment',
 'comparison_scope':'Mexico all strict N20 admissions; Brazil intervention models use urgent N20, with different code systems and ascertainment'}
(OUT/'mexico_specific_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
print('MEXICO_SPECIFIC_COMPLETE',len(models),'models',flush=True)
