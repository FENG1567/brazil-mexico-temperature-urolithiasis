"""Predetermined observational Poisson analyses, never modifying frozen results."""
import os
os.environ['OMP_NUM_THREADS']='8'
os.environ['NUMBA_NUM_THREADS']='8'
os.environ['OPENBLAS_NUM_THREADS']='8'
import sys
import json
import time
from pathlib import Path
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy.stats import norm, chi2

BASE = Path(__file__).resolve().parents[2]
WORK = BASE/'work/upgrade_20261112'
OUT = BASE/'outputs/upgrade_20261112'
OLD = Path('F:/analysis/climate_biliary_20261102/work')
db=duckdb.connect()
db.execute('SET threads=8')
def read(path):
    return db.sql("SELECT * FROM read_parquet('"+Path(path).as_posix()+"')").df()

results=[]
def fit(label, outcome, terms, frame, fixed='stratum', vcov='hrid', contrasts=None):
    start=time.time()
    formula=outcome+' ~ '+' + '.join(terms)+' | '+fixed
    print('START',label,len(frame),formula,flush=True)
    m=pf.fepois(formula,data=frame,vcov={'CRV1':vcov})
    coef=m.coef()
    names=list(coef.index)
    covariance=np.asarray(m._vcov)
    if contrasts is None:
        contrasts={n:[n] for n in names}
    for name, selected in contrasts.items():
        weights=np.array([1.0 if n in selected else 0.0 for n in names])
        if not any(weights):
            raise ValueError('Contrast absent '+name)
        b=float(weights@coef.to_numpy())
        se=float(np.sqrt(weights@covariance@weights))
        row={'model':label,'outcome':outcome,'contrast':name,'RR':np.exp(b),
             'CI_low':np.exp(b-1.96*se),'CI_high':np.exp(b+1.96*se),
             'p':2*norm.sf(abs(b/se)),'beta':b,'se':se,'input_rows':len(frame),
             'model_rows':int(m._N),'events_input':float(frame[outcome].sum()),
             'regions_input':int(frame.hrid.nunique()),'formula':formula,'cluster':vcov}
        results.append(row)
        print('RESULT',row,flush=True)
    pd.DataFrame(results).to_csv(OUT/'new_model_results.csv',index=False,encoding='utf-8-sig')
    (OUT/(label+'_fit.json')).write_text(json.dumps({'terms':names,'coef':coef.tolist(),
        'vcov':covariance.tolist(),'N':int(m._N),'seconds':time.time()-start},indent=2),encoding='utf-8')
    return m

p=read(OLD/'panel_composite_humidity.parquet')
p=p[p.date.dt.year.between(2018,2022)].copy()
p['hrid']=p.hrid.astype(str)
p=p.sort_values(['hrid','date']).reset_index(drop=True)
if p.duplicated(['hrid','date']).any():
    raise ValueError('Duplicate panel days')
gaps=p.groupby('hrid').date.diff().dropna().dt.days
if not gaps.eq(1).all():
    raise ValueError('Panel date gaps before lag creation')
p['fire10']=p.fire_pm25/10
p['h_only']=((p.heat==1)&(p.smoke==0)).astype(int)
p['s_only']=((p.heat==0)&(p.smoke==1)).astype(int)
p['both']=((p.heat==1)&(p.smoke==1)).astype(int)
p['ym']=p.date.dt.strftime('%Y%m')
p['day_id']=p.date.dt.strftime('%Y%m%d')
p['stratum']=p.hrid+'_'+p.ym+'_'+p.date.dt.dayofweek.astype(str)
categorical=['h_only','s_only','both']
lag_terms=[]
for lag in range(8):
    for prefix,col in [('t','tmean'),('f','fire10')]:
        name=f'{prefix}{lag}'
        p[name]=p.groupby('hrid')[col].shift(lag)
        lag_terms.append(name)
contrasts={'heat_per_C_lag0_7':[f't{i}' for i in range(8)],
           'fire_per10_lag0_7':[f'f{i}' for i in range(8)]}

# Occupancy specifications reuse the exact existing aggregate for implementation checking.
o=read(OLD/'n20_occupied_beds_health_region_day.parquet')
o.hrid=o.hrid.astype(str)
o.date=pd.to_datetime(o.date)
q=p.merge(o,on=['hrid','date'],how='left',validate='one_to_one')
q['occupied']=q.occupied.fillna(0)
q=q.dropna(subset=['tmean','fire10'])
fit('occupancy_original_reproduction','occupied',['tmean','fire10'],q,fixed='hrid + ym + dow')
fit('occupancy_time_stratified','occupied',['tmean','fire10'],q)
fit('occupancy_two_way_cluster','occupied',['tmean','fire10'],q,vcov='hrid + day_id')

# Clinical code groups are fixed from the public monthly dictionary before association fitting.
code_groups={
 'drainage': ['0409010103','0409010162','0409010170','0409010278','0409010286','0409010294'],
 'stone_removal': ['0309030102','0309030110','0309030129','0309030137',
                  '0409010146','0409010154','0409010189','0409010227',
                  '0409010235','0409010316','0409010391','0409010561','0409010596'],
 'generic_surgery': ['0415010012','0415020034'],
}
dictionary=pd.read_csv(OUT/'sigtap_monthly_dictionary.csv',dtype=str)
whitelist=[]
for group,codes in code_groups.items():
    for code in codes:
        names=dictionary.loc[dictionary.code==code,'name'].unique().tolist()
        if not names:
            raise ValueError('Unverified procedure code '+code)
        whitelist.append({'group':group,'code':code,'names':' / '.join(names),
                          'months_present':dictionary.loc[dictionary.code==code,'month'].nunique()})
pd.DataFrame(whitelist).to_csv(OUT/'procedure_endpoint_whitelist.csv',index=False,encoding='utf-8-sig')
raw=read(WORK/'clinical_code_day_aggregate.parquet')
raw.date=pd.to_datetime(raw.date)
raw=raw[(raw.diagnosis=='N20')&(raw.CAR_INT=='02')&(raw.IDENT=='1')].copy()
raw['group']='other_recorded'
for group,codes in code_groups.items():
    raw.loc[raw.PROC_REA.isin(codes),'group']=group
valid_pairs=set(zip(dictionary.month,dictionary.code))
raw['dictionary_month_match']=[(month,code) in valid_pairs for month,code in zip(raw.competence,raw.PROC_REA)]
raw.groupby(['competence','group','dictionary_month_match'],dropna=False).admissions.sum().reset_index().to_csv(
 OUT/'procedure_mapping_coverage.csv',index=False,encoding='utf-8-sig')
flow={'urgent_n20_records':int(raw.admissions.sum()),
      'unmapped_region':int(raw.loc[raw.hrid.isna(),'admissions'].sum()),
      'unparsed_date':int(raw.loc[raw.date.isna(),'admissions'].sum()),
      'unmatched_dictionary_month':int(raw.loc[~raw.dictionary_month_match,'admissions'].sum())}
(OUT/'procedure_flow.json').write_text(json.dumps(flow,indent=2),encoding='utf-8')
print('Clinical flow',flow,flush=True)
if flow['unmapped_region'] or flow['unparsed_date']:
    raise ValueError('Resolve clinical region/date missingness first')
if flow['unmatched_dictionary_month']:
    raise ValueError('Resolve monthly procedure dictionary missingness before novel associations')
daily=raw.groupby(['hrid','date','group']).admissions.sum().unstack('group',fill_value=0).reset_index()
q=p.merge(daily,on=['hrid','date'],how='left',validate='one_to_one')
groups=list(code_groups)+['other_recorded']
for group in groups:
    q[group]=q[group].fillna(0)
q['urgent_n20']=q[groups].sum(axis=1)
counts=[]
for group in groups+['urgent_n20']:
    counts.append({'endpoint':group,'admissions':int(q[group].sum()),
                   'regions_with_events':int(q.loc[q[group]>0,'hrid'].nunique()),
                   'days_with_events':int((q[group]>0).sum())})
pd.DataFrame(counts).to_csv(OUT/'procedure_endpoint_counts.csv',index=False,encoding='utf-8-sig')
for outcome in ['urgent_n20','drainage','stone_removal','generic_surgery','other_recorded']:
    fit(outcome+'_categorical',outcome,categorical,q.dropna(subset=['fire_pm25','tmean']))
    fit(outcome+'_continuous',outcome,lag_terms,q.dropna(subset=lag_terms),contrasts=contrasts)

# Formal endpoint heterogeneity, allowing an independent baseline for every endpoint/stratum.
stack=[]
for group in groups:
    part=q[['hrid','date','stratum']+lag_terms].copy()
    part['count']=q[group].to_numpy()
    part['endpoint']=group
    part['endpoint_stratum']=part.stratum+'_'+group
    stack.append(part)
stack=pd.concat(stack,ignore_index=True).dropna(subset=lag_terms)
all_terms=[]
endpoint_contrasts={}
for group in groups:
    for term in lag_terms:
        name=group+'_'+term
        stack[name]=stack[term]*(stack.endpoint==group)
        all_terms.append(name)
    endpoint_contrasts[group+'_heat']=[group+f'_t{i}' for i in range(8)]
m=fit('procedure_joint_continuous','count',all_terms,stack,
      fixed='endpoint_stratum',contrasts=endpoint_contrasts)
names=list(m.coef().index)
b=m.coef().to_numpy()
V=np.asarray(m._vcov)
tests=[]
for first,second in [('drainage','other_recorded'),('stone_removal','other_recorded'),('drainage','stone_removal')]:
    z=np.array([float(n in endpoint_contrasts[first+'_heat'])-float(n in endpoint_contrasts[second+'_heat']) for n in names])
    delta=float(z@b)
    se=float(np.sqrt(z@V@z))
    tests.append({'comparison':first+' / '+second,'ratio_of_RRs':np.exp(delta),
                  'low':np.exp(delta-1.96*se),'high':np.exp(delta+1.96*se),'p':2*norm.sf(abs(delta/se))})
pd.DataFrame(tests).to_csv(OUT/'procedure_heterogeneity.csv',index=False,encoding='utf-8-sig')
print('COMPLETE',flush=True)
