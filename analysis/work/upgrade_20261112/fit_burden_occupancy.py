import os
for variable in ['OMP_NUM_THREADS','NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS']:
    os.environ[variable]='8'
import sys
from pathlib import Path
import json
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
import pandas as pd
import numpy as np
import pyfixest as pf
from scipy.stats import norm

base=Path(__file__).resolve().parents[2]
out=base/'outputs/upgrade_20261112'
work=base/'work/upgrade_20261112'
old=Path('F:/analysis/climate_biliary_20261102/work')
db=duckdb.connect()
db.execute('SET threads=8')
def read(path):
    return db.sql("SELECT * FROM read_parquet('"+Path(path).as_posix()+"')").df()
p=read(old/'panel_composite_humidity.parquet')
p=p[p.date.dt.year.between(2018,2022)].copy()
p.hrid=p.hrid.astype(str)
p['ym']=p.date.dt.strftime('%Y%m')
p['day_id']=p.date.dt.strftime('%Y%m%d')
p['stratum']=p.hrid+'_'+p.ym+'_'+p.date.dt.dayofweek.astype(str)
p['fire10']=p.fire_pm25/10
p['h_only']=((p.heat==1)&(p.smoke==0)).astype(int)
p['s_only']=((p.heat==0)&(p.smoke==1)).astype(int)
p['both']=((p.heat==1)&(p.smoke==1)).astype(int)
terms=['h_only','s_only','both']
results=[]
burdens=[]
def fit(label,outcome,frame,regressors,cluster='hrid'):
    print('START',label,flush=True)
    formula=outcome+' ~ '+' + '.join(regressors)+' | stratum'
    m=pf.fepois(formula,data=frame,vcov={'CRV1':cluster})
    c=m.coef()
    V=np.asarray(m._vcov)
    for name in c.index:
        b=float(c[name]); se=float(m.se()[name])
        row={'model':label,'outcome':outcome,'term':name,'RR':np.exp(b),
             'low':np.exp(b-1.96*se),'high':np.exp(b+1.96*se),'p':2*norm.sf(abs(b/se)),
             'input_rows':len(frame),'model_rows':int(m._N),'outcome_sum':float(frame[outcome].sum()),
             'formula':formula,'cluster':cluster}
        results.append(row)
        print('RESULT',row,flush=True)
    pd.DataFrame(results).to_csv(out/'burden_occupancy_model_results.csv',index=False,encoding='utf-8-sig')
    (out/(label+'_fit.json')).write_text(json.dumps({'terms':list(c.index),'coef':c.tolist(),'vcov':V.tolist(),'N':int(m._N)},indent=2),encoding='utf-8')
    return m

def absolute_contrast(label,outcome,frame,model):
    """Model-fitted heat-only day contrast conditional on observed stratum totals.

    Keep estimated stratum intercepts fixed in the exposure contrast. Recalculate
    the conditional-Poisson intercept at each coefficient draw for its uncertainty.
    Baseline totals are held fixed. These are parameter-uncertainty intervals.
    """
    frame=frame.copy()
    frame['category']=np.select([frame.h_only==1,frame.s_only==1,frame.both==1],[1,2,3],default=0)
    n=frame.groupby(['stratum','category']).size().unstack(fill_value=0).reindex(columns=range(4),fill_value=0)
    y=frame.groupby('stratum')[outcome].sum().reindex(n.index).to_numpy()
    counts=n.to_numpy(dtype=float)
    c=model.coef().reindex(terms).to_numpy()
    names=list(model.coef().index)
    idx=[names.index(t) for t in terms]
    V=np.asarray(model._vcov)[np.ix_(idx,idx)]
    rng=np.random.default_rng(20261112)
    draws=rng.multivariate_normal(c,V,size=4096)
    samples=np.vstack([c,draws])
    estimates=[]
    fitted_on_heat=[]
    for start in range(0,len(samples),64):
        b=samples[start:start+64]
        weights=np.exp(np.column_stack([np.zeros(len(b)),b])).T
        denominators=counts@weights
        expected_h=((y[:,None]*counts[:,1,None])*weights[1]/denominators).sum(axis=0)
        estimates.extend((expected_h*(1-np.exp(-b[:,0]))).tolist())
        fitted_on_heat.extend(expected_h.tolist())
    interval=np.quantile(estimates[1:],[.025,.975])
    row={'outcome':outcome,'model':label,'contrast':'heat-only days reassigned to neither, fitted stratum intercept held fixed',
         'observed_total':float(frame[outcome].sum()),'observed_heat_only':float(frame.loc[frame.h_only==1,outcome].sum()),
         'fitted_heat_only':fitted_on_heat[0],'model_associated_difference':estimates[0],
         'parameter_interval_low':interval[0],'parameter_interval_high':interval[1],
         'uncertainty':'4096 joint normal coefficient draws, regional cluster covariance, conditional on observed stratum totals',
         'interpretation':'associational model contrast, noncausal; excludes baseline outcome-count uncertainty'}
    burdens.append(row)
    pd.DataFrame(burdens).to_csv(out/'absolute_service_contrasts.csv',index=False,encoding='utf-8-sig')
    print('ABSOLUTE',row,flush=True)

clinical=read(work/'clinical_code_day_aggregate.parquet')
clinical.date=pd.to_datetime(clinical.date)
clinical=clinical[(clinical.IDENT=='1')&(clinical.diagnosis=='N20')&(clinical.CAR_INT=='02')].copy()
prices=pd.read_csv(out/'ipca_monthly.csv',dtype={'month':str})
factors=dict(zip(prices.month,prices.factor_to_2022_12))
clinical['price_factor']=clinical.competence.map(factors)
if clinical.price_factor.isna().any():
    raise ValueError('Price index missing for claim month')
clinical['reimbursement_2022brl']=clinical.nominal_cost*clinical.price_factor
daily=clinical.groupby(['hrid','date']).agg(urgent=('admissions','sum'),beddays=('beddays','sum'),
 nominal_reimbursement=('nominal_cost','sum'),reimbursement_2022brl=('reimbursement_2022brl','sum')).reset_index()
q=p.merge(daily,on=['hrid','date'],how='left',validate='one_to_one')
for col in ['urgent','beddays','nominal_reimbursement','reimbursement_2022brl']:
    q[col]=q[col].fillna(0)
    m=fit(col+'_service',col,q,terms)
    absolute_contrast(col+'_service',col,q,m)

for definition in ['inclusive','exclusive']:
    o=read(work/f'occupied_boundary_{definition}.parquet')
    o.date=pd.to_datetime(o.date)
    q=p.merge(o,on=['hrid','date'],how='left',validate='one_to_one')
    if q.occupied.isna().any() or (q.occupied<0).any():
        raise ValueError('Reconstructed occupancy invalid')
    fit('occupancy_boundary_'+definition,'occupied',q,['tmean','fire10'],cluster='hrid + day_id')
    print('OCCUPANCY_SUM',definition,q.occupied.sum(),flush=True)
print('BURDEN_COMPLETE',flush=True)
