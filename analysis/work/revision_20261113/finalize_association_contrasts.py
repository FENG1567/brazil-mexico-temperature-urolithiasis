from pathlib import Path
import sys,json,csv
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/revision_20261113/covariance'
fit=json.loads((ROOT/'outputs/revision_20261113/nonlinear/fit_brazil_primary_same_day_compound_two_way.json').read_text(encoding='utf-8'))
names=fit['coefficient_names']
b=np.asarray(fit['coefficients']); V=np.asarray(fit['covariance'])
c=np.asarray([-1 if n in ['heat_only','smoke_only'] else 1 if n=='both' else 0 for n in names],float)
beta=float(c@b); variance=float(c@V@c); se=np.sqrt(variance)
df=449; critical=t.ppf(.975,df)
row=dict(model='brazil_primary_same_day_compound_two_way',estimand='multiplicative_heat_smoke_interaction',
    beta=beta,variance=variance,RR=float(np.exp(beta)),CI_low=float(np.exp(beta-critical*se)),CI_high=float(np.exp(beta+critical*se)),
    p=float(2*t.sf(abs(beta/se),df)),df=df,contrast='compound - heat_only - smoke_only on log scale',
    covariance='two-way health-region/date CRV1; full joint covariance',model_rows=1484094)
pd.DataFrame([row]).to_csv(OUT/'brazil_compound_interaction_two_way.csv',index=False,encoding='utf-8-sig')
old=ROOT/'work/frozen_snapshot/tables/biliary_comparators.csv'
comparators=pd.read_csv(old)
mask=comparators['group'].eq('composite'); ids=comparators.index[mask]
p=comparators.loc[ids,'p'].to_numpy(); order=np.argsort(p)
qsort=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
q=np.empty(len(p)); q[order]=np.minimum(qsort,1)
comparators['q_BH_three_compound_diagnoses']=np.nan
comparators.loc[ids,'q_BH_three_compound_diagnoses']=q
comparators['interpretation']='Exploratory broad administrative diagnoses; compound-day contrasts do not represent heat alone or uniform acute abdominal events.'
comparators.to_csv(OUT/'biliary_complete_categories_and_compound_multiplicity.csv',index=False,encoding='utf-8-sig')
print(json.dumps(row,indent=2)); print(comparators.loc[mask,['outcome','RR','p','q_BH_three_compound_diagnoses']].to_string(index=False))
