"""Diagnose archived Mexico procedure covariance without re-fitting the means."""
from pathlib import Path
import json
import sys
sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/revision_20261113/covariance'
OUT.mkdir(parents=True, exist_ok=True)
source = ROOT / 'outputs/upgrade_20261112/mexico_boundary_full_covariance.json'
models = json.loads(source.read_text(encoding='utf-8'))['models']
diagnostics, contrasts, corrected = [], [], {}
for key, item in models.items():
    if 'temperature_fire_lag07_holiday' not in key:
        continue
    names = item['terms']
    b = np.asarray(item['coefficients'], dtype=float)
    original = np.asarray(item['vcov'], dtype=float)
    V = (original + original.T) / 2
    eigen, Q = np.linalg.eigh(V)
    psd = (Q * np.maximum(eigen, 0.0)) @ Q.T
    diagnostics.append(dict(model=key, endpoint=item['outcome'], dimension=len(b),
        minimum_eigenvalue=float(eigen.min()), maximum_eigenvalue=float(eigen.max()),
        relative_negative_eigenvalue=float(abs(min(eigen.min(), 0)) / max(abs(eigen).max(), 1e-300)),
        negative_eigenvalues=int((eigen < 0).sum()),
        relative_frobenius_correction=float(np.linalg.norm(psd-V)/np.linalg.norm(V)),
        state_clusters=item['state_clusters'], date_clusters=item['date_clusters'],
        t_df=item['t_df'], model_rows=item['model_rows'], events_after_lag=item['events_after_lag_window']))
    for exposure, prefix, unit in [('temperature','t','1 degree C'), ('fire_pm25','f','10 ug/m3')]:
        c = np.zeros(len(names))
        for lag in range(8):
            c[names.index(f'{prefix}{lag}')] = 1
        beta = float(c @ b)
        rawvar, psdvar = float(c @ V @ c), float(c @ psd @ c)
        for method, variance in [('original_CGM',rawvar),('eigenvalue_PSD_CGM',psdvar)]:
            se = np.sqrt(variance) if variance > 0 else float('nan')
            df = int(item['t_df'])
            crit = t.ppf(.975, df)
            contrasts.append(dict(model=key, endpoint=item['outcome'], exposure=exposure, unit=unit,
                covariance_method=method, beta=beta, variance=variance, se=se, RR=float(np.exp(beta)),
                CI_low=float(np.exp(beta-crit*se)), CI_high=float(np.exp(beta+crit*se)),
                p=float(2*t.sf(abs(beta/se),df)), t_df=df,
                variance_ratio_PSD_to_original=psdvar/rawvar if rawvar>0 else None))
    corrected[key] = dict(terms=names,coefficients=b.tolist(), original_vcov=V.tolist(), PSD_vcov=psd.tolist(),
        eigenvalues=eigen.tolist(), method='Symmetrize archived CGM covariance, replace all negative eigenvalues by zero; unchanged coefficient vector.',
        model_metadata={k:v for k,v in item.items() if k not in ['terms','coefficients','vcov']})
df = pd.DataFrame(contrasts)
for method in df.covariance_method.unique():
    idx = df.index[df.covariance_method.eq(method)]
    p = df.loc[idx,'p'].to_numpy()
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    q = np.empty(len(p)); q[order] = np.minimum(adjusted,1)
    df.loc[idx,'q_BH_six_contrasts'] = q
df.to_csv(OUT/'mexico_procedure_covariance_contrasts.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(diagnostics).to_csv(OUT/'mexico_procedure_covariance_diagnostics.csv',index=False,encoding='utf-8-sig')
(OUT/'mexico_procedure_PSD_covariance.json').write_text(json.dumps(corrected,indent=2),encoding='utf-8')
print(pd.DataFrame(diagnostics).to_string(index=False))
print(df[['endpoint','exposure','covariance_method','RR','CI_low','CI_high','p','q_BH_six_contrasts']].to_string(index=False))
