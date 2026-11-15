"""Scan local SIH-SP with bounded memory. Persist only deidentified daily counts.

Authorization numbers are used for local in-memory linkage and never emitted
as linkage outputs. Full source copies and decompressed DBFs are staged locally
and deleted after each file. Professional identifier fields are never extracted
into arrays or reported.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys
from pathlib import Path
import json
import struct
import mmap
import shutil
import time
import hashlib
import argparse
import pickle
import subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
import threading
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import numpy as np
from pyreaddbc import dbc2dbf

BASE=Path(__file__).resolve().parents[2]
WORK=BASE/'work/upgrade_20261112'
OUT=BASE/'outputs/upgrade_20261112'
DRAIN=['0409010103','0409010162','0409010170','0409010278','0409010286','0409010294']
REMOVE=['0309030102','0309030110','0309030129','0309030137','0409010146','0409010154',
        '0409010189','0409010227','0409010235','0409010316','0409010391','0409010561','0409010596']
CODE_BITS={c.encode():1 for c in DRAIN}|{c.encode():2 for c in REMOVE}
ELIGIBLE=set()
DECOMPRESS_LOCK=threading.Lock()
def init(keys):
    global ELIGIBLE
    ELIGIBLE=keys

def scan(item):
    index,filename=item
    start=time.time()
    f=Path(filename)
    staging=WORK/'sp_staging'
    staging.mkdir(exist_ok=True)
    source=staging/f'source_{index}.dbc'
    target=staging/f'source_{index}.dbf'
    entry={'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
    linked={}
    try:
        shutil.copyfile(f,source)
        # Serialize the native decompressor. Array scanning and file I/O can overlap.
        with DECOMPRESS_LOCK:
            dbc2dbf(str(source),str(target))
        with target.open('rb') as file:
            header=file.read(32)
            n=struct.unpack_from('<I',header,4)[0]
            hlen,rlen=struct.unpack_from('<HH',header,8)
            layout={}
            offset=1
            while True:
                descriptor=file.read(32)
                if descriptor[0]==13:
                    break
                name=descriptor[:11].split(b'\0')[0].decode('ascii')
                length=descriptor[16]
                layout[name]=(offset,length)
                offset+=length
            if offset!=rlen:
                raise ValueError('DBF layout/record-length mismatch')
            for col in ['SP_NAIH','SP_ATOPROF','SP_CIDPRI']:
                if col not in layout:
                    raise ValueError('SP required field missing '+col)
            if target.stat().st_size<hlen+n*rlen:
                raise ValueError('Truncated DBF')
            with mmap.mmap(file.fileno(),0,access=mmap.ACCESS_READ) as data:
                co,cl=layout['SP_ATOPROF']; ao,al=layout['SP_NAIH']; do,dl=layout['SP_CIDPRI']
                codes=np.ndarray((n,),dtype=f'S{cl}',buffer=data,offset=hlen+co,strides=(rlen,))
                diag=np.ndarray((n,),dtype='S3',buffer=data,offset=hlen+do,strides=(rlen,))
                deleted=np.ndarray((n,),dtype='S1',buffer=data,offset=hlen,strides=(rlen,))
                selected=np.flatnonzero(((diag==b'N20')|np.isin(codes,[c.ljust(cl) for c in CODE_BITS]))&(deleted!=b'*'))
                authors=np.ndarray((n,),dtype=f'S{al}',buffer=data,offset=hlen+ao,strides=(rlen,))
                for a,c in zip(authors[selected],codes[selected]):
                    key=bytes(a).strip()
                    if key in ELIGIBLE:
                        linked[key]=linked.get(key,0)|CODE_BITS.get(bytes(c).strip(),0)
                del authors,codes,diag,deleted
        entry.update(status='ok',records=n,candidate_rows=len(selected),linked_authorizations=len(linked),seconds=time.time()-start)
    finally:
        source.unlink(missing_ok=True)
        target.unlink(missing_ok=True)
    return linked,entry

def main():
    import pandas as pd
    import duckdb
    parser=argparse.ArgumentParser()
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--backend',choices=['threads','subprocess'],default='subprocess')
    args=parser.parse_args()
    print('SP_START','pilot' if args.pilot else 'full','workers',args.workers,flush=True)
    db=duckdb.connect()
    db.execute('SET threads=8')
    raw_root=Path('D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/private/health_partitions/rd')
    raw_files=[p.as_posix() for directory in raw_root.iterdir()
               if directory.is_dir() and '201801'<=directory.name<='202312'
               for p in directory.glob('*/*.parquet')]
    xw='D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/reference_step03/municipality_sih6_ibge7_health_region_2023_v2023_1.csv'
    rd=db.execute("""SELECT r.N_AIH,x.health_region_2023::VARCHAR hrid,r.DT_INTER date,r.PROC_REA
      FROM read_parquet(?) r JOIN read_csv_auto(?) x
      ON r.MUNIC_RES=x.municipality_sih6::VARCHAR
      WHERE substr(r.DIAG_PRINC,1,3)='N20' AND r.CAR_INT='02' AND r.IDENT='1'
      AND r.DT_INTER BETWEEN '20180101' AND '20221231'""",[raw_files,xw]).df()
    assert len(rd)==226402 and rd.N_AIH.nunique()==226402
    keys={s.encode() for s in rd.N_AIH}
    print('RD_LINK_KEYS_READY',len(keys),flush=True)
    root=Path('D:/数据库集合/DATASUS数据库/01_raw/SIHSUS/SP')
    files=[]
    for year in range(2018,2024):
        files.extend(sorted((root/str(year)).rglob('*.dbc')))
    if len(files)!=1944:
        raise ValueError('Expected all27states×72months SP')
    if args.pilot:
        files=[f for f in files if f.parent.name in ['AC','SP'] and f.stem[-4:] in ['1801','2201']]
    manifest=[]
    linked={}
    start=time.time()
    init(keys)
    if args.backend=='subprocess':
        buckets=[[] for _ in range(args.workers)]
        loads=[0]*args.workers
        for i,f in sorted(enumerate(files),key=lambda item:item[1].stat().st_size,reverse=True):
            worker=min(range(args.workers),key=lambda j:loads[j])
            buckets[worker].append((i,str(f)))
            loads[worker]+=f.stat().st_size
        def run_worker(worker,items):
            payload=pickle.dumps((keys,items,worker,args.pilot),protocol=5)
            log=WORK/f'sp_worker_{worker}.log'
            with log.open('wb') as stderr:
                process=subprocess.Popen([sys.executable,'-u',str(WORK/'sp_worker.py')],
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=stderr)
                output,_=process.communicate(payload)
            if process.returncode:
                raise RuntimeError(f'SP worker{worker} failed, see {log.name}')
            return pickle.loads(output)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        if args.backend=='subprocess':
            jobs={pool.submit(run_worker,i,items):i for i,items in enumerate(buckets) if items}
        else:
            jobs={pool.submit(scan,(i,str(f))):f for i,f in enumerate(files)}
        for job in as_completed(jobs):
            hits,entry=job.result()
            for key,bit in hits.items():
                linked[key]=linked.get(key,0)|bit
            manifest.extend(entry if isinstance(entry,list) else [entry])
            if args.backend=='subprocess' or len(manifest)%10==0 or len(manifest)==len(files):
                print('SP_PROGRESS',len(manifest),'/',len(files),'seconds',round(time.time()-start),flush=True)
                pd.DataFrame(manifest).to_csv(OUT/('sp_pilot_manifest.csv' if args.pilot else 'sp_source_manifest.csv'),index=False,encoding='utf-8-sig')
    summary={'files':len(files),'bytes':sum(r['bytes'] for r in manifest),'seconds':time.time()-start,
             'linked_authorizations_via_N20_or_target_code_rows':len(linked),'eligible_RD_admissions':len(rd),
             'pilot_only':args.pilot}
    (OUT/('sp_pilot_summary.json' if args.pilot else 'sp_linkage_summary.json')).write_text(json.dumps(summary,indent=2),encoding='utf-8')
    if args.pilot:
        print('SP_PILOT',summary,flush=True)
        return
    rd['sp_seen_in_screened_rows']=[s.encode() in linked for s in rd.N_AIH]
    rd['sp_drainage']=[bool(linked.get(s.encode(),0)&1) for s in rd.N_AIH]
    rd['sp_removal']=[bool(linked.get(s.encode(),0)&2) for s in rd.N_AIH]
    rd['rd_drainage']=rd.PROC_REA.isin(DRAIN)
    rd['rd_removal']=rd.PROC_REA.isin(REMOVE)
    rd['any_drainage']=rd.sp_drainage|rd.rd_drainage
    rd['any_removal']=rd.sp_removal|rd.rd_removal
    rd['generic_principal']=rd.PROC_REA.isin(['0415010012','0415020034'])
    rd['generic_resolved']=rd.generic_principal&(rd.sp_drainage|rd.sp_removal)
    flags=['sp_seen_in_screened_rows','sp_drainage','sp_removal','rd_drainage','rd_removal','any_drainage','any_removal','generic_principal','generic_resolved']
    rd.groupby(['hrid','date'])[flags].sum().reset_index().to_csv(OUT/'sp_linked_daily_outcomes.csv',index=False,encoding='utf-8-sig')
    summary.update({c:int(rd[c].sum()) for c in flags})
    (OUT/'sp_linkage_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('SP_COMPLETE',summary,flush=True)

if __name__=='__main__':
    main()
