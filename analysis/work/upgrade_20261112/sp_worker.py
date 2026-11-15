"""Local binary-pipe worker. Linkage keys stay in memory and pipes."""
import sys
import pickle
import json
import time
from scan_secondary_procedures import scan,init,OUT

keys,items,worker,pilot=pickle.loads(sys.stdin.buffer.read())
init(keys)
linked={}
manifest=[]
start=time.time()
for item in items:
    hits,entry=scan(item)
    for key,bit in hits.items():
        linked[key]=linked.get(key,0)|bit
    manifest.append(entry)
    if len(manifest)%10==0 or len(manifest)==len(items):
        progress={'worker':worker,'completed':len(manifest),'total':len(items),
                  'bytes':sum(r['bytes'] for r in manifest),'seconds':time.time()-start,'pilot':pilot}
        target=OUT/f'sp_worker_{worker}_progress.json'
        target.write_text(json.dumps(progress,indent=2),encoding='utf-8')
sys.stdout.buffer.write(pickle.dumps((linked,manifest),protocol=5))
sys.stdout.buffer.flush()
