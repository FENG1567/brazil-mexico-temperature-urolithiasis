"""Verify and losslessly reassemble large geographic-day analysis panels."""
from pathlib import Path
import argparse
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[2]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def safe_relative(value):
    path = (ROOT / value).resolve()
    path.relative_to(ROOT)
    return path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true')
    args=parser.parse_args()
    records=json.loads((ROOT/'data/processed_inputs/inventory.json').read_text(encoding='utf-8'))
    for record in records:
        if 'parts' not in record:
            path=safe_relative(record['path'])
            if path.stat().st_size!=record['bytes'] or digest(path)!=record['sha256']:
                raise ValueError('Processed input integrity failed: '+record['name'])
            continue
        name=record['name']
        if Path(name).name!=name:
            raise ValueError('Invalid reassembled filename')
        target=ROOT/'artifacts/processed_inputs'/name
        whole=hashlib.sha256()
        total=0
        temporary=target.with_suffix(target.suffix+'.tmp')
        if not args.verify_only:
            target.parent.mkdir(parents=True,exist_ok=True)
        output=temporary.open('wb') if not args.verify_only else None
        try:
            for part in record['parts']:
                path=safe_relative(part['path'])
                if path.stat().st_size!=part['bytes'] or digest(path)!=part['sha256']:
                    raise ValueError('Panel chunk integrity failed: '+part['path'])
                with path.open('rb') as source:
                    for block in iter(lambda:source.read(1048576),b''):
                        whole.update(block);total+=len(block)
                        if output:output.write(block)
        finally:
            if output:output.close()
        if whole.hexdigest()!=record['sha256'] or total!=record['bytes']:
            temporary.unlink(missing_ok=True)
            raise ValueError('Whole-panel integrity failed: '+name)
        if not args.verify_only:os.replace(temporary,target)
    print(json.dumps({'processed_inputs':len(records),'mode':'verify' if args.verify_only else 'reassemble','status':'passed'}))

if __name__=='__main__':
    main()
