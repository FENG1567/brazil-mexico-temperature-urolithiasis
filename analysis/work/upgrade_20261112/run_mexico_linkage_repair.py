"""Detached sequential execution of the corrected Mexico result builders."""
from pathlib import Path
from datetime import datetime, timezone
import json
import subprocess
import sys

BASE = Path(__file__).resolve().parents[2]
WORK = BASE / 'work/upgrade_20261112'
STATUS = WORK / 'mexico_linkage_repair_status.json'
steps = ['mexico_specific_procedure_models.py', 'repair_mexico_surgery_models.py']
state = {'status':'running', 'steps_completed':[], 'python':sys.executable,
         'started_UTC':datetime.now(timezone.utc).isoformat()}
for script in steps:
    state['current_step'] = script
    STATUS.write_text(json.dumps(state, indent=2), encoding='utf-8')
    print('START_REPAIR_STEP', script, flush=True)
    completed = subprocess.run([sys.executable, '-u', str(WORK/script)], cwd=BASE)
    if completed.returncode:
        state.update(status='failed', exit_code=completed.returncode)
        STATUS.write_text(json.dumps(state, indent=2), encoding='utf-8')
        sys.exit(completed.returncode)
    state['steps_completed'].append(script)
state.update(status='complete', finished_UTC=datetime.now(timezone.utc).isoformat())
STATUS.write_text(json.dumps(state, indent=2), encoding='utf-8')
print('ALL_MEXICO_LINKAGE_REPAIR_STEPS_COMPLETE', flush=True)
