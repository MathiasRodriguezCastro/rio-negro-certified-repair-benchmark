#!/usr/bin/env python3
"""Replay the retained v1.1.0 examples plus v1.2.0 global-bound example."""
import json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(ROOT/'reproduce_full_state.py')],check=True)
subprocess.run([sys.executable,str(ROOT/'reproduce_global_bounds.py')],check=True)
old=json.loads((ROOT/'FULL_STATE_REPRODUCTION_CHECKS.json').read_text())
new=json.loads((ROOT/'GLOBAL_BOUNDS_REPRODUCTION_CHECKS.json').read_text())
result=dict(status='passed',version='1.2.0',retained_full_state_replay=old,global_bounds_example=new,
            total_LP_replays=old['total_LP_replays']+new['replayed_LP_solves'])
(ROOT/'ALL_REPRODUCTION_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
