#!/usr/bin/env python3
"""Replay v1.0.1 and the v1.1.0 full-volume extension using only open solvers."""
import argparse,csv,json,platform,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
import benchmark as b
import full_state as f
from reproduce import independent_point,close
from test_full_state import interior_example

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--audit-only',action='store_true');args=parser.parse_args()
    p=json.loads((ROOT/'protocol_full_state.json').read_text())
    assert p['implementation_sha256']==b.digest(ROOT/'full_state.py') and p['configuration_sha256']==b.digest(f.CONFIG)
    command=[sys.executable,str(ROOT/'reproduce.py')]
    if args.audit_only:command.append('--audit-only')
    subprocess.run(command,check=True)
    subprocess.run([sys.executable,str(ROOT/'test_full_state.py'),'-v'],check=True)
    assert json.loads((ROOT/'interior_volume_witness.json').read_text())==interior_example()
    checks=points=interior_rows=0
    for path in sorted((ROOT/'data/test').glob('*/scenario.json')):
        instance=json.loads(path.read_text());case=instance['id']
        fixed=json.loads((ROOT/'results/test'/case/'candidate/certificates.json').read_text())
        for method in f.CONFIG['methods']:
            folder=ROOT/'results/full_state'/case/method
            rows=json.loads((folder/'certificates.json').read_text());new,_=f.certificates(fixed,method)
            assert b.digest(rows)==b.digest(new)
            point=json.loads((folder/'point.json').read_text());rec=json.loads((folder/'record.json').read_text())
            assert rec['solver']['success'] and point is not None
            audit=independent_point(instance,point)
            for row,old in zip(rows,fixed):
                assert row['box']==old['box'] and row['reference']==old['reference']
                values=(point['x'][row['plant']][row['time']],point['y'][row['plant']][row['time']],point['v'][row['plant']][row['time']])
                for value,lo,hi in zip(values,row['box'][::2],row['box'][1::2]):assert lo-b.CONFIG['point_tolerance']<=value<=hi+b.CONFIG['point_tolerance']
                limit=row['px']*values[0]+row['py']*values[1]+row['pv']*values[2]+row['intercept']
                assert point['g'][row['plant']][row['time']]-limit<=b.CONFIG['point_tolerance']
                coeff=b.CONFIG['coefficients'][row['plant']]
                maximum=f.oracle(coeff,(row['px'],row['py'],row['pv'],row['intercept']),row['box'])
                assert maximum[0]<=0
                interior_rows+=row['box'][4]<maximum[3]<row['box'][5]
                checks+=1
            assert close(audit['scheduled_cost'],rec['point_audit']['scheduled_cost']) and close(audit['final_cost'],rec['point_audit']['final_cost'])
            points+=1
    solves=0
    if not args.audit_only:
        with tempfile.TemporaryDirectory(prefix='full-state-replay-') as temporary:
            output=Path(temporary);f.run(output)
            old=list(csv.DictReader((ROOT/'results/full_state/metrics.csv').open()));new=list(csv.DictReader((output/'metrics.csv').open()))
            old={(r['instance'],r['method']):r for r in old};new={(r['instance'],r['method']):r for r in new};assert old.keys()==new.keys()
            for key,a in old.items():
                c=new[key];assert a['admissible']==c['admissible']
                for field in ('scheduled_cost','final_cost'):assert close(float(a[field]),float(c[field]))
                assert (ROOT/'results/full_state'/key[0]/key[1]/'certificates.json').read_bytes()==(output/key[0]/key[1]/'certificates.json').read_bytes()
            solves=len(new)
    old=json.loads((ROOT/'REPRODUCTION_CHECKS.json').read_text());assert checks==2592 and points==36
    result=dict(status='passed',version='1.1.0',original_replay=old,full_state_LP_replays=solves,total_LP_replays=old['replayed_LP_solves']+solves,
                full_state_exact_rows=checks,total_exact_rows=old['independent_full_box_rational_checks']+checks,full_state_points=points,
                state_interior_exported_rows=interior_rows,essential_interior_volume_witness_verified=True,
                maximum_candidates_per_plant=[8,12,12],fixed_and_full_domains_identical=True,
                proprietary_solver_required=False,inherited_data_or_code_required=False,network_access_required=False,
                python=platform.python_version(),numpy=np.__version__,scipy=b.scipy.__version__)
    b.write_json(ROOT/'FULL_STATE_REPRODUCTION_CHECKS.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
