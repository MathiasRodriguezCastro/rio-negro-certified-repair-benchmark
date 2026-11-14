#!/usr/bin/env python3
"""Verify the supplied exact example and optionally repeat its four open LPs."""
import argparse,hashlib,json,subprocess,sys,tempfile
from fractions import Fraction as F
from pathlib import Path
from verify_global_bounds_demo import verify
ROOT=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--audit-only',action='store_true');args=ap.parse_args()
    protocol=json.loads((ROOT/'protocol_global_bounds.json').read_text())
    for name,expected in protocol['sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected
    data=json.loads((ROOT/'global_bounds_example/input.json').read_text());stored=json.loads((ROOT/'global_bounds_example/certificate.json').read_text())
    subprocess.run([sys.executable,'-S',str(ROOT/'verify_global_bounds_demo.py')],check=True,stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable,str(ROOT/'test_global_bounds_demo.py'),'-v'],check=True)
    result=verify(data,stored);solves=0
    if not args.audit_only:
        from global_bounds_demo import generate
        with tempfile.TemporaryDirectory(prefix='global-bounds-replay-') as temporary:
            evidence=generate(ROOT/'global_bounds_example/input.json',Path(temporary)/'certificate.json')
            fresh=verify(data,evidence)
            for name in ('lower','scheduled_upper','final_upper'):
                assert abs(float(F(fresh[name])-F(result[name])))<=1e-7*(1+abs(float(F(result[name]))))
            solves=evidence['LP_solves']
    result.update(version='1.2.0',replayed_LP_solves=solves,adversarial_tests_passed=6,standard_library_only_verification_passed=True,
                  private_data_required=False,claim='Complete synthetic verification example, not reproduction of restricted real-bank cost caps')
    (ROOT/'GLOBAL_BOUNDS_REPRODUCTION_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
