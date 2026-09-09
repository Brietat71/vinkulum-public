"""Coût de k_c_m_z à état imposé identique, indépendamment du solveur statique."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import time

import numpy as np
from vinkulum import _vinkulum
from confronte_mbdyn import princeton


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie',type=Path,required=True)
    p.add_argument('--formulation', choices=('milieu', 'integree'), default='milieu')
    args=p.parse_args()
    out={'extension_sha256':hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
         'formulation':args.formulation,'cases':{}}
    for ne in (10,30,60,120):
        n,ids=princeton(ne, formulation=args.formulation)
        t,pos,rot,v,w,vi=n.etat()
        for i in range(len(pos)):
            pos[i][1]=.005*pos[i][0]**2
            w[i]=[.1,.2,.3]
        n.pose_etat(t,pos,rot,v,w,vi)
        n.k_c_m_z()
        times=[]
        for _ in range(3):
            start=time.perf_counter(); result=n.k_c_m_z(); times.append(time.perf_counter()-start)
        rng=np.random.default_rng(731)
        vectors=rng.normal(size=(len(pos)*6,4))
        signatures={name:(np.asarray(result[i])@vectors).tolist() for i,name in enumerate(('K','C'))}
        out['cases'][ne]={'seconds':times,'median_seconds':statistics.median(times),'signatures':signatures}
        print(ne,statistics.median(times),flush=True)
        args.sortie.write_text(json.dumps(out,separators=(',',':'))+'\n')


if __name__=='__main__':main()
