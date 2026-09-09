"""Tangentes de plusieurs familles et coût statique dans des processus isolés.

Comparer deux extensions avec le même programme et les mêmes modèles.
Les échecs restent des échecs ; les temps de k_c_m_z incluent les copies Python.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import expm
from vinkulum import Noyau, _vinkulum
from vinkulum.test_noyau import (_pale_inflow_spatial, _tangentes_champ_fd,
                                 _roulement_charge, _couple_precharge)

ROOT = Path(__file__).resolve().parents[1]


def modele(famille, nb):
    n = Noyau([0., 0., -9.81])
    for i in range(nb):
        n.corps(str(i), 1., np.diag([.2, .3, .4]).ravel().tolist(),
                [.4*i, .03*np.sin(i), .46], v=[.2, -.3, .1], w=[.1, .2, .3])
    if famille == 'couples':
        for i in range(nb):
            n.couple(str(i), i-1 if i else None, i, [0., 0., 1.],
                     ('ressort', [30., .7, .2]))
    elif famille == 'contacts':
        for i in range(nb):
            n.contact('sol'+str(i), i, [.02, .03, 0.], .5,
                      k=110., c=.2, mu=.3, v_eps=.05)
            if i:
                n.contact('paire'+str(i), i, [.03, .01, 0.], .3,
                          b=i-1, pb=[0., .02, 0.], rayon_b=.2,
                          k=170., c=.3, mu=.2, v_eps=.05)
    elif famille == 'supers':
        for i in range(nb-2):
            k = np.diag(np.arange(1., 19.))
            n.superelement(str(i), [i+1, i, i+2], k.ravel().tolist(), beta=.07)
    else:
        raise ValueError(famille)
    s = n.etat_precis()
    for i in range(nb):
        a = .01*np.sin(i+.3)
        s[1][i][1] += .01*np.cos(i)
        s[2][i] = expm(np.array([[0., -.3*a, a], [.3*a, 0., -.2*a],
                                 [-a, .2*a, 0.]])).ravel().tolist()
    n.pose_etat(*s)
    return n


def controle(famille):
    n = modele(famille, 4)
    matrices = n.k_c_m_z()[:2]
    ref = _tangentes_champ_fd(n)
    errors = {}
    for name, m, r in zip(('K', 'C'), matrices, ref):
        m = np.asarray(m)
        errors[name] = dict(max_abs=float(np.max(np.abs(m-r))),
                            relatif=float(np.linalg.norm(m-r)/max(np.linalg.norm(r), 1e-300)))
        np.testing.assert_allclose(m, r, rtol=3e-6, atol=2e-6)
    return errors


def mesure(famille, nb):
    n = modele(famille, nb)
    n.k_c_m_z()
    times = []
    for _ in range(3):
        start = time.perf_counter()
        value = n.k_c_m_z()
        times.append(time.perf_counter()-start)
        matrices = value[:2]
        del value
    vectors = np.random.default_rng(732).normal(size=(nb*6, 4))
    return dict(famille=famille, corps=nb, secondes=times,
                mediane_s=statistics.median(times),
                projections={name:(np.asarray(m)@vectors).tolist()
                             for name, m in zip(('K', 'C'), matrices)})


def statique(ne):
    from mesure_statique import run
    def high_water():
        # Sous Linux, getrusage peut garder le pic du parent après exec.
        # VmHWM concerne l'espace mémoire courant ; ne pas confondre les deux.
        for line in Path('/proc/self/status').read_text().splitlines():
            if line.startswith('VmHWM:'):
                return int(line.split()[1])
        raise RuntimeError('VmHWM indisponible : mesure mémoire non qualifiée')
    before = high_water()
    start = time.perf_counter()
    try:
        r = run(ne, strict=True)
        r['erreur'] = None
    except Exception as e:
        r = dict(erreur=str(e))
    r.update(intervalles=ne, total_s=time.perf_counter()-start,
             rss_initial_kib=before, rss_max_kib=high_water(),
             rusage_max_kib_non_utilise=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path)
    p.add_argument('--tailles', type=int, nargs='+', default=[8, 32, 128])
    p.add_argument('--statique', type=int, nargs='*', default=[])
    p.add_argument('--repetitions-statique', type=int, default=3)
    p.add_argument('--worker-statique', type=int)
    args = p.parse_args()
    if args.worker_statique is not None:
        print(json.dumps(statique(args.worker_statique)))
        return
    if args.sortie is None:
        p.error('--sortie requis')
    if args.repetitions_statique < 1:
        p.error('--repetitions-statique doit être positif')
    out = dict(python=platform.python_version(), plateforme=platform.platform(),
               extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
               programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               modeles_et_controles_sha256=hashlib.sha256((ROOT/'python/vinkulum/test_noyau.py').read_bytes()).hexdigest(),
               fils={k:os.environ.get(k) for k in ('RAYON_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')},
               controles={}, inflows={}, mesures=[], statique=[],
               roulement=_roulement_charge().audit_jacobien(.05),
               amortissement_sous_precharge=dict(attendu=1., observe=_couple_precharge().k_c_m_z()[1][5][5]))
    def save():
        args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')
    for kind in ('harmoniques', 'profil', 'carte'):
        n = _pale_inflow_spatial(kind)
        out['inflows'][kind] = n.audit_jacobien(.05)[0]
    for famille in ('couples', 'contacts', 'supers'):
        out['controles'][famille] = controle(famille)
        for nb in args.tailles:
            r = mesure(famille, nb)
            out['mesures'].append(r)
            print(famille, nb, r['mediane_s'], flush=True)
            save()
    for ne, repetition in ((ne, r) for ne in args.statique for r in range(args.repetitions_statique)):
        cmd = [sys.executable, str(Path(__file__).resolve()), '--worker-statique', str(ne)]
        try:
            process = subprocess.run(cmd, text=True, capture_output=True, timeout=120)
            r = dict(json.loads(process.stdout)) if process.returncode == 0 else dict(erreur=process.stderr)
            r.update(code_sortie=process.returncode, stderr=process.stderr)
        except subprocess.TimeoutExpired:
            r = dict(erreur='délai de 120 s dépassé')
        r['intervalles'] = ne
        r['repetition'] = repetition
        out['statique'].append(r)
        print('statique', ne, r.get('seconds'), r.get('rss_max_kib'), r['erreur'], flush=True)
        save()
    save()


if __name__ == '__main__':
    main()
