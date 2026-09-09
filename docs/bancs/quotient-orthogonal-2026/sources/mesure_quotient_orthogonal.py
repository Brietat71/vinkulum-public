"""Sondage reproductible du quotient local face à un QR dense des normales.

Chaque essai s'exécute dans un processus frais. Les durées excluent les
imports, la construction du modèle et les oracles ; elles incluent la
normalisation, la factorisation, une projection et un calcul de réactions
sur trois vecteurs. Ce n'est pas une confrontation de solveurs multicorps.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CASES = [('chaine', n) for n in (257, 1025, 2049)] + [
    (f, n) for f in ('cascade', 'cascade_tournee') for n in (4, 20, 60)] + [
    ('couplages', n) for n in (129, 513)]
VARIANTS = ('local', 'qr_dense')
THREADS = {k: '1' for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS')}
BUDGET = 2_000_000


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def modele(family, size):
    import numpy as np
    from scipy.sparse import csc_matrix, csr_matrix, diags, vstack
    if family in ('chaine', 'couplages'):
        g = diags((np.ones(size-1), -np.ones(size-1)), (0, 1), shape=(size-1, size), format='csr')
        if family == 'chaine':
            g = vstack((g, g[::2]+g[1::2]), format='csr')
        else:
            rng = np.random.default_rng(473)
            cols = np.array([rng.choice(size, 3, replace=False) for _ in range(size)])
            extra = csr_matrix((np.tile([1., 1., -2.], size), (np.repeat(np.arange(size), 3), cols.ravel())), shape=(size, size))
            g = vstack((g, extra), format='csr')
        return g, None
    from scipy.spatial.transform import Rotation
    from vinkulum.test_noyau import _cascade_initiale
    rotation = Rotation.from_rotvec([.3, -.7, .4]).as_matrix() if family.endswith('tournee') else np.eye(3)
    n = _cascade_initiale(size, rotation)
    nr, nc, ptr, idx, data = n.k_c_m_g_creux()[-1]
    return csc_matrix((data, idx, ptr), shape=(nr, nc)).tocsr(), rotation


def reference(family, size, rotation, v):
    import numpy as np
    if family in ('chaine', 'couplages'):
        return np.broadcast_to(v.mean(axis=0), v.shape), size-1
    # Cinématique indépendante du jacobien assemblé : un angle par cellule.
    # Les deux bielles tournent ensemble ; la traverse conserve sa direction.
    # Un angle amont translate tous les corps en aval de deux bras dérivés.
    arm_prime = rotation@np.array([-np.sin(np.pi/3), np.cos(np.pi/3), 0.])/(2*size)
    axis = rotation@np.array([0., 0., 1.])
    z = np.zeros((18*size, size))
    for cell in range(size):
        for body in range(3):
            at = 18*cell+6*body
            z[at:at+3, :cell] = 2*arm_prime[:, None]
            z[at:at+3, cell] = (2 if body == 1 else 1)*arm_prime
            if body != 1:
                z[at+3:at+6, cell] = axis
    orth, _ = np.linalg.qr(z, mode='reduced')
    return orth@(orth.T@v), 17*size


class QRDense:
    def __init__(self, g):
        from scipy.linalg import qr
        import numpy as np
        from contraintes_orthogonales import normalise, RangAmbigu
        self.g, self.echelles, self.normes = normalise(g)
        self.seuil = 64*np.finfo(float).eps*max(g.shape)
        q, r, indices = qr(self.g.T.toarray(), mode='economic', pivoting=True)
        diag = abs(np.diag(r))
        if np.any((diag >= self.seuil/4) & (diag <= 4*self.seuil)):
            raise RangAmbigu('rang QR dense ambigu')
        self.rang = int(np.count_nonzero(diag > self.seuil))
        self.q = q[:, :self.rang]
        self.r = r[:self.rang, :self.rang]
        self.indices = indices[:self.rang]
        self.diagnostic = dict(rang=self.rang, coefficients_base_normales=int(self.q.size), base_dense=True)

    def projette(self, v):
        return v-self.q@(self.q.T@v)

    def reactions(self, v):
        import numpy as np
        from scipy.linalg import solve_triangular
        out = np.zeros((self.g.shape[0], v.shape[1]))
        out[self.indices] = solve_triangular(self.r, self.q.T@v)
        return out/self.normes[:, None]/self.echelles[:, None]


def worker(family, size, variant, path):
    import numpy as np
    import scipy
    import resource
    from contraintes_orthogonales import QuotientOrthogonal, RangAmbigu, BudgetDepasse
    g, rotation = modele(family, size)
    v = np.random.default_rng(946).normal(size=(g.shape[1], 3))
    result = dict(famille=family, taille=size, variante=variant, forme=list(g.shape), nnz=g.nnz,
                  python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                  cpu=sorted(os.sched_getaffinity(0)), fils={k: os.environ.get(k) for k in THREADS})
    if family.startswith('cascade'):
        from vinkulum import _vinkulum
        import vinkulum.test_noyau
        result['extension_sha256'] = sha(_vinkulum.__file__)
        result['modele_sha256'] = sha(vinkulum.test_noyau.__file__)
    start = time.perf_counter()
    try:
        q = QuotientOrthogonal(g, budget=BUDGET) if variant == 'local' else QRDense(g)
        factorized = time.perf_counter()
        p = q.projette(v)
        projected = time.perf_counter()
        eta = q.reactions(v)
        reacted = time.perf_counter()
        result.update(statut='calcule', factorisation_s=factorized-start,
                      projection_s=projected-factorized, reactions_s=reacted-projected,
                      total_s=reacted-start, diagnostic=q.diagnostic)
    except (RangAmbigu, BudgetDepasse) as exc:
        result.update(statut=type(exc).__name__, message=str(exc), total_s=time.perf_counter()-start)
    result['rss_avant_oracle_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if result['statut'] == 'calcule':
        ref, rank = reference(family, size, rotation, v)
        checks = dict(erreur_projection=float(np.max(abs(p-ref))),
                      erreur_force=float(np.max(abs(g.T@eta-(v-ref)))),
                      fermeture=float(np.max(abs(q.g@p), initial=0.)),
                      idempotence=float(np.max(abs(q.projette(p)-p))),
                      rang_attendu=rank, rang_correct=q.rang == rank)
        result['controles'] = checks
        result['accepte'] = checks['rang_correct'] and max(checks[k] for k in ('erreur_projection', 'erreur_force', 'fermeture', 'idempotence')) <= 1e-9
    Path(path).write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def campaign(path, cpu):
    if path.exists():
        raise SystemExit('sortie déjà présente : choisir un nouveau chemin')
    if cpu not in os.sched_getaffinity(0):
        raise SystemExit('CPU indisponible')
    path.mkdir(parents=True)
    sources = path/'sources'
    sources.mkdir()
    for name in ('mesure_quotient_orthogonal.py', 'contraintes_orthogonales.py', 'test_contraintes_orthogonales.py'):
        (sources/name).write_bytes((ROOT/'ci'/name).read_bytes())
    env = dict(os.environ, **THREADS)
    records = []
    for family, size in CASES:
        for repetition in range(4):
            variants = VARIANTS if repetition % 2 == 0 else VARIANTS[::-1]
            for variant in variants:
                name = f'{family}-{size}-{variant}-{repetition}.json'
                subprocess.run(['taskset', '-c', str(cpu), sys.executable, __file__, '--worker', family, str(size), variant, str(path/name)],
                               env=env, check=True, timeout=120)
                result = json.loads((path/name).read_text())
                result.update(fichier=name, repetition=repetition, echauffement=repetition == 0, sha256=sha(path/name))
                records.append(result)
                print(f'{family} {size} {variant} {repetition}: {result["statut"]} {result["total_s"]:.4g}s', flush=True)
    out = dict(schema=1, date=datetime.now(timezone.utc).isoformat(), cas=CASES, variantes=VARIANTS,
               programme_sha256=sha(__file__), prototype_sha256=sha(ROOT/'ci/contraintes_orthogonales.py'),
               cpu_info=Path('/proc/cpuinfo').read_text().split('model name\t: ')[1].split('\n')[0],
               budget_coefficients_transformes=BUDGET, repetitions=4, essais=records)
    (path/'manifest.json').write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    verify(path)


def verify(path):
    data = json.loads((path/'manifest.json').read_text())
    assert sha(path/'sources/mesure_quotient_orthogonal.py') == data['programme_sha256']
    assert sha(path/'sources/contraintes_orthogonales.py') == data['prototype_sha256']
    expected = {(f, n, v, r) for f, n in data['cas'] for v in data['variantes'] for r in range(data['repetitions'])}
    found = set()
    for trial in data['essais']:
        key = (trial['famille'], trial['taille'], trial['variante'], trial['repetition'])
        assert key in expected and key not in found, key
        found.add(key)
        assert sha(path/trial['fichier']) == trial['sha256']
        raw = json.loads((path/trial['fichier']).read_text())
        assert all(trial[k] == v for k, v in raw.items())
        assert trial['fils'] == THREADS and len(trial['cpu']) == 1
        assert trial['echauffement'] == (trial['repetition'] == 0)
        assert trial['statut'] in ('calcule', 'RangAmbigu', 'BudgetDepasse')
        assert math.isfinite(trial['total_s']) and trial['total_s'] > 0
        if trial['statut'] == 'calcule':
            assert trial['accepte'], key
            checks = trial['controles']
            assert checks['rang_correct'] and trial['diagnostic']['rang'] == checks['rang_attendu']
            for metric in ('erreur_projection', 'erreur_force', 'fermeture', 'idempotence'):
                assert math.isfinite(checks[metric]) and 0 <= checks[metric] <= 1e-9, (key, metric)
            for metric in ('factorisation_s', 'projection_s', 'reactions_s'):
                assert math.isfinite(trial[metric]) and trial[metric] >= 0
            assert math.isclose(sum(trial[k] for k in ('factorisation_s', 'projection_s', 'reactions_s')), trial['total_s'], rel_tol=1e-12)
    assert found == expected
    for family, size in data['cas']:
        line = []
        for variant in data['variantes']:
            trials = [r for r in data['essais'] if (r['famille'], r['taille'], r['variante']) == (family, size, variant) and not r['echauffement']]
            statuses = {r['statut'] for r in trials}
            assert len(statuses) == 1
            line.append((variant, next(iter(statuses)), statistics.median(r['total_s'] for r in trials)))
        print(family, size, line)
    refus = sum(r['statut'] != 'calcule' for r in data['essais'])
    print(f'Archive vérifiée : {len(found)} essais, dont {refus} refus et les échauffements.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--worker', nargs=4)
    group.add_argument('--sortie', type=Path)
    group.add_argument('--verifier', type=Path)
    parser.add_argument('--cpu', type=int, default=8)
    args = parser.parse_args()
    if args.worker:
        family, size, variant, output = args.worker
        worker(family, int(size), variant, output)
    elif args.verifier:
        verify(args.verifier)
    else:
        campaign(args.sortie, args.cpu)
