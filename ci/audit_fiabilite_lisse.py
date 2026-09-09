"""Corpus lisse de qualification : mesures, références et limites explicites.

Aucune de ces comparaisons flottantes n'est un certificat de trajectoire.
Les critères sont fixés ici ; les durées ne sont pas des benchmarks de vitesse.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import vinkulum
from vinkulum import Noyau, _vinkulum, bancs, domaines
from audit_assemblage_physique import campagne as assemblage

ROOT = Path(__file__).resolve().parents[1]
CRITERES = dict(toupie_moment_relatif_max=1e-11, toupie_energie_relative_max=1e-11,
               toupie_erreur_omega_max=1e-6, toupie_ordre_min=1.8,
               double_pendule_ordre_min=1.8, quatre_barres_ecart_max=1e-8,
               chute_erreur_max=1e-11)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rotation():
    # Même modèle/critères exécutés contre chacune des roues, y compris celle
    # qui ne contient pas encore les nouvelles contre-épreuves.
    p = ROOT/'python/vinkulum/test_fiabilite_lisse.py'
    spec = importlib.util.spec_from_file_location('modele_toupie', p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    rhs = lambda t, w: np.linalg.solve(m.J, np.cross(m.J@w, w))
    refs = [solve_ivp(rhs, (0., .1), m.W, method='DOP853', rtol=t, atol=t/100)
            for t in (1e-11, 1e-13)]
    if not all(s.success and s.t[-1] == .1 for s in refs):
        raise RuntimeError('référence de toupie incomplète')
    ref = refs[-1].y[:, -1]
    ecart_ref = float(np.linalg.norm(refs[0].y[:, -1]-ref))
    if ecart_ref > 1e-11:
        raise RuntimeError('raffinement de référence insuffisant')
    e0, l0 = .5*m.W@m.J@m.W, m.ROT@m.J@m.W
    rows = []
    for scale in m.ECHELLES:
        for i, cadre in enumerate(m.CADRES):
            for pas in (.002, .001):
                r, w = m.toupie(scale, cadre, pas)
                wb = r.T@w
                rows.append(dict(echelle=scale, cadre=i, pas=pas,
                                 rotation=r.tolist(), omega=w.tolist(),
                                 erreur_omega=float(np.linalg.norm(wb-ref)),
                                 moment_relatif=float(np.linalg.norm(r@m.J@wb-l0)/np.linalg.norm(l0)),
                                 energie_relative=float(abs(.5*wb@m.J@wb-e0)/e0)))
    return dict(reference=ref.tolist(), ecart_reference=ecart_ref, cas=rows,
                modele_sha256=sha(p))


def chute():
    cas = []
    for scale in (2.**-20, 1., 2.**20):
        for schema in ('simule', 'simule_em'):
            g = np.array([2., -3., -9.81])*scale
            v = np.array([.1, -.2, .3])*scale
            n = Noyau(g.tolist())
            n.corps('libre', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], v=v.tolist())
            for t in (.25, .5):
                getattr(n, schema)(t, .002, tous=10**9)
                etat = n.etat()
                erreur = max(np.max(np.abs(np.array(etat[1][0])-(t*v+.5*t*t*g))),
                             np.max(np.abs(np.array(etat[3][0])-(v+t*g))))/scale
                cas.append(dict(echelle=scale, schema=schema, temps=t, erreur=float(erreur)))
    return cas


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    a = p.parse_args()
    result = dict(schema='vinkulum.fiabilite_lisse.1', version=vinkulum.__version__,
                  extension_sha256=sha(_vinkulum.__file__), programme_sha256=sha(__file__),
                  criteres=CRITERES, certification_trajectoire=False,
                  sources_sha256={str(f.relative_to(ROOT)): sha(f) for f in (
                      ROOT/'ci/audit_assemblage_physique.py', ROOT/'python/vinkulum/test_fiabilite_lisse.py')})
    with contextlib.redirect_stdout(io.StringIO()):
        result['toupie'] = rotation()
        result['chute'] = chute()
        result['assemblage'] = assemblage()
        result['stewart'] = domaines.hexapode(bavard=False)
        result['double_pendule'] = [dict(h=h, **bancs.double_pendule(h)) for h in (.002, .001, .0005)]
        result['quatre_barres'] = [dict(h=h, **bancs.quatre_barres(h)) for h in (.002, .001, .0005)]
    failures = []
    for row in result['toupie']['cas']:
        for variable, maximum in [('moment_relatif', CRITERES['toupie_moment_relatif_max']),
                                  ('energie_relative', CRITERES['toupie_energie_relative_max'])]:
            if row[variable] > maximum:
                failures.append(dict(famille='toupie', variable=variable, cas=row))
        if row['pas'] == .001 and row['erreur_omega'] > CRITERES['toupie_erreur_omega_max']:
            failures.append(dict(famille='toupie', variable='erreur_omega', cas=row))
    for coarse, fine in zip(result['toupie']['cas'][::2], result['toupie']['cas'][1::2], strict=True):
        ordre = float(np.log2(coarse['erreur_omega']/fine['erreur_omega']))
        if ordre < CRITERES['toupie_ordre_min']:
            failures.append(dict(famille='toupie', variable='ordre', echelle=fine['echelle'], ordre=ordre))
    for row in result['chute']:
        if row['erreur'] > CRITERES['chute_erreur_max']:
            failures.append(dict(famille='chute', cas=row))
    for coarse, fine in zip(result['double_pendule'][:-1], result['double_pendule'][1:]):
        ordre = float(np.log2(coarse['erreur']/fine['erreur']))
        if ordre < CRITERES['double_pendule_ordre_min']:
            failures.append(dict(famille='double_pendule', ordre=ordre))
    for row in result['quatre_barres']:
        if row['erreur'] > CRITERES['quatre_barres_ecart_max']:
            failures.append(dict(famille='quatre_barres', cas=row))
    if not result['assemblage']['admis']:
        failures.append(dict(famille='assemblage', echecs=result['assemblage']['echecs']))
    result.update(echecs=failures, admis=not failures)
    a.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(vinkulum.__version__, 'corpus lisse :', 'admis' if not failures else f'{len(failures)} échecs')
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
