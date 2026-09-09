"""Qualification numérique reproductible, distincte d'un certificat spectral.

Les archives contiennent les opérateurs natifs et les modes : la relecture
ne rappelle ni le noyau ni son solveur modal. Les références sont assemblées
indépendamment dans modeles_modes. Les temps ne servent pas de benchmark.
"""
import argparse
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np
from scipy.sparse import csc_matrix
import scipy

import modeles_modes

SCHEMA = 'vinkulum.qualification.modale.1'
SEUIL = 1e-8
TAILLES = (4, 16, 64, 256, 512, 1024)
FAMILLES = ('console', 'articulee', 'gravite')


def empreinte(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference(famille, nb, count):
    if famille == 'console':
        return modeles_modes.reference_console(nb, count)
    return modeles_modes.reference_chaine(nb, famille == 'gravite')[:count]


def matrice(export):
    nr, nc, ptr, idx, data = export
    return csc_matrix((data, idx, ptr), shape=(nr, nc))


def norme(a):
    return float(np.max(np.asarray(abs(a).sum(axis=1)), initial=0.))


def ratio(num, den):
    return num / den if den else (0. if num == 0. else float('inf'))


def verifier_cas(cas):
    """Recalcul binary64 indépendant ; aucune borne d'arrondi validée."""
    if cas['schema'] != SCHEMA:
        raise ValueError('schéma inconnu')
    famille, nb = cas['famille'], cas['taille']
    if famille not in FAMILLES or nb not in TAILLES:
        raise ValueError('cas hors grille')
    result = cas['resultat']
    modes = result['modes']
    if len(modes) != (6 if famille == 'console' else min(6, nb)):
        raise ValueError('nombre de modes incorrect')
    k, _, m, g = map(matrice, cas['operateurs'])
    k = .5*k + .5*k.T
    x = np.array([a['forme'] for a in modes]).T
    eta = np.array([a['reactions'] for a in modes]).T
    values = np.array([a['valeur_propre'] for a in modes])
    ref = reference(famille, nb, len(values))
    if not np.isfinite(x).all() or not np.isfinite(eta).all() or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('modes non finis ou non positifs')
    freq_error = float(np.max(abs(np.sqrt(values/ref)-1)))
    orth = float(np.max(abs(x.T @ (m @ x)-np.eye(len(values)))))
    residual = k @ x + g.T @ eta - (m @ x)*values
    errors, closures = [], []
    for j, lam in enumerate(values):
        nx, ne = float(max(abs(x[:, j]))), float(np.max(abs(eta[:, j]), initial=0.))
        den = (norme(k)+abs(lam)*norme(m))*nx+norme(g.T)*ne
        errors.append(ratio(float(max(abs(residual[:, j]))), den))
        closures.append(ratio(float(np.max(abs(g @ x[:, j]), initial=0.)), norme(g)*nx))
    reported = max(max(a['erreur_arriere'], a['erreur_arriere_originale'], a['residu_contrainte']) for a in modes)
    metrics = dict(erreur_relative_frequence=freq_error, orthogonalite_massique=orth,
                   residu_physique_relatif=max(errors), contrainte_relative=max(closures),
                   controle_solveur=reported)
    if not cas['etat_inchange'] or result['certification_spectrale'] or result['rang_certifie']:
        raise ValueError('contrat de qualification incorrect')
    if any(not np.isfinite(v) or v > SEUIL for v in metrics.values()):
        raise ValueError(f'qualification refusée : {metrics}')
    return metrics


def verifier(dossier):
    dossier = Path(dossier)
    manifest = json.loads((dossier/'qualification.json').read_text())
    if manifest['schema'] != SCHEMA or manifest['seuil'] != SEUIL:
        raise ValueError('contrat de manifeste incorrect')
    expected = {(f, n) for f in FAMILLES for n in TAILLES}
    seen = set()
    results = []
    for entry in manifest['cas']:
        name = entry['fichier']
        if Path(name).name != name:
            raise ValueError('chemin incorrect')
        path = dossier/name
        if empreinte(path) != entry['sha256']:
            raise ValueError('empreinte incorrecte')
        cas = json.loads(gzip.decompress(path.read_bytes()))
        key = (cas['famille'], cas['taille'])
        if key in seen:
            raise ValueError('cas dupliqué')
        seen.add(key)
        results.append(verifier_cas(cas))
    if seen != expected:
        raise ValueError('grille incomplète')
    return results


def produire(dossier, roue):
    import vinkulum
    import vinkulum.modes_creux as modal
    import modes_creux_avant_qualification as avant
    dossier = Path(dossier)
    dossier.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema=SCHEMA, seuil=SEUIL, certification_spectrale=False,
                    python=platform.python_version(), plateforme=platform.platform(),
                    numpy=np.__version__, scipy=scipy.__version__,
                    vinkulum=importlib.metadata.version('vinkulum'),
                    roue=dict(nom=Path(roue).name, sha256=empreinte(roue)),
                    extensions={p.name: empreinte(p) for p in Path(vinkulum.__file__).parent.glob('*.so')},
                    sources={str(Path(p).name): empreinte(p) for p in
                             (modal.__file__, avant.__file__, modeles_modes.__file__, __file__)},
                    cas=[])
    for famille in FAMILLES:
        for nb in TAILLES:
            n = modeles_modes.vinkulum(famille, nb)
            etat = n.etat_precis()
            result = n.modes_creux(6)
            old = avant.analyse(n, 6)
            values = np.array([a['valeur_propre'] for a in old['modes']])
            ref = reference(famille, nb, len(values))
            cas = dict(schema=SCHEMA, famille=famille, taille=nb, resultat=result,
                       operateurs=n.k_c_m_g_creux(),
                       avant_qualification=dict(valeurs=values.tolist(),
                           erreur_relative_frequence=float(np.max(abs(np.sqrt(values/ref)-1)))),
                       etat_inchange=(etat == n.etat_precis()))
            metrics = verifier_cas(cas)
            name = f'{famille}-{nb}.json.gz'
            (dossier/name).write_bytes(gzip.compress(json.dumps(cas, allow_nan=False).encode(), mtime=0))
            manifest['cas'].append(dict(fichier=name, sha256=empreinte(dossier/name), **metrics))
            (dossier/'qualification.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
            print(f'{famille} {nb}: {metrics}', flush=True)
    verifier(dossier)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sortie')
    parser.add_argument('--roue')
    parser.add_argument('--verifier')
    args = parser.parse_args()
    if args.verifier:
        print(json.dumps(verifier(args.verifier), indent=2))
    elif args.sortie and args.roue:
        produire(args.sortie, args.roue)
    else:
        parser.error('--verifier ou --sortie et --roue requis')
