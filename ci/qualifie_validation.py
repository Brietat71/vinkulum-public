"""Rejoue la validation d'une roue identifiée, sans modifier le rapport historique.

Les familles tournent dans des processus isolés. Un écart physique ne devient
jamais une erreur d'exécution, et un cas absent ne devient jamais un succès.
Les versions et les empreintes identifient l'expérience, pas une certification.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.metadata as md
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import zipfile

FAMILLES = ('cas_wright', 'cas_johnson', 'cas_lock', 'cas_maryland', 'cas_coleman',
            'cas_elliott', 'cas_s809', 'cas_jantzen', 'cas_kim', 'cas_references')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ecrit(p, obj):
    Path(p).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def identite(ligne):
    return tuple(ligne[k] for k in ('cas', 'source', 'grandeur'))


def historique(p):
    rows = []
    for line in Path(p).read_text().splitlines():
        if not line.startswith('| ') or '| OK |' not in line and '| ÉCART |' not in line:
            continue
        c = [v.strip() for v in line.split('|')[1:-1]]
        rows.append(dict(cas=c[0], source=c[1], grandeur=c[2], reference_affichee=c[3],
                         valeur_affichee=c[4], verdict=c[6]))
    if not rows or len({identite(r) for r in rows}) != len(rows):
        raise ValueError('inventaire historique vide ou ambigu')
    return rows


def rapproche(ancien, courant):
    index = {}
    for row in courant:
        key = identite(row)
        if key in index:
            raise ValueError(f'identité de résultat dupliquée : {key}')
        index[key] = row
    return [dict(historique=old, actuel=index.get(identite(old)),
                 suivi='absent' if identite(old) not in index else index[identite(old)]['statut'])
            for old in ancien]


def paquet(roue):
    import vinkulum
    from vinkulum import _vinkulum
    root = Path(vinkulum.__file__).parent
    with zipfile.ZipFile(roue) as z:
        fichiers = {n: hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist()
                    if n.startswith('vinkulum/') and not n.endswith('/')}
    for name, empreinte in fichiers.items():
        if sha(root / Path(name).relative_to('vinkulum')) != empreinte:
            raise ValueError(f'paquet installé différent de la roue : {name}')
    return dict(version=vinkulum.__version__, roue=Path(roue).name, roue_sha256=sha(roue),
                extension_sha256=sha(_vinkulum.__file__), fichiers_sha256=fichiers,
                python=sys.version, plateforme=platform.platform(),
                dependances={d.metadata['Name']: d.version for d in sorted(md.distributions(),
                              key=lambda d: d.metadata['Name'].lower())})


def travail(famille, sortie):
    from vinkulum import validation as v
    v.LIGNES.clear()
    tolerances = {}
    original = v._ligne

    def mesure(cas, source, grandeur, ref, val, tol, note='', re=None):
        tolerances[len(v.LIGNES)] = float(tol)
        return original(cas, source, grandeur, ref, val, tol, note, re)

    v._ligne = mesure
    erreur = None
    debut = time.monotonic()
    try:
        getattr(v, famille)(False)
    except Exception as ex:
        traceback.print_exc()
        erreur = dict(type=type(ex).__name__, message=str(ex))
    rows = []
    for i, (cas, src, gr, ref, val, ecart, ok, note, re) in enumerate(v.LIGNES):
        banc = gr in ('banc', 'banc `verification`')
        invalide = (not banc and not all(math.isfinite(x) for x in (ref, val, ecart)))
        execution = note.startswith('ERREUR') or invalide or banc and not ok and not math.isfinite(ecart)
        # L'ancienne branche Andrews encode l'exception dans une ligne NaN
        # sans préfixe ERREUR ; son écart 1 ne permet pas de la qualifier.
        if banc and not ok and not math.isfinite(val):
            execution = True
        statut = 'erreur_execution' if execution else ('ok' if ok else 'ecart')
        propre = lambda x: float(x) if x is not None and math.isfinite(x) else None
        rows.append(dict(cas=cas, source=src, grandeur=gr, reference=propre(ref), valeur=propre(val),
                         ecart_relatif=propre(ecart), tolerance_relative=tolerances.get(i),
                         statut=statut, note=note, reynolds=propre(re),
                         nature='banc_assertions' if banc else 'comparaison_numerique'))
    ecrit(sortie, dict(famille=famille, secondes=time.monotonic()-debut,
                      erreur=erreur, lignes=rows))
    return int(erreur is not None or any(r['statut'] == 'erreur_execution' for r in rows))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--roue', type=Path)
    p.add_argument('--historique', type=Path)
    p.add_argument('--sortie', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--delai', type=float, default=1800.)
    p.add_argument('--famille', choices=FAMILLES, help=argparse.SUPPRESS)
    a = p.parse_args()
    if a.famille:
        return travail(a.famille, a.sortie)
    if not a.roue or not a.historique or a.jobs < 1 or not math.isfinite(a.delai) or a.delai <= 0:
        p.error('--roue, --historique, jobs et délai positifs requis')
    preuve = paquet(a.roue)
    anciens = historique(a.historique)
    a.sortie.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MPLCONFIGDIR=str(a.sortie.resolve() / 'matplotlib'))

    def execute(f):
        out, log = a.sortie / (f + '.json'), a.sortie / (f + '.log')
        code, erreur = None, None
        try:
            with log.open('w') as stream:
                code = subprocess.run([sys.executable, '-I', str(Path(__file__).resolve()),
                                       '--famille', f, '--sortie', str(out.resolve())],
                                      stdout=stream, stderr=subprocess.STDOUT, env=env,
                                      cwd=a.sortie, timeout=a.delai).returncode
        except subprocess.TimeoutExpired:
            erreur = dict(type='TimeoutExpired', message=f'délai {a.delai} s dépassé')
        if out.exists():
            data = json.loads(out.read_text())
        else:
            data = dict(famille=f, lignes=[], erreur=erreur or dict(type='Processus', message=f'code {code}'))
        data.update(code_sortie=code, log_sha256=sha(log))
        print(f'{f}: {len(data["lignes"])} lignes, code {code}', flush=True)
        return data

    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        familles = list(pool.map(execute, FAMILLES))
    if paquet(a.roue) != preuve:
        raise ValueError('le paquet ou son environnement a changé pendant la campagne')
    rows = [r for f in familles for r in f['lignes']]
    suivi = rapproche(anciens, rows)
    complet = all(f['code_sortie'] == 0 and f['erreur'] is None for f in familles)
    bilan = dict(schema='vinkulum.qualification_physique.1', environnement=preuve,
                 programme_sha256=sha(__file__), historique_sha256=sha(a.historique),
                 mode='complet', tolerances_modifiees=False, familles=familles,
                 suivi_historique=suivi, nouvelles_lignes=[r for r in rows
                     if identite(r) not in {identite(h) for h in anciens}],
                 execution_complete=complet,
                 historique_entier_rejoue=all(r['suivi'] != 'absent' for r in suivi),
                 comptes={s: sum(r['statut'] == s for r in rows)
                          for s in ('ok', 'ecart', 'erreur_execution')},
                 certification_mathematique=False)
    ecrit(a.sortie / 'rapport.json', bilan)
    print(json.dumps({k: bilan[k] for k in ('execution_complete', 'historique_entier_rejoue', 'comptes')}))
    return int(not complet or not bilan['historique_entier_rejoue'])


if __name__ == '__main__':
    sys.exit(main())
