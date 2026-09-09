"""Confrontation par interfaces publiques, précision commune avant classement.

Aucun code d'implémentation des solveurs concurrents n'est lu ou importé
comme source. Leurs interfaces publiques exécutent les mêmes modèles.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import resource
import signal
import subprocess
import sys
import time
import threading

SCHEMA = 'vinkulum.confrontation.modale.1'
TAILLES = (4, 16, 64, 256, 1024)
MOTEURS = ('vinkulum', 'exudyn_cartesien_dense', 'exudyn_arbre_dense',
           'exudyn_arbre_creux', 'mbdyn_lapack', 'mbdyn_arpack', 'mbdyn_arpack_large')
SEUIL = 1e-8
REPETITIONS = 3
DELAI = 60
MEMOIRE = 4*1024**3
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'ci'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ecrire(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def oracle(nb):
    import numpy as np
    from modeles_modes import reference_chaine
    valeurs = reference_chaine(nb)
    k = min(6, nb-1)
    freq = np.sqrt(valeurs)/(2*np.pi)
    return dict(frequences=freq[:k].tolist(),
                bande_hz=[float(freq[0]/2), float((freq[k-1]+freq[k])/2)])


def modele_mbdyn(nb, moteur, bande):
    """Chaîne à pivots Y et ressorts, syntaxe du manuel public MBDyn."""
    l = m = 1/nb
    j = m*l*l/12
    k = min(6, nb-1)
    dimension = 17*nb+18  # 12 par nœud, 5 par pivot, 6 pour la base bloquée.
    methode = 'use lapack, balance, permute'
    if moteur == 'mbdyn_arpack':
        methode = f'use arpack, {2*k+4}, {max(32, 4*k+8)}, 1e-12'
    elif moteur == 'mbdyn_arpack_large':
        methode = f'use arpack, {dimension-4}, {dimension}, 1e-12'
    elif moteur != 'mbdyn_lapack':
        raise ValueError('moteur MBDyn inconnu')
    s = ['begin: data; problem: initial value; end: data;',
         'begin: initial value; initial time: 0.; final time: 0.001; '
         'time step: 0.001; tolerance: 1e-10; max iterations: 20; '
         'linear solver: umfpack; threads: disable; '
         'eigenanalysis: 0., output eigenvectors, results output precision, 17, '
         f'parameter, 0.02, lower frequency limit, {bande[0]:.17g}, '
         f'upper frequency limit, {bande[1]:.17g}, {methode}; end: initial value;',
         f'begin: control data; structural nodes: {nb+1}; rigid bodies: {nb+1}; '
         f'joints: {2*nb+1}; default output: none; end: control data;',
         'begin: nodes;']
    for i in range(nb+1):
        z = 0. if i == 0 else -(i-.5)*l
        s.append(f'structural: {i}, dynamic, 0., 0., {z:.17g}, eye, null, null;')
    s.extend(['end: nodes;', 'begin: elements;', 'body: 0, 0, 1., null, eye;',
              'joint: 0, clamp, 0, node, node;'])
    for i in range(1, nb+1):
        z = -(i-1)*l
        s.append(f'body: {i}, {i}, {m:.17g}, null, diag, {j:.17g}, {j:.17g}, {j:.17g};')
        position = (f'position, reference, global, 0.,0.,{z:.17g}, '
                    'position orientation, reference, global, eye, '
                    'rotation orientation, reference, global, eye')
        s.append(f'joint: {i}, total joint, {i-1}, {position}, {i}, {position}, '
                 'position constraint, 1,1,1, null, orientation constraint, 1,0,1, null;')
        s.append(f'joint: {nb+i}, deformable hinge, {i-1}, {i}, '
                 f'linear elastic generic, diag, 0., {10*nb}, 0.;')
    s.append('end: elements;')
    return '\n'.join(s)+'\n'


def lire_mbdyn(path):
    s = Path(path).read_text()
    match = re.search(r'dCoef\s*=\s*([^;]+);', s)
    alpha = re.search(r'alpha\s*=\s*\[(.*?)\];', s, re.S)
    if match is None or (alpha is None and 'alpha' in s):
        raise ValueError('sortie spectrale MBDyn absente ou mal formée')
    c = float(match[1])
    if not math.isfinite(c) or c <= 0:
        raise ValueError('coefficient de transformation incorrect')
    rows = ([[float(x) for x in row.split()] for row in alpha[1].split(';') if row.strip()]
            if alpha is not None else [])
    return dict(dCoef=c, alpha=rows)


def frequences_mbdyn(data):
    """Transformation homogène de Cayley, formule publique du manuel."""
    c = data['dCoef']
    values = []
    for row in data['alpha']:
        if len(row) != 3 or not all(math.isfinite(x) for x in row):
            raise ValueError('triplet spectral non fini')
        scale = max(map(abs, row))
        if scale == 0:
            continue  # triplet nul d'une équation algébrique, pas un mode fini.
        a, b, beta = (x/scale for x in row)
        den = (a+beta)**2+b*b
        if den == 0:
            continue  # valeur infinie, associée aux contraintes.
        real = (a*a+b*b-beta*beta)/den/c
        imag = 2*b*beta/den/c
        if imag > 0:
            values.append((imag/(2*math.pi), real))
    return sorted(values)


def juger(moteur, resultat, ref):
    """Ne classe que le spectre demandé et qualifié, sans apparier au plus proche."""
    if moteur.startswith('mbdyn'):
        if not math.isfinite(resultat['dCoef']) or resultat['dCoef'] <= 0:
            raise ValueError('coefficient de transformation incorrect')
        candidats = frequences_mbdyn(resultat)
        lo, hi = ref['bande_hz']
        selected = [(f, r) for f, r in candidats if lo < f < hi]
        freq = [x[0] for x in selected]
        real = [abs(x[1])/(2*math.pi*x[0]) for x in selected]
        hors_bande = len(candidats)-len(selected)
    else:
        values = resultat['valeurs_propres']
        if not all(math.isfinite(v) and v > 0 for v in values):
            return dict(decision='spectre_invalide')
        freq = sorted(math.sqrt(v)/(2*math.pi) for v in values)
        real = [0.]*len(freq)
        hors_bande = 0
    if len(freq) != len(ref['frequences']):
        return dict(decision='spectre_incomplet_ou_surnumeraire', frequences=freq,
                    attendu=len(ref['frequences']), hors_bande=hors_bande)
    error = max(abs(f/r-1) for f, r in zip(freq, ref['frequences']))
    damping = max(real)
    return dict(decision=('qualifie' if error <= SEUIL and damping <= SEUIL else 'imprecis'),
                frequences=freq, erreur_relative_frequence=error,
                partie_reelle_relative=damping, hors_bande=hors_bande)


def worker(moteur, nb, sortie):
    import importlib.metadata
    import numpy as np
    import scipy
    from modeles_modes import vinkulum, exudyn
    k = min(6, nb-1)
    start = time.perf_counter()
    if moteur == 'vinkulum':
        n = vinkulum('articulee', nb)
        built = time.perf_counter()
        result = n.modes_creux(k)
        values = [m['valeur_propre'] for m in result['modes']]
    else:
        sc, n = exudyn(nb, arbre='arbre' in moteur)
        built = time.perf_counter()
        values, _ = n.ComputeODE2Eigenvalues(numberOfEigenvalues=k,
                    useSparseSolver=moteur.endswith('creux'), convert2Frequencies=False)
        values = np.asarray(values).tolist()
    end = time.perf_counter()
    package = 'vinkulum' if moteur == 'vinkulum' else 'exudyn'
    ecrire(sortie, dict(valeurs_propres=values, construction_s=built-start,
                       analyse_s=end-built, version=importlib.metadata.version(package),
                       numpy=np.__version__, scipy=scipy.__version__))


def limites():
    resource.setrlimit(resource.RLIMIT_AS, (MEMOIRE, MEMOIRE))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def commande_mbdyn(binaire):
    # Un chemin absolu avec des points est tronqué par la sortie spectrale
    # du binaire qualifié. Une extension explicite réserve le dernier point
    # au nom de fichier, même après conversion interne en chemin absolu.
    return [str(Path(binaire).resolve()), '-f', 'modele.mbd', '-o', 'natif.out']


def executer(moteur, nb, repetition, ref, args, dossier):
    trial = dossier/f'{moteur}-{nb}-{repetition}'
    trial.mkdir()
    raw = trial/'resultat.json'
    if moteur.startswith('mbdyn'):
        deck = trial/'modele.mbd'
        deck.write_text(modele_mbdyn(nb, moteur, ref['bande_hz']))
        command = commande_mbdyn(args.mbdyn)
    else:
        python = args.python_vinkulum if moteur == 'vinkulum' else args.python_exudyn
        command = [str(Path(python).absolute()), '-I', '-B', str(Path(__file__).resolve()),
                   '--worker', moteur, '--taille', str(nb), '--sortie', str(raw)]
    env = dict(os.environ)
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    env['PYTHONHASHSEED'] = '0'
    started = time.perf_counter()
    with (trial/'stdout.log').open('w') as log:
        p = subprocess.Popen(['/usr/bin/time', '-f', '%M', '-o', str(trial/'rss.txt'), *command],
                             cwd=trial, env=env, stdout=log, stderr=subprocess.STDOUT,
                             start_new_session=True, preexec_fn=limites)
        expired = threading.Event()
        def interrompre():
            if p.poll() is None:
                expired.set()
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        timer = threading.Timer(DELAI, interrompre)
        timer.start()
        try:
            code = p.wait()
        except BaseException:
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
            raise
        finally:
            timer.cancel()
            timer.join()
        timeout = expired.is_set()
    elapsed = time.perf_counter()-started
    rss_text = (trial/'rss.txt').read_text() if (trial/'rss.txt').exists() else ''
    rss_rows = [int(x) for x in rss_text.splitlines() if x.isdigit()]
    line = dict(moteur=moteur, taille=nb, repetition=repetition,
                echauffement=(repetition == 0), secondes_processus=elapsed,
                rss_kio=max(rss_rows) if rss_rows else None, code_retour=code,
                charge_systeme=list(os.getloadavg()))
    if timeout:
        line.update(decision='delai_depasse')
    elif code != 0:
        line.update(decision='echec_execution')
    else:
        try:
            result = lire_mbdyn(trial/'natif_0.m') if moteur.startswith('mbdyn') else json.loads(raw.read_text())
            line.update(resultat=result, **juger(moteur, result, ref))
        except (ValueError, OSError, KeyError) as exc:
            line.update(decision='sortie_invalide', raison=str(exc))
    # Seuls les résultats numériques sont publiés. Les journaux natifs MBDyn
    # peuvent inclure l'environnement du processus : ils restent hors archive.
    return line


def identite_python(python, package):
    script = """
import importlib, importlib.metadata, hashlib, json, sys
from pathlib import Path
import numpy, scipy
m=importlib.import_module(sys.argv[1])
p=Path(m.__file__).parent
print(json.dumps(dict(version=importlib.metadata.version(sys.argv[1]),
    python=sys.version, numpy=numpy.__version__, scipy=scipy.__version__,
    extensions={str(q.relative_to(p)):hashlib.sha256(q.read_bytes()).hexdigest()
                for q in p.rglob('*.so')})))
"""
    return json.loads(subprocess.check_output([python, '-I', '-c', script, package], text=True))


def identite_mbdyn(binaire):
    binaire = Path(binaire).resolve()
    racine = binaire.parent.parent
    result = dict(sha256=sha(binaire),
        version=subprocess.check_output([str(binaire), '--version'], text=True).strip())
    makefile = racine/'Makefile'
    if makefile.exists():
        # Configuration de compilation uniquement, sans code du solveur.
        champs = ('CFLAGS', 'CXXFLAGS', 'FFLAGS', 'FCFLAGS', 'ARPACK_LIBS', 'UMFPACK_LIBS', 'BLAS_LIBS', 'LAPACK_LIBS')
        result['compilation'] = {k: match[1] for k in champs
            if (match := re.search(r'^'+k+r' = (.*)$', makefile.read_text(), re.M))}
    git = subprocess.run(['git','rev-parse','HEAD'],cwd=racine,text=True,
                         stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    result['revision_sources'] = git.stdout.strip() if git.returncode == 0 else None
    return result


def campagne(args):
    import importlib.metadata
    dossier = Path(args.sortie).resolve()
    dossier.mkdir(parents=True, exist_ok=True)
    destination = dossier/'campagne.json'
    if destination.exists():
        raise ValueError('campagne existante : ne pas écraser des observations')
    cpus = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpus[0]})
    references = {str(n): oracle(n) for n in TAILLES}
    manifest = dict(schema=SCHEMA, seuil=SEUIL, tailles=list(TAILLES), moteurs=list(MOTEURS),
                    repetitions=REPETITIONS, delai_s=DELAI, memoire_octets=MEMOIRE,
                    cpu=cpus[0], plateforme=platform.platform(), python=platform.python_version(),
                    sources={p.name: sha(p) for p in (Path(__file__), ROOT/'ci/modeles_modes.py')},
                    mbdyn=identite_mbdyn(args.mbdyn),
                    cpu_modele=next((l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),None),
                    identites=dict(vinkulum=identite_python(args.python_vinkulum,'vinkulum'),
                                   exudyn=identite_python(args.python_exudyn,'exudyn')),
                    references=references, essais=[])
    ecrire(destination, manifest)
    for nb in TAILLES:
        for repetition in range(REPETITIONS+1):
            # Rotation déterministe des moteurs pour limiter un biais d'ordre.
            shift = repetition % len(MOTEURS)
            for moteur in MOTEURS[shift:]+MOTEURS[:shift]:
                line = executer(moteur, nb, repetition, references[str(nb)], args, dossier)
                manifest['essais'].append(line)
                ecrire(destination, manifest)
                print(json.dumps({k: v for k, v in line.items() if k not in ('resultat', 'frequences')}, allow_nan=False), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker', choices=MOTEURS)
    p.add_argument('--taille', type=int)
    p.add_argument('--sortie', required=True)
    p.add_argument('--mbdyn')
    p.add_argument('--python-vinkulum')
    p.add_argument('--python-exudyn')
    a = p.parse_args()
    if a.worker:
        worker(a.worker, a.taille, a.sortie)
    elif a.mbdyn and a.python_vinkulum and a.python_exudyn:
        campagne(a)
    else:
        p.error('les trois exécutables sont requis')
