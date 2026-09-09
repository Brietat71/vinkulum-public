"""Ordonnanceur interne des validations ; processus isolés et budget CPU partagé.

VINKULUM_BANCS_JOBS : nombre de processus (défaut min(8, CPU disponibles)).
VINKULUM_BANCS_CPUS : budget total (défaut affinité du processus).
VINKULUM_BANCS_TIMEOUT : délai maximal par cas en secondes (défaut 900).
VINKULUM_BANCS_LOGS : répertoire parent des journaux (défaut temporaire).
Les mesures de temps par cas sous concurrence ne sont pas des benchmarks isolés.
"""
import concurrent.futures
import importlib
import json
import math
import os
from pathlib import Path
import pickle
import queue
import signal
import subprocess
import sys
import tempfile
import threading
import time
import traceback


def tache(nom, module, fonction, *args, garder=False, exclusif=False, **kwargs):
    return dict(nom=nom, module=module, fonction=fonction, args=args,
                kwargs=kwargs, garder=garder, exclusif=exclusif)


def _appel(t):
    fonction = getattr(importlib.import_module(t['module']), t['fonction'])
    return fonction(*t['args'], **t['kwargs'])


def _entier(nom, defaut):
    valeur = int(os.environ.get(nom, defaut))
    if valeur < 1:
        raise ValueError(f'{nom} doit être >= 1')
    return valeur


def _configuration(nombre, defaut_jobs=None):
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else list(range(os.process_cpu_count() or 1))
    budget = min(len(cpus), _entier('VINKULUM_BANCS_CPUS', len(cpus)))
    jobs = min(nombre, budget, _entier('VINKULUM_BANCS_JOBS', defaut_jobs or min(8, budget)))
    timeout = float(os.environ.get('VINKULUM_BANCS_TIMEOUT', '900'))
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('VINKULUM_BANCS_TIMEOUT doit être fini et > 0')
    if int(os.environ.get('RAYON_NUM_THREADS', '0') or '0') < 0:
        raise ValueError('RAYON_NUM_THREADS doit être >= 0')
    return [cpus[i:budget:jobs] for i in range(jobs)], timeout


def _arrete(p):
    if p.poll() is not None:
        return
    try:
        if os.name == 'posix':
            os.killpg(p.pid, signal.SIGTERM)
        else:
            p.terminate()
        p.wait(timeout=2)
    except subprocess.TimeoutExpired:
        if os.name == 'posix':
            os.killpg(p.pid, signal.SIGKILL)
        else:
            p.kill()
        p.wait()
    except ProcessLookupError:
        pass


def executer(taches, titre, defaut_jobs=None):
    if len({t['nom'] for t in taches}) != len(taches):
        raise ValueError('noms de cas dupliqués')
    if not taches:
        return {}
    # Un cas composite déjà placé dans un worker garde ses sous-cas dans ce
    # processus : aucun deuxième niveau de pools ou de réservation de CPU.
    if os.environ.get('VINKULUM_BANCS_WORKER') == '1':
        return {t['nom']: _appel(t) for t in taches}
    groupes, timeout = _configuration(len(taches), defaut_jobs)
    base = os.environ.get('VINKULUM_BANCS_LOGS')
    if base:
        Path(base).mkdir(parents=True, exist_ok=True)
    dossier = Path(tempfile.mkdtemp(prefix='vinkulum-campagne-', dir=base))
    ressources = queue.Queue()
    for cpus in groupes:
        ressources.put(cpus)
    actifs, verrou = set(), threading.Lock()
    annule = threading.Event()
    debut = time.perf_counter()
    print(f'╔═ {titre} : {len(taches)} cas, {len(groupes)} processus, '
          f'{sum(map(len, groupes))} CPU au total', flush=True)
    print(f'║ journaux : {dossier}', flush=True)
    if len(groupes) > 1:
        print('║ durées sous concurrence ; JOBS=1 pour mesurer des cas isolés', flush=True)

    def lance(index, t, reservation=None):
        cpus = ressources.get() if reservation is None else reservation
        threads = min(len(cpus), int(os.environ.get('RAYON_NUM_THREADS', '0') or '0') or len(cpus))
        fichier = dossier / f'{index:03d}.log'
        resultat = dossier / f'{index:03d}.pickle'
        env = dict(os.environ, VINKULUM_BANCS_WORKER='1',
                   RAYON_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1',
                   OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
                   VECLIB_MAXIMUM_THREADS='1', PYTHONUNBUFFERED='1')
        d = time.perf_counter()
        valeur, statut, erreur = None, 'ok', None
        try:
            if annule.is_set():
                raise RuntimeError('campagne interrompue')
            with fichier.open('w') as log:
                p = subprocess.Popen([sys.executable, '-m', 'vinkulum._campagnes',
                                      '--worker', str(resultat), ','.join(map(str, cpus))],
                                     stdin=subprocess.PIPE, stdout=log, stderr=subprocess.STDOUT,
                                     text=True, env=env, start_new_session=os.name == 'posix')
                with verrou:
                    actifs.add(p)
                try:
                    if annule.is_set():
                        _arrete(p)
                        raise RuntimeError('campagne interrompue')
                    p.communicate(json.dumps(t), timeout=timeout)
                    if p.returncode:
                        raise RuntimeError(f'code de sortie {p.returncode}')
                except BaseException:
                    _arrete(p)
                    raise
                finally:
                    with verrou:
                        actifs.discard(p)
            # Fichier produit uniquement par notre processus enfant, dans
            # un répertoire temporaire privé ; jamais une entrée utilisateur.
            with resultat.open('rb') as source:
                valeur = pickle.load(source)
        except Exception as e:
            statut, erreur = 'echec', str(e)
        finally:
            if reservation is None:
                ressources.put(cpus)
        return valeur, dict(nom=t['nom'], statut=statut, secondes=time.perf_counter()-d,
                            cpus=cpus, threads=threads, log=str(fichier), erreur=erreur,
                            exclusif=t['exclusif'])

    sorties, bilans = {}, {}
    def recoit(t, resultat):
        valeur, bilan = resultat
        sorties[t['nom']], bilans[t['nom']] = valeur, bilan
        print(f"║ {bilan['statut'].upper():5s} {t['nom']} — {bilan['secondes']:.2f} s", flush=True)
        if bilan['statut'] != 'ok':
            print(bilan['erreur'], flush=True)
            print(Path(bilan['log']).read_text() if Path(bilan['log']).exists() else '', flush=True)

    pool = concurrent.futures.ThreadPoolExecutor(max_workers=len(groupes))
    try:
        futures = {pool.submit(lance, i, t): t for i, t in enumerate(taches) if not t['exclusif']}
        for future in concurrent.futures.as_completed(futures):
            recoit(futures[future], future.result())
        for i, t in enumerate(taches):
            if t['exclusif']:
                recoit(t, lance(i, t, sorted(c for groupe in groupes for c in groupe)))
    except BaseException:
        annule.set()
        with verrou:
            for p in list(actifs):
                _arrete(p)
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
        rapport = dict(titre=titre, jobs=len(groupes), cpu_budget=sum(map(len, groupes)),
                       secondes=time.perf_counter()-debut,
                       cas=[bilans[t['nom']] for t in taches if t['nom'] in bilans],
                       prevus=[t['nom'] for t in taches])
        (dossier / 'rapport.json').write_text(json.dumps(rapport, indent=2)+'\n')
    echecs = [b['nom'] for b in bilans.values() if b['statut'] != 'ok']
    if echecs:
        raise RuntimeError(f"{titre} : {len(echecs)} échec(s) : {', '.join(echecs)} ; journaux {dossier}")
    print(f"╚ {titre} OK — {len(bilans)}/{len(taches)} cas, {rapport['secondes']:.2f} s", flush=True)
    return {t['nom']: sorties[t['nom']] for t in taches}


def _worker():
    destination, cpus = sys.argv[2:]
    if hasattr(os, 'sched_setaffinity'):
        os.sched_setaffinity(0, {int(c) for c in cpus.split(',')})
    t = json.load(sys.stdin)
    valeur = _appel(t)
    with open(destination, 'wb') as sortie:
        pickle.dump(valeur if t['garder'] else None, sortie)


if __name__ == '__main__':
    try:
        if len(sys.argv) != 4 or sys.argv[1] != '--worker':
            raise ValueError('point d’entrée réservé aux campagnes')
        _worker()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
