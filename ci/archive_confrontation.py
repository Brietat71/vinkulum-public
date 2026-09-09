"""Archive ou vérifie une confrontation, ses journaux et ses sources de mesure.

Les fichiers .mov volumineux restent dans le répertoire de calcul : leurs
empreintes sont vérifiées à l'archivage, les observables échantillonnés sont
conservés sans arrondi dans l'archive des trajectoires. Le juge peut ainsi
être rejoué hors de la machine qui possède les deux exécutables.
"""
import argparse
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

from bilan_confrontation import analyse, load_report

ROOT=Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def payload(value):
    return json.dumps(value,ensure_ascii=False,separators=(',',':')).encode()


def compact(report):
    result=deepcopy(report)
    trajectories=result['trajectoires']={}
    for case in result['cases'].values():
        for run in case['runs']:
            if not run['ok']:
                continue
            trajectory={key:run.pop(key) for key in ('times','values')}
            # This canonical encoding is also used by load_report().
            digest=sha(json.dumps(trajectory,separators=(',',':')).encode())
            trajectories[digest]=trajectory
            run['trajectory_sha256']=digest
    return result


def archive(args):
    report=json.loads((args.calcul/'resultats.json').read_text())
    if args.probes:
        report['exploratory_probes']=json.loads(args.probes.read_text())
    metadata=report['metadata']
    version=metadata['vinkulum_version']
    release=json.loads((ROOT/f'docs/bancs/version-{version}.json').read_text())
    if sha(args.roue.read_bytes())!=release['roue']['sha256']:
        raise ValueError('la roue ne correspond pas à la livraison archivée')
    with zipfile.ZipFile(args.roue) as wheel:
        extensions=[n for n in wheel.namelist() if n.endswith('.so')]
        if len(extensions)!=1 or sha(wheel.read(extensions[0]))!=metadata['extension_sha256']:
            raise ValueError('le binaire mesuré ne correspond pas à la roue')
        if sha(wheel.read('vinkulum/andrews.py'))!=metadata['andrews_model_sha256']:
            raise ValueError('modèle Andrews différent de celui de la roue')

    journal=dict(definition=__doc__,contenus={},essais=[],sources={})
    def keep(data):
        digest=sha(data)
        journal['contenus'][digest]=data.decode('utf8')
        return digest
    for name,key in [('ci/confronte_mbdyn.py','script_sha256'),('ci/modeles_confrontation.py','models_sha256')]:
        data=(ROOT/name).read_bytes()
        if sha(data)!=metadata[key]:
            raise ValueError(f'source de mesure changée : {name}')
    for name in ('ci/confronte_mbdyn.py','ci/modeles_confrontation.py','ci/bilan_confrontation.py',
                 'ci/test_confrontation.py','ci/archive_confrontation.py','ci/trace_confrontation.py'):
        journal['sources'][name]=keep((ROOT/name).read_bytes())
    runs=[r for c in report['cases'].values() for r in c['runs']]+report.get('exploratory_probes',[])
    for run in runs:
        entry=dict(directory=run['directory'],fichiers={},non_embarques={})
        for name,digest in run['artifact_sha256'].items():
            data=(Path(run['directory'])/name).read_bytes()
            if sha(data)!=digest:
                raise ValueError(f'artefact de calcul changé : {run["directory"]}/{name}')
            if name in ('model','stdout.txt','run.out') or Path(name).suffix in ('.ref','.set','.nod','.elm'):
                entry['fichiers'][name]=keep(data)
            else:
                entry['non_embarques'][name]=digest
        journal['essais'].append(entry)
    products={'.json':(json.dumps(analyse(report),indent=2,ensure_ascii=False)+'\n').encode(),
              '-trajectoires.json.gz':gzip.compress(payload(compact(report)),mtime=0),
              '-journaux.json.gz':gzip.compress(payload(journal),mtime=0)}
    manifest=dict(version=version,roue=release['roue'],extension_sha256=metadata['extension_sha256'],
                  release_commit=subprocess.check_output(['git','rev-list','-n','1',f'v{version}'],cwd=ROOT,text=True).strip(),
                  scripts_sha256={name:sha((ROOT/name).read_bytes()) for name in journal['sources']},
                  fichiers={Path(str(args.sortie)+suffix).name:dict(sha256=sha(data),octets=len(data)) for suffix,data in products.items()})
    products['-manifest.json']=(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n').encode()
    paths=[Path(str(args.sortie)+suffix) for suffix in products]
    if any(path.exists() for path in paths):
        raise FileExistsError('archive existante : choisir un nouveau préfixe de sortie')
    args.sortie.parent.mkdir(parents=True,exist_ok=True)
    for suffix,data in products.items():
        Path(str(args.sortie)+suffix).write_bytes(data)
    verify(args.sortie)


def verify(prefix):
    manifest=json.loads(Path(str(prefix)+'-manifest.json').read_text())
    for name,info in manifest['fichiers'].items():
        data=(prefix.parent/name).read_bytes()
        if sha(data)!=info['sha256'] or len(data)!=info['octets']:
            raise ValueError(f'empreinte d’archive invalide : {name}')
    report=load_report(Path(str(prefix)+'-trajectoires.json.gz'))
    if analyse(report)!=json.loads(Path(str(prefix)+'.json').read_text()):
        raise ValueError('bilan différent du recalcul des trajectoires')
    journal=json.loads(gzip.decompress(Path(str(prefix)+'-journaux.json.gz').read_bytes()))
    for digest,data in journal['contenus'].items():
        if sha(data.encode())!=digest:
            raise ValueError('contenu de journal altéré')
    if journal['sources']!=manifest['scripts_sha256']:
        raise ValueError('sources de mesure non liées au manifeste')
    if any(digest not in journal['contenus'] for digest in journal['sources'].values()):
        raise ValueError('source de mesure manquante')
    runs=[r for c in report['cases'].values() for r in c['runs']]+report.get('exploratory_probes',[])
    if len(runs)!=len(journal['essais']):
        raise ValueError('journaux manquants')
    for run,entry in zip(runs,journal['essais']):
        if entry['directory']!=run['directory'] or dict(entry['fichiers'],**entry['non_embarques'])!=run['artifact_sha256']:
            raise ValueError('journal non lié à son essai')
        if any(digest not in journal['contenus'] for digest in entry['fichiers'].values()):
            raise ValueError('contenu de journal manquant')
    if report['metadata']['extension_sha256']!=manifest['extension_sha256']:
        raise ValueError('extension non liée au manifeste')
    print(f'{len(runs)} essais, trajectoires, bilans, journaux et empreintes vérifiés.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verifier',type=Path,help='Préfixe d’archive à vérifier, sans les exécutables')
    p.add_argument('--calcul',type=Path)
    p.add_argument('--probes',type=Path)
    p.add_argument('--roue',type=Path)
    p.add_argument('--sortie',type=Path)
    args=p.parse_args()
    if args.verifier:
        verify(args.verifier)
    elif all((args.calcul,args.roue,args.sortie)):
        archive(args)
    else:
        p.error('--calcul, --roue et --sortie sont requis pour archiver')


if __name__=='__main__':
    main()
