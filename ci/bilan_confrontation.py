"""Réduit les trajectoires de confronte_mbdyn en mesures auditables.

Erreur : maximum absolu sur les composantes et les instants échantillonnés.
Chaque candidat est comparé aux DEUX références fines ; aucune n'est tenue
pour une vérité exacte. Le classement exige deux raffinements réussis par
moteur et une variation de référence inférieure à 10 % du seuil demandé.
"""
import argparse
import hashlib
import gzip
import json
from pathlib import Path
import statistics

import numpy as np


def load_report(path):
    data = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
    report = json.loads(data)
    if 'trajectoires' in report:
        for case in report['cases'].values():
            for run in case['runs']:
                if 'trajectory_sha256' in run:
                    trajectory = report['trajectoires'][run['trajectory_sha256']]
                    payload = json.dumps(trajectory,separators=(',',':')).encode()
                    if hashlib.sha256(payload).hexdigest() != run['trajectory_sha256']:
                        raise ValueError('empreinte de trajectoire invalide')
                    run.update(trajectory)
    return report


def difference(a, b):
    ta,tb=np.asarray(a['times']),np.asarray(b['times'])
    if ta.ndim!=1 or ta.shape!=tb.shape or ta.size==0 or not np.isfinite(ta).all() or not np.isfinite(tb).all() or not np.all(np.diff(ta)>0) or not np.all(np.diff(tb)>0) or not np.allclose(ta,tb,rtol=0.,atol=1e-12):
        raise ValueError('grilles temporelles différentes')
    x, y = np.asarray(a['values']), np.asarray(b['values'])
    if x.ndim!=2 or x.shape != y.shape or x.shape[0]!=ta.size or x.size==0 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('trajectoires invalides')
    return float(np.max(np.abs(x-y)))


def analyse(report):
    out = dict(metadata=report['metadata'], definition=__doc__, cases={})
    repetitions=report['metadata'].get('repetitions')
    if report['metadata'].get('protocol_version',1)>=2 and (type(repetitions) is not int or repetitions<1):
        raise ValueError('nombre de répétitions manquant ou invalide')
    for name, case in report['cases'].items():
        spec, runs = case['spec'], case['runs']
        if report['metadata'].get('protocol_version',1)>=2 and 'formulation' in spec:
            for r in runs:
                if r['engine']!='vinkulum' or not r['ok']:
                    continue
                ds=r.get('static_diagnostics',[])
                if len(ds)!=50 or [d['palier_charge'] for d in ds]!=list(range(1,51)) or any(
                    d['statut']!='tolerance' or d['strict'] is not True or not d['tolerance_finale_atteinte']
                    or d['paliers_stagnation']!=0 or d['tol']>1e-8
                    or not 0<=d['residu_relatif']<=d['tol'] or not 0<=d['contraintes']<=d['tol'] for d in ds):
                    raise ValueError('Princeton : cinquante paliers stricts non vérifiés')
        identities=[(r['engine'],r['step'],r['rep']) for r in runs]
        if report['metadata'].get('protocol_version',1)>=2 and len(set(identities))!=len(identities):
            raise ValueError('répétition dupliquée')
        result = dict(spec=spec, configurations=[], failures=[{k:r[k] for k in ('engine','step','rep','error','directory')} for r in runs if not r['ok'] and r['rep']>=0],
                      warmup_failures=[{k:r[k] for k in ('engine','step','rep','error','directory')} for r in runs if not r['ok'] and r['rep']<0])
        if 'reference_policy' in case:
            result['reference_policy']=case['reference_policy']
        if 'omitted_references' in case:
            result['omitted_references']=case['omitted_references']
        out['cases'][name] = result
        refs = {}
        for engine in ('vinkulum', 'mbdyn'):
            ref_steps=spec.get('ref_by_engine',{}).get(engine,spec['ref'])
            selected = [next((r for r in runs if r['engine']==engine and r['step']==step and r['rep']==0 and r['ok']), None) for step in ref_steps]
            if all(r is not None for r in selected):
                refs[engine] = selected[-1]
                result[engine+'_reference_change'] = difference(*selected)
        result['reference_cross_gap'] = difference(*refs.values()) if len(refs)==2 else None
        result['references_qualified'] = len(refs)==2 and all(result[e+'_reference_change'] <= .1*spec['target'] for e in refs) and result['reference_cross_gap'] <= .2*spec['target']
        for engine in ('vinkulum', 'mbdyn'):
            ref_steps=spec.get('ref_by_engine',{}).get(engine,spec['ref'])
            for step in dict.fromkeys(spec['steps']+ref_steps):
                samples = [r for r in runs if r['engine']==engine and r['step']==step and r['rep']>=0]
                good = [r for r in samples if r['ok']]
                c = dict(engine=engine,step=step,successes=len(good),attempts=len(samples),candidate=step in spec['steps'])
                expected=repetitions if c['candidate'] else 1
                c['repetitions_complete']=expected is None or sorted(r['rep'] for r in samples)==list(range(expected))
                if good:
                    c.update(wall_seconds=[r['wall_seconds'] for r in good], median_seconds=statistics.median(r['wall_seconds'] for r in good),
                             rss_kib=[r['rss_kib'] for r in good], repeat_difference=max(difference(good[0],r) for r in good),
                             errors_to_references={e:max(difference(r,ref) for r in good) for e,ref in refs.items()},
                             trajectory_sha256=hashlib.sha256(json.dumps(good[0]['values']).encode()).hexdigest(),
                             final_values=good[0]['values'][-1])
                    c['error_envelope'] = max(c['errors_to_references'].values(),default=None)
                    c['reference_variation_margin'] = max((result[e+'_reference_change'] for e in refs),default=0.)
                    c['estimated_error_with_margin'] = c['error_envelope']+c['reference_variation_margin'] if c['error_envelope'] is not None else None
                    c['eligible'] = result['references_qualified'] and len(good)==len(samples) and c['repetitions_complete'] and c['candidate'] and c['estimated_error_with_margin']<=spec['target']
                result['configurations'].append(c)
            choices=[c for c in result['configurations'] if c['engine']==engine and c.get('eligible')]
            result[engine+'_best'] = min(choices,key=lambda c:c['median_seconds']) if choices else None
        if all(result[e+'_best'] for e in ('mbdyn','vinkulum')):
            result['mbdyn_over_vinkulum_time'] = result['mbdyn_best']['median_seconds']/result['vinkulum_best']['median_seconds']
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path); p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(analyse(load_report(args.input)),indent=2,ensure_ascii=False)+'\n')
