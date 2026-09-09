"""Capture les moments portés et confronte les sorties à la roue publiée.

Les traces ne sont pas des certificats ; le calcul rationnel des défauts
et de leur propagation est une étape distincte.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from vinkulum import Noyau, __version__, _vinkulum
import diagnostic_rotation as d

PAS = (.001, .0005, .00025)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def modele(cadre):
    monde, matiere = d.CADRES[cadre]
    n=Noyau([0.,0.,0.])
    n.corps('corps',.7,(matiere.T@d.J@matiere).ravel().tolist(),[0.,0.,0.],
            rot=(monde@d.R0@matiere).ravel().tolist(),w=(monde@d.W0).tolist())
    return n


def reguliere(n,h):
    initial=n.etat_precis()
    rows=[(initial[0],initial[2],initial[4])]
    rows.extend((r[0],r[2],r[4]) for r in n.simule_em(20.,h,tous=1))
    return rows


def matrice(rows):
    return np.array([[t,*r[0],*w[0]] for t,r,w in rows],dtype=np.float64)


def comparer(a,b):
    a,b=matrice(a),matrice(b)
    if a.shape!=b.shape:
        raise ValueError('grilles différentes')
    bits=int(np.count_nonzero(a.view(np.uint64)!=b.view(np.uint64)))
    return dict(composantes_differentes_bits=bits,ecart_absolu_max=float(np.max(abs(a-b),initial=0.)))


def produire(out,baseline=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'manifest.json').exists():raise ValueError('manifest existant')
    manifest=dict(schema='vinkulum.diagnostic.trace_em.1',version=__version__,
                  extension_sha256=sha(_vinkulum.__file__),programme_sha256=sha(__file__),
                  modele_sha256=sha(d.__file__),propagation_certifiee=False,cas=[])
    for h in PAS:
        for cadre in d.CADRES:
            n=modele(cadre)
            before=n.etat_precis()
            name=f'{cadre}-{h:.5f}.json.gz'
            record=dict(cadre=cadre,pas_demande=h,duree=20.,
                        monde=d.CADRES[cadre][0].tolist(),matiere=d.CADRES[cadre][1].tolist())
            if baseline:
                trace=n._trace_em(20.,h)
                if n.etat_precis()!=before:raise ValueError('instrumentation mutante')
                rows=[(r[0],r[2],r[3]) for r in trace['echantillons']]
                regular=reguliere(n,h)
                local=comparer(rows,regular)
                old=json.loads(gzip.decompress((Path(baseline)/name).read_bytes()))
                previous=comparer(rows,old['sorties'])
                record.update(trace=trace,comparaison_meme_roue=local,comparaison_roue_publiee=previous)
                if local['composantes_differentes_bits'] or previous['composantes_differentes_bits']:
                    print('DIFFÉRENCE',cadre,h,local,previous,flush=True)
            else:
                rows=reguliere(n,h)
                record.update(sorties=rows)
            if rows[-1][0]!=20.:raise ValueError('trajectoire incomplète')
            data=gzip.compress(json.dumps(record,allow_nan=False,separators=(',',':')).encode(),mtime=0)
            (out/name).write_bytes(data)
            manifest['cas'].append(dict(fichier=name,sha256=sha(out/name),points=len(rows),
                **({k:record[k] for k in ('comparaison_meme_roue','comparaison_roue_publiee')} if baseline else {})))
            (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            print(cadre,h,len(rows),len(data),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie',required=True);p.add_argument('--baseline')
    a=p.parse_args();produire(a.sortie,a.baseline)
