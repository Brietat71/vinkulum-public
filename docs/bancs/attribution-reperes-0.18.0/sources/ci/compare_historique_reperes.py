"""Rattache les traces instrumentées aux maxima publiés dans la 0.14.1."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


def comparer(traces,archive):
    archive=Path(archive);payload=archive.read_bytes()
    files=json.loads(gzip.decompress(payload))['fichiers']
    historical=json.loads(files['reperes-candidat.json'])
    out=[]
    for experiment in historical['essais']:
        if experiment['duree']!=20.:continue
        h=experiment['pas'];stride=max(1,round(20./h)//200)
        def read(frame):
            path=Path(traces)/f'{frame}-{h:.5f}.json.gz'
            return json.loads(gzip.decompress(path.read_bytes()))
        base=read('initial')['trace']['echantillons']
        for record in experiment['mesures']:
            frame=record['cadre'];data=read(frame)
            world,material=np.asarray(data['monde']),np.asarray(data['matiere'])
            values=[]
            for i in range(stride,len(base),stride):
                b=data['trace']['echantillons'][i];a=base[i]
                dr=world.T@np.asarray(b[2][0]).reshape(3,3)@material.T-np.asarray(a[2][0]).reshape(3,3)
                dw=world.T@np.asarray(b[3][0])-np.asarray(a[3][0])
                values.append([float(np.linalg.norm(dr)),float(np.linalg.norm(dw))])
            maxima=np.max(values,axis=0).tolist()
            old=record['ecarts_max']
            # Les maxima historiques sont des valeurs observées, pas des seuils
            # d'acceptation physique. Conserver toute différence, sans tolérance.
            out.append(dict(cadre=frame,pas=h,historique=old,trace=maxima,
                            identiques=maxima==old,differences=[a-b for a,b in zip(maxima,old)]))
    return dict(schema='vinkulum.diagnostic.rattachement_reperes.1',
                archive_sha256=hashlib.sha256(payload).hexdigest(),
                programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                version_historique=historical['version'],comparaisons=out,
                toutes_identiques=all(x['identiques'] for x in out))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('traces');p.add_argument('archive');p.add_argument('--sortie',required=True)
    a=p.parse_args();d=comparer(a.traces,a.archive)
    Path(a.sortie).write_text(json.dumps(d,indent=2)+'\n')
    print('maxima historiques tous identiques:',d['toutes_identiques'])
