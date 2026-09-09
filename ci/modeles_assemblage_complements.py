"""Trois branches flexibles distinctes reliées à une jonction rigide linéarisée."""
import json
from pathlib import Path
import numpy as np
from scipy.sparse import block_diag,csr_matrix
from confronte_ports_exudyn import lire,sparse_pack
from experience_controle_facteurs import sha


def geometrie():
    rotations=[np.eye(3),np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]]),
               np.array([[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]])]
    offsets=[np.array([0.,0.,0.]),np.array([.125,0.,.0625]),np.array([0.,-.125,.0625])]
    apps=[]
    for r,v in zip(rotations,offsets,strict=True):
        x,y,z=v;cross=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
        a=np.zeros((6,6));a[:3,:3]=r;a[:3,3:]=-r@cross;a[3:,3:]=r;apps.append(a)
    return apps


def creer_cas(source,destination,n):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=False)
    entry=Path(source)/f'n{n}-f40.npz';manifest=json.loads((Path(source)/'manifest.json').read_text())
    if sha(entry)!=manifest[entry.name]:raise ValueError('modèle à masse consistante altéré')
    d,m,metric,forces,omega=lire(entry);apps=geometrie();h=np.zeros((6,6))
    for j,(sd,sm,a) in enumerate(zip((1.,1.125,1.25),(1.,.875,1.125),apps,strict=True)):
        dj,mj=(sd*d).tocsr(),(sm*m).tocsr();hj=sd**2*metric
        np.savez_compressed(destination/f'piece{j}.npz',**sparse_pack('d',dj),**sparse_pack('m',mj),
                            metrique=hj,forces=forces,omega=omega)
        h+=a.T@hj@a
    h=(h+h.T)*.5
    de=np.diag([2.,3.,4.,.125,.25,.375]);me=np.diag([.125]*3+[.000125]*3)
    np.savez_compressed(destination/'jonction.npz',applications=np.asarray(apps),metrique=h,
                        facteur_externe=de,masse_externe=me,forces=forces,omega=omega)
    info=dict(n=n,source_sha256=sha(entry),source=str(entry),
        modele='trois branches linéaires encastrées, jonction rigide avec bras de levier',
        facteurs_D=[1.,1.125,1.25],facteurs_M=[1.,.875,1.125],
        comptabilite='énergie de chaque branche plus énergie de jonction ajoutée une fois',
        metrique='somme A_j.T H_j A_j, évaluée puis stockée en binary64',
        fichiers_sha256={p.name:sha(p) for p in destination.glob('*.npz')})
    (destination/'modele.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))


def charger_cas(path):
    path=Path(path);info=json.loads((path/'modele.json').read_text())
    for name,expected in info['fichiers_sha256'].items():
        if sha(path/name)!=expected:raise ValueError('entrée de branche ou jonction altérée')
    pieces=[lire(path/f'piece{j}.npz')[:3] for j in range(3)]
    with np.load(path/'jonction.npz',allow_pickle=False) as z:
        jonction={k:z[k] for k in z.files}
    return pieces,jonction


def espace_juge(pieces,j):
    # Copie de chaque champ de branche et champ physique de jonction terminal.
    # La somme d'énergies est exactement la définition du modèle assemblé.
    return (block_diag([p[0] for p in pieces]+[csr_matrix(j['facteur_externe'])],format='csr'),
            block_diag([p[1] for p in pieces]+[csr_matrix(j['masse_externe'])],format='csr'))


def empiler(rep):return np.vstack(rep['champs']+[rep['champ_ports']])
