"""Défauts de conformité des champs conservés, hors mesures et sans certificat."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky,qr,solve_triangular,svdvals
from modeles_assemblage_complements import charger_cas
from experience_controle_facteurs import sha,ecrire


def auditer(campagne,references,out):
    out.mkdir(parents=True,exist_ok=False);shutil.copy2(__file__,out/Path(__file__).name);results=[]
    for r in json.loads((campagne/'bilan.json').read_text()):
        path=campagne/r['dossier']/'champs.npy'
        if sha(path)!=r['champs_sha256']:raise ValueError('champs altérés')
        fields=np.load(path,mmap_mode='r',allow_pickle=False);pieces,j=charger_cas(references/f'n{r["n"]}')
        roots=[cholesky(h,lower=True).T for d,m,h in pieces];end=np.cumsum([d.shape[1] for d,m,h in pieces])
        columns=[];operators=[]
        for x in fields:
            target=np.vstack([root@(a@x[-6:]) for root,a in zip(roots,j['applications'],strict=True)])
            error=np.vstack([root@(x[stop-6:stop]-a@x[-6:]) for root,a,stop in zip(roots,j['applications'],end,strict=True)])
            scales=np.linalg.norm(target,axis=0)
            if not np.all(np.isfinite(target)) or not np.all(np.isfinite(error)) or np.any(scales<=0):raise ArithmeticError('champ non représentable')
            _,rr=qr(target/scales,mode='economic');values=svdvals(rr)
            if values[-1]<=np.finfo(float).eps*max(target.shape)*values[0]:raise ArithmeticError('rang numérique de conformité insuffisant')
            op=svdvals(solve_triangular(rr.T,(error/scales).T,lower=True).T)[0]
            columns.append((np.linalg.norm(error,axis=0)/scales).tolist());operators.append(float(op))
        results.append(dict(dossier=r['dossier'],champs_sha256=r['champs_sha256'],colonnes=columns,operateurs=operators,
            maximum_colonnes=max(v for row in columns for v in row),maximum_operateur=max(operators)))
    ecrire(out/'audit.json',dict(portee='écart u_port_j - A_j u_global, norme directe somme H_j ; QR/SVD numériques',
        certification_machine=False,essais=results))
    ecrire(out/'manifest.json',{p.name:sha(p) for p in out.iterdir() if p.is_file()})
    print('Maximum conformité :',max(r['maximum_operateur'] for r in results))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campagne',type=Path);p.add_argument('references',type=Path);p.add_argument('sortie',type=Path)
    a=p.parse_args();auditer(a.campagne,a.references,a.sortie)
