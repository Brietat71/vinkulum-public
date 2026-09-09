"""Reproduit deux sous-majorations sur le snapshot local, jamais le module courant."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS'):
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.sched_setaffinity(0,{0})
from pathlib import Path
import hashlib
import json
import sys
from fractions import Fraction as F

import numpy as np
import scipy

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/tmp/vinkulum-energie-native-travail/ci')
sys.path.insert(0,str(ROOT))
from controle_facteurs import ControleFacteurs, ImageFaibleRang

image=np.full((128,2),2.**-27);image[0]=1.
gauche=np.zeros((128,1));droite=np.zeros((1,2));v=np.array([[1.],[0.]])
dx=np.array([[2.**-600,0.],[0.,2.**-600],[0.,0.]])
np.savez_compressed(ROOT/'donnees.npz',image=image,gauche=gauche,droite=droite,v=v,dx=dx)

h=ImageFaibleRang(image,gauche,droite)
borne=float(h.normes(v)[0]);exact2=F(1)+127*F(2)**-54
ecart=exact2-F(borne)**2
image_result=dict(borne=borne,borne_hex=borne.hex(),
    norme_exacte_carre=dict(numerateur=str(exact2.numerator),denominateur=str(exact2.denominator)),
    ecart_exact_carre=dict(numerateur=str(ecart.numerator),denominateur=str(ecart.denominator)),
    sous_majore=ecart>0,diagnostic=h.diagnostic,
    conditions='Entrées et produits normaux, deux colonnes,128 lignes ; aucun sous-flux requis.')

# Sonde isolée de la routine de norme : ses deux dimensions sont celles
# d'un contrôleur à3coordonnées physiques et2coordonnées conservées.
c=object.__new__(ControleFacteurs);c.extension='gram';u=np.finfo(float).eps/2
c._gram_gamma=3*u/(1-3*u);c._somme_gamma=2*u/(1-2*u)
g=dx.T@dx
borne_gram=float(c._norme_extension_d(None,dx,g))
gram_result=dict(borne=borne_gram,norme_exacte=2.**-600,norme_exacte_hex=(2.**-600).hex(),
    gram=g.tolist(),sous_majore=F(borne_gram)<F(2.**-600),
    conditions='Routine de norme isolée ; Gram sous-flué, pas une réponse mécanique complète.')
result=dict(statuts='Deux réfutations reproductibles, aucun changement du module courant',
    python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
    snapshot_sha256=hashlib.sha256((ROOT/'controle_facteurs.py').read_bytes()).hexdigest(),
    cpu=sorted(os.sched_getaffinity(0)),threads={k:os.environ[k] for k in
        ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS')},
    image=image_result,gram=gram_result)
(ROOT/'resultats.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
assert image_result['sous_majore'] and gram_result['sous_majore']
print(json.dumps(result,ensure_ascii=False),flush=True)
