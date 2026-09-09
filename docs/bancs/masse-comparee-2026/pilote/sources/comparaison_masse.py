"""Comparaison de Loewner vérifiée : alpha diag(ell) < M < beta diag(ell).

Les opérateurs d'équilibre gardent M complet ; la diagonale sert uniquement
à transporter des majorants de normes primales et duales. Hors API publique.
"""
import ctypes as ct
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from inertie_binaire_compile import Pivot,Summary,_preuve
from inertie_complement_dirigee import InertieImpossible,_empreinte
from matrices_certificat import matrice_entree
from vinkulum._ports.inverse_selectionnee import _entier_positif


def compiler_comparaison(destination):
    destination=Path(destination).resolve();destination.mkdir(parents=True,exist_ok=False)
    captures={}
    for name in ('inertie_binaire_native.cpp','comparaison_masse_native.cpp'):
        data=Path(__file__).with_name(name).read_bytes();(destination/name).write_bytes(data)
        captures[name]=hashlib.sha256(data).hexdigest()
    lib=destination/'comparaison_masse.so'
    command=['c++','-std=c++17','-O3','-fPIC','-shared','-fno-fast-math','-ffp-contract=off',
        '-frounding-math','-Wall','-Wextra','-Werror',str(destination/'comparaison_masse_native.cpp'),'-o',str(lib)]
    result=subprocess.run(command,capture_output=True,text=True)
    (destination/'compilation.log').write_text(result.stdout+result.stderr);result.check_returncode()
    compiler=subprocess.run(['c++','--version'],check=True,capture_output=True,text=True).stdout
    identity=dict(sources_sha256=captures,commande=command,compilateur=compiler,
                  bibliotheque_sha256=hashlib.sha256(lib.read_bytes()).hexdigest())
    (destination/'identite.json').write_text(json.dumps(identity,ensure_ascii=False,indent=2)+'\n')
    return lib


class BibliothequeComparaison:
    def __init__(self,path):
        path=Path(path).resolve();self.identite=dict(chemin=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        self.lib=ct.CDLL(str(path));self.appeler=self.lib.vinkulum_compare_mass
        pi,pd=ct.POINTER(ct.c_int64),ct.POINTER(ct.c_double)
        self.appeler.argtypes=[ct.c_int64,pi,pi,pd,pd,ct.c_double,ct.c_double,pi,ct.c_int64,ct.c_int64,
                              ct.POINTER(Pivot),ct.POINTER(Pivot),ct.POINTER(Summary),pi,pd]
        self.appeler.restype=ct.c_int


def certifier_comparaison_masse(m,alpha,beta,*,bibliotheque,diagonale=None,permutation=None,
                               budget_operations=100_000_000,budget_coefficients=2_000_000):
    start=time.perf_counter()
    budgets=[_entier_positif(x,k) for x,k in ((budget_operations,'budget_operations'),(budget_coefficients,'budget_coefficients'))]
    if max(budgets)>np.iinfo(np.int64).max:raise ValueError('budgets int64 requis')
    if any(np.iscomplexobj(v) or np.ndim(v)!=0 for v in (alpha,beta)):
        raise ValueError('alpha et beta réels scalaires requis')
    alpha,beta=float(alpha),float(beta)
    if not 0<alpha<beta or not math.isfinite(beta):raise ValueError('0 < alpha < beta finis requis')
    m=matrice_entree(m,'M');n=m.shape[0]
    if n==0 or m.shape!=(n,n) or np.any((m-m.T).data!=0):raise ValueError('M carrée exactement symétrique requise')
    if n>budgets[1]:raise InertieImpossible('budget de dimension dépassé',dict(dimension=n,budget=budgets[1]))
    if diagonale is not None and np.iscomplexobj(diagonale):raise ValueError('diagonale réelle requise')
    ell=np.array(m.diagonal() if diagonale is None else diagonale,dtype=np.float64,order='C',copy=True)
    if ell.shape!=(n,) or not np.all(np.isfinite(ell)) or np.any(ell<=0):raise ValueError('diagonale finie strictement positive requise')
    if permutation is None:perm=np.arange(n,dtype=np.int64)
    else:
        perm=np.asarray(permutation)
        if perm.dtype.kind not in 'iu' or perm.shape!=(n,) or not np.array_equal(np.sort(perm),np.arange(n)):
            raise ValueError('permutation entière complète requise')
        perm=np.ascontiguousarray(perm,dtype=np.int64)
    if not isinstance(bibliotheque,BibliothequeComparaison):raise TypeError('bibliothèque explicitement chargée requise')
    arrays=[np.ascontiguousarray(x,dtype=t) for x,t in ((m.indptr,np.int64),(m.indices,np.int64),(m.data,np.float64))]
    def ptr(a):return a.ctypes.data_as(ct.POINTER(ct.c_int64 if a.dtype==np.int64 else ct.c_double))
    lower,upper=(Pivot*n)(),(Pivot*n)();summaries=(Summary*2)();diagnostic=(ct.c_int64*3)();interval=(ct.c_double*2)()
    status=bibliotheque.appeler(n,*[ptr(a) for a in arrays],ptr(ell),alpha,beta,ptr(perm),*budgets,lower,upper,summaries,diagnostic,interval)
    if status:
        raise InertieImpossible('comparaison massique non démontrée',dict(code=status,
            phase='borne_inferieure' if diagnostic[1]==0 else 'borne_superieure',
            pivot=diagnostic[2],intervalle=[v.hex() for v in interval],operations_binary64=diagnostic[0]))
    result=dict(alpha=alpha,beta=beta,dimension=n,diagonale=ell.tolist(),
        masse_sha256=_empreinte(m),diagonale_sha256=hashlib.sha256(np.asarray(ell,dtype='<f8').tobytes()).hexdigest(),
        certification_machine=True,encodage_intervalles='hex_binary64',
        preuve_inferieure=_preuve(summaries[0],lower),preuve_superieure=_preuve(summaries[1],upper),
        permutation_physique=perm.tolist(),operations_binary64=diagnostic[0],bibliotheque=bibliotheque.identite,
        portee='alpha diag(ell) < M < beta diag(ell), coefficients stockés exacts ; aucune réponse certifiée')
    result['preparation_s']=time.perf_counter()-start
    return result
