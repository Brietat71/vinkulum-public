"""Pont privé vers le prototype C++ d'inertie, compilation explicite séparée."""
import ctypes as ct
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from inertie_complement_dirigee import InertieImpossible, _empreinte, _matrice_entree
from vinkulum._ports.inverse_selectionnee import _entier_positif


class Pivot(ct.Structure):
    _fields_ = [(x, ct.c_int64) for x in ('i', 'j', 'positive', 'negative')] + [
        (x, ct.c_double) for x in ('alo','ahi','blo','bhi','clo','chi','dlo','dhi')]


class Summary(ct.Structure):
    _fields_ = [(x, ct.c_int64) for x in
               ('positive','negative','zero','count','width','peak','created','initial')]


def compiler(destination):
    """Un dossier neuf conserve source, commande, journal et identité du binaire."""
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).with_name('inertie_binaire_native.cpp').read_bytes()
    cpp, lib = destination/'inertie_binaire_native.cpp', destination/'inertie_binaire.so'
    cpp.write_bytes(source)
    command = ['c++','-std=c++17','-O3','-fPIC','-shared','-fno-fast-math',
               '-ffp-contract=off','-frounding-math','-Wall','-Wextra','-Werror',str(cpp),'-o',str(lib)]
    result = subprocess.run(command, capture_output=True, text=True)
    (destination/'compilation.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    version = subprocess.run(['c++','--version'], check=True, capture_output=True, text=True).stdout
    identity = dict(source_sha256=hashlib.sha256(source).hexdigest(),
                    bibliotheque_sha256=hashlib.sha256(lib.read_bytes()).hexdigest(),
                    commande=command, compilateur=version)
    (destination/'identite.json').write_text(json.dumps(identity, ensure_ascii=False, indent=2))
    return lib


class BibliothequeInertie:
    def __init__(self, chemin):
        chemin = Path(chemin).resolve()
        self.identite = dict(chemin=str(chemin), sha256=hashlib.sha256(chemin.read_bytes()).hexdigest())
        self.lib = ct.CDLL(str(chemin))
        self.appeler = self.lib.vinkulum_inertia_binary64
        pi, pd = ct.POINTER(ct.c_int64), ct.POINTER(ct.c_double)
        self.appeler.argtypes = [ct.c_int64]*3 + [pi,pi,pd,pi,pi,pd,pd,ct.c_double,pi,
            ct.c_int64,ct.c_int64,ct.c_bool,ct.POINTER(Pivot),ct.POINTER(Pivot),ct.POINTER(Summary),pi,pd,pd]
        self.appeler.restype = ct.c_int


def _preuve(summary, pivots):
    resultat = dict(signature=[summary.positive, summary.negative, summary.zero],
        largeur_max=summary.width, pic_coefficients=summary.peak,
        coefficients_crees=summary.created, pivots=[])
    for p in pivots[:summary.count]:
        ligne = dict(indices=[p.i] if p.j < 0 else [p.i,p.j], signature=[p.positive,p.negative])
        if p.j < 0:
            ligne['diagonal'] = [p.alo.hex(),p.ahi.hex()]
        else:
            for name,lo,hi in [('a',p.alo,p.ahi),('b',p.blo,p.bhi),('c',p.clo,p.chi),('determinant',p.dlo,p.dhi)]:
                ligne[name] = [lo.hex(),hi.hex()]
        resultat['pivots'].append(ligne)
    return resultat


def certifier_inertie_compilee(d_i, mii, b, gamma, *, bibliotheque,
        budget_operations=100_000_000, budget_coefficients=2_000_000,
        budget_rectangulaire=2_000_000, permutation=None):
    """Même KKT exact que le témoin Decimal ; aucun repli ou seuil modifié."""
    debut = time.perf_counter()
    budgets = [_entier_positif(v,k) for v,k in [(budget_operations,'budget_operations'),
        (budget_coefficients,'budget_coefficients'),(budget_rectangulaire,'budget_rectangulaire')]]
    budget_operations,budget_coefficients,budget_rectangulaire = budgets
    if max(budgets) > np.iinfo(np.int64).max:
        raise ValueError('budgets représentables en int64 requis')
    if np.iscomplexobj(gamma) or np.ndim(gamma) != 0:
        raise ValueError('gamma scalaire réel requis')
    gamma = float(gamma)
    if not math.isfinite(gamma) or gamma <= 0:
        raise ValueError('gamma fini strictement positif requis')
    d,m = _matrice_entree(d_i,'D_i'),_matrice_entree(mii,'M_ii')
    n = d.shape[1]
    if m.shape != (n,n) or n == 0 or np.any((m-m.T).data != 0):
        raise ValueError('M_ii exactement symétrique et compatible requise')
    if np.iscomplexobj(b):
        raise ValueError('B réel requis')
    b = np.asarray(b)
    if b.ndim != 2 or b.shape[0] != n or not 0 < b.shape[1] < n or not np.all(np.isfinite(b)):
        raise ValueError('B fini de forme (n,s), 0<s<n requis')
    s = b.shape[1]
    if b.size > budget_rectangulaire or n+s > budget_coefficients:
        raise InertieImpossible('budget de dimension ou de contraintes dépassé',
            dict(dimension=n+s, coefficients_contraintes=b.size,
                 budget_coefficients=budget_coefficients, budget_rectangulaire=budget_rectangulaire))
    b = np.array(b,dtype=np.float64,order='C',copy=True)
    if not np.all(np.isfinite(b)):
        raise ValueError('B converti fini requis')
    if permutation is None:
        perm = np.arange(n,dtype=np.int64)
    else:
        perm = np.asarray(permutation)
        if perm.dtype.kind not in 'iu' or perm.shape != (n,) or not np.array_equal(np.sort(perm),np.arange(n)):
            raise ValueError('permutation physique entière complète requise')
        perm = np.ascontiguousarray(perm,dtype=np.int64)
    diagonal = (m.nnz == n and np.array_equal(m.indptr,np.arange(n+1)) and np.array_equal(m.indices,np.arange(n)))
    if diagonal and np.any(m.data <= 0):
        raise InertieImpossible('masse diagonale non définie positive',dict(phase='masse'))
    if not isinstance(bibliotheque,BibliothequeInertie):
        raise TypeError('bibliothèque explicitement chargée requise')
    arrays = [np.ascontiguousarray(x,dtype=dtype) for x,dtype in
        [(d.indptr,np.int64),(d.indices,np.int64),(d.data,np.float64),
         (m.indptr,np.int64),(m.indices,np.int64),(m.data,np.float64)]]
    def ptr(x):
        return x.ctypes.data_as(ct.POINTER(ct.c_int64 if x.dtype == np.int64 else ct.c_double))
    mass_pivots,kkt_pivots = (Pivot*n)(),(Pivot*(n+s))()
    summaries,diagnostic,interval,phases = (Summary*2)(),(ct.c_int64*3)(),(ct.c_double*2)(),(ct.c_double*3)()
    status = bibliotheque.appeler(n,s,d.shape[0],*[ptr(x) for x in arrays],ptr(b),gamma,ptr(perm),
        budget_operations,budget_coefficients,diagonal,mass_pivots,kkt_pivots,summaries,diagnostic,interval,phases)
    if status:
        messages = {1:'budget de coefficients dépassé',2:"budget d'opérations binary64 dépassé",
            3:'encadrement binary64 fini indisponible',4:'diviseur non séparé de zéro',
            5:'masse non définie positive',6:'pivot non séparé de zéro',
            7:'complément non coercif au seuil demandé',8:'environnement IEEE au plus proche et sous-flux graduel requis',
            9:'échec interne ou allocation indisponible'}
        raise InertieImpossible(messages.get(status,'échec natif inconnu'),
            dict(code=status,phase=('masse','assemblage','elimination')[diagnostic[1]],
                operations_binary64=diagnostic[0],pivot=diagnostic[2],
                intervalle=[x.hex() for x in interval],encodage_intervalles='hex_binary64',
                preuve_masse=_preuve(summaries[0],mass_pivots),preuve_kkt=_preuve(summaries[1],kkt_pivots)))
    mass = (dict(methode='diagonale_binary64_positive',signature=[n,0,0],minimum_binary64=float(m.data.min()))
            if diagonal else dict(_preuve(summaries[0],mass_pivots),methode='congruences_binary64'))
    result = dict(lambda_min=gamma,certification_machine=True,arithmetique='binary64_nextafter_cpp',
        encodage_intervalles='hex_binary64',dimension_interieure=n,rang_contraintes=s,dimension_complement=n-s,
        contraintes_sha256=hashlib.sha256(np.asarray(b,dtype='<f8').tobytes()).hexdigest(),
        d_sha256=_empreinte(d),masse_sha256=_empreinte(m),permutation_physique=perm.tolist(),
        preuve_masse=mass,preuve_kkt=_preuve(summaries[1],kkt_pivots),
        coefficients_initiaux=summaries[1].initial,operations_binary64=diagnostic[0],coefficients_contraintes=b.size,
        bibliotheque=bibliotheque.identite,
        phases_natives_s=dict(zip(('masse','assemblage','elimination'),phases)),
        portee='coercivité stricte de D_i.T D_i-gamma M_ii sur ker(B.T) ; aucune réponse certifiée')
    result['preparation_s'] = time.perf_counter()-debut
    return result
