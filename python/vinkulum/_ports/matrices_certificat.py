"""Validation CSR par frontières de lignes, sans drapeau canonique en cache."""
import numpy as np
from scipy.sparse import csr_matrix

from .validation_certificat import _matrice_entree


def matrice_entree(x, nom):
    if getattr(x,'format',None) != 'csr':
        return _matrice_entree(x,nom)
    if np.iscomplexobj(x):
        raise ValueError(nom+' réel requis')
    ptr,ids,data=x.indptr,x.indices,x.data
    rows,columns=x.shape
    if (ptr.ndim!=1 or ids.ndim!=1 or data.ndim!=1
            or ptr.dtype.kind not in 'iu' or ids.dtype.kind not in 'iu'
            or len(ptr)!=rows+1 or ptr[0]!=0 or ptr[-1]!=len(ids)
            or len(data)!=len(ids) or np.any(ptr[1:]<ptr[:-1])
            or np.any(ids<0) or np.any(ids>=columns)):
        raise ValueError(nom+' : stockage CSR invalide')
    # Chaque paire adjacente d'indices appartient à la même ligne, sauf
    # aux frontières ptr[i]-1. Les lignes vides peuvent répéter ces frontières.
    descentes=ids[1:]<=ids[:-1]
    frontieres=ptr[1:-1]
    frontieres=frontieres[(frontieres>0)&(frontieres<len(ids))]-1
    descentes[frontieres]=False
    if np.any(descentes):
        raise ValueError(nom+' : doublons ou indices non triés')
    result=csr_matrix(x,dtype=np.float64,copy=True)
    if not np.all(np.isfinite(result.data)):
        raise ValueError(nom+' converti fini requis')
    result.eliminate_zeros()
    return result
