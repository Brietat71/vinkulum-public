"""Validation et identité des coefficients exacts visés par les certificats."""
import hashlib
import numpy as np
from scipy.sparse import csr_matrix
from .inverse_selectionnee import _verifier_canonique

class InertieImpossible(RuntimeError):
    def __init__(self, message, bilan):
        super().__init__(message)
        self.bilan = dict(bilan)
        self.diagnostic = self.bilan

def _matrice_entree(x, nom):
    _verifier_canonique(x, nom)
    a = csr_matrix(x, dtype=float, copy=True)
    _verifier_canonique(a, nom+" converti")
    if not np.all(np.isfinite(a.data)):
        raise ValueError(nom+" fini requis")
    a.eliminate_zeros()
    return a

def _empreinte(a):
    h = hashlib.sha256()
    h.update(np.asarray(a.shape, dtype="<i8").tobytes())
    for x, dtype in ((a.indptr, "<i8"), (a.indices, "<i8"), (a.data, "<f8")):
        h.update(np.asarray(x, dtype=dtype).tobytes())
    return h.hexdigest()

