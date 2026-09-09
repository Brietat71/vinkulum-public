"""Réduction expérimentale du complément par résolutions statiques contraintes.

Une base physique fixe conserve ports, directions retenues et Krylov du
complément. Les fréquences ne factorisent que les petits blocs projetés.
Le modèle reste D.T D et M ; les couplages sont tous conservés. La borne
spectrale est fournie, et aucune certification machine des réponses n'est
revendiquée par ce prototype.
"""
import operator
import time
import warnings

import numpy as np
from scipy.linalg import eigh, solve, lu_factor, lu_solve, LinAlgWarning
from scipy.sparse import csc_matrix

from .inverse_contrainte_energie import InverseContrainteEnergie


def _sym(a):
    return (a+a.T)*.5


def _entier(x, nom):
    if isinstance(x, (bool, np.bool_)):
        raise ValueError(nom+" entier strictement positif requis")
    x = operator.index(x)
    if x < 1:
        raise ValueError(nom+" entier strictement positif requis")
    return x


class KrylovContraint:
    """Krylov contraint autour d'un QR fourni, sans KKT fréquentiel.

    Hypothèses : K_ii et M_ii SPD, masse complète PSD, B de rang plein,
    B.T Phi inversible, et lambda_complement minorant sur ker(B.T).
    Les tests numériques ne démontrent pas seuls ces propriétés.
    """

    def __init__(self, qr, m, b, phi, lambda_complement, omega_max, *,
                 blocs=4, max_directions=128, seuil_dependance=1e-11,
                 seuil_semence=1e-12):
        debut = time.perf_counter()
        if any(np.iscomplexobj(v) for v in (m, b, phi, lambda_complement, omega_max)):
            raise ValueError("données réelles requises")
        self.qr = qr
        self.d, self.i, self.s = qr.d.copy(), qr.i.copy(), qr.s.copy()
        self.di = self.d[:, self.i].tocsc()
        self.m = csc_matrix(m, dtype=float, copy=True)
        n, ni, p = self.d.shape[1], len(self.i), len(self.s)
        if (self.m.shape != (n, n) or not np.all(np.isfinite(self.m.data))
                or np.any((self.m-self.m.T).data != 0)):
            raise ValueError("masse finie exactement symétrique requise")
        self.mi = self.m[self.i][:, self.i].tocsc()
        self.mi.sort_indices()
        self.b, self.phi = np.array(b, dtype=float, copy=True), np.array(phi, dtype=float, copy=True)
        if (self.b.ndim != 2 or self.b.shape[0] != ni or not 0 < self.b.shape[1] < ni
                or self.phi.shape != self.b.shape or not np.all(np.isfinite(self.b))
                or not np.all(np.isfinite(self.phi))):
            raise ValueError("B et Phi compatibles finis requis")
        self.r, self.p = self.b.shape[1], p
        if np.linalg.matrix_rank(self.b.T@self.phi) != self.r:
            raise ValueError("B.T Phi inversible requis")
        self.lambda_complement, self.omega_max = float(lambda_complement), float(omega_max)
        if (not np.isfinite(self.lambda_complement) or self.lambda_complement <= 0
                or not np.isfinite(self.omega_max) or self.omega_max < 0
                or self.omega_max**2 >= self.lambda_complement):
            raise ValueError("bande strictement sous le minorant complémentaire requise")
        if (not np.isfinite(seuil_dependance) or not 0 < seuil_dependance < 1
                or not np.isfinite(seuil_semence) or seuil_semence < 0):
            raise ValueError("seuils de construction invalides")
        self.seuil_dependance = seuil_dependance
        self.max_directions = min(_entier(max_directions, "max_directions"), ni-self.r)
        blocs = _entier(blocs, "blocs")
        self.inverse = InverseContrainteEnergie(qr.r, qr.echelles, self.b)
        self.c = self.inverse.canonique()
        # Les coordonnées conservées sont normalisées en énergie. Les
        # covecteurs B restent exactement ceux du certificat fourni.
        normes = np.linalg.norm(self.di@self.phi, axis=0)
        if np.any(normes <= 0) or not np.all(np.isfinite(normes)):
            raise ValueError("directions retenues d'énergie non positive")
        self.phi /= normes
        self.t = np.zeros((n, p+self.r))
        self.t[self.i, :p], self.t[self.s, :p] = qr.psi, qr.w
        self.t[self.i, p:] = self.phi
        self.dt, self.mt = self.d@self.t, self.m@self.t
        self.g0, self.g1 = self.di.T@self.dt, self.mt[self.i]
        self.v = np.empty((ni, 0))
        self.dv = np.empty((self.d.shape[0], 0))
        self.historique = []
        self.resolutions_statiques = 0
        semence = self._inverse(np.column_stack((self.g0, self.lambda_complement*self.g1)))
        # Écarter une semence minuscule ne signifie pas que son résidu
        # est nul : G0 et G1 complets restent disponibles pour l'audit.
        semence = semence[:, np.linalg.norm(self.di@semence, axis=0) > seuil_semence]
        self.dernier_bloc = np.empty((ni, 0))
        self.statut = "budget_blocs"
        for etape in range(blocs):
            nouveau = self._ajouter(semence)
            if not nouveau.shape[1]:
                self.statut = "stagnation"
                break
            if self.v.shape[1] >= self.max_directions:
                self.statut = "budget_directions"
                break
            if etape+1 == blocs:
                break
            semence = self._inverse(self.lambda_complement*(self.mi@nouveau))
        self._projeter()
        self.preparation_s = time.perf_counter()-debut

    def _inverse(self, rhs):
        self.resolutions_statiques += 1
        return self.inverse.appliquer(rhs)

    def _ajouter(self, candidat):
        if not candidat.shape[1]:
            return candidat
        w = candidat.copy()
        initiales = np.linalg.norm(self.di@w, axis=0)
        # Réparation flottante : son défaut résiduel sera rapporté ; il
        # n'est jamais identifié à une contrainte satisfaite exactement.
        for _ in range(2):
            w -= self.c@(self.b.T@w)
        if self.v.shape[1]:
            gram = _sym(self.dv.T@self.dv)
            for _ in range(2):
                w -= self.v@solve(gram, self.dv.T@(self.di@w), assume_a="pos")
        dw = self.di@w
        normes = np.linalg.norm(dw, axis=0)
        actifs = normes > self.seuil_dependance*np.maximum(initiales, np.finfo(float).tiny)
        w, dw = w[:, actifs], dw[:, actifs]
        normes = normes[actifs]
        if not len(normes):
            return w
        w, dw = w/normes, dw/normes
        valeurs, vecteurs = eigh(_sym(dw.T@dw), check_finite=False)
        indices = np.flatnonzero(valeurs > self.seuil_dependance*max(valeurs[-1], np.finfo(float).tiny))[::-1]
        indices = indices[:max(0, self.max_directions-self.v.shape[1])]
        nouveau = (w@vecteurs[:, indices])/np.sqrt(valeurs[indices])
        # La normalisation d'une direction presque dépendante peut
        # amplifier son défaut de contrainte : réparer après cette étape.
        for _ in range(2):
            nouveau -= self.c@(self.b.T@nouveau)
        self.v = np.column_stack((self.v, nouveau))
        self.dv = self.di@self.v
        self.dernier_bloc = nouveau
        self.historique.append(dict(directions=self.v.shape[1],
            defaut_contrainte=float(np.linalg.norm(self.b.T@self.v, 2))))
        return nouveau

    def enrichir(self):
        """Ajoute un bloc avec le même facteur statique, puis reprojette."""
        if self.v.shape[1] >= self.max_directions or not self.dernier_bloc.shape[1]:
            return False
        nouveau = self._ajouter(self._inverse(self.lambda_complement*(self.mi@self.dernier_bloc)))
        if not nouveau.shape[1]:
            self.statut = "stagnation"
            return False
        self._projeter()
        return True

    def _projeter(self):
        ev = np.zeros((self.d.shape[1], self.v.shape[1]))
        ev[self.i] = self.v
        self.base = np.column_stack((self.t, ev))
        self.db = self.d@self.base
        self.mb = self.m@self.base
        self.kr = _sym(self.db.T@self.db)
        self.mr = _sym(self.base.T@self.mb)
        self.taille_reduite = self.base.shape[1]
        self.taille_complement_reduit = self.v.shape[1]
        self.defaut_contrainte = float(np.linalg.norm(self.b.T@self.v, 2)) if self.v.shape[1] else 0.
        self.defaut_canonique = float(np.linalg.norm(self.b.T@self.c-np.eye(self.r), 2))

    def _pulsation(self, omega):
        if np.iscomplexobj(omega) or np.ndim(omega) != 0:
            raise ValueError("pulsation scalaire réelle requise")
        omega = float(omega)
        if not np.isfinite(omega) or not 0 <= omega <= self.omega_max:
            raise ValueError("pulsation hors bande construite")
        return omega

    def reponses(self, omega, forces):
        """Réponse groupée dans la base fixe ; pas encore de borne de champ."""
        omega = self._pulsation(omega)
        if np.iscomplexobj(forces):
            raise ValueError("forces réelles requises")
        forces = np.asarray(forces, dtype=float)
        if (forces.ndim != 2 or forces.shape[0] != self.p or not forces.shape[1]
                or not np.all(np.isfinite(forces))):
            raise ValueError("forces finies de forme (ports, charges) requises")
        ar = self.kr-omega**2*self.mr
        rhs = np.zeros((self.taille_reduite, forces.shape[1]))
        rhs[:self.p] = self.qr.w.T@forces
        with warnings.catch_warnings():
            warnings.simplefilter("error", LinAlgWarning)
            try:
                lu = lu_factor(ar, check_finite=False)
            except LinAlgWarning as exc:
                raise np.linalg.LinAlgError("système réduit singulier") from exc
        q = lu_solve(lu, rhs, check_finite=False)
        champ = self.base@q
        if not np.all(np.isfinite(champ)):
            raise ArithmeticError("champ réduit non fini")
        return dict(champ=champ, coordonnees=q, residu_reduit=rhs-ar@q,
                    certification_machine=False, taille_reduite=self.taille_reduite,
                    taille_complement_reduit=self.taille_complement_reduit,
                    defaut_contrainte=self.defaut_contrainte)
