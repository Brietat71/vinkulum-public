"""Témoin de condensation à coordonnées intérieures retenues, sans base Z.

Ce programme de recherche refactorise un KKT creux de taille n+r à chaque
fréquence. Il contrôle la composition mathématique, sans gain de temps
revendiqué et sans modifier le noyau. Modèle réel conservatif K=D.T D ;
forces de port seulement. La coercivité sur ker(B.T) est une hypothèse
fournie, pas un certificat établi par ce témoin.

B désigne les valeurs binary64 stockées, même si son origine est M Phi.
Tous les couplages du modèle D sont conservés. Les calculs, résidus,
valeurs singulières et contraintes restent évalués sans encadrement.
"""
import warnings

import numpy as np
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from scipy.sparse import bmat, csc_matrix, diags
from scipy.sparse.linalg import splu


def _reel(a, nom):
    if np.iscomplexobj(a):
        raise ValueError(nom + " réel requis")
    a = np.asarray(a, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(nom + " fini requis")
    return a.copy()


def _indices(a, n, nom):
    a = np.asarray(a)
    if (a.ndim != 1 or a.dtype.kind not in "iu" or not len(a)
            or len(np.unique(a)) != len(a) or np.any(a < 0) or np.any(a >= n)):
        raise ValueError(nom + " : indices entiers distincts requis")
    return a.astype(np.intp)


class TemoinRetention:
    """Ports et r coordonnées intérieures conservées ; complément éliminé.

    Hypothèses : D_i.T D_i et M_ii SPD, masse complète PSD, B de rang r,
    B.T Phi inversible, et coercivité fournie sur le noyau du B stocké.
    Les tests numériques de rang ne certifient aucune de ces hypothèses.

    Le KKT utilise le produit creux D_i.T D_i calculé. Deux corrections
    résiduelles évaluées via D_i rapprochent la résolution du modèle
    énergétique d'entrée ; ce n'est pas une preuve d'arrondi.
    """
    def __init__(self, d, m, interieur, interface, b, phi, lambda_complement,
                 omega_max, *, psi=None, normalisation=None):
        if np.iscomplexobj(d) or np.iscomplexobj(m):
            raise ValueError("D et M réels requis")
        self.d, self.m = csc_matrix(d, dtype=float, copy=True), csc_matrix(m, dtype=float, copy=True)
        n = self.d.shape[1]
        if (self.m.shape != (n, n) or not np.all(np.isfinite(self.d.data))
                or not np.all(np.isfinite(self.m.data))):
            raise ValueError("D et M : dimensions et valeurs finies requises")
        if np.any((self.m-self.m.T).data != 0):
            raise ValueError("M exactement symétrique dans les données stockées requis")
        self.i = _indices(interieur, n, "intérieur")
        self.s = _indices(interface, n, "interface")
        if not np.array_equal(np.sort(np.r_[self.i, self.s]), np.arange(n)):
            raise ValueError("partition complète disjointe requise")
        self.b, self.phi = _reel(b, "B"), _reel(phi, "Phi")
        ni, self.p = len(self.i), len(self.s)
        if (self.b.ndim != 2 or self.b.shape[0] != ni
                or not 0 < self.b.shape[1] < ni or self.phi.shape != self.b.shape):
            raise ValueError("B et Phi de forme (intérieur, retenues), 0 < retenues < intérieur")
        self.r = self.b.shape[1]
        if (np.linalg.matrix_rank(self.b) != self.r
                or np.linalg.matrix_rank(self.b.T@self.phi) != self.r):
            raise ValueError("B de rang plein et B.T Phi inversible requis")
        if (np.iscomplexobj(lambda_complement) or np.ndim(lambda_complement) != 0
                or np.iscomplexobj(omega_max) or np.ndim(omega_max) != 0):
            raise ValueError("bande et minorant scalaires réels requis")
        self.lambda_complement = float(lambda_complement)
        self.omega_max = float(omega_max)
        if (not np.isfinite(self.lambda_complement) or self.lambda_complement <= 0
                or not np.isfinite(self.omega_max) or self.omega_max < 0
                or self.omega_max**2 >= self.lambda_complement):
            raise ValueError("omega_max² < lambda_complement positif requis")
        self.psi = np.zeros((ni, self.p)) if psi is None else _reel(psi, "Psi")
        self.w = np.eye(self.p) if normalisation is None else _reel(normalisation, "W")
        if (self.psi.shape != (ni, self.p) or self.w.shape != (self.p, self.p)
                or np.linalg.matrix_rank(self.w) != self.p):
            raise ValueError("relèvement compatible et W inversible requis")
        self.t = np.zeros((n, self.p+self.r))
        self.t[self.i, :self.p], self.t[self.s, :self.p] = self.psi, self.w
        self.t[self.i, self.p:] = self.phi
        self.di = self.d[:, self.i]
        self.ki = (self.di.T@self.di).tocsc()
        self.mi = self.m[self.i][:, self.i].tocsc()
        self.dt, self.mt = self.d@self.t, self.m@self.t
        if np.any(self.ki.diagonal() <= 0):
            raise ValueError("diagonale intérieure positive requise")
        equilibre = 1/np.sqrt(self.ki.diagonal())
        poids_b = np.linalg.norm(equilibre[:, None]*self.b, axis=0)
        self.equilibre_kkt = np.r_[equilibre, 1/poids_b]
        if not np.all(np.isfinite(self.equilibre_kkt)):
            raise ArithmeticError("équilibration KKT non finie")

    def _pulsation(self, omega):
        if np.iscomplexobj(omega) or np.ndim(omega) != 0:
            raise ValueError("pulsation scalaire réelle requise")
        omega = float(omega)
        if not np.isfinite(omega) or not 0 <= omega <= self.omega_max:
            raise ValueError("pulsation hors bande du complément")
        return omega

    def condensation(self, omega):
        """Reconstruit X et le Schur conservé ; aucune factorisation globale dense."""
        omega = self._pulsation(omega)
        z = omega**2
        ni = len(self.i)
        g = self.di.T@self.dt-z*self.mt[self.i]
        kkt = bmat([[self.ki-z*self.mi, csc_matrix(self.b)],
                    [csc_matrix(self.b.T), None]], format="csc")
        equilibrage = diags(self.equilibre_kkt, format="csc")
        try:
            facteur = splu((equilibrage@kkt@equilibrage).tocsc())
        except RuntimeError as exc:
            raise np.linalg.LinAlgError("KKT du complément singulier") from exc

        def resoudre(rhs):
            return self.equilibre_kkt[:, None]*facteur.solve(self.equilibre_kkt[:, None]*rhs)

        rhs = np.vstack((g, np.zeros((self.r, self.p+self.r))))
        solution = resoudre(rhs)

        def residu_kkt(sol):
            v, mu = sol[:ni], sol[ni:]
            return np.vstack((g-self.di.T@(self.di@v)+z*(self.mi@v)-self.b@mu,
                              -self.b.T@v))

        for _ in range(2):
            solution += resoudre(residu_kkt(solution))
        if not np.all(np.isfinite(solution)):
            raise ArithmeticError("résolution KKT non finie")
        correction = solution[:ni]
        x = self.t.copy()
        x[self.i] -= correction
        dx, mx = self.d@x, self.m@x
        energetique = dx.T@dx-z*(x.T@mx)
        schur = (energetique+energetique.T)*0.5
        algebrique = self.dt.T@self.dt-z*(self.t.T@self.mt)-g.T@correction
        residu = residu_kkt(solution)
        return dict(schur=schur, reconstruction=x,
                    sigma_min_bloc_conserve=float(np.linalg.svd(schur, compute_uv=False)[-1]),
                    defaut_contrainte=float(np.linalg.norm(self.b.T@(x[self.i]-self.t[self.i]), 2)),
                    residu_kkt=float(np.linalg.norm(residu, 2)),
                    ecart_fonctionnel=float(np.linalg.norm(schur-algebrique, 2)),
                    taille_conservee=self.p+self.r,
                    taille_kkt=ni+self.r,
                    lambda_complement_hypothese=self.lambda_complement,
                    convention_b="valeurs binary64 stockées",
                    certification_machine=False)

    def reponses(self, omega, forces_ports):
        """Colonnes de forces physiques ; zéro force directe des coordonnées retenues."""
        omega = self._pulsation(omega)
        forces = _reel(forces_ports, "forces")
        if forces.ndim != 2 or forces.shape[0] != self.p or not forces.shape[1]:
            raise ValueError("forces de forme (ports, charges) requises")
        etat = self.condensation(omega)
        rhs = np.vstack((self.w.T@forces, np.zeros((self.r, forces.shape[1]))))
        with warnings.catch_warnings():
            warnings.simplefilter("error", LinAlgWarning)
            try:
                facteur = lu_factor(etat["schur"], check_finite=False)
            except LinAlgWarning as exc:
                raise np.linalg.LinAlgError("bloc conservé singulier : réponse non unique") from exc
        q = lu_solve(facteur, rhs, check_finite=False)
        champ = etat["reconstruction"]@q
        if not np.all(np.isfinite(champ)):
            raise ArithmeticError("réponse non finie")
        f = np.zeros_like(champ)
        f[self.s] = forces
        residu = f-self.d.T@(self.d@champ)+omega**2*(self.m@champ)
        return dict(etat, coordonnees_conservees=q, champ=champ,
                    residu_physique=residu,
                    residu_bloc_conserve=rhs-etat["schur"]@q)
