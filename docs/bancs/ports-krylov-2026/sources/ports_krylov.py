"""Prototype de réduction par ports, avec enveloppe résiduelle sur une bande.

Hypothèse mathématique fournie par l'appelant : K_II >= lambda_min M_II > 0.
Les contrôles numériques détectent certaines contradictions, sans certifier
la coercivité ni les arrondis. ``borne_uniforme`` est une évaluation flottante
d'un majorant démontré en arithmétique exacte, pas un certificat machine.
Le prototype reste séparé du noyau généraliste (linéaire, conservatif, M_IS=0).
"""
import operator
import time

import numpy as np
from scipy.linalg import cholesky, eigh, solve, solve_triangular
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import splu


def _sym(a):
    return (a+a.T)*0.5


def _entier(x, nom):
    if isinstance(x, (bool, np.bool_)):
        raise ValueError(nom+" doit être un entier positif")
    x = operator.index(x)
    if x < 1:
        raise ValueError(nom+" doit être un entier positif")
    return x


def _matrice(a, nom):
    a = csc_matrix(a, dtype=float)
    if a.shape[0] != a.shape[1] or not np.all(np.isfinite(a.data)):
        raise ValueError(nom+" doit être carrée et finie")
    diff = a-a.T
    echelle = max(np.max(np.abs(a.data), initial=0.), np.finfo(float).tiny)
    if np.max(np.abs(diff.data), initial=0.) > 64*np.finfo(float).eps*echelle:
        raise ValueError(nom+" doit être symétrique")
    return (a+a.T)*0.5


def _norme_gram(gram):
    """Norme induite par un Gram ; refuse une perte de positivité substantielle."""
    if not gram.size:
        return 0.
    valeurs = eigh(_sym(gram), eigvals_only=True, check_finite=False)
    maximum = max(float(valeurs[-1]), 0.)
    if valeurs[0] < -1e-8*max(maximum, np.finfo(float).tiny):
        raise ArithmeticError("Gram résiduel numériquement indéfini")
    return np.sqrt(maximum)


class InterieurKrylov:
    """Construit un espace n×r par résolutions statiques et orthogonalisation K.

    ``b`` porte des forces normalisées dans la métrique des ports. Aucun
    spectre intérieur complet ni matrice n×n dense n'est construit. La LU
    creuse peut toutefois subir du remplissage ; son stockage est rapporté.
    Une déflation numérique ne vaut jamais acceptation de la tolérance.
    """

    def __init__(self, k, m, b, lambda_min, omega_max, tolerance=1e-6,
                 max_blocs=12, max_directions=256, seuil_dependance=1e-11):
        debut = time.perf_counter()
        self.k, self.m = _matrice(k, "K"), _matrice(m, "M")
        self.b = np.asarray(b, dtype=float)
        if (self.k.shape != self.m.shape or self.k.shape[0] == 0
                or self.b.ndim != 2 or self.b.shape[0] != self.k.shape[0]
                or self.b.shape[1] == 0 or not np.all(np.isfinite(self.b))):
            raise ValueError("dimensions K, M, B incompatibles")
        self.lambda_min, self.omega_max = float(lambda_min), float(omega_max)
        self.tolerance = float(tolerance)
        valeurs = [self.lambda_min, self.omega_max, self.tolerance, seuil_dependance]
        if (not np.all(np.isfinite(valeurs)) or self.lambda_min <= 0
                or self.omega_max < 0 or self.tolerance <= 0
                or not 0 < seuil_dependance < 1):
            raise ValueError("paramètres de bande ou tolérance invalides")
        self.rho = self.omega_max**2/self.lambda_min
        if not self.rho < 1:
            raise ValueError("bande hors coercivité établie : omega_max² < lambda_min requis")
        max_blocs = _entier(max_blocs, "max_blocs")
        max_directions = _entier(max_directions, "max_directions")
        if np.any(self.k.diagonal() <= 0) or np.any(self.m.diagonal() <= 0):
            raise ValueError("diagonale non positive : hypothèse SPD contredite")
        self._equilibre = 1/np.sqrt(self.k.diagonal())
        e = diags(self._equilibre, format="csc")
        self._lu = splu((e@self.k@e).tocsc(), permc_spec="MMD_AT_PLUS_A")
        self.nnz_facteurs = self._lu.L.nnz+self._lu.U.nnz
        self.resolutions = 0
        self.seconds_membres = 0
        self.residu_resolution_max = 0.
        self._knorme = float(np.max(np.asarray(abs(self.k).sum(axis=1))))
        self.v = np.zeros((self.k.shape[0], 0))
        self.historique = []
        self.statut = "budget_blocs"
        self.borne_uniforme = np.inf
        candidat = self._resoudre(self.b)
        for q in range(1, max_blocs+1):
            nouveau = self._orthogonaliser(candidat, seuil_dependance)
            reste = max_directions-self.v.shape[1]
            nouveau = nouveau[:, :max(0, reste)]
            if nouveau.shape[1] == 0:
                if not self.v.shape[1]:
                    # B nul : espace vide exact, y compris plusieurs ports nuls.
                    self._preparer(q)
                self.statut = ("tolerance_estimee" if self.borne_uniforme <= self.tolerance
                               else "stagnation")
                break
            self.v = np.column_stack((self.v, nouveau))
            self._preparer(q)
            self.historique.append({"blocs": q, "directions": self.v.shape[1],
                                    "borne_uniforme": self.borne_uniforme,
                                    "defaut_gram": self.defaut_gram})
            if self.borne_uniforme <= self.tolerance:
                self.statut = "tolerance_estimee"
                break
            if self.v.shape[1] >= max_directions:
                self.statut = "budget_directions"
                break
            candidat = self._resoudre(self.lambda_min*(self.m@nouveau))
        self.taille_interieure = self.v.shape[1]
        self.preparation_s = time.perf_counter()-debut

    def _resoudre(self, rhs):
        self.resolutions += 1
        self.seconds_membres += rhs.shape[1]
        x = self._equilibre[:, None]*self._lu.solve(self._equilibre[:, None]*rhs)
        # Deux corrections sur la matrice d'origine, jamais sur une matrice dense.
        for _ in range(2):
            r = rhs-self.k@x
            denom = self._knorme*np.max(np.abs(x), axis=0)+np.max(np.abs(rhs), axis=0)
            rel = np.max(np.abs(r), axis=0)/np.maximum(denom, np.finfo(float).tiny)
            if np.max(rel, initial=0.) < 8*np.finfo(float).eps:
                break
            x += self._equilibre[:, None]*self._lu.solve(self._equilibre[:, None]*r)
        r = rhs-self.k@x
        denom = self._knorme*np.max(np.abs(x), axis=0)+np.max(np.abs(rhs), axis=0)
        rel = np.max(np.abs(r), axis=0)/np.maximum(denom, np.finfo(float).tiny)
        self.residu_resolution_max = max(self.residu_resolution_max, np.max(rel, initial=0.))
        return x

    def _orthogonaliser(self, candidat, seuil):
        w = candidat.copy()
        if self.v.shape[1]:
            g = _sym(self.v.T@(self.k@self.v))
            for _ in range(2):
                w -= self.v@solve(g, self.v.T@(self.k@w), assume_a="pos")
        gram = _sym(w.T@(self.k@w))
        normes = np.sqrt(np.maximum(np.diag(gram), 0.))
        initiales = np.sqrt(np.maximum(np.sum(candidat*(self.k@candidat), axis=0), 0.))
        actifs = normes > seuil*np.maximum(initiales, np.finfo(float).tiny)
        if not np.any(actifs):
            return w[:, :0]
        w = w[:, actifs]/normes[actifs]
        valeurs, vecteurs = eigh(_sym(w.T@(self.k@w)))
        garder = valeurs > seuil*max(float(valeurs[-1]), np.finfo(float).tiny)
        # Les directions dominantes passent en premier si le budget coupe un bloc.
        indices = np.flatnonzero(garder)[::-1]
        return (w@vecteurs[:, indices])/np.sqrt(valeurs[indices])

    def _preparer(self, profondeur):
        r = self.v.shape[1]
        if r:
            g = _sym(self.v.T@(self.k@self.v))
            l = cholesky(g, lower=True)
            self.w = solve_triangular(l, self.v.T, lower=True).T
        else:
            self.w = self.v.copy()
        kw = self.k@self.w
        mw = self.lambda_min*(self.m@self.w)
        self.g = _sym(self.w.T@kw)
        self.t = _sym(self.w.T@mw)
        self.d = self.w.T@self.b
        self.defaut_gram = float(np.linalg.norm(self.g-np.eye(r), 2)) if r else 0.
        self._tvaleurs, self._tvecteurs = eigh(self.t)
        self._dmodal = self._tvecteurs.T@self.d
        tmax = float(self._tvaleurs[-1]) if r else 0.
        if tmax > (1+self.defaut_gram)*(1+1e-8):
            raise ValueError("borne lambda_min contredite par un quotient de Rayleigh")
        if self.rho*tmax >= 1:
            raise ArithmeticError("résolvante réduite hors marge positive")
        self.e = self.b-kw@self.d
        self.f = mw-kw@self.t
        ke = self._resoudre(self.e)
        kf = self._resoudre(self.f) if r else self.f.copy()
        delta = _norme_gram(self.e.T@ke)
        p = np.eye(r)
        dnorme = float(np.linalg.norm(self.d, 2)) if r else 0.
        meilleur = np.inf
        self.termes_borne = []
        for q in range(1, profondeur+1):
            # Former les résidus avant les Gram évite J-T^j J T^j par cancellation.
            fp, kfp = self.f@p, kf@p
            queue = (self.rho**q*_norme_gram(fp.T@kfp)*dnorme
                     /(1-self.rho*tmax))
            borne = (delta+queue)**2/(1-self.rho)
            self.termes_borne.append(float(borne))
            meilleur = min(meilleur, borne)
            fd, kfd = fp@self.d, kfp@self.d
            delta += self.rho**q*_norme_gram(fd.T@kfd)
            p = p@self.t
        self.borne_uniforme = float(meilleur)
        if not np.isfinite(self.borne_uniforme):
            raise ArithmeticError("majorant hors domaine flottant")

    def _mu(self, omega):
        omega = float(omega)
        if not np.isfinite(omega) or not 0 <= omega <= self.omega_max:
            raise ValueError("fréquence hors bande construite")
        return omega**2/self.lambda_min

    def coefficients(self, omega):
        mu = self._mu(omega)
        return self._tvecteurs@(self._dmodal/(1-mu*self._tvaleurs[:, None]))

    def transfert(self, omega):
        """Fonctionnel énergétique corrigé, égal à Galerkin en calcul exact."""
        mu = self._mu(omega)
        y = self.coefficients(omega)
        return _sym(self.d.T@y+y.T@self.d-y.T@(self.g-mu*self.t)@y)

    def borne_ponctuelle(self, omega):
        """Audit O(nr), par résidu explicite, indépendant de l'expansion uniforme."""
        mu = self._mu(omega)
        x = self.w@self.coefficients(omega)
        residu = self.b-self.k@x+mu*self.lambda_min*(self.m@x)
        gram = _sym(residu.T@self._resoudre(residu))
        return gram/(1-mu)


class ReductionPorts:
    def __init__(self, k, m, interieur, interface, metrique, **options):
        debut = time.perf_counter()
        k, m = _matrice(k, "K"), _matrice(m, "M")
        i, s = np.asarray(interieur), np.asarray(interface)
        if (i.ndim != 1 or s.ndim != 1 or not len(i) or not len(s)
                or i.dtype.kind not in "iu" or s.dtype.kind not in "iu"
                or k.shape != m.shape
                or not np.array_equal(np.sort(np.concatenate((i, s))), np.arange(k.shape[0]))):
            raise ValueError("intérieur et interface doivent partitionner les indices physiques")
        if m[i][:, s].nnz and np.any(m[i][:, s].data != 0):
            raise ValueError("prototype limité au couplage de masse M_IS nul")
        metrique = np.asarray(metrique, dtype=float)
        if metrique.shape != (len(s), len(s)) or not np.all(np.isfinite(metrique)):
            raise ValueError("métrique des ports incompatible")
        _matrice(metrique, "métrique")
        l = cholesky(metrique, lower=True)
        self.normalisation = solve_triangular(l.T, np.eye(len(s)), lower=False)
        w = self.normalisation
        self.kss = _sym(w.T@(k[s][:, s]@w))
        self.mss = _sym(w.T@(m[s][:, s]@w))
        b = k[i][:, s]@w
        self.interieur = InterieurKrylov(k[i][:, i], m[i][:, i], b, **options)
        self.indices_i, self.indices_s = i, s
        self.taille_interieure = self.interieur.taille_interieure
        self.preparation_s = time.perf_counter()-debut

    def schur(self, omega):
        return self.kss-omega**2*self.mss-self.interieur.transfert(omega)

    def reconstruire(self, omega, port_normalise):
        """Champ candidat pour un port imposé ; charges intérieures nulles."""
        port = np.asarray(port_normalise, dtype=float)
        champ = np.empty((len(self.indices_i)+len(self.indices_s),)+port.shape[1:])
        champ[self.indices_s] = self.normalisation@port
        champ[self.indices_i] = -self.interieur.w@self.interieur.coefficients(omega)@port
        return champ


def reduire(k, m, interieur, interface, metrique, **options):
    return ReductionPorts(k, m, interieur, interface, metrique, **options)
