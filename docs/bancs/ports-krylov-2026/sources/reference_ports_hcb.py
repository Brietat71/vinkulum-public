"""Référence HCB creuse indépendante, construite à partir de K et M physiques.

L'interface conserve ses déplacements physiques. La normalisation W issue
de la métrique de port intervient ensuite : u_S = W y, W=L_metric^-T.
Les modes fixes complètent la déformation statique Psi u_S. Les masses
statiques et les couplages massiques sont intégralement projetés.

La même factorisation de K_II sert à Psi et à eigsh par shift-invert en 0.
La petite paire modale projetée est ensuite diagonalisée : une interrogation
fréquentielle n'effectue que des divisions modales et des petits produits,
sans nouvelle factorisation. Aucun code de solveur concurrent n'est utilisé.
"""
from dataclasses import dataclass
import operator
import time

import numpy as np
from scipy.linalg import cholesky, eigh, solve, solve_triangular
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import LinearOperator, eigsh, splu


def _symetrique(a, nom):
    a = csc_matrix(a, dtype=float, copy=True)
    a.sum_duplicates()
    a.eliminate_zeros()
    if a.shape[0] != a.shape[1] or not np.all(np.isfinite(a.data)):
        raise ValueError(f"{nom} doit être carrée et finie")
    delta = a-a.T
    maximum = np.max(np.abs(a.data), initial=0.)
    if np.max(np.abs(delta.data), initial=0.) > 64*np.finfo(float).eps*maximum:
        raise ValueError(f"{nom} doit être symétrique")
    # Seules les dissymétries d'arrondi ci-dessus sont absorbées.
    return ((a+a.T)*0.5).tocsc()


def _indices(indices, n, nom):
    indices = np.asarray(indices)
    if indices.ndim != 1 or indices.dtype.kind not in "iu":
        raise ValueError(f"{nom} doit contenir des indices entiers")
    if len(np.unique(indices)) != len(indices) or np.any(indices < 0) or np.any(indices >= n):
        raise ValueError(f"{nom} contient des indices invalides ou répétés")
    return indices.astype(np.intp, copy=True)


def _factorise_spd(a, nom):
    """LU à permutation symétrique, contrôle des pivots positifs de LDL^T.

    Pour une matrice symétrique, sans permutation de lignes supplémentaire,
    les diagonales de U sont les pivots D. Un pivot non positif ou une
    permutation supplémentaire fait refuser le problème. Ce contrôle
    flottant ne fournit pas de borne certifiée sur lambda_min.
    """
    try:
        facteur = splu(a, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.,
                       options={"SymmetricMode": True})
    except RuntimeError as exc:
        raise ValueError(f"{nom} n'admet pas la factorisation SPD requise") from exc
    pivots = facteur.U.diagonal()
    if (not np.array_equal(facteur.perm_r, facteur.perm_c)
            or not np.all(np.isfinite(pivots)) or np.any(pivots <= 0)):
        raise ValueError(f"{nom} doit être définie positive")
    return facteur


def _norme_gram(gram):
    if not gram.size:
        return 0.
    valeurs = eigh(0.5*(gram+gram.T), eigvals_only=True, check_finite=False)
    maximum = max(float(valeurs[-1]), 0.)
    if valeurs[0] < -1e-8*max(maximum, np.finfo(float).tiny):
        raise ArithmeticError("Gram résiduel HCB numériquement indéfini")
    return float(np.sqrt(maximum))


@dataclass
class ReferenceHCB:
    """Réduction HCB ; les coordonnées réduites sont (q_modal, y_port).

    ``taille_interieure`` est le nombre de modes retenus, et
    ``taille_interieure_initiale`` le nombre de DDL avant réduction.
    ``k_reduit`` et ``m_reduit`` gardent les projections complètes, y compris
    les écarts d'arrondi à la diagonalité/orthogonalité modale.
    """
    interieur: np.ndarray
    interface: np.ndarray
    psi: np.ndarray
    phi: np.ndarray
    w: np.ndarray
    k_reduit: np.ndarray
    m_reduit: np.ndarray
    valeurs_modales: np.ndarray
    changement_modal: np.ndarray
    couplage_k: np.ndarray
    couplage_m: np.ndarray
    preparation_s: float
    compteurs: dict
    _n: int
    _kii: csc_matrix
    _mii: csc_matrix
    _residu_statique: np.ndarray
    _facteur: object
    temps_borne_s: float = 0.
    bilan_borne: dict | None = None

    @property
    def taille_interieure(self):
        return self.phi.shape[1]

    @property
    def taille_interieure_initiale(self):
        return len(self.interieur)

    @property
    def temps_preparation(self):
        return self.preparation_s

    def _frequence(self, omega):
        omega = float(omega)
        if not np.isfinite(omega) or omega < 0 or omega > np.sqrt(np.finfo(float).max):
            raise ValueError("omega doit être fini, positif ou nul et représentable au carré")
        z = omega*omega
        diviseurs = self.valeurs_modales-z
        if len(diviseurs):
            echelle = np.maximum(np.abs(self.valeurs_modales), z)
            if np.any(np.abs(diviseurs) <= 32*np.finfo(float).eps*echelle):
                raise ValueError("fréquence au voisinage numérique d'un pôle intérieur HCB")
        return z, diviseurs

    def schur(self, omega):
        """W.T @ S_HCB(omega) @ W, sans factorisation par fréquence."""
        z, diviseurs = self._frequence(omega)
        r = self.taille_interieure
        s = self.k_reduit[r:, r:]-z*self.m_reduit[r:, r:]
        if r:
            couplage = self.couplage_k-z*self.couplage_m
            s = s-couplage.T@(couplage/diviseurs[:, None])
            self.compteurs["eliminations_modales"] += 1
            self.compteurs["divisions_modales_rhs"] += r*len(self.interface)
        self.compteurs["appels_schur"] += 1
        return 0.5*(s+s.T)

    def reconstruit(self, omega, y):
        """Champ physique, sans charge intérieure, pour des coordonnées de port y.

        Les ports physiques valent W y. y peut être un vecteur ou une
        matrice de plusieurs colonnes. La reconstruction utilise les mêmes
        matrices projetées que ``schur``.
        """
        z, diviseurs = self._frequence(omega)
        y = np.asarray(y, dtype=float)
        if y.ndim not in (1, 2) or y.shape[0] != len(self.interface) or not np.all(np.isfinite(y)):
            raise ValueError("y doit contenir des coordonnées finies dans l'espace du port normalisé")
        vectoriel = y.ndim == 1
        y = y[:, None] if vectoriel else y
        us = self.w@y
        ui = self.psi@us
        if self.taille_interieure:
            modal = -((self.couplage_k-z*self.couplage_m)@y)/diviseurs[:, None]
            ui = ui+self.phi@(self.changement_modal@modal)
            self.compteurs["eliminations_modales"] += 1
            self.compteurs["divisions_modales_rhs"] += self.taille_interieure*y.shape[1]
        u = np.zeros((self._n, y.shape[1]))
        u[self.interface], u[self.interieur] = us, ui
        self.compteurs["appels_reconstruction"] += 1
        return u[:, 0] if vectoriel else u

    def borne_uniforme(self, lambda_min, omega_max):
        """Majorant conditionnel de ||S_HCB-S_exact|| sur [0, omega_max].

        Hypothèse externe : K_II >= lambda_min M_II > 0. On pose
        rho=omega_max²/lambda_min<1 et Mhat=lambda_min M_II. Aucun Ritz
        intérieur n'est utilisé comme borne inférieure du spectre complet.
        La normalisation statique porte Psi=W au port et psi W à l'intérieur.

        Pour Q=Phi L_K^-T, les produits réels définissent G=Q.T K Q,
        T=Q.T Mhat Q, C=Mhat psi W, D=Q.T C, J=Q.T R0. Avec
        E0=R0-KQ G^-1 J, E=C-KQ G^-1 D, F=MhatQ-KQ G^-1 T et
        Y=(G-mu T)^-1(mu D-J), le résidu vaut exactement

            R = E0 - mu E - mu F Y.

        Les normes K^-1 sont calculées par Gram, avec la LU déjà construite.
        La formule garde donc les défauts statiques et le Gram réel G.
        Le carré de la majoration du résidu, divisé par 1-rho, majore l'erreur
        du Schur énergétique projeté. En arithmétique exacte, R0=0 et G=I
        retrouvent la formule résiduelle HCB habituelle.

        Les valeurs et contrôles flottants ne constituent pas un certificat
        machine : ni coercivité externe ni arrondis dirigés ne sont prouvés.
        Le coût s'ajoute à la préparation, dans ``temps_borne_s`` ; les
        détails et hypothèses sont conservés dans ``bilan_borne``.
        """
        debut = time.perf_counter()
        lambda_min, omega_max = float(lambda_min), float(omega_max)
        if (not np.isfinite(lambda_min) or lambda_min <= 0
                or not np.isfinite(omega_max) or omega_max < 0
                or omega_max > np.sqrt(np.finfo(float).max)):
            raise ValueError("lambda_min doit être positif, omega_max positif ou nul, tous deux finis")
        rho = omega_max**2/lambda_min
        if not np.isfinite(rho) or rho >= 1:
            raise ValueError("omega_max² < lambda_min est requis")
        r, p = self.taille_interieure, len(self.interface)
        c = lambda_min*(self._mii@(self.psi@self.w))
        if r:
            lk = cholesky(self.k_reduit[:r, :r], lower=True)
            q = solve_triangular(lk, self.phi.T, lower=True).T
            kq, mq = self._kii@q, lambda_min*(self._mii@q)
            g, t = q.T@kq, q.T@mq
            g, t = 0.5*(g+g.T), 0.5*(t+t.T)
            d, j = q.T@c, q.T@self._residu_statique
            valeurs_gram = eigh(g, eigvals_only=True, check_finite=False)
            gmin = float(valeurs_gram[0])
            tmax = float(eigh(t, eigvals_only=True, check_finite=False)[-1])
            quotient_max = float(eigh(t, g, eigvals_only=True, check_finite=False)[-1])
            if quotient_max > 1+1e-8:
                raise ValueError("lambda_min contredite par un quotient de Rayleigh intérieur")
            marge = gmin-rho*tmax
            if marge <= 0:
                raise ArithmeticError("résolvante modale HCB sans marge positive")
            coefficients = solve(g, np.column_stack((j, d, t)), assume_a="pos", check_finite=False)
            self.compteurs["factorisations_Gram_borne"] += 1
            self.compteurs["colonnes_rhs_Gram_borne"] += 2*p+r
            e0 = self._residu_statique-kq@coefficients[:, :p]
            e = c-kq@coefficients[:, p:2*p]
            f = mq-kq@coefficients[:, 2*p:]
            defaut_gram = float(np.max(np.abs(valeurs_gram-1.)))
            dnorme, jnorme = float(np.linalg.norm(d, 2)), float(np.linalg.norm(j, 2))
        else:
            e0, e, f = self._residu_statique, c, np.zeros((len(self.interieur), 0))
            gmin, tmax, marge = 1., 0., 1.
            defaut_gram = dnorme = jnorme = 0.
        rhs = np.column_stack((e0, e, f))
        resolu = self._facteur.solve(rhs)
        self.compteurs["resolutions_K_II"] += 1
        self.compteurs["colonnes_rhs_K_II"] += rhs.shape[1]
        self.compteurs["resolutions_borne_K_II"] += 1
        self.compteurs["colonnes_rhs_borne_K_II"] += rhs.shape[1]
        delta_0 = _norme_gram(e0.T@resolu[:, :p])
        delta_e = _norme_gram(e.T@resolu[:, p:2*p])
        delta_f = _norme_gram(f.T@resolu[:, 2*p:])
        residu = delta_0+rho*delta_e+rho*delta_f*(rho*dnorme+jnorme)/marge
        borne = residu**2/(1-rho)
        if not np.isfinite(borne):
            raise ArithmeticError("majorant HCB hors domaine flottant")
        temps = time.perf_counter()-debut
        self.temps_borne_s += temps
        self.bilan_borne = {
            "lambda_min": lambda_min, "omega_max": omega_max, "rho": rho,
            "borne_uniforme": float(borne), "temps_s": temps,
            "delta_statique": delta_0, "delta_E": delta_e, "delta_F": delta_f,
            "defaut_gram": defaut_gram,
            "norme_D": dnorme, "norme_J": jnorme,
            "lambda_min_Gram": gmin, "norme_T": tmax, "marge_resolvante": marge,
            "statut": "majorant_conditionnel_flottant_non_certifie",
        }
        return float(borne)


def reduire(k, m, interieur, interface, n_modes, metrique):
    """Prépare une référence HCB creuse à interfaces fixes.

    Domaine : K_II et M_II définies positives, K/M symétriques, M_IS=0.
    Les coordonnées d'intérieur et d'interface doivent partitionner tous
    les DDL. La métrique SPD définit uniquement la normalisation du port ;
    sa valeur ne remplace pas la raideur statique obtenue par projection.

    eigsh(sigma=0, tol=1e-11) est utilisé si n_II>32 et n_modes<n_II-1.
    Le calcul dense est limité aux petits intérieurs ou aux demandes de
    quasi-totalité/tous les modes, explicitement indiquées dans les compteurs.
    Tous les coûts de préparation sont inclus dans ``preparation_s`` ; les
    appels et colonnes de résolution par K_II sont aussi dénombrés.
    """
    debut = time.perf_counter()
    k, m = _symetrique(k, "K"), _symetrique(m, "M")
    n = k.shape[0]
    if m.shape != k.shape:
        raise ValueError("K et M doivent avoir les mêmes dimensions")
    ii, ss = _indices(interieur, n, "interieur"), _indices(interface, n, "interface")
    if not len(ii) or not len(ss) or not np.array_equal(np.sort(np.r_[ii, ss]), np.arange(n)):
        raise ValueError("intérieur et interface non vides doivent partitionner tous les DDL")
    if isinstance(n_modes, (bool, np.bool_)):
        raise ValueError("n_modes doit être un entier")
    try:
        r = operator.index(n_modes)
    except TypeError as exc:
        raise ValueError("n_modes doit être un entier") from exc
    if not 0 <= r <= len(ii):
        raise ValueError("n_modes doit appartenir à [0, taille intérieure]")
    metrique = _symetrique(metrique, "metrique").toarray()
    if metrique.shape != (len(ss), len(ss)):
        raise ValueError("la métrique doit avoir la dimension du port")
    w = solve_triangular(cholesky(metrique, lower=True).T, np.eye(len(ss)), lower=False)
    kii, kis, kss = k[ii, :][:, ii].tocsc(), k[ii, :][:, ss].tocsc(), k[ss, :][:, ss].toarray()
    mii, mis, mss = m[ii, :][:, ii].tocsc(), m[ii, :][:, ss].tocsc(), m[ss, :][:, ss].toarray()
    if mis.nnz:
        raise ValueError("la référence initiale exige M_IS=0")
    compteurs = {
        "factorisations_K_II": 0, "factorisations_M_II_controle": 0,
        "resolutions_K_II": 0, "colonnes_rhs_K_II": 0,
        "resolutions_statiques_K_II": 0, "resolutions_modes_K_II": 0,
        "resolutions_borne_K_II": 0, "colonnes_rhs_borne_K_II": 0,
        "factorisations_Gram_borne": 0, "colonnes_rhs_Gram_borne": 0,
        "eigensolves_denses_interieurs": 0, "eigensolves_creux_interieurs": 0,
        "eigensolves_petits_projetes": 0, "taille_dense_interieure": 0,
        "appels_schur": 0, "appels_reconstruction": 0,
        "factorisations_frequence": 0, "eliminations_modales": 0,
        "divisions_modales_rhs": 0,
    }
    # Les masses diagonales des fixtures n'exigent aucune factorisation.
    off = mii-diags(mii.diagonal(), format="csc")
    off.eliminate_zeros()
    if not off.nnz:
        if np.any(mii.diagonal() <= 0):
            raise ValueError("M_II doit être définie positive")
    else:
        _factorise_spd(mii, "M_II")
        compteurs["factorisations_M_II_controle"] += 1
    facteur = _factorise_spd(kii, "K_II")
    compteurs["factorisations_K_II"] += 1
    compteurs["nnz_facteurs_K_II"] = facteur.L.nnz+facteur.U.nnz

    def resout(rhs, mode=False):
        compteurs["resolutions_K_II"] += 1
        compteurs["colonnes_rhs_K_II"] += 1 if rhs.ndim == 1 else rhs.shape[1]
        compteurs["resolutions_modes_K_II" if mode else "resolutions_statiques_K_II"] += 1
        return facteur.solve(rhs)

    psi = resout(-kis.toarray())
    if r == 0:
        phi = np.zeros((len(ii), 0))
        compteurs["methode_modes"] = "aucun"
    elif len(ii) <= 32 or r >= len(ii)-1:
        _, phi = eigh(kii.toarray(), mii.toarray(), subset_by_index=(0, r-1), check_finite=False)
        compteurs["eigensolves_denses_interieurs"] += 1
        compteurs["taille_dense_interieure"] = len(ii)
        compteurs["methode_modes"] = "dense_petit" if len(ii) <= 32 else "dense_quasi_complet_demande"
    else:
        inverse = LinearOperator(kii.shape, dtype=float,
                                 matvec=lambda v: resout(v, True),
                                 matmat=lambda v: resout(v, True))
        valeurs, phi = eigsh(kii, k=r, M=mii, sigma=0., which="LM", OPinv=inverse,
                             tol=1e-11, ncv=min(len(ii), max(2*r+1, 20)),
                             v0=np.random.default_rng(7321).standard_normal(len(ii)))
        phi = phi[:, np.argsort(valeurs)]
        compteurs["eigensolves_creux_interieurs"] += 1
        compteurs["methode_modes"] = "eigsh_shift_inverse_zero"
    # T=[Phi, Psi W ; 0, W] dans l'ordre intérieur/interface.
    # Le petit résidu K_II Psi+K_IS est conservé dans tous les blocs de K_r.
    residu_statique = kii@psi+kis.toarray()
    kqq, mqq = phi.T@(kii@phi), phi.T@(mii@phi)
    kqp = phi.T@residu_statique@w
    mqp = phi.T@(mii@psi)@w
    kpp = w.T@(kss+kis.T@psi+psi.T@residu_statique)@w
    mpp = w.T@(mss+psi.T@(mii@psi))@w
    kr = np.block([[kqq, kqp], [kqp.T, kpp]])
    mr = np.block([[mqq, mqp], [mqp.T, mpp]])
    kr, mr = 0.5*(kr+kr.T), 0.5*(mr+mr.T)
    if r:
        valeurs, changement = eigh(kr[:r, :r], mr[:r, :r], check_finite=False)
        if np.any(valeurs <= 0):
            raise ValueError("la paire modale projetée doit rester définie positive")
        couplage_k, couplage_m = changement.T@kqp, changement.T@mqp
        compteurs["eigensolves_petits_projetes"] += 1
    else:
        valeurs, changement = np.empty(0), np.empty((0, 0))
        couplage_k, couplage_m = kqp, mqp
    return ReferenceHCB(ii, ss, psi, phi, w, kr, mr, valeurs, changement,
                        couplage_k, couplage_m, time.perf_counter()-debut, compteurs, n,
                        kii, mii, residu_statique@w, facteur)
