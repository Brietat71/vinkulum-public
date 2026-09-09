"""Ports relevés en énergie et assemblage de sous-structures linéaires.

Le QR définit un modèle énergétique calculé, distinct des matrices K issues
de produits flottants D.T D. L'enveloppe finale inclut son écart au modèle D
d'entrée. Elle ne certifie ni les arrondis, ni la constante spectrale fournie.
"""
import time

import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import splu

from condensation_energie import CondensationEnergie
from ports_krylov import InterieurKrylov, _matrice, _sym, _norme_gram


class PortsReleves:
    def __init__(self, d, m, interieur, interface, metrique, lambda_min, omega_max,
                 tolerance=1e-8, max_blocs=12, max_directions=128, budget_qr=10_000_000,
                 methode="qr"):
        debut = time.perf_counter()
        self.methode = methode
        if methode == "qr":
            self.qr = CondensationEnergie(d, interieur, interface, metrique, budget_qr)
        elif methode == "lu_energie":
            from condensation_energie_lu import CondensationEnergieLU
            self.qr = CondensationEnergieLU(d, interieur, interface, metrique)
        else:
            raise ValueError("méthode inconnue : qr ou lu_energie requis")
        self.m = _matrice(m, "M")
        if self.m.shape[0] != self.qr.d.shape[1]:
            raise ValueError("masse et énergie de dimensions différentes")
        i, s, w = self.qr.i, self.qr.s, self.qr.w
        mii, mis, mss = self.m[i][:, i], self.m[i][:, s], self.m[s][:, s]
        psi = self.qr.psi
        c = mii@psi+mis@w
        self.m0 = _sym(psi.T@(mii@psi)+psi.T@(mis@w)+w.T@(mis.T@psi)+w.T@(mss@w))
        self.lambda_min, self.omega_max = float(lambda_min), float(omega_max)
        self.tolerance = float(tolerance)
        rho = self.omega_max**2/self.lambda_min
        if not 0 <= rho < 1:
            raise ValueError("bande hors coercivité établie")
        self.p = len(s)
        b = np.column_stack((self.qr.r0, self.lambda_min*c))
        # Le défaut statique minuscule peut être absent de la semence retenue :
        # B reste entier dans les résidus, donc aucune erreur n'est effacée.
        seuil = 1e-3*np.sqrt(self.tolerance*(1-rho)/(1+rho*rho))
        self.interieur = InterieurKrylov(
            self.qr.matrice(), mii, b, self.lambda_min, self.omega_max,
            tolerance=self.tolerance/(1+rho*rho), max_blocs=max_blocs,
            max_directions=max_directions, facteur_energie=self.qr,
            seuil_semence=seuil)
        self.borne_qr = (1+rho*rho)*self.interieur.borne_uniforme
        self.taille_interieure = self.interieur.taille_interieure
        self.audit_bande = self.audit_uniforme()
        self.borne_uniforme = float(self.audit_bande["majorant_total"])
        self.statut = "tolerance_estimee" if self.borne_uniforme <= self.tolerance else "tolerance_non_atteinte"
        self.preparation_s = time.perf_counter()-debut

    def _charge(self, omega):
        mu = self.interieur._mu(omega)
        return np.vstack((np.eye(self.p), -mu*np.eye(self.p)))

    def enrichir(self):
        """Enrichit l'intérieur et réévalue l'enveloppe dans D d'entrée."""
        debut = time.perf_counter()
        change = self.interieur.enrichir()
        if change:
            self.borne_qr = (1+self.interieur.rho**2)*self.interieur.borne_uniforme
            self.taille_interieure = self.interieur.taille_interieure
            self.audit_bande = self.audit_uniforme()
            self.borne_uniforme = float(self.audit_bande["majorant_total"])
            self.statut = ("tolerance_estimee" if self.borne_uniforme <= self.tolerance
                           else "tolerance_non_atteinte")
        self.preparation_s += time.perf_counter()-debut
        return change

    def schur(self, omega):
        charge = self._charge(omega)
        h = self.interieur.transfert(omega)
        return _sym(self.qr.k0-omega**2*self.m0-charge.T@h@charge)

    def reconstruire(self, omega, port_normalise):
        port = np.asarray(port_normalise, dtype=float)
        if port.ndim not in (1, 2) or port.shape[0] != self.p or not np.all(np.isfinite(port)):
            raise ValueError("coordonnées de port finies requises")
        champ = np.empty((self.m.shape[0],)+port.shape[1:])
        champ[self.qr.s] = self.qr.w@port
        y = self.interieur.coefficients(omega)@self._charge(omega)
        champ[self.qr.i] = self.qr.psi@port-self.interieur.w@y@port
        return champ

    def audit_modele(self):
        return self.qr.audit()

    def audit_energie(self, omega):
        """Audit dans D d'entrée, par énergie et coercivité massique fournie.

        Exactement : erreur de l'énergie candidate <= R.T M_II^-1 R /
        (lambda_min-omega²). Ce diagnostic comprend l'écart entre énergie
        candidate et Schur renvoyé. Ses opérations flottantes ne sont pas
        encadrées ; la borne lambda_min doit viser le modèle D d'entrée.
        Coût hors préparation : reconstruction, produits D et résolution M.
        """
        x = self.reconstruire(omega, np.eye(self.p))
        dx, mx = self.qr.d@x, self.m@x
        candidat = _sym(dx.T@dx-omega**2*(x.T@mx))
        residu = self.qr.di.T@dx-omega**2*mx[self.qr.i]
        mii = self.m[self.qr.i][:, self.qr.i]
        dual = splu(mii).solve(residu)
        erreur = _sym(residu.T@dual)/(self.lambda_min-omega**2)
        ecart = float(np.linalg.norm(candidat-self.schur(omega), 2))
        return dict(schur_energetique=candidat, residu_massique_borne=erreur,
                    ecart_evaluation=ecart,
                    majorant_total=ecart+float(np.linalg.norm(erreur, 2)),
                    certification_machine=False)

    def audit_uniforme(self):
        """Enveloppe dans D d'entrée, incluant l'écart du modèle QR.

        L'expansion résiduelle utilise la norme duale de M_II ; l'écart
        entre les deux fonctionnels est borné sur leurs petites projections.
        Les arrondis et la constante spectrale restent non certifiés.
        """
        v, t = self.interieur.w, self.interieur.t
        r, p = v.shape[1], self.p
        rho = self.interieur.rho
        a0, a1 = self.interieur.d[:, :p], -self.interieur.d[:, p:]
        x0 = np.zeros((self.m.shape[0], p))
        x0[self.qr.i], x0[self.qr.s] = self.qr.psi, self.qr.w
        dx0, dv = self.qr.d@x0, self.qr.di@v
        mx0 = self.m@x0
        mii = self.m[self.qr.i][:, self.qr.i]
        mv = self.lambda_min*(mii@v)
        kv = self.qr.di.T@dv
        c = self.lambda_min*mx0[self.qr.i]
        r0 = self.qr.di.T@dx0
        e0, e1 = r0-kv@a0, -c-kv@a1
        f = mv-kv@t
        c1 = e1+f@a0
        dd = t@a0+a1
        mf = splu(mii)

        def norme(x):
            return _norme_gram(x.T@mf.solve(x)) if x.shape[1] else 0.

        delta = norme(e0)+rho*norme(c1)
        tf = float(np.linalg.norm(t, 2)) if r else 0.
        dn = float(np.linalg.norm(dd, 2)) if r else 0.
        y_max = ((float(np.linalg.norm(a0, 2))+rho*float(np.linalg.norm(a1, 2)))
                 /(1-rho*tf)) if r else 0.
        puissance = np.eye(r)
        meilleure = np.inf
        enveloppes = []
        profondeur = max(2, len(self.interieur.historique)+2)
        for q in range(2, profondeur+1):
            fp = f@puissance
            queue = rho**q*norme(fp)*dn/(1-rho*tf)
            borne = (delta+queue)**2/(self.lambda_min*(1-rho))
            enveloppes.append(float(borne))
            meilleure = min(meilleure, borne)
            delta += rho**q*norme(fp@dd)
            puissance = puissance@t
        k0d = _sym(dx0.T@dx0)
        m0d = _sym(x0.T@mx0)
        ecart0 = np.linalg.norm(k0d-self.qr.k0, 2)+self.omega_max**2*np.linalg.norm(m0d-self.m0, 2)
        if r:
            ecart_cross = (np.linalg.norm(dv.T@dx0-a0, 2)
                           +rho*np.linalg.norm(v.T@c+a1, 2))
            ecart_interieur = (np.linalg.norm(dv.T@dv-self.interieur.g, 2)
                               +rho*np.linalg.norm(v.T@mv-t, 2))
        else:
            ecart_cross = ecart_interieur = 0.
        modele = float(ecart0+2*ecart_cross*y_max+ecart_interieur*y_max*y_max)
        return dict(majorant_total=modele+meilleure, majorant_residu=meilleure,
                    ecart_modele=modele, enveloppes=enveloppes,
                    norme_champ_dynamique=y_max, certification_machine=False)


class AssemblagePorts:
    """Assemble des ports physiques par applications linéaires E_s.

    Les applications données expriment les déplacements de port physiques
    à partir des coordonnées globales choisies. Leur échelle définit la
    norme d'erreur globale ; les unités doivent être choisies par l'appelant.
    La borne additionne les majorantes bilatérales en ordre de Loewner.
    """
    def __init__(self, sous_structures, applications, k_externe=None, m_externe=None):
        if not sous_structures or len(sous_structures) != len(applications):
            raise ValueError("une application par sous-structure requise")
        self.sous_structures = list(sous_structures)
        self.applications = []
        taille = None
        for r, e in zip(sous_structures, applications):
            e = np.asarray(e, dtype=float)
            if e.ndim != 2 or e.shape[0] != r.p or not np.all(np.isfinite(e)):
                raise ValueError("application de port invalide")
            if taille is None:
                taille = e.shape[1]
            if e.shape[1] != taille or not 0 < taille <= 1024:
                raise ValueError("assemblage réduit limité à 1024 coordonnées globales")
            self.applications.append(solve(r.qr.w, e))
        self.k_externe = np.zeros((taille, taille)) if k_externe is None else np.asarray(k_externe, dtype=float)
        self.m_externe = np.zeros((taille, taille)) if m_externe is None else np.asarray(m_externe, dtype=float)
        for a in (self.k_externe, self.m_externe):
            if a.shape != (taille, taille):
                raise ValueError("matrice externe de taille incorrecte")
            _matrice(a, "matrice externe")
        self.omega_max = min(r.omega_max for r in sous_structures)
        self.erreur_matrice = sum(r.borne_uniforme*(e.T@e)
                                   for r, e in zip(sous_structures, self.applications))
        self.borne_uniforme = float(np.linalg.norm(self.erreur_matrice, 2))

    def schur(self, omega):
        if not np.isfinite(omega) or not 0 <= omega <= self.omega_max:
            raise ValueError("fréquence hors bande commune")
        s = self.k_externe-omega**2*self.m_externe
        for r, e in zip(self.sous_structures, self.applications):
            s = s+e.T@r.schur(omega)@e
        return _sym(s)

    def reponse(self, omega, force):
        s = self.schur(omega)
        u = solve(s, force, assume_a="sym")
        sigma = float(np.linalg.svd(s, compute_uv=False)[-1])
        marge = sigma-self.borne_uniforme
        borne = (self.borne_uniforme/marge*np.linalg.norm(u, 2)) if marge > 0 else None
        return dict(deplacement=u, borne_erreur_norme=borne,
                    sigma_min_reduit=sigma, marge=marge, certification_machine=False)
