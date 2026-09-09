"""Chaîne construite par sous-structures indépendantes, avec ports physiques.

Chaque segment a ses deux extrémités comme ports, une masse demi-nodale
à chaque bout et un intérieur coercif. Les segments sont tous construits ;
aucune réutilisation par périodicité n'est utilisée dans ce témoin.
"""
import operator
import time

import numpy as np
from scipy.sparse import diags

from ports_releves import PortsReleves, AssemblagePorts


class ChaineAssemblee:
    def __init__(self, n, segments, omega_max, tolerance=1e-10, raideur=1., masse=1.):
        debut = time.perf_counter()
        n, segments = operator.index(n), operator.index(segments)
        if n < 2 or segments < 1 or n % segments or n//segments < 2:
            raise ValueError("segments égaux de longueur au moins deux requis")
        if not np.isfinite(raideur) or not np.isfinite(masse) or min(raideur, masse) <= 0:
            raise ValueError("raideur et masse positives requises")
        longueur = n//segments
        racine = np.sqrt(raideur)
        d = diags((-np.full(longueur, racine), np.full(longueur, racine)),
                  (0, 1), shape=(longueur, longueur+1), format="csc")
        masses = np.full(longueur+1, masse)
        masses[[0, -1]] *= .5
        m = diags(masses, format="csc")
        i, s = np.arange(1, longueur), np.array([0, longueur])
        metric = raideur/longueur*np.eye(2)
        lam = .999*4*raideur/masse*np.sin(np.pi/(2*longueur))**2
        self.echelle = np.sqrt(n/raideur)
        self.structures, applications = [], []
        for bloc in range(segments):
            r = PortsReleves(d, m, i, s, metric, lam, omega_max,
                            tolerance=tolerance/(2*segments))
            e = np.zeros((2, segments))
            if bloc:
                e[0, bloc-1] = self.echelle
            e[1, bloc] = self.echelle
            self.structures.append(r)
            applications.append(e)
        supplement = np.zeros((segments, segments))
        supplement[-1, -1] = .5*masse*self.echelle**2
        self.reduction = AssemblagePorts(self.structures, applications, m_externe=supplement)
        self.n, self.segments, self.longueur = n, segments, longueur
        self.preparation_s = time.perf_counter()-debut
        self.taille_interieure = sum(r.taille_interieure for r in self.structures)
        self.borne_uniforme = self.reduction.borne_uniforme

    def reponse(self, omega):
        force = np.zeros(self.segments)
        force[-1] = 1.
        return self.reduction.reponse(omega, force)

    def reconstruire(self, omega, deplacement_global):
        champ = np.zeros(self.n+1)
        for j, (r, e) in enumerate(zip(self.structures, self.reduction.applications)):
            local = r.reconstruire(omega, e@deplacement_global)
            champ[j*self.longueur:(j+1)*self.longueur+1] = local
        return champ[1:]
