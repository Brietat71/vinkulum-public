"""Énergies élémentaires des deux fixtures homogènes de ``modeles_ports``.

La convention est E(u)=1/2 ||D u||², donc K=D.T D en arithmétique exacte.
D est assemblé directement depuis les différences et déformations locales ;
aucune matrice de raideur globale ni factorisation n'est nécessaire ici.
Les masses et coordonnées restent physiques. Ces représentations concernent
exclusivement la chaîne scalaire et la console droite, homogène, intégrée.

Le produit D.T D et l'assemblage natif de K suivent des chemins d'arrondi
différents. Leur proximité est testable ; leur égalité bit à bit n'est ni
supposée ni requise. Ce module n'est pas une factorisation universelle du
noyau, ni un certificat d'arrondi.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, diags

from modeles_ports import _nombre, _positif


@dataclass(frozen=True)
class EnergiePorts:
    d: csc_matrix
    m: csc_matrix
    interieur: np.ndarray
    interface: np.ndarray
    metadata: dict


def energie_chaine(n, raideur=1.0, masse=1.0):
    """n ressorts et n masses, u_0=0 ; (D u)_j=sqrt(k)(u_j-u_{j-1})."""
    n = _nombre(n)
    raideur, masse = _positif("raideur", raideur), _positif("masse", masse)
    racine = np.sqrt(raideur)
    d = diags((-np.full(n-1, racine), np.full(n, racine)), (-1, 0), format="csc")
    m = diags(np.full(n, masse), format="csc")
    return EnergiePorts(d, m, np.arange(n-1), np.array([n-1]), {
        "famille": "chaine", "n": n, "ddl": n, "ports": [n-1],
        "raideur_N_par_m": raideur, "masse_nodale_kg": masse,
        "masse_mobile_kg": n*masse,
        "representation_energie": "differences_elementaires",
        "convention_energie": "E = 1/2 ||D u||^2",
    })


def _parametres_console(metadata):
    if (metadata.get("famille") != "console"
            or metadata.get("formulation_poutre") != "integree"
            or metadata.get("coefficient_cisaillement") != 1.0
            or metadata.get("poids_extremites_masse_et_inertie") != 0.5
            or metadata.get("ordre_ddl_par_corps") != ["tx", "ty", "tz", "rx", "ry", "rz"]):
        raise ValueError("représentation réservée à la fixture console homogène intégrée")
    n = _nombre(metadata["n"])
    noms = ("longueur_m", "EA_N", "GA_N", "EI_N_m2", "GJ_N_m2",
            "rayon_m", "masse_volumique_kg_par_m3")
    valeurs = [_positif(nom, metadata[nom]) for nom in noms]
    if metadata.get("ports") != list(range(6*(n-1), 6*n)):
        raise ValueError("les ports doivent être les six coordonnées du dernier corps")
    return (n, *valeurs)


def energie_console(metadata):
    """Énergie nodale exacte de la console linéaire de Timoshenko intégrée.

    Pour un élément de longueur l, Ceff=1/(1/GA+l²/(12 EI)). L'énergie
    est l/2 [EA epsilon² + GJ kappa_x² + Ceff (gamma_y²+gamma_z²)
    + EI (kappa_y²+kappa_z²)], avec

      gamma_y=(ty_b-ty_a)/l-(rz_a+rz_b)/2,
      gamma_z=(tz_b-tz_a)/l+(ry_a+ry_b)/2,
      kappa_j=(rj_b-rj_a)/l, epsilon=(tx_b-tx_a)/l.

    La correction Ceff condense la flexion interne de l'élément intégré.
    Les lignes sont (epsilon,gamma_y,gamma_z,kappa_x,kappa_y,kappa_z),
    multipliées par les racines des rigidités et de l. Les six colonnes de
    racine sont omises dès l'assemblage. D a la forme (6n,6n), avec 16n-8
    coefficients non nuls. Aucun appel au noyau n'est effectué.
    """
    n, longueur, ea, ga, ei, gj, rayon, rho = _parametres_console(metadata)
    pas = longueur/n
    ceff = 1/(1/ga+pas**2/(12*ei))
    axial, torsion, flexion = np.sqrt(ea/pas), np.sqrt(gj/pas), np.sqrt(ei/pas)
    shear = np.sqrt(pas*ceff)
    # Triplets locaux (ligne de déformation, DDL de a ou b, coefficient).
    # Les DDL de b portent les indices 6,...,11.
    lignes_locales = np.array([0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5])
    colonnes_locales = np.array([0, 6, 1, 7, 5, 11, 2, 8, 4, 10, 3, 9, 4, 10, 5, 11])
    valeurs_locales = np.array([-axial, axial, -shear/pas, shear/pas, -shear/2, -shear/2,
                                -shear/pas, shear/pas, shear/2, shear/2,
                                -torsion, torsion, -flexion, flexion, -flexion, flexion])
    decalages = (6*np.arange(n))[:, None]
    lignes = (decalages+lignes_locales).ravel()
    colonnes = (decalages+colonnes_locales-6).ravel()
    valeurs = np.broadcast_to(valeurs_locales, (n, len(valeurs_locales))).ravel()
    libres = colonnes >= 0
    d = coo_matrix((valeurs[libres], (lignes[libres], colonnes[libres])),
                   shape=(6*n, 6*n)).tocsc()
    aire, moment = np.pi*rayon**2, np.pi*rayon**4/4
    masse_par_noeud = rho*pas*np.array([aire, aire, aire, 2*moment, moment, moment])
    masses = np.broadcast_to(masse_par_noeud, (n, 6)).copy()
    masses[-1] *= 0.5
    m = diags(masses.ravel(), format="csc")
    description = dict(metadata)
    description.update({"representation_energie": "deformations_elementaires_integrees",
                        "convention_energie": "E = 1/2 ||D u||^2",
                        "ordre_lignes_element": ["epsilon", "gamma_y", "gamma_z",
                                                 "kappa_x", "kappa_y", "kappa_z"]})
    return EnergiePorts(d, m, np.arange(6*(n-1)), np.arange(6*(n-1), 6*n), description)


def lifting_chaine(n):
    """Déplacement intérieur physique sous déplacement unitaire du port : i/n."""
    n = _nombre(n)
    return (np.arange(1, n)/n)[:, None]


def lifting_console(metadata):
    """Relèvement statique intérieur physique Psi, de forme (6(n-1),6).

    Le port contient (tx,ty,tz,rx,ry,rz) en m/rad, sans normalisation.
    Pour t=x/L, z=t(1-t), eta=EI/(GA L²), delta=1+12eta, chaque plan
    (v,theta) est donné par le produit exact F_xL F_LL^-1, factorisé en

      h_vv     = t[12eta+t(3-2t)]/delta,
      h_vtheta = -L z(6eta+t)/delta,
      h_thetav = 6z/(L delta),
      h_tt     = t(12eta+3t-2)/delta.

    Ces polynômes n'exigent ni résolution ni soustraction de flexibilités
    presque égales. Les plans physiques sont (ty,rz) et (tz,-ry). L'axial
    et la torsion sont linéaires en t. Les déplacements de racine et de port
    sont respectivement 0 et I ; ils ne sont pas inclus dans ce retour.
    """
    n, longueur, _, ga, ei, _, _, _ = _parametres_console(metadata)
    t = np.arange(1, n)/n
    z = t*(1-t)
    eta = ei/(ga*longueur**2)
    delta = 1+12*eta
    vv = t*(12*eta+t*(3-2*t))/delta
    vr = -longueur*z*(6*eta+t)/delta
    rv = 6*z/(longueur*delta)
    rr = t*(12*eta+3*t-2)/delta
    psi = np.zeros((n-1, 6, 6))
    psi[:, 0, 0], psi[:, 3, 3] = t, t
    psi[:, 1, 1], psi[:, 5, 5] = vv, rr
    psi[:, 1, 5], psi[:, 5, 1] = vr, rv
    psi[:, 2, 2], psi[:, 4, 4] = vv, rr
    psi[:, 2, 4], psi[:, 4, 2] = -vr, -rv
    return psi.reshape(6*(n-1), 6)
