"""Sous-structures linéaires reproductibles, en coordonnées et masses physiques.

Les matrices CSC suivent ``(K - omega**2 M) u = f``, omega en rad/s. Les
indices d'interface portent des déplacements physiques, sans base modale ni
normalisation préalable par la masse. L'intérieur précède l'interface dans
les deux familles ; le réducteur ne doit pas dépendre de cet ordre particulier.

Ce fichier construit les données et un oracle analytique de chaîne. Il ne
chronomètre pas les solveurs et ne dépend d'aucun prototype modal.
"""
from dataclasses import dataclass
import operator

import numpy as np
from scipy.sparse import csc_matrix, diags


@dataclass(frozen=True)
class ModelePorts:
    k: csc_matrix
    m: csc_matrix
    interieur: np.ndarray
    interface: np.ndarray
    metadata: dict


def _nombre(n):
    if isinstance(n, (bool, np.bool_)):
        raise ValueError("n doit être un entier supérieur ou égal à 2")
    try:
        n = operator.index(n)
    except TypeError as exc:
        raise ValueError("n doit être un entier supérieur ou égal à 2") from exc
    if n < 2:
        raise ValueError("n doit être un entier supérieur ou égal à 2")
    return n


def _positif(nom, valeur):
    valeur = float(valeur)
    if not np.isfinite(valeur) or valeur <= 0:
        raise ValueError(f"{nom} doit être fini et strictement positif")
    return valeur


def _parametres_chaine(n, raideur, masse):
    return _nombre(n), _positif("raideur", raideur), _positif("masse", masse)


def chaine(n, raideur=1.0, masse=1.0):
    """n masses égales, n ressorts égaux, racine fixe et dernière masse au port.

    Pour u_0=0, l'énergie potentielle vaut k/2 somme (u_i-u_{i-1})²
    et l'énergie cinétique m/2 somme v_i², i=1,...,n. Le dernier
    ressort se termine sur la masse n : aucun ressort ne la relie au sol.
    Lorsque le port est bloqué, les n-1 masses intérieures ont des conditions
    de Dirichlet aux deux bouts et mu_min=4k/m sin²(pi/(2n)).
    """
    n, raideur, masse = _parametres_chaine(n, raideur, masse)
    diagonale = np.full(n, 2.0 * raideur)
    diagonale[-1] = raideur
    k = diags((-np.full(n-1, raideur), diagonale,
               -np.full(n-1, raideur)), (-1, 0, 1), format="csc")
    m = diags(np.full(n, masse), format="csc")
    mu_min = 4.0 * (raideur / masse) * np.sin(np.pi / (2*n))**2
    return ModelePorts(k, m, np.arange(n-1), np.array([n-1]), {
        "famille": "chaine", "n": n, "ddl": n, "ports": [n-1],
        "raideur_N_par_m": raideur, "masse_nodale_kg": masse,
        "masse_mobile_kg": n * masse,
        "mu_min_interieur_s_moins_2": float(mu_min),
        "omega_min_interieur_rad_s": float(np.sqrt(mu_min)),
        "ddl_unites": "m", "force_interface_unite": "N",
        "convention": "(K - omega^2 M) u = f; u_0 = 0",
    })


def reference_chaine(n, omega, raideur=1.0, masse=1.0, force=1.0):
    """Schur et réponse harmonique de la chaîne, sans résolution matricielle.

    Dans la bande 0 <= omega*sqrt(m/k) <= 2, les déplacements vérifient
    u_i/u_n = sin(i theta)/sin(n theta), cos(theta)=1-m omega²/(2k).
    Le Schur est k[sin((n+1)theta)-sin(n theta)]/sin(n theta). Les
    formes en sinc évitent la perte de précision à omega=0 et au bord de
    bande. Au-delà, des rapports d'exponentielles décroissantes remplacent
    les sinus hyperboliques, sans débordement exponentiel.

    L'oracle renvoie aussi le champ complet sous une force appliquée au port
    seulement. Son coût est O(n), dont O(1) pour le Schur lui-même. Il refuse
    les pôles intérieurs et résonances globales numériquement indiscernables
    à la précision machine ; aucune borne d'erreur certifiée n'est revendiquée.
    """
    n, raideur, masse = _parametres_chaine(n, raideur, masse)
    omega, force = float(omega), float(force)
    if not np.isfinite(omega) or omega < 0 or not np.isfinite(force):
        raise ValueError("omega doit être fini et positif ou nul, force finie")
    rho = omega * np.sqrt(masse / raideur)
    if not np.isfinite(rho) or rho > np.sqrt(np.finfo(float).max):
        raise ValueError("fréquence adimensionnée hors domaine flottant")
    i = np.arange(1, n+1)
    eps = np.finfo(float).eps
    if rho <= 2:
        basse = rho <= np.sqrt(2.0)
        angle = 2*np.arcsin(rho/2) if basse else 2*np.arccos(rho/2)
        denominateur = np.sinc(n*angle/np.pi)
        if n*angle > 1e-3 and abs(np.sin(n*angle)) <= 32*eps*max(1., n*angle):
            raise ValueError("fréquence au voisinage numérique d'un pôle intérieur")
        rapports = (i/n) * np.sinc(i*angle/np.pi) / denominateur
        if basse:
            schur = (raideur/n * np.cos((n+0.5)*angle)
                     * np.sinc(angle/(2*np.pi)) / denominateur)
        else:
            rapports *= np.where((n-i) % 2, -1.0, 1.0)
            schur = (-raideur*(2*n+1)/n * np.sinc((n+0.5)*angle/np.pi)
                     * np.cos(angle/2) / denominateur)
    else:
        alpha = 2*np.arccosh(rho/2)
        rapports = (np.exp((i-n)*alpha) * np.expm1(-2*i*alpha)
                    / np.expm1(-2*n*alpha))
        rapports *= np.where((n-i) % 2, -1.0, 1.0)
        schur = raideur * (1-rho*rho-rapports[-2])
    # La force nulle ne rend pas une résonance exploitable comme oracle :
    # le problème homogène y est lui aussi non unique.
    echelle = raideur * max(1.0/n, rho*rho)
    if abs(schur) <= 32*eps*echelle:
        raise ValueError("fréquence au voisinage numérique d'une résonance globale")
    port = force/schur
    return {
        "schur": float(schur),
        "deplacement_interface": float(port),
        "deplacement": port * rapports,
        "rapports_deplacement": rapports,
        "omega_rad_s": omega,
        "mu_min_interieur_s_moins_2": float(
            4*(raideur/masse)*np.sin(np.pi/(2*n))**2),
    }


def console(n, longueur=1.0, rayon=0.005, module_young=70e9,
            masse_volumique=2700.0, poisson=0.3):
    """Console native droite, section circulaire, 6 DDL physiques au dernier nœud.

    La masse et l'inertie nodales sont obtenues par lumping de la longueur :
    poids 1/2 aux deux extrémités et 1 ailleurs. La section donne A=pi r²,
    I=pi r⁴/4 et J=2I. Les rigidités sont EA, GA, GJ, EI ; GA utilise un
    coefficient de cisaillement égal à 1, explicitement fixé pour ce modèle.
    Il ne s'agit pas d'une référence analytique de poutre continue.

    La formulation ``integree`` reproduit la flexibilité nodale linéaire de
    Timoshenko pour les rigidités constantes ci-dessus. Le noyau assemble K
    et M à la pose droite non précontrainte. Les six
    coordonnées du nœud racine sont supprimées explicitement. Les indices
    restants sont (tx,ty,tz,rx,ry,rz) par corps, translations en m et petites
    rotations en rad. Les efforts conjugués sont des N et des N.m.
    """
    from vinkulum import Noyau

    n = _nombre(n)
    longueur = _positif("longueur", longueur)
    rayon = _positif("rayon", rayon)
    module_young = _positif("module_young", module_young)
    masse_volumique = _positif("masse_volumique", masse_volumique)
    poisson = float(poisson)
    if not np.isfinite(poisson) or not -1 < poisson < 0.5:
        raise ValueError("poisson doit appartenir à ]-1, 0.5[")
    aire = np.pi*rayon**2
    moment = np.pi*rayon**4/4
    polaire = 2*moment
    cisaillement = module_young/(2*(1+poisson))
    pas = longueur/n
    noyau = Noyau([0., 0., 0.])
    corps = []
    for j in range(n+1):
        poids = 0.5 if j in (0, n) else 1.0
        masse = masse_volumique*aire*pas*poids
        inertie = masse_volumique*pas*poids*np.diag([polaire, moment, moment])
        corps.append(noyau.corps(f"section{j}", masse, inertie.ravel().tolist(),
                                 [j*pas, 0., 0.]))
    noyau.liaison("encastrement", None, corps[0])
    for j in range(n):
        noyau.poutre(f"poutre{j}", corps[j], corps[j+1],
                     module_young*aire, cisaillement*aire,
                     cisaillement*polaire, module_young*moment,
                     formulation="integree")
    k, c, m, g = [csc_matrix((data, indices, pointers), shape=(nr, nc))
                  for nr, nc, pointers, indices, data in noyau.k_c_m_g_creux()]
    racine = np.arange(6*corps[0], 6*corps[0]+6)
    libres = np.setdiff1d(np.arange(k.shape[0]), racine)
    contraintes_restantes = g[:, libres].copy()
    contraintes_restantes.eliminate_zeros()
    if contraintes_restantes.nnz or g.shape[0] != 6:
        raise ValueError("l'encastrement ne se réduit pas aux six DDL de racine")
    if c.nnz and np.any(c.data):
        raise ValueError("la console de référence doit être sans amortissement")
    k, m = k[libres, :][:, libres].tocsc(), m[libres, :][:, libres].tocsc()
    port_physique = np.arange(6*corps[-1], 6*corps[-1]+6)
    interface = np.searchsorted(libres, port_physique)
    interieur = np.setdiff1d(np.arange(k.shape[0]), interface)
    return ModelePorts(k, m, interieur, interface, {
        "famille": "console", "n": n, "ddl": int(k.shape[0]),
        "ports": interface.tolist(), "ports_avant_encastrement": port_physique.tolist(),
        "racine_supprimee": racine.tolist(), "longueur_m": longueur,
        "rayon_m": rayon, "module_young_Pa": module_young,
        "masse_volumique_kg_par_m3": masse_volumique, "poisson": poisson,
        "masse_totale_kg": masse_volumique*aire*longueur,
        "masse_mobile_kg": masse_volumique*aire*pas*(n-0.5),
        "poids_extremites_masse_et_inertie": 0.5,
        "coefficient_cisaillement": 1.0, "formulation_poutre": "integree",
        "EA_N": module_young*aire, "GA_N": cisaillement*aire,
        "EI_N_m2": module_young*moment, "GJ_N_m2": cisaillement*polaire,
        "ordre_ddl_par_corps": ["tx", "ty", "tz", "rx", "ry", "rz"],
        "ddl_unites_par_corps": ["m", "m", "m", "rad", "rad", "rad"],
        "convention": "(K - omega^2 M) u = f; racine éliminée en coordonnées physiques",
        "assemblage": "Noyau.k_c_m_g_creux",
    })


def flexibilites_conditionnees_console(metadata):
    """Diagonale analytique de K_II^-1, par nœud intérieur et DDL physique.

    La référence est le noyau de Green statique de la poutre de Timoshenko,
    conditionné aux déplacements et rotations nuls aux deux extrémités.
    Avec t=x/L, z=t(1-t), eta=EI/(GA L²), les deux diagonales de flexion
    se factorisent en

        Fvv = (L³/EI) z [eta + z(eta+z/3)/(1+12 eta)],
        Frr = (L/EI) z [1 - 3z/(1+12 eta)].

    L'axial vaut Lz/EA et la torsion Lz/GJ. Ces expressions évitent la
    soustraction Fxx-FxL FLL^-1 FLx près du bout. La matrice renvoyée a
    forme (n-1, 6), dans l'ordre tx,ty,tz,rx,ry,rz. Cette identité concerne
    la console droite, non précontrainte, intégrée et homogène de ce fichier.
    """
    if metadata.get("famille") != "console" or metadata.get("formulation_poutre") != "integree":
        raise ValueError("référence réservée à la console homogène intégrée")
    n = _nombre(metadata["n"])
    longueur = _positif("longueur_m", metadata["longueur_m"])
    ea, ga, ei, gj = [_positif(nom, metadata[nom])
                      for nom in ("EA_N", "GA_N", "EI_N_m2", "GJ_N_m2")]
    t = np.arange(1, n)/n
    z = t*(1-t)
    eta = ei/(ga*longueur**2)
    v = longueur**3/ei*z*(eta+z*(eta+z/3)/(1+12*eta))
    r = longueur/ei*z*(1-3*z/(1+12*eta))
    return np.column_stack((longueur/ea*z, v, v, longueur/gj*z, r, r))


def borne_lambda_trace_console(metadata):
    """Borne inférieure de lambda_min(K_II,M_II), en s^-2, calculée en O(n).

    En arithmétique exacte, 1/trace(M_II K_II^-1) <= lambda_min car les
    valeurs propres de M_II^(1/2) K_II^-1 M_II^(1/2) sont positives et
    leur maximum est inférieur à leur somme. Les diagonales analytiques de
    flexibilité et les masses physiques suffisent : aucun eigensolve.

    Ce calcul flottant n'utilise pas d'arrondis dirigés. Sa valeur n'est donc
    pas une borne certifiée de la matrice effectivement assemblée en machine.
    La borne vise uniquement les paramètres de la fixture, sans modification.
    """
    souplesse = flexibilites_conditionnees_console(metadata)
    n, longueur = metadata["n"], metadata["longueur_m"]
    rayon, rho = metadata["rayon_m"], metadata["masse_volumique_kg_par_m3"]
    aire, moment = np.pi*rayon**2, np.pi*rayon**4/4
    masses = rho*longueur/n*np.array([aire, aire, aire, 2*moment, moment, moment])
    trace = float(np.sum(souplesse*masses))
    if not np.isfinite(trace) or trace <= 0:
        raise ValueError("trace de flexibilité massique hors domaine flottant")
    return 1.0/trace


def metrique_ports(modele):
    """Raideur statique analytique au port, dans les unités physiques de la fixture.

    Pour la chaîne, la raideur équivalente des n ressorts en série vaut k/n.
    Pour la console, on inverse analytiquement sa flexibilité au bout L,
    sans former la différence K_SS-K_SI K_II^-1 K_IS. Chaque plan possède
    la flexibilité [[L/(GA)+L³/(3EI), L²/(2EI)], [L²/(2EI), L/(EI)]].
    Les couples conjugués sont (ty,rz) et (tz,-ry), ce qui fixe les signes
    des couplages. L'axial vaut EA/L et la torsion GJ/L.

    La matrice est définie positive ; v.T @ metrique @ v mesure une énergie
    (deux fois l'énergie élastique statique), même si v mélange m et rad.
    Elle dépend des paramètres non modifiés de la fixture. Une permutation
    des seuls indices d'interface est respectée ; une transformation des
    coordonnées physiques exige de transformer aussi cette métrique.
    """
    metadata = modele.metadata
    if metadata.get("famille") == "chaine":
        n = _nombre(metadata["n"])
        ressort = _positif("raideur_N_par_m", metadata["raideur_N_par_m"])
        metrique = np.array([[ressort/n]])
    elif metadata.get("famille") == "console" and metadata.get("formulation_poutre") == "integree":
        longueur = _positif("longueur_m", metadata["longueur_m"])
        ea, ga, ei, gj = [_positif(nom, metadata[nom])
                          for nom in ("EA_N", "GA_N", "EI_N_m2", "GJ_N_m2")]
        eta = ei/(ga*longueur**2)
        denominateur = 1+12*eta
        plan = np.array([[12*ei/(longueur**3*denominateur),
                          -6*ei/(longueur**2*denominateur)],
                         [-6*ei/(longueur**2*denominateur),
                          ei/longueur*(1+3/denominateur)]])
        metrique = np.zeros((6, 6))
        metrique[0, 0], metrique[3, 3] = ea/longueur, gj/longueur
        metrique[np.ix_([1, 5], [1, 5])] = plan
        signes = np.array([1., -1.])
        metrique[np.ix_([2, 4], [2, 4])] = signes[:, None]*plan*signes[None, :]
    else:
        raise ValueError("métrique réservée aux deux fixtures analytiques de ce module")
    origine = list(metadata["ports"])
    interface = list(modele.interface)
    if len(interface) != len(origine) or set(interface) != set(origine):
        raise ValueError("indices d'interface différents des ports physiques de la fixture")
    ordre = [origine.index(j) for j in interface]
    return metrique[np.ix_(ordre, ordre)]
