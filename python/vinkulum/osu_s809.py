"""vinkulum.osu_s809 — le profil S809 en tangage oscillant, MESURÉ (OSU / NREL 1995).

Source : R. Reuss Ramsay, M. J. Hoffman, G. M. Gregorek, *Effects of Grit
Roughness and Pitch Oscillations on the S809 Airfoil*, NREL/TP-442-7817
(Ohio State University, 1995 ; OSTI 205563). Le rapport est un SCAN sans
couche texte : la table statique B3 est TRANSCRITE de l'image de la page
B-10 (valeur par valeur), et les repères des boucles dynamiques sont LUS sur
la figure C25 (page C-13) — un tracé à grille 0,5 en CL et 5° en α, donc
±0,05 en CL et ±1° en α, déclarés. Aucune courbe n'est « digitalisée » :
seuls les extrêmes des boucles, lisibles sans ambiguïté, sont pris.

Statique (table B3) : S809 propre, Re 1,25·10⁶ ; colonnes AOA (deg), Cl,
Cdp (traînée de pression), Cm¼.
Dynamique (fig. C25) : S809 propre, Re 1,01·10⁶, tangage sinusoïdal ±5,5°
à f = 1,85 Hz, k = ωc/(2V) = 0,077, trois angles moyens 8°, 14°, 20°.
"""
import numpy as np

# table B3 — (AOA, Cl, Cdp, Cm)
STATIQUE_B3 = np.array([
    [-20.2, -0.56, 0.3016, 0.0599], [-18.1, -0.60, 0.2783, 0.0747], [-16.0, -0.81, 0.1250, -0.0042],
    [-14.2, -0.76, 0.0675, -0.0058], [-12.0, -0.70, 0.0555, -0.0009], [-10.1, -0.62, 0.0371, -0.0066],
    [-8.1, -0.57, 0.0293, 0.0008], [-6.1, -0.59, 0.0142, -0.0021], [-4.1, -0.39, 0.0011, -0.0251],
    [-2.1, -0.15, 0.0011, -0.0301], [0.0, 0.09, 0.0012, -0.0354], [2.1, 0.33, 0.0030, -0.0423],
    [4.1, 0.59, 0.0037, -0.0490], [6.1, 0.80, 0.0033, -0.0495], [8.1, 0.90, 0.0112, -0.0387],
    [10.2, 0.93, 0.0203, -0.0362], [11.2, 0.94, 0.0256, -0.0310], [12.2, 1.00, 0.0423, -0.0380],
    [13.2, 1.02, 0.0544, -0.0392], [14.2, 1.04, 0.0613, -0.0379], [15.2, 1.06, 0.0750, -0.0419],
    [16.2, 1.02, 0.0894, -0.0435], [17.2, 0.96, 0.1060, -0.0491], [18.2, 0.91, 0.1373, -0.0619],
    [19.0, 0.72, 0.3302, -0.1395], [20.1, 0.70, 0.3394, -0.1317], [22.1, 0.70, 0.3725, -0.1274],
    [24.0, 0.76, 0.4299, -0.1416], [26.0, 0.87, 0.5161, -0.1688],
])

K_RED = 0.077
AMPLITUDE_DEG = 5.5
# repères lus sur la figure C25 : angle moyen → CL max de la boucle, α du max,
# CL MIN de la boucle sur la descente (pour 8° c'est le retour à α_min = 2,5°,
# le premier repère transcrit — 0,95 vers 10° — n'était pas un minimum)
DYNAMIQUE_C25 = {
    8.0: dict(cl_max=1.22, alpha_max=13.0, cl_min=0.50),   # min de boucle = retour à 2,5°
    14.0: dict(cl_max=1.20, alpha_max=15.0, cl_min=0.65),
    20.0: dict(cl_max=1.63, alpha_max=24.0, cl_min=0.65),
}


def polaire_c81():
    """La table B3 prolongée en plateau vers ±180° (le noyau exige la plage)."""
    a, cl, cd, cm = STATIQUE_B3.T
    cd = cd + 0.008                       # traînée de frottement absente de Cdp : forfait, déclaré
    a_ext = np.concatenate([[-180, -90, -40], a, [40, 90, 180]])
    cl_ext = np.concatenate([[0, 0, -0.6], cl, [0.9, 0, 0]])
    cd_ext = np.concatenate([[0.02, 1.6, 1.0], cd, [1.0, 1.6, 0.02]])
    cm_ext = np.concatenate([[0, 0, 0.1], cm, [-0.2, 0, 0]])
    return list(a_ext), list(cl_ext), list(cd_ext), list(cm_ext)
