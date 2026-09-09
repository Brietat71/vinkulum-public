"""Jantzen, Taira, Granlund & Ol 2014 — plaque plane en RAMPE de tangage 0 → 45°, Re 20 000.

« Vortex dynamics around pitching plates », Phys. Fluids 26 (2014) / AFRL, DTIC
ADA603996. Figure 13, panneau 2D (plaque de paroi à paroi, tunnel hydraulique
AFRL, Re 20 000, courbes en TIRETS) : C_L(t) mesuré pour les rampes C1
(k = π/8, la plaque tourne sur UNE corde de course) et C6 (k = π/48, six
cordes). Pivot au BORD D'ATTAQUE. Loi de rampe (Eldredge, cas canonique
NATO AVT-202) :

    α(t) = (Ω₀/2a) · ln[cosh(a(t − t₁)) / cosh(a(t − t₂))] + α_max/2
    Ω₀ = α_max / t_p,  t_p = t₂ − t₁ = α_max c / (2 U K),  a = 11 (lissage)

PROVENANCE DES CHIFFRES : le PDF est vectoriel ; les deux courbes en tirets
du panneau 2D-C_L ont été extraites du SVG (pdftocairo), calibrées sur les
GRADUATIONS des axes (t 0…8 par pas de 1, C_L 0…7 par pas de 1), puis
rééchantillonnées à Δt 0,05 (642 et 676 points d'origine). Erreur de
lecture : l'épaisseur du trait, ~0,03 en C_L. Le pic NON CIRCULATOIRE de C1
(C_L 6,06 à t 0,85 — la masse ajoutée, ∝ α̈, pendant l'accélération) est un
fait mesuré que le noyau ne modélise PAS (Wagner est circulatoire seul) : il
est exclu du jugement, déclaré. Après la rampe la plaque à 45° LÂCHE des
tourbillons (C_L oscille entre 1,5 et 2,9 sur C1) : on juge des MOYENNES de
fenêtre, pas des instants.

Le Reynolds (2e4) est dans la décade de FRELON ; le Mach est nul (eau).
"""
import numpy as np

RE = 20_000.0
ALPHA_MAX_DEG = 45.0
A_LISSAGE = 11.0
K_RED = {"C1": np.pi / 8, "C6": np.pi / 48}     # K = Ω₀ c / (2 U)

T_CONV = np.array([-0.250, -0.200, -0.150, -0.100, -0.050, -0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.050, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.450, 1.500, 1.550, 1.600, 1.650, 1.700, 1.750, 1.800, 1.850, 1.900, 1.950, 2.000, 2.050, 2.100, 2.150, 2.200, 2.250, 2.300, 2.350, 2.400, 2.450, 2.500, 2.550, 2.600, 2.650, 2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000, 3.050, 3.100, 3.150, 3.200, 3.250, 3.300, 3.350, 3.400, 3.450, 3.500, 3.550, 3.600, 3.650, 3.700, 3.750, 3.800, 3.850, 3.900, 3.950, 4.000, 4.050, 4.100, 4.150, 4.200, 4.250, 4.300, 4.350, 4.400, 4.450, 4.500, 4.550, 4.600, 4.650, 4.700, 4.750, 4.800, 4.850, 4.900, 4.950, 5.000, 5.050, 5.100, 5.150, 5.200, 5.250, 5.300, 5.350, 5.400, 5.450, 5.500, 5.550, 5.600, 5.650, 5.700, 5.750, 5.800, 5.850, 5.900, 5.950, 6.000, 6.050, 6.100, 6.150, 6.200, 6.250, 6.300, 6.350, 6.400, 6.450, 6.500, 6.550, 6.600, 6.650, 6.700, 6.750, 6.800, 6.850, 6.900, 6.950, 7.000, 7.050, 7.100, 7.150, 7.200, 7.250, 7.300, 7.350, 7.400, 7.450, 7.500, 7.550, 7.600, 7.650, 7.700, 7.750, 7.800, 7.850, 7.900, 7.950, 8.000])
CL_C1 = np.array([0.083, 0.103, 0.197, 1.523, 3.766, 5.014, 5.350, 4.636, 3.751, 4.262, 4.481, 4.531, 4.960, 4.896, 5.246, 5.470, 5.515, 5.802, 5.702, 5.939, 5.979, 5.828, 6.045, 4.809, 2.790, 1.641, 1.350, 2.052, 3.079, 2.873, 2.693, 3.077, 2.991, 2.962, 2.979, 2.828, 2.875, 2.849, 2.760, 2.750, 2.745, 2.719, 2.664, 2.629, 2.594, 2.578, 2.518, 2.446, 2.433, 2.372, 2.313, 2.268, 2.198, 2.172, 2.118, 2.035, 1.984, 1.907, 1.837, 1.781, 1.711, 1.680, 1.644, 1.605, 1.595, 1.577, 1.575, 1.565, 1.547, 1.539, 1.525, 1.511, 1.495, 1.475, 1.468, 1.468, 1.468, 1.476, 1.475, 1.473, 1.479, 1.472, 1.472, 1.470, 1.468, 1.468, 1.472, 1.472, 1.468, 1.468, 1.468, 1.476, 1.484, 1.491, 1.502, 1.513, 1.527, 1.540, 1.560, 1.580, 1.604, 1.627, 1.657, 1.687, 1.711, 1.727, 1.743, 1.762, 1.782, 1.802, 1.822, 1.842, 1.858, 1.865, 1.867, 1.867, 1.878, 1.877, 1.867, 1.873, 1.879, 1.885, 1.905, 1.898, 1.901, 1.915, 1.915, 1.925, 1.927, 1.923, 1.927, 1.925, 1.921, 1.915, 1.907, 1.903, 1.895, 1.887, 1.879, 1.869, 1.857, 1.841, 1.826, 1.810, 1.790, 1.770, 1.751, 1.735, 1.719, 1.699, 1.680, 1.658, 1.628, 1.604, 1.580, 1.560, 1.541, 1.514, 1.486, 1.463, 1.437, 1.410, 1.386, 1.366, 1.346, 1.325])
CL_C6 = np.array([0.249, 0.296, 0.351, 0.407, 0.458, 0.510, 0.558, 0.597, 0.629, 0.649, 0.673, 0.697, 0.716, 0.737, 0.764, 0.788, 0.812, 0.848, 0.875, 0.892, 0.932, 0.958, 0.980, 1.022, 1.037, 1.064, 1.116, 1.136, 1.160, 1.199, 1.227, 1.263, 1.291, 1.322, 1.354, 1.381, 1.414, 1.449, 1.473, 1.502, 1.533, 1.560, 1.597, 1.631, 1.647, 1.682, 1.728, 1.750, 1.773, 1.809, 1.837, 1.870, 1.924, 1.947, 1.981, 2.043, 2.070, 2.116, 2.158, 2.191, 2.236, 2.245, 2.296, 2.346, 2.348, 2.394, 2.417, 2.438, 2.504, 2.512, 2.530, 2.569, 2.585, 2.612, 2.631, 2.644, 2.663, 2.660, 2.677, 2.698, 2.693, 2.694, 2.698, 2.691, 2.700, 2.697, 2.669, 2.667, 2.667, 2.647, 2.629, 2.619, 2.604, 2.574, 2.572, 2.552, 2.503, 2.516, 2.469, 2.419, 2.453, 2.374, 2.337, 2.361, 2.240, 2.246, 2.263, 2.151, 2.177, 2.155, 2.081, 2.137, 2.091, 2.032, 2.062, 1.997, 1.963, 1.958, 1.881, 1.852, 1.817, 1.731, 1.703, 1.643, 1.573, 1.556, 1.491, 1.460, 1.462, 1.422, 1.434, 1.438, 1.416, 1.438, 1.434, 1.424, 1.432, 1.418, 1.413, 1.409, 1.394, 1.389, 1.389, 1.381, 1.381, 1.385, 1.381, 1.381, 1.379, 1.377, 1.377, 1.375, 1.373, 1.381, 1.389, 1.397, 1.411, 1.420, 1.435, 1.451, 1.463, 1.475, 1.490, 1.502, 1.515, 1.527])
CL = {"C1": CL_C1, "C6": CL_C6}


def t_pitch(cas):
    """Durée de la rampe en temps convectif : t_p = α_max / (2 K) (c/U = 1)."""
    return np.radians(ALPHA_MAX_DEG) / (2.0 * K_RED[cas])


def alpha_deg(t_conv, cas):
    """Rampe lissée d'Eldredge, t en cordes parcourues, t₁ = 0."""
    tp = t_pitch(cas)
    am = np.radians(ALPHA_MAX_DEG)
    om0 = am / tp
    a = A_LISSAGE
    x = om0 / (2 * a) * (np.logaddexp(a * t_conv, -a * t_conv) - np.logaddexp(a * (t_conv - tp), -a * (t_conv - tp))) + am / 2
    return np.degrees(x)


def moyenne(cas, t0, t1):
    m = (T_CONV >= t0) & (T_CONV <= t1)
    return float(CL[cas][m].mean())


def polaire_plaque_c81():
    """Polaire STATIQUE de plaque plane à Re 2e4, CONSTRUITE (HYPOTHESE, déclarée) :
    attachée 4,6 /rad jusqu'à 8° (Pelletier & Mueller 2000, plaque 1,93 % à Re
    6e4 : pente ~0,08/°, C_L max ~0,75 à 10°), C_L max 0,75 à 10°, puis plaque
    décollée CN = 1,9·sin α (2D, Hoerner ; C_L(45°) = 0,95), raccord de 10 à 16°.
    Le papier ne publie pas de statique ; le niveau post-décrochage est l'entrée
    faible, et `validation` publie la sensibilité ±15 % de CN."""
    al = np.arange(-180.0, 180.01, 1.0)
    cl = np.zeros_like(al)
    for i, a in enumerate(al):
        s = 1.0 if a >= 0 else -1.0
        aa = abs(a)
        if aa <= 8.0:
            v = 4.6 * np.radians(aa)
        elif aa <= 10.0:
            v = 4.6 * np.radians(8.0) + (0.75 - 4.6 * np.radians(8.0)) * (aa - 8.0) / 2.0
        elif aa <= 16.0:
            w = (aa - 10.0) / 6.0
            v = (1 - w) * 0.75 + w * 1.9 * np.sin(np.radians(16)) * np.cos(np.radians(16))
        elif aa <= 90.0:
            v = 1.9 * np.sin(np.radians(aa)) * np.cos(np.radians(aa))
        else:
            v = -1.9 * np.sin(np.radians(180 - aa)) * np.cos(np.radians(180 - aa))
        cl[i] = s * v
    cd = 0.02 + 1.9 * np.sin(np.radians(al)) ** 2 * np.where(np.abs(al) > 10, 1.0, 0.3)
    return list(al), list(cl), list(cd), [0.0] * len(al)


def demo():
    """Contrôles : la rampe atteint 45° en t_p (1 et 6 cordes), ses bornes sont
    0 et 45, la pente médiane vaut Ω₀ ; les données couvrent [−0,25 ; 8] ; le
    pic non circulatoire de C1 est bien là (C_L > 4 avant t 1)."""
    for cas, tp in (("C1", 1.0), ("C6", 6.0)):
        assert abs(t_pitch(cas) - tp) < 1e-12, (cas, t_pitch(cas))
        a = alpha_deg(T_CONV, cas)
        assert a.min() > -0.01 and abs(a.max() - 45.0) < 0.05, (a.min(), a.max())
        assert abs(alpha_deg(np.array([tp / 2]), cas)[0] - 22.5) < 0.01
        pente = (alpha_deg(np.array([tp / 2 + 1e-3]), cas) - alpha_deg(np.array([tp / 2 - 1e-3]), cas))[0] / 2e-3
        assert abs(pente - 45.0 / tp) / (45.0 / tp) < 1e-3, pente
    assert T_CONV[0] == -0.25 and T_CONV[-1] == 8.0 and len(T_CONV) == len(CL_C1) == len(CL_C6)
    assert CL_C1.max() > 4.0 and T_CONV[np.argmax(CL_C1)] < 1.0, "le pic non circulatoire de C1 doit être là"
    assert CL_C6.max() < 3.0 and 3.5 < T_CONV[np.argmax(CL_C6)] < 4.5, "le pic de C6 est pendant la rampe, vers t 4"
    print(f"╚═ Jantzen 2014 OK — C1 pic {CL_C1.max():.2f} à t {T_CONV[np.argmax(CL_C1)]:.2f}, moyenne [2;8] {moyenne('C1', 2, 8):.3f} ; "
          f"C6 pic {CL_C6.max():.2f} à t {T_CONV[np.argmax(CL_C6)]:.2f}, moyenne [3;6] {moyenne('C6', 3, 6):.3f}")


if __name__ == "__main__":
    demo()
