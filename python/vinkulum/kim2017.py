"""Kim, Chang & Sohn 2017 — NACA 0012 en tangage sinusoïdal à Re 2,3e4, quatre fréquences réduites.

« Unsteady aerodynamic characteristics depending on reduced frequency for a
pitching NACA0012 airfoil at Re_c = 2.3×10⁴ », Int'l J. Aeronautical & Space
Sci. 18(1) 8–16 (2017), DOI 10.5139/IJASS.2017.18.1.8 — libre d'accès
(koreascience). Mêmes essais que Kim & Chang 2014 (AST 32) à k 0,1. Soufflerie,
corde 0,18 m, pivot au QUART de corde, α_m = 0°, α_a = 6°, pression instationnaire
intégrée en C_L ; visualisation par fil de fumée.

REPÈRES lus sur la figure 8 (imprimés EN CLAIR sur chaque cadre, pas
digitalisés) : C_L max de la boucle et pente moyenne dC_L/dα (/rad), plus le
SENS de rotation de la boucle lu à la flèche — anti-horaire (retard de phase)
à k 0,1 · 0,2 · 0,4, HORAIRE à k 0,76. Theodorsen (potentiel) place ce
renversement entre 0,1 et 0,2 : le décalage est l'effet visqueux, et c'est ce
qu'un modèle de section à bas Re doit rendre. À k 0,1 le C_L saute de 0,4 à 0,6
entre α 5° et 6° en montée (tourbillons de bord de fuite, bulle) et la
descente porte PLUS que la montée (réattachement retardé) — la boucle
tourne « à l'envers » d'un décrochage dynamique classique.

Le Reynolds (2,3e4) est dans la décade de FRELON — juste en dessous ; Mach nul.
"""
RE = 23_000.0
ALPHA_M_DEG = 0.0
ALPHA_A_DEG = 6.0
PIVOT_X_C = 0.25
FIG8 = {   # k : (C_L max, dC_L/dα /rad, sens : +1 anti-horaire, −1 horaire)
    0.10: (0.65, 4.5, +1),
    0.20: (0.52, 4.1, +1),
    0.40: (0.47, 3.8, +1),
    0.76: (0.41, 3.7, -1),
}
LECTURE = 0.03   # C_L max imprimé à deux décimales ; ±0,03 de jugement


def demo():
    ks = sorted(FIG8)
    assert all(FIG8[a][0] > FIG8[b][0] for a, b in zip(ks, ks[1:])), "C_L max décroît avec k (fig. 9 du papier)"
    assert [FIG8[k][2] for k in ks] == [1, 1, 1, -1], "le sens ne se renverse qu'à k 0,76"
    print("╚═ Kim 2017 OK — C_L max 0,65 → 0,41 de k 0,1 à 0,76, boucle anti-horaire jusqu'à 0,4, horaire à 0,76")


if __name__ == "__main__":
    demo()
