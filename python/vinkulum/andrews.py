"""LE MÉCANISME D'ANDREWS — un modèle que l'auteur du solveur n'a pas conçu.

Tous les autres bancs de ce noyau ont été écrits par celui qui a écrit le
solveur. Bennett et Kapitza lèvent une part du biais — problèmes posés
ailleurs, réponses publiées — mais la mise en équation restait de moi.

Ici elle ne l'est pas. Le mécanisme d'Andrews (« squeezer ») vient du
*Multibody Systems Handbook* de W. Schiehlen (Springer, 1990) et fait partie
de l'**effort de benchmarking IFToMM**. Ses données, sa topologie et le choix
composante-par-composante de ses contraintes sont lus tels quels dans le jeu
de tests livré avec MBDyn (Masarati & Mantegazza, Politecnico di Milano) :
sept corps, onze liaisons, un ressort, un couple moteur.

Ce sont des DONNÉES de problème, pas du code : la distinction est celle que
le projet applique déjà à sa règle clean-room. Aucune ligne de MBDyn n'est
reprise ; ce sont les paramètres d'un mécanisme publié dans un livre.

CE QUE LE MÉCANISME A DE PARTICULIER, et pourquoi c'est un bon juge :
· il est PLAN monté en 3D, donc un solveur qui contraint naïvement se retrouve
  sur-contraint. Les auteurs choisissent 41 contraintes pour 42 ddl —
  composante par composante — et obtiennent exactement la mobilité 1 ;
· il est RAIDE : le ressort vaut 4530 N/m sur des corps de quelques grammes,
  ce qui en fait un cas d'école des intégrateurs DAE (Hairer & Wanner) ;
· sa réponse est celle que MBDyn calcule sur le MÊME modèle — deux solveurs
  indépendants, une seule définition.
"""
import numpy as np

from . import Noyau

__all__ = ["monte", "demo"]

# ── données du benchmark (Schiehlen 1990 ; jeu de tests MBDyn) ───────────────
O = np.array([0.0, 0.0, 0.0])
A = np.array([-0.06934, -0.00227, 0.0])
B = np.array([-0.03635, 0.03273, 0.0])
C = np.array([0.01400, 0.07200, 0.0])
# configuration assemblée initiale
E = np.array([-2.0960022568369124e-02, 1.2951690850118839e-03, 0.0])
F = np.array([6.9865503092124805e-03, -4.3372200410285764e-04, 0.0])
G = np.array([-3.3997203915751013e-02, 1.6461971731432712e-02, 0.0])
H = np.array([-3.1633134579402290e-02, -1.5618868871665989e-02, 0.0])

K_RESSORT, L0, TAU = 4530.0, 0.07785, 0.033
E_B = np.array([0.0, 0.035, 0.0])
E_D = np.array([0.02, 0.017, 0.0])

# (nom, origine du repère, direction de l'axe x, masse, CdM local, Izz)
LIENS = [
    ("O_F", O, F - O, 0.04325, (0.00092, 0.0), 2.194e-6),
    ("E_F", E, F - E, 0.00365, (0.0165, 0.0), 4.410e-7),
    ("H_E", H, E - H, 0.00706, (0.00579, 0.0), 5.667e-7),
    ("G_E", G, E - G, 0.00706, (0.00579, 0.0), 5.667e-7),
    ("A_G", A, G - A, 0.07050, (0.02308, 0.00916), 1.169e-5),
    ("A_H", A, H - A, 0.05498, (0.01228, -0.00449), 1.912e-5),
    ("E_B_D", E, None, 0.02373, (0.01043, 0.01626), 5.255e-6),   # axe y vers B
]

# (nom, corps a, corps b, point, bloque_t, bloque_r) — LE CHOIX DES AUTEURS,
# composante par composante : c'est lui qui donne mobilité 1 sans redondance
JOINTS = [
    ("O", None, "O_F", O, [0, 1, 2], [0, 1]),
    ("F", "O_F", "E_F", F, [0, 1], []),
    ("E1", "E_B_D", "E_F", E, [0, 1, 2], [0, 1]),
    ("E2", "E_B_D", "G_E", E, [0, 1, 2], [0, 1]),
    ("E3", "E_B_D", "H_E", E, [0, 1, 2], [0, 1]),
    ("G", "A_G", "G_E", G, [0, 1, 2], []),
    ("H", "A_H", "H_E", H, [0, 1, 2], []),
    ("A1", None, "A_G", A, [0, 1], [0, 1]),
    ("A2", None, "A_H", A, [0, 1], [0, 1]),
    ("Bp", None, "E_B_D", B, [0, 1, 2], [0, 1]),
]


def _repere(ex_dir, ey_dir=None):
    """Repère dont l'axe x (ou y) suit une direction, z hors du plan."""
    z = np.array([0.0, 0.0, 1.0])
    if ey_dir is not None:
        y = ey_dir / np.linalg.norm(ey_dir)
        x = np.cross(y, z)
        return np.column_stack([x, y, z])
    x = ex_dir / np.linalg.norm(ex_dir)
    y = np.cross(z, x)
    return np.column_stack([x, y, z])


def monte(couple=TAU):
    """Le mécanisme dans vinkulum. Rend (Noyau, index des corps)."""
    n = Noyau([0.0, 0.0, 0.0])                 # pas de gravité : le benchmark n'en a pas
    idx, poses = {}, {}
    for nom, org, dirx, m, cm, izz in LIENS:
        r = _repere(dirx, ey_dir=(B - E) if nom == "E_B_D" else None)
        cdm = org + r @ np.array([cm[0], cm[1], 0.0])
        # inertie AU CdM, axes corps : le benchmark ne donne que Izz (plan)
        j = [1e-9, 0, 0, 0, 1e-9, 0, 0, 0, izz]
        idx[nom] = n.corps(nom, m, j, list(cdm))
        poses[nom] = (cdm, r)
    t, _, _, v, w, vi = n.etat()
    n.pose_etat(t, [list(poses[k][0]) for k, *_ in [(l[0],) for l in LIENS]],
                [list(poses[k][1].ravel()) for k, *_ in [(l[0],) for l in LIENS]], v, w, vi)
    for nom, a, b, pt, bt, br in JOINTS:
        if a is None:
            pa, ra = list(pt), list(np.eye(3).ravel())
        else:
            cdm, r = poses[a]
            pa, ra = list(r.T @ (np.asarray(pt) - cdm)), list(np.eye(3).ravel())
        n.liaison(nom, None if a is None else idx[a], idx[b],
                  pa=pa, ra=ra, bloque_t=bt, bloque_r=br)
    # RESSORT D–C. Il est COMPRIMÉ dans tout le benchmark : la distance C–D
    # vaut ~0,053 m pour une longueur au repos de 0,078. Un « câble » (qui ne
    # tire que si d > L0) ne ferait donc RIEN — c'est la première version, et
    # c'est le juge externe qui l'a montrée fausse : vinkulum tournait
    # « plausiblement » mais 2,4× trop vite, et rien d'interne ne pouvait le
    # dire. Le bon élément est le contact NON inversé, F = K(L0 − d) quand
    # d < L0, et un assert vérifie qu'on ne sort pas de ce régime.
    cdm, r = poses["E_B_D"]
    d_loc = r.T @ ((E + r @ E_D) - cdm)
    anc = n.corps("ancre_C", 1.0, [1e-9, 0, 0, 0, 1e-9, 0, 0, 0, 1e-9], list(C))
    n.liaison("fixC", None, anc, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    n.contact("ressort", idx["E_B_D"], list(d_loc), L0, b=anc, pb=[0.0] * 3,
              rayon_b=1e-12, k=K_RESSORT, expo=1.0, c=0.0)
    if couple:
        n.couple("moteur", None, idx["O_F"], [0.0, 0.0, 1.0], ("constant", [couple]))
    return n, idx


# Référence MBDyn sur LE MÊME modèle : angle cumulé de la manivelle (deg).
# Obtenue en lançant le jeu de tests livré avec MBDyn — `mbdyn -f
# andrewssqueezer` — et en lisant rz du nœud LINK_O_F. Stockée ici pour que le
# banc tourne sans MBDyn installé ; le script qui la régénère est dans le
# docstring de `reference_mbdyn`.
REF_MBDYN = {0.005: 15.5582, 0.010: 126.9110, 0.020: 472.1024,
             0.035: 1135.3837, 0.050: 1940.9948}


def reference_mbdyn(chemin=None):
    """Relit la référence depuis un `.mov` de MBDyn, si on en a un.

        cp -r <mbdyn>/tests/benchmarks/andrewssqueezer /tmp/and && cd /tmp/and
        sed -i 's/default output: none/default output: all/' andrewssqueezer
        mbdyn -f andrewssqueezer -o run

    Rend {t: angle cumulé en degrés} aux mêmes instants que `REF_MBDYN`.
    """
    import os
    if chemin is None or not os.path.exists(chemin):
        return None
    m = []
    for l in open(chemin):
        c = l.split()
        if c and c[0] == "1600":
            m.append(float(c[5]))
    if not m:
        return None
    th = np.unwrap(np.radians(np.array(m)))
    dt = 0.05 / (len(th) - 1)
    return {t: float(np.degrees(th[min(int(round(t / dt)), len(th) - 1)] - th[0]))
            for t in REF_MBDYN}


def demo(t_end=0.05, h=2e-6, mov=None):
    """Andrews : vinkulum contre MBDyn, sur une définition qu'aucun n'a choisie.

    Le juge n'est ni une formule que j'aurais écrite, ni un raffinement de
    maillage : c'est un SECOND SOLVEUR, sur un modèle publié par des tiers.
    C'est le plus près de « éprouvé par quelqu'un d'autre » qu'on puisse
    atteindre sans quelqu'un d'autre.

    ET IL A TROUVÉ UNE FAUTE. La première version modélisait le ressort D–C
    par un câble ; or il est COMPRIMÉ dans tout le benchmark (d ≈ 0,053 m pour
    L0 = 0,078), donc il ne faisait rien. vinkulum tournait « plausiblement »
    — mais 2,4 fois trop vite, et **rien d'interne ne pouvait le dire** : les
    contraintes étaient tenues à 1e-17, le mécanisme bougeait, la mobilité
    était bonne. Seul un modèle dont la réponse vient d'ailleurs le montre.
    """
    n, idx = monte()
    print("╔═ vinkulum — mécanisme d'ANDREWS (Schiehlen 1990 · benchmark IFToMM)")
    n.assemble()
    print(f"║ 7 corps (42 ddl) · 10 liaisons (41 contraintes) · mobilité 1 · "
          f"|Φ| {max(abs(x) for x in n.phi()):.1e}")
    tr = n.simule(t_end, h, tous=250)
    tv = np.array([f[0] for f in tr])
    th = np.unwrap([np.arctan2(f[2][idx["O_F"]][3], f[2][idx["O_F"]][0]) for f in tr])
    ref = reference_mbdyn(mov) or REF_MBDYN
    ecarts = []
    print("║   t (s)     MBDyn (deg)   vinkulum (deg)   écart")
    for t, a in sorted(ref.items()):
        iv = int(np.argmin(abs(tv - t)))
        b = float(np.degrees(th[iv] - th[0]))
        ecarts.append(abs(a - b))
        print(f"║   {t:.3f}    {a:+11.4f}    {b:+11.4f}    {abs(a - b):8.4f}")
    phif = max(abs(x) for x in n.phi())
    enf = n.contacts()[0][1]
    print(f"║ |Φ| final {phif:.1e} · ressort comprimé de {enf * 1e3:.2f} mm "
          f"(il POUSSE : un câble ne ferait rien)")
    print(f"╚═ Andrews OK — {max(ecarts):.3f}° d'écart après {max(ref.values()):.0f}° parcourus, "
          f"soit {100 * max(ecarts) / max(ref.values()):.4f} %")
    assert phif < 1e-9, ("les contraintes du benchmark ne sont pas tenues", phif)
    assert enf > 0.0, ("le ressort n'est plus comprimé : le modèle sort de son domaine", enf)
    assert max(ecarts) < 0.5, ("vinkulum s'écarte de MBDyn sur le benchmark", ecarts)
    return dict(ecart_max=float(max(ecarts)), phi=float(phif), ref=ref)


if __name__ == "__main__":
    demo()
