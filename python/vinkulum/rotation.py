"""vinkulum.rotation — le raidissement centrifuge d'une poutre qui tourne.

    python -m vinkulum.rotation        # courbe de Southwell et cas limite

C'est la brique qui manque entre « poutre géométriquement exacte » et « pale
flexible » : une poutre en rotation est BEAUCOUP plus raide qu'à l'arrêt, et
ce surcroît ne vient d'aucune raideur de matériau — il vient de la TENSION
centrifuge, c'est-à-dire de la précontrainte. Un code qui l'oublie sous-estime
les fréquences de pale de dizaines de pour-cent, et c'est précisément ce que
perd une réduction modale faite à l'arrêt (cf. `reduction`, et le manuel
ABAQUS : « a linear perturbation about the state »).

    ω² = ω₀² + S·Ω²        (Southwell)

Le contrôle décisif ne demande AUCUNE donnée extérieure : quand la raideur
élastique tend vers zéro, la poutre devient une chaîne, c'est-à-dire une pale
ARTICULÉE à charnière centrale, dont la fréquence de battement vaut
EXACTEMENT Ω — soit ν = ω/Ω = 1. C'est une identité, pas une corrélation :
pour une charnière centrale, l'inertie et la raideur centrifuge s'annulent
terme à terme.
"""
import sys

import numpy as np

from vinkulum import Noyau


def poutre_tournante(omega, n=6, L=1.0, b=0.02, hh=0.004, rho=2700.0, e_mod=70e9,
                     ei_fac=1.0, amorti=5e-2):
    """Poutre encastrée sur un moyeu en rotation imposée à `omega`.

    `amorti` : amortisseur de relaxation entre moyeu et nœuds — il sert à
    ATTEINDRE le régime établi, pas à modifier la physique.
    """
    a = b * hh
    i2 = b * hh ** 3 / 12.0
    i3 = hh * b ** 3 / 12.0
    jt = i2 + i3
    g = e_mod / (2 * 1.3)
    dl = L / n
    mel = rho * a * dl
    N = Noyau([0.0, 0.0, 0.0])
    moyeu = N.corps("moyeu", 1.0, list(np.diag([1e-3] * 3).ravel()), [0, 0, 0],
                    w=[0, 0, omega])
    N.liaison("arbre", None, moyeu, cible_r=([0, 0, 1.0], ("lineaire", [0.0, omega])))
    idx = []
    for i in range(n + 1):
        mm = mel * (0.5 if i in (0, n) else 1.0)
        j = np.diag([rho * jt * dl, rho * i2 * dl, rho * i3 * dl])
        r = [i * dl, 0.0, 0.0]
        idx.append(N.corps(f"n{i}", max(mm, 1e-9), list(j.ravel()), r,
                           w=[0, 0, omega], v=list(np.cross([0, 0, omega], r))))
    N.liaison("pied", moyeu, idx[0])
    for i in range(n):
        N.poutre(f"b{i}", idx[i], idx[i + 1], e_mod * a, g * a, g * jt,
                 ei_fac * e_mod * i2)
    if amorti:
        for i in idx[1:]:
            N.couple(f"a{i}", moyeu, i, [0.0, 1.0, 0.0], ("ressort", [0.0, amorti, 0.0]))
    return N, idx


def frequence(omega, t_end=1.5, h=2e-4, **kw):
    """Première fréquence de flexion à `omega`, une fois le régime ÉTABLI.

    ⚠ LE RÉGIME ÉTABLI N'EST PAS UN DÉTAIL. La raideur centrifuge est portée
    par la PRÉCONTRAINTE, donc par les tensions. Posée à l'état initial, la
    poutre n'est pas encore étirée : ses tensions oscillent, et la fréquence
    lue vaut 0,70·Ω au lieu de 1,00·Ω — 30 % d'erreur, silencieuse. On relaxe
    avant de lire.
    """
    N, _ = poutre_tournante(omega, **kw)
    if omega:
        N.simule(t_end, h, tous=10 ** 9)
    m = N.modes(3)
    return m[0][0], N


def demo(rapide=False):
    print("╔═ vinkulum — raidissement centrifuge : Southwell, et le cas limite exact")
    omegas = (0.0, 20.0, 50.0, 100.0) if not rapide else (0.0, 50.0, 100.0)

    # 1. LA COURBE : f monte avec Ω, et f² est AFFINE en Ω²
    fs = []
    for om in omegas:
        f, _ = frequence(om)
        fs.append(f)
        nu = f / (om / (2 * np.pi)) if om else float("nan")
        print(f"║ Ω {om:6.1f} rad/s   f {f:8.3f} Hz" + (f"   ν {nu:6.4f}" if om else ""))
    x = np.array([o ** 2 for o in omegas])
    y = np.array([f ** 2 for f in fs])
    a, b = np.polyfit(x, y, 1)
    r2 = 1.0 - np.sum((y - (a * x + b)) ** 2) / np.sum((y - y.mean()) ** 2)
    print(f"║ Southwell : f² = {b:.3f} + {a:.5f}·Ω²   ·   R² {r2:.6f}   "
          f"(f₀ mesurée {fs[0]:.3f} Hz, ajustée {np.sqrt(b):.3f})")
    assert fs[-1] > 3 * fs[0], (fs, "le raidissement centrifuge ne se voit pas")
    assert r2 > 0.999, ("f² n'est pas affine en Ω² : ce n'est pas du Southwell", r2)
    # SOUTHWELL EST UNE APPROXIMATION, pas une identité : le coefficient S
    # dépend faiblement de Ω, donc l'ordonnée à l'origine du fit ne retombe pas
    # exactement sur f₀. Mesuré ici : 9,6 %. On l'asserte large et on PUBLIE
    # l'écart, plutôt que de serrer un seuil sur une loi approchée.
    ec_f0 = abs(np.sqrt(b) / fs[0] - 1)
    print(f"║   l'ordonnée du fit s'écarte de f₀ de {100 * ec_f0:.1f} % — c'est "
          f"l'erreur de l'APPROXIMATION de Southwell, pas du calcul")
    assert ec_f0 < 0.15, (b, fs[0])

    # 2. LE CAS LIMITE, et il ne demande aucune donnée extérieure : raideur
    #    élastique → 0 fait de la poutre une CHAÎNE, c'est-à-dire une pale
    #    articulée à charnière centrale, dont ν vaut EXACTEMENT 1
    om = 100.0
    print("║ cas limite — la poutre devient une chaîne, donc une pale articulée :")
    nus = []
    for ei in (1.0, 1e-2, 1e-3, 1e-5):
        f, _ = frequence(om, ei_fac=ei)
        nu = f / (om / (2 * np.pi))
        nus.append(nu)
        print(f"║   EI ×{ei:<8g} ν {nu:7.4f}" + ("   ← ν = 1 exactement attendu" if ei <= 1e-3 else ""))
    assert abs(nus[-1] - 1.0) < 0.01, (nus, "la chaîne tournante ne rend pas ν = 1")
    assert nus[0] > nus[-1], nus

    # 3. CE QUE LE RÉGIME NON ÉTABLI COÛTE — mesuré, pour que le piège soit su
    f_brut, _ = frequence(om, t_end=0.4, ei_fac=1e-5, amorti=0.0)
    nu_brut = f_brut / (om / (2 * np.pi))
    print(f"║ sans relaxation (précontrainte non établie) : ν {nu_brut:.4f} "
          f"au lieu de {nus[-1]:.4f} — {100 * abs(nu_brut / nus[-1] - 1):.0f} % d'erreur, "
          f"et rien ne le signale")
    assert abs(nu_brut - 1.0) > 0.1, nu_brut
    print("╚═ rotation OK")
    return dict(southwell=(float(b), float(a), float(r2)), nu_limite=float(nus[-1]))


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
