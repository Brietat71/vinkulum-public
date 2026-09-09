"""Peters–He : l'inflow dynamique à ÉTATS FINIS d'un disque rotor (6–7 septembre).

Le niveau qui manquait entre Pitt–Peters (3 états : λ₀, λ1c, λ1s — dans le
noyau) et le sillage libre (`vinkulum.sillage`) : le potentiel d'accélération
d'un disque, développé en harmoniques d'azimut et en polynômes radiaux
(He 1989, Peters & He 1995). Le champ induit au disque s'écrit

    w(r̄, ψ) = Σ_r Σ_j φ_j^r(r̄) [α_j^r cos rψ + β_j^r sin rψ],   j = r+1, r+3, …

et les états obéissent à

    [K]{α̇} + V [L̃^c]⁻¹{α} = ½{τ^c}      (idem sinus, L̃^s)

K = diag(2H_n^m/π) (masse apparente), L̃ en forme fermée (fonction Γ, X = tan(χ/2),
χ l'angle d'obliquité du sillage), τ les coefficients de la charge de pale
projetée sur les mêmes fonctions. Temps en Ω t, vitesses en ΩR.

SOURCES ET CE QU'ELLES ONT TRANCHÉ. Le rapport de S. Li (UC Davis 2020, transmis
par Paul) donne φ, H, L̃, Γ, V, τ (ses éq. 28–41) ; Ferreira et al. 2021 (TEMA
22(2), libre) donne le « 2 » du Γ impair, K = 2H/π, et le ½ devant τ. Sans ce
½, la convergence à N états rend exactement 2× le momentum — mesuré ici avant
de le trouver dans le texte (le contrôle `demo` le garde : 8 états m = 0,
λ̄ = C_T/(2V) à 0,6 %). Trois formes fermées SORTENT du modèle et sont
assertées : un état = 9/8 du momentum (la troncature, pas une erreur) ;
trois états, charge uniforme, k_x = (2π/3)·X — ni Coleman (X) ni Pitt–Peters
(15π/32·X) : le modèle a sa propre valeur, publiée ; X → 0 rend le champ
axisymétrique.

NON COUVERT, déclaré : le couplage AU NOYAU (imposer un champ w(r̄, ψ) à une
pale) n'existe pas encore — `Inflow` n'accepte qu'un profil radial et une
harmonique 1/rev ; ici le modèle est jugé avec sa propre théorie des tranches
(rotor rigide, pas de battement) contre l'inflow mesuré d'Elliott. Les V de
He (V_T pour m = 0, V pour m > 0, avec λ_m = √3·α_1^0) sont ceux du rapport.
"""
from math import pi, sqrt, tan, atan2, sin, cos, radians
import numpy as np


def _dfact(n):
    return 1.0 if n <= 0 else n * _dfact(n - 2)


def H(n, m):
    return _dfact(n + m - 1) * _dfact(n - m - 1) / (_dfact(n + m) * _dfact(n - m))


def phi(j, r, rb):
    """Fonction de forme radiale φ_j^r(r̄) (He, éq. 28 du rapport)."""
    rb = np.asarray(rb, float)
    s = np.zeros_like(rb)
    for q in range(r, j, 2):
        s += rb ** q * (-1) ** ((q - r) // 2) * _dfact(j + q) / (_dfact(q - r) * _dfact(q + r) * _dfact(j - q - 1))
    return sqrt((2 * j + 1) * H(j, r)) * s


def gamma(r, j, m, n):
    if (r + m) % 2 == 0:
        return ((-1) ** ((n + j - 2 * r) // 2) * 2 * sqrt((2 * n + 1) * (2 * j + 1))
                / (sqrt(H(n, m) * H(j, r)) * (j + n) * (j + n + 2) * ((j - n) ** 2 - 1)))
    if abs(j - n) == 1:
        return pi * np.sign(r - m) / (2 * sqrt(H(n, m) * H(j, r) * (2 * n + 1) * (2 * j + 1)))
    return 0.0


def etats(M, Q):
    """Liste (m, n) par la méthode des tables : n = m+1, m+3, … ≤ Q+1."""
    return [(m, n) for m in range(M + 1) for n in range(m + 1, Q + 2, 2)]


class PetersHe:
    def __init__(self, M=4, Q=6):
        self.idx = etats(M, Q)                      # états cosinus (m ≥ 0)
        self.idx_s = [(m, n) for m, n in self.idx if m > 0]   # états sinus
        self.a = np.zeros(len(self.idx))
        self.b = np.zeros(len(self.idx_s))
        self.K = np.array([2 * H(n, m) / pi for m, n in self.idx])
        self.Ks = np.array([2 * H(n, m) / pi for m, n in self.idx_s])

    def matrices(self, X):
        Lc = np.zeros((len(self.idx), len(self.idx)))
        for i, (r, j) in enumerate(self.idx):
            for k, (m, n) in enumerate(self.idx):
                g = gamma(r, j, m, n)
                l = min(r, m)
                Lc[i, k] = X ** m * g if r == 0 else (X ** abs(m - r) + (-1) ** l * X ** (m + r)) * g
        Ls = np.zeros((len(self.idx_s), len(self.idx_s)))
        for i, (r, j) in enumerate(self.idx_s):
            for k, (m, n) in enumerate(self.idx_s):
                l = min(r, m)
                Ls[i, k] = (X ** abs(m - r) - (-1) ** l * X ** (m + r)) * gamma(r, j, m, n)
        return Lc, Ls

    def w(self, rb, psi):
        """Inflow adimensionnel (>0 vers le bas) aux points (r̄, ψ)."""
        rb, psi = np.asarray(rb, float), np.asarray(psi, float)
        out = np.zeros(np.broadcast(rb, psi).shape)
        for k, (m, n) in enumerate(self.idx):
            out = out + phi(n, m, rb) * self.a[k] * np.cos(m * psi)
        for k, (m, n) in enumerate(self.idx_s):
            out = out + phi(n, m, rb) * self.b[k] * np.sin(m * psi)
        return out

    def tau(self, rb, psi, dL, rho, omega, R):
        """Coefficients de charge (éq. 41) : dL = force normale par élément (N,
        déjà intégrée sur Δr), aux (r̄, ψ) donnés (tous les éléments de toutes les
        pales d'un tour, en moyenne d'azimut : on divise par le nombre de pas
        d'azimut n_psi représentés)."""
        rb, psi, dL = map(lambda x: np.asarray(x, float).ravel(), (rb, psi, dL))
        n_psi = len(np.unique(np.round(psi, 9)))
        s = rho * omega ** 2 * R ** 4
        tc = np.array([np.sum(dL * phi(n, m, rb) * np.cos(m * psi)) / (s * (2 * pi if m == 0 else pi)) / n_psi
                       for m, n in self.idx])
        ts = np.array([np.sum(dL * phi(n, m, rb) * np.sin(m * psi)) / (s * pi) / n_psi for m, n in self.idx_s])
        return tc, ts

    def vitesses(self, mu, lam_f):
        """(V_T, V) de He : λ = λ_f + λ_m, λ_m = √3·α_1^0."""
        lam_m = sqrt(3) * self.a[0]
        lam = lam_f + lam_m
        # plancher : en stationnaire pur, V_T = λ et λ = 0 au départ — le point
        # fixe λ ← c/λ oscille sans sous-relaxation, d'où relax = 0,5 par défaut
        vt = max(sqrt(mu * mu + lam * lam), 1e-3)
        v = max((mu * mu + (lam + lam_m) * lam) / vt, 1e-3)
        return vt, v, lam

    def stationnaire(self, tc, ts, mu, lam_f, relax=0.5):
        """α = L̃ τ /(2Ṽ) — Ṽ dépend de λ_m, donc de α : point fixe (sous-relaxé)."""
        if not np.any(self.a):
            self.a[0] = sqrt(max(abs(tc[0]), 1e-9)) / sqrt(3)     # ordre de grandeur du momentum
        for _ in range(200):
            vt, v, lam = self.vitesses(mu, lam_f)
            X = tan(0.5 * atan2(mu, max(abs(lam), 1e-9)))
            Lc, Ls = self.matrices(X)
            vv = np.array([vt if m == 0 else v for m, n in self.idx])
            vs = np.array([v for m, n in self.idx_s])
            a_new = Lc @ (0.5 * tc) / vv
            b_new = Ls @ (0.5 * ts) / vs if len(ts) else self.b
            da = np.max(np.abs(a_new - self.a)) if self.a.size else 0.0
            self.a = self.a + relax * (a_new - self.a)
            self.b = self.b + relax * (b_new - self.b)
            if da < 1e-10:
                break
        return self

    def avance(self, tc, ts, mu, lam_f, dpsi):
        """Un pas en temps Ω·t : Euler implicite sur [K]α̇ + Ṽ L̃⁻¹ α = ½τ."""
        vt, v, lam = self.vitesses(mu, lam_f)
        X = tan(0.5 * atan2(mu, max(abs(lam), 1e-9)))
        Lc, Ls = self.matrices(X)
        vv = np.array([vt if m == 0 else v for m, n in self.idx])
        A = np.diag(self.K / dpsi) + np.diag(vv) @ np.linalg.inv(Lc)
        self.a = np.linalg.solve(A, self.K * self.a / dpsi + 0.5 * tc)
        if len(ts):
            vs = np.array([v for m, n in self.idx_s])
            B = np.diag(self.Ks / dpsi) + np.diag(vs) @ np.linalg.inv(Ls)
            self.b = np.linalg.solve(B, self.Ks * self.b / dpsi + 0.5 * ts)
        return self


def bet_charges(ph, R, n_b, omega, corde, mu, alpha_arbre, theta0, theta1c, theta1s, vrillage,
                r_pied=0.0, a0=2 * pi, n_r=24, n_psi=36, rho=1.225):
    """Théorie des tranches d'un rotor RIGIDE (pas de battement) : rend (r̄, ψ, dL)
    pour toutes les stations d'un tour, avec l'inflow courant de `ph`. θ = θ₀ +
    θ_tw(r̄ − 0,75) + θ1c cos ψ + θ1s sin ψ ; U_T = Ω r + μΩR sin ψ ;
    U_P = ΩR(λ_f + w) ; α = θ − atan(U_P/U_T) ; dL = ½ρ U² c a₀ α Δr."""
    rb = (np.arange(n_r) + 0.5) / n_r * (1 - r_pied) + r_pied
    dr = (1 - r_pied) / n_r * R
    psi = (np.arange(n_psi) + 0.5) / n_psi * 2 * pi
    RB, PSI = np.meshgrid(rb, psi, indexing="ij")
    lam_f = mu * tan(alpha_arbre)
    theta = theta0 + vrillage * (RB - 0.75) + theta1c * np.cos(PSI) + theta1s * np.sin(PSI)
    ut = omega * R * (RB + mu * np.sin(PSI))
    up = omega * R * (lam_f + ph.w(RB, PSI))
    al = theta - np.arctan2(up, ut)
    dL = 0.5 * rho * (ut * ut + up * up) * corde * a0 * al * dr
    return RB, PSI, dL, lam_f


def trim_ct(ph, ct_cible, R, n_b, omega, corde, mu, alpha_arbre, theta1c, theta1s, vrillage, r_pied=0.0,
            rho=1.225, iters=40, sur=0.5):
    """Point fixe charge ↔ inflow (sous-relaxation SUR, éq. 49 du rapport) avec θ₀
    ajusté à chaque passage pour rendre C_T : rend (θ₀, C_T obtenu)."""
    theta0 = 0.1
    for _ in range(iters):
        for _k in range(3):
            RB, PSI, dL, lam_f = bet_charges(ph, R, n_b, omega, corde, mu, alpha_arbre, theta0, theta1c, theta1s,
                                              vrillage, r_pied=r_pied, rho=rho)
            T = n_b * dL.mean(axis=1).sum()            # moyenne d'azimut, somme radiale
            ct = T / (rho * pi * R ** 2 * (omega * R) ** 2)
            theta0 += 0.5 * (ct_cible - ct) / max(ct, 1e-6) * theta0
        tc, ts = ph.tau(np.repeat(RB, n_b, axis=0), np.repeat(PSI, n_b, axis=0), np.repeat(dL, n_b, axis=0),
                        rho, omega, R)
        a_old, b_old = ph.a.copy(), ph.b.copy()
        ph.stationnaire(tc, ts, mu, lam_f)
        ph.a = a_old + sur * (ph.a - a_old)
        ph.b = b_old + sur * (ph.b - b_old)
    return theta0, ct


def demo():
    rb = np.linspace(0, 1, 2001)
    # 1. un état = 9/8 du momentum ; huit états m = 0 → C_T/(2V) — le ½ de τ, mesuré
    for nq, att, tol in ((0, 9 / 8 * 0.5, 1e-9), (14, 0.5, 0.01)):
        ph = PetersHe(M=0, Q=nq)
        tc = np.array([float(np.trapezoid(phi(n, 0, rb) * rb, rb)) for m, n in ph.idx])   # P = C_T = 1
        ph.a = ph.matrices(0.0)[0] @ (0.5 * tc)                                            # V = 1
        lam_bar = 2 * float(np.trapezoid(ph.w(rb, 0.0) * rb, rb))
        assert abs(lam_bar - att) < tol, (nq, lam_bar, att)
    print(f"║ m = 0 : 1 état λ̄ = 9/8·C_T/2V (troncature) · 8 états λ̄ = {lam_bar:.4f}·C_T/V → momentum C_T/(2V) à {100 * abs(lam_bar / 0.5 - 1):.1f} %")
    # 2. trois états, charge uniforme, X = 0,3 : k_x = (2π/3)·X — la forme fermée du modèle
    ph = PetersHe(M=1, Q=1)
    Lc, _ = ph.matrices(0.3)
    tc = np.array([float(np.trapezoid(phi(1, 0, rb) * rb, rb)), 0.0])
    ph.a = Lc @ (0.5 * tc)
    lam0 = sqrt(3) * ph.a[0]
    lam1c = phi(2, 1, 1.0) * ph.a[1]
    kx = lam1c / lam0 / 0.3
    assert abs(kx - 2 * pi / 3) < 1e-6, kx
    print(f"║ 3 états, charge uniforme : k_x/X = {kx:.4f} = 2π/3 (Pitt–Peters 15π/32 = {15 * pi / 32:.4f}, Coleman 1)")
    # 3. X → 0 : champ axisymétrique ; X > 0 : plus d'inflow à l'arrière (ψ = 0)
    ph = PetersHe(M=3, Q=5)
    tc = np.array([float(np.trapezoid(phi(n, m, rb) * rb, rb)) if m == 0 else 0.0 for m, n in ph.idx])
    ph.a = ph.matrices(0.0)[0] @ (0.5 * tc)
    w0 = ph.w(0.7, np.linspace(0, 2 * pi, 8))
    assert np.ptp(w0) < 1e-12, "X = 0 doit être axisymétrique"
    ph.a = ph.matrices(0.5)[0] @ (0.5 * tc)
    assert ph.w(0.7, 0.0) > ph.w(0.7, pi), "à X > 0 l'arrière du disque doit être plus chargé en inflow"
    # 4. marche en temps : le stationnaire est la limite, τ_1^0 = K/(V·L⁻¹) ≈ 0,85/V (8/(3π)·… l'échelle de Pitt–Peters)
    ph = PetersHe(M=0, Q=0)
    tc = np.array([float(np.trapezoid(phi(1, 0, rb) * rb, rb))])
    ph.a[:] = 0.0
    hist = []
    for _ in range(400):
        ph.avance(tc, np.zeros(0), 0.0, 0.0, 0.02)          # V_T = λ (hover), non linéaire
        hist.append(sqrt(3) * ph.a[0])
    ph2 = PetersHe(M=0, Q=0)
    ph2.stationnaire(tc, np.zeros(0), 0.0, 0.0)
    assert abs(hist[-1] - sqrt(3) * ph2.a[0]) < 1e-3, (hist[-1], sqrt(3) * ph2.a[0])
    print(f"║ marche en temps : λ(8 rad) {hist[-1]:.4f} = stationnaire {sqrt(3) * ph2.a[0]:.4f} — masse apparente K = 2H/π")
    print("╚═ Peters–He OK")


if __name__ == "__main__":
    demo()
