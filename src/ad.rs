//! ad — différentiation automatique en mode AVANT par nombres duaux, à N
//! tangentes, EMBOÎTABLE (Dual<Dual<f64, N>, N> donne les dérivées secondes).
//!
//! Clean-room : c'est l'arithmétique des nombres duaux (Clifford 1873), rien
//! d'autre. Écrite ici plutôt qu'importée parce que (1) le noyau doit pouvoir
//! emboîter deux niveaux pour la matrice tangente des efforts de liaison
//! (Hessienne de λ·Φ), (2) `std::autodiff`/Enzyme est nightly, (3) 150 lignes
//! qu'on possède valent mieux qu'une API qu'on subit.
//!
//! `Scalar` est le trait que les fonctions d'ÉLÉMENT (Φ, efforts, inertie)
//! sont écrites contre : f64 le satisfait (production), Dual aussi (dérivées).
//! Une fonction écrite une fois rend sa valeur, son gradient exact et sa
//! Hessienne exacte, à l'arrondi près — plus de différences finies.

use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait Scalar:
    Copy
    + Clone
    + std::fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + PartialOrd
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(x: f64) -> Self;
    fn sqrt(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn atan2(self, x: Self) -> Self;
    fn tanh(self) -> Self;
    /// puissance réelle, dérivable en 0⁺ (plancher sur la valeur)
    fn powf(self, e: f64) -> Self;
    /// logarithme népérien, planché pour rester défini
    fn ln(self) -> Self;
    /// max/min par BRANCHE (la dérivée suit la branche choisie — saturations)
    fn max(self, o: Self) -> Self;
    fn min(self, o: Self) -> Self;
    /// la partie réelle (valeur), tous niveaux d'emboîtement retirés
    fn val(self) -> f64;
    /// produit par une CONSTANTE : sur un dual, N produits au lieu de 2N + 1
    fn scale(self, s: f64) -> Self;
    /// profondeur d'emboîtement : 0 pour f64, 1 pour Dual<f64, N>, 2 pour Dual<Dual<…>>
    const ORDRE: usize;
}

impl Scalar for f64 {
    const ORDRE: usize = 0;
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_f64(x: f64) -> Self {
        x
    }
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    fn sin(self) -> Self {
        f64::sin(self)
    }
    fn cos(self) -> Self {
        f64::cos(self)
    }
    fn atan2(self, x: Self) -> Self {
        f64::atan2(self, x)
    }
    fn tanh(self) -> Self {
        f64::tanh(self)
    }
    fn powf(self, e: f64) -> Self {
        f64::powf(self.max(0.0), e)
    }
    fn ln(self) -> Self {
        f64::ln(self.max(1e-300))
    }
    fn max(self, o: Self) -> Self {
        if self >= o {
            self
        } else {
            o
        }
    }
    fn min(self, o: Self) -> Self {
        if self <= o {
            self
        } else {
            o
        }
    }
    fn val(self) -> f64 {
        self
    }
    fn scale(self, s: f64) -> Self {
        self * s
    }
}

/// Nombre dual à N tangentes : valeur `v`, dérivées partielles `d[k]`.
#[derive(Copy, Clone, Debug)]
pub struct Dual<T: Scalar, const N: usize> {
    pub v: T,
    pub d: [T; N],
}

impl<T: Scalar, const N: usize> Dual<T, N> {
    pub fn cst(v: T) -> Self {
        Dual {
            v,
            d: [T::zero(); N],
        }
    }
    /// variable k : valeur v, tangente unité sur la composante k
    pub fn var(v: T, k: usize) -> Self {
        let mut d = [T::zero(); N];
        d[k] = T::one();
        Dual { v, d }
    }
    fn map2(self, o: Self, fv: T, fa: T, fb: T) -> Self {
        // d(f(a, b)) = fa·da + fb·db
        let mut d = [T::zero(); N];
        for (k, dk) in d.iter_mut().enumerate() {
            *dk = fa * self.d[k] + fb * o.d[k];
        }
        Dual { v: fv, d }
    }
    fn map1(self, fv: T, fa: T) -> Self {
        let mut d = [T::zero(); N];
        for (k, dk) in d.iter_mut().enumerate() {
            *dk = fa * self.d[k];
        }
        Dual { v: fv, d }
    }
}

impl<T: Scalar, const N: usize> Add for Dual<T, N> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        self.map2(o, self.v + o.v, T::one(), T::one())
    }
}
impl<T: Scalar, const N: usize> Sub for Dual<T, N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self.map2(o, self.v - o.v, T::one(), -T::one())
    }
}
impl<T: Scalar, const N: usize> Mul for Dual<T, N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.map2(o, self.v * o.v, o.v, self.v)
    }
}
impl<T: Scalar, const N: usize> Div for Dual<T, N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = T::one() / o.v;
        self.map2(o, self.v * inv, inv, -self.v * inv * inv)
    }
}
impl<T: Scalar, const N: usize> Neg for Dual<T, N> {
    type Output = Self;
    fn neg(self) -> Self {
        self.map1(-self.v, -T::one())
    }
}
impl<T: Scalar, const N: usize> PartialEq for Dual<T, N> {
    fn eq(&self, o: &Self) -> bool {
        self.v == o.v
    }
}
impl<T: Scalar, const N: usize> PartialOrd for Dual<T, N> {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        self.v.partial_cmp(&o.v)
    }
}

impl<T: Scalar, const N: usize> Scalar for Dual<T, N> {
    const ORDRE: usize = T::ORDRE + 1;
    fn zero() -> Self {
        Self::cst(T::zero())
    }
    fn one() -> Self {
        Self::cst(T::one())
    }
    fn from_f64(x: f64) -> Self {
        Self::cst(T::from_f64(x))
    }
    fn scale(self, s: f64) -> Self {
        self.map1(self.v.scale(s), T::from_f64(s))
    }
    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        self.map1(s, T::from_f64(0.5) / s)
    }
    fn sin(self) -> Self {
        self.map1(self.v.sin(), self.v.cos())
    }
    fn cos(self) -> Self {
        self.map1(self.v.cos(), -self.v.sin())
    }
    fn atan2(self, x: Self) -> Self {
        // d atan2(y, x) = (x·dy − y·dx) / (x² + y²)
        let den = x.v * x.v + self.v * self.v;
        self.map2(x, self.v.atan2(x.v), x.v / den, -self.v / den)
    }
    fn tanh(self) -> Self {
        let t = self.v.tanh();
        self.map1(t, T::one() - t * t)
    }
    /// d(x^e) = e·x^(e−1)·dx — plancher sur la valeur pour que x = 0 ne rende
    /// pas une dérivée infinie quand e < 1 (le contact de Hertz a e = 1,5,
    /// donc la dérivée y est nulle, mais l'exposant est un paramètre)
    fn powf(self, e: f64) -> Self {
        let x = self.v.max(T::from_f64(1e-300));
        self.map1(x.powf(e), T::from_f64(e) * x.powf(e - 1.0))
    }
    /// d(ln x) = dx/x
    fn ln(self) -> Self {
        let x = self.v.max(T::from_f64(1e-300));
        self.map1(x.ln(), T::one() / x)
    }
    fn max(self, o: Self) -> Self {
        if self.v >= o.v {
            self
        } else {
            o
        }
    }
    fn min(self, o: Self) -> Self {
        if self.v <= o.v {
            self
        } else {
            o
        }
    }
    fn val(self) -> f64 {
        self.v.val()
    }
}

// ── petite algèbre générique : vecteurs et matrices 3×3 sur tout Scalar ─────
pub type V<T> = [T; 3];
pub type M<T> = [[T; 3]; 3];

pub fn v_add<T: Scalar>(a: V<T>, b: V<T>) -> V<T> {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn v_sub<T: Scalar>(a: V<T>, b: V<T>) -> V<T> {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn v_scale<T: Scalar>(s: T, a: V<T>) -> V<T> {
    [s * a[0], s * a[1], s * a[2]]
}
pub fn v_dot<T: Scalar>(a: V<T>, b: V<T>) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn v_cross<T: Scalar>(a: V<T>, b: V<T>) -> V<T> {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn v_norm<T: Scalar>(a: V<T>) -> T {
    v_dot(a, a).sqrt()
}
pub fn v_from<T: Scalar>(a: [f64; 3]) -> V<T> {
    [T::from_f64(a[0]), T::from_f64(a[1]), T::from_f64(a[2])]
}
pub fn m_from<T: Scalar>(m: &nalgebra::Matrix3<f64>) -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = T::from_f64(m[(i, j)]);
        }
    }
    o
}
pub fn m_ident<T: Scalar>() -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for (i, oi) in o.iter_mut().enumerate() {
        oi[i] = T::one();
    }
    o
}
pub fn m_mul<T: Scalar>(a: &M<T>, b: &M<T>) -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                o[i][j] = o[i][j] + a[i][k] * b[k][j];
            }
        }
    }
    o
}
/// a · B avec B CONSTANTE (repère matériel, rotation de référence) : 27
/// produits scalaires-par-constante au lieu de 27 produits de duaux
pub fn m_mul_cst<T: Scalar>(a: &M<T>, b: &nalgebra::Matrix3<f64>) -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                o[i][j] = o[i][j] + a[i][k].scale(b[(k, j)]);
            }
        }
    }
    o
}
pub fn m_t<T: Scalar>(a: &M<T>) -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = a[j][i];
        }
    }
    o
}
pub fn m_v<T: Scalar>(a: &M<T>, x: V<T>) -> V<T> {
    [v_dot(a[0], x), v_dot(a[1], x), v_dot(a[2], x)]
}
pub fn m_tv<T: Scalar>(a: &M<T>, x: V<T>) -> V<T> {
    [
        a[0][0] * x[0] + a[1][0] * x[1] + a[2][0] * x[2],
        a[0][1] * x[0] + a[1][1] * x[1] + a[2][1] * x[2],
        a[0][2] * x[0] + a[1][2] * x[1] + a[2][2] * x[2],
    ]
}
pub fn m_trace<T: Scalar>(a: &M<T>) -> T {
    a[0][0] + a[1][1] + a[2][2]
}
pub fn skew<T: Scalar>(w: V<T>) -> M<T> {
    let z = T::zero();
    [[z, -w[2], w[1]], [w[2], z, -w[0]], [-w[1], w[0], z]]
}
pub fn vee<T: Scalar>(s: &M<T>) -> V<T> {
    [s[2][1], s[0][2], s[1][0]]
}

/// exp([w]×) — Rodrigues, avec la série en dessous de 1e-4 pour que les
/// dérivées (duales) restent exactes à w → 0 (sin t / t n'est pas dérivable
/// numériquement à 0 sans ça).
pub fn expm<T: Scalar>(w: V<T>) -> M<T> {
    let k = skew(w);
    // POSE NON PERTURBÉE (duaux du premier ordre à valeur NULLE — le cas de
    // toutes les tangentes du noyau, `pose_pert`) : exp([w]×) = I + [w]×
    // EXACTEMENT, parce que [w]×² a valeur 0 et tangente 0·e + e·0 = 0. La
    // série rendait le même résultat en payant la racine, k² (27 produits
    // de duaux) et 18 produits de plus ; mesuré 14 % du pas sur Princeton.
    // Pas pour les duaux EMBOÎTÉS (ordre 2) : le terme mixte ½s[u, δ] y vit.
    if T::ORDRE == 1 && w.iter().all(|c| c.val() == 0.0) {
        let mut o = m_ident();
        for i in 0..3 {
            for j in 0..3 {
                o[i][j] = o[i][j] + k[i][j];
            }
        }
        return o;
    }
    let t2 = v_dot(w, w);
    let k2 = m_mul(&k, &k);
    let (a, b) = if t2.val() < 1e-8 {
        // sin t/t = 1 − t²/6 + t⁴/120 ; (1−cos t)/t² = 1/2 − t²/24 + t⁴/720
        (
            T::one() - t2 / T::from_f64(6.0) + t2 * t2 / T::from_f64(120.0),
            T::from_f64(0.5) - t2 / T::from_f64(24.0) + t2 * t2 / T::from_f64(720.0),
        )
    } else {
        let t = t2.sqrt();
        (t.sin() / t, (T::one() - t.cos()) / t2)
    };
    let mut o = m_ident();
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = o[i][j] + a * k[i][j] + b * k2[i][j];
        }
    }
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_et_hessienne_exacts() {
        // f(x, y) = x·sin(y) + √x / y ; gradient et Hessienne analytiques
        type D1 = Dual<f64, 2>;
        type D2 = Dual<Dual<f64, 2>, 2>;
        let (x, y) = (1.7_f64, 0.6_f64);
        let f = |x: D2, y: D2| x * y.sin() + x.sqrt() / y;
        let xd = D2::var(D1::var(x, 0), 0);
        let yd = D2::var(D1::var(y, 1), 1);
        let r = f(xd, yd);
        let fx = y.sin() + 0.5 / (x.sqrt() * y);
        let fy = x * y.cos() - x.sqrt() / (y * y);
        let fxx = -0.25 / (x.powf(1.5) * y);
        let fxy = y.cos() - 0.5 / (x.sqrt() * y * y);
        let fyy = -x * y.sin() + 2.0 * x.sqrt() / (y * y * y);
        assert!((r.v.v - (x * y.sin() + x.sqrt() / y)).abs() < 1e-15);
        assert!((r.v.d[0] - fx).abs() < 1e-14 && (r.v.d[1] - fy).abs() < 1e-14);
        assert!(
            (r.d[0].d[0] - fxx).abs() < 1e-14,
            "fxx {} {}",
            r.d[0].d[0],
            fxx
        );
        assert!((r.d[0].d[1] - fxy).abs() < 1e-14 && (r.d[1].d[0] - fxy).abs() < 1e-14);
        assert!((r.d[1].d[1] - fyy).abs() < 1e-14);
    }

    #[test]
    fn expm_derivee_a_zero_et_loin() {
        // d/dε exp([ε e_k]×) à ε = 0 vaut [e_k]× ; et exp reste orthonormale
        type D = Dual<f64, 3>;
        let w = [D::var(0.0, 0), D::var(0.0, 1), D::var(0.0, 2)];
        let r = expm(w);
        assert!((r[0][1].d[2] + 1.0).abs() < 1e-15 && (r[1][0].d[2] - 1.0).abs() < 1e-15);
        let w2 = [D::var(0.3, 0), D::var(-0.2, 1), D::var(0.5, 2)];
        let r2 = expm(w2);
        let rt = m_t(&r2);
        let i = m_mul(&rt, &r2);
        for a in 0..3 {
            for b in 0..3 {
                let e = if a == b { 1.0 } else { 0.0 };
                assert!((i[a][b].v - e).abs() < 1e-14);
            }
        }
    }
}

#[cfg(test)]
mod tests_emboite {
    use super::*;

    /// Deux niveaux de duaux à perturbations SUCCESSIVES sur SO(3) :
    /// U(R) = tr(AᵀR) avec R = exp(ε)·exp(δ)·R₀. On veut g_i = ∂U/∂ε_i (le
    /// gradient au point perturbé) puis K_ij = ∂g_i/∂δ_j. Contrôle : DF du
    /// gradient. C'est EXACTEMENT le patron de `Poutre::tangente`, isolé.
    #[test]
    fn double_niveau_perturbations_successives() {
        type D1 = Dual<f64, 3>;
        type D2 = Dual<Dual<f64, 3>, 3>;
        let r0 = nalgebra::Matrix3::new(
            0.936, -0.289, 0.199, 0.302, 0.951, -0.062, -0.178, 0.118, 0.977,
        );
        let a = nalgebra::Matrix3::new(0.3, -1.1, 0.7, 2.0, 0.5, -0.9, 1.3, 0.2, -1.7);
        // gradient seul (un niveau), pour la DF
        let grad = |rot: &nalgebra::Matrix3<f64>| -> [f64; 3] {
            type D = Dual<f64, 3>;
            let e = expm([D::var(0.0, 0), D::var(0.0, 1), D::var(0.0, 2)]);
            let r = m_mul(&e, &m_from::<D>(rot));
            let at = m_from::<D>(&a);
            let mut u = D::zero();
            for i in 0..3 {
                for j in 0..3 {
                    u = u + at[i][j] * r[i][j];
                }
            }
            [u.d[0], u.d[1], u.d[2]]
        };
        // double niveau : ext = δ (extérieur), int = ε (intérieur, appliqué à GAUCHE)
        let ext = |k: usize| -> D2 {
            let mut d = D2::cst(D1::cst(0.0));
            d.d[k] = D1::cst(1.0);
            d
        };
        let int = |k: usize| -> D2 { D2::cst(D1::var(0.0, k)) };
        let e_ext = expm([ext(0), ext(1), ext(2)]);
        let e_int = expm([int(0), int(1), int(2)]);
        let rot = m_mul(&e_int, &m_mul(&e_ext, &m_from::<D2>(&r0)));
        let at = m_from::<D2>(&a);
        let mut u = D2::zero();
        for i in 0..3 {
            for j in 0..3 {
                u = u + at[i][j] * rot[i][j];
            }
        }
        let eps = 1e-6;
        for j in 0..3 {
            let mut dp = [0.0; 3];
            dp[j] = eps;
            let ep = expm([dp[0], dp[1], dp[2]]);
            let em = expm([-dp[0], -dp[1], -dp[2]]);
            let (gp, gm) = (grad(&(mat(&ep) * r0)), grad(&(mat(&em) * r0)));
            for i in 0..3 {
                let fd = (gp[i] - gm[i]) / (2.0 * eps);
                let ad_ = u.d[j].d[i];
                assert!(
                    (ad_ - fd).abs() < 1e-6 * (1.0 + fd.abs()),
                    "K[{i},{j}] : AD {ad_} vs DF {fd}"
                );
            }
        }
    }

    fn mat(m: &M<f64>) -> nalgebra::Matrix3<f64> {
        nalgebra::Matrix3::new(
            m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2],
        )
    }
}
