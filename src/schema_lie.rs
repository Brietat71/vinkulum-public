//! Cinématique sigma de generalized-alpha en rotations SPATIALES.
//!
//! Dérivation propre des équations (12)-(16) de Holzinger, Arnold et Gerstmayr,
//! MMT 217 (2025), 106236, DOI 10.1016/j.mechmachtheory.2025.106236.
//! Le papier écrit une composition à droite ; ici R1 = exp(theta) R0 et
//! omega = J_l(theta) theta_dot. Avec b l'incrément Newmark classique :
//!
//! theta = b + s (J_l(theta)^-1 - I) omega1,  s = sigma h beta/gamma.
//!
//! On résout cette relation, plutôt que les mises à jour approximées du
//! correcteur de l'article. La correction delta = theta-b est l'inconnue
//! locale pour ne pas la perdre par soustraction de deux grands termes.
//! La tangente implicite conserve la dérivée du Jacobien inverse.

use crate::ad::{self, Dual, Scalar};
use crate::{M3, V3};

/// (J_l(theta)^-1-I) w, sans soustraire deux matrices proches.
fn inverse_moins_ident<T: Scalar>(theta: [T; 3], w: [T; 3]) -> [T; 3] {
    let q = ad::v_dot(theta, theta);
    let c = if q.val() < 1e-3 {
        // (1 - (t/2) cot(t/2))/t^2, série paire en t.
        T::from_f64(1.0 / 12.0)
            + q * (T::from_f64(1.0 / 720.0)
                + q * (T::from_f64(1.0 / 30240.0) + q.scale(1.0 / 1209600.0)))
    } else {
        let t = q.sqrt().scale(0.5);
        (T::one() - t * t.cos() / t.sin()) / q
    };
    let croix = ad::v_cross(theta, w);
    ad::v_add(
        ad::v_scale(T::from_f64(-0.5), croix),
        ad::v_scale(c, ad::v_cross(theta, croix)),
    )
}

fn correction(theta: V3, w: V3) -> V3 {
    V3::from(inverse_moins_ident::<f64>(theta.into(), w.into()))
}

fn derivee_correction(theta: V3, w: V3) -> M3 {
    type D = Dual<f64, 3>;
    let r = inverse_moins_ident(
        std::array::from_fn(|i| D::var(theta[i], i)),
        std::array::from_fn(|i| D::cst(w[i])),
    );
    M3::from_fn(|i, j| r[i].d[j])
}

fn dans_carte(theta: &V3) -> bool {
    theta.iter().all(|x| x.is_finite()) && theta.norm() < std::f64::consts::PI
}

/// Résolution locale de la relation cinématique dans la carte ||theta|| < pi.
/// Le résidu est recalculé ; aucun succès sur le seul petit incrément Newton.
pub(crate) fn rotation(b: V3, w: V3, s: f64) -> Result<V3, String> {
    if !b.iter().chain(w.iter()).all(|x| x.is_finite()) || !s.is_finite() {
        return Err("cinématique sigma : données non finies".into());
    }
    if s == 0.0 {
        return Ok(b);
    }
    if !dans_carte(&b) {
        return Err("cinématique sigma : incrément initial hors de la carte ||theta|| < pi".into());
    }
    let mut delta = V3::zeros();
    for _ in 0..12 {
        let theta = b + delta;
        let r = delta - s * correction(theta, w);
        let echelle = (delta.norm() + s.abs() * w.norm() * theta.norm()).max(1e-300);
        if !r.iter().all(|x| x.is_finite()) || !echelle.is_finite() {
            return Err("cinématique sigma : résidu non fini".into());
        }
        if r.norm() <= 16.0 * f64::EPSILON * echelle {
            return Ok(theta);
        }
        let a = M3::identity() - s * derivee_correction(theta, w);
        let d = a
            .lu()
            .solve(&r)
            .ok_or("cinématique sigma : tangente singulière")?;
        let mut accepte = false;
        let mut alpha = 1.0;
        for _ in 0..10 {
            let candidat = delta - alpha * d;
            let angle = b + candidat;
            if dans_carte(&angle) {
                let essai = candidat - s * correction(angle, w);
                if essai.norm() < r.norm() {
                    delta = candidat;
                    accepte = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if !accepte {
            return Err("cinématique sigma : correction locale sans descente".into());
        }
    }
    Err("cinématique sigma : plafond d'itérations locales atteint".into())
}

/// Dérivées spatiales de exp(theta) :
/// K = J_l A^-1 B pour l'accélération (avant beta h² c),
/// A = I-s d[(J_l^-1-I)omega]/dtheta, B = I+sigma(J_l^-1-I).
/// Pour un incrément indépendant de b, comme GGL : K B^-1.
pub(crate) fn tangentes(theta: V3, w: V3, s: f64, sigma: f64) -> Result<(M3, M3), String> {
    let jl = crate::tangent::jac_gauche(&theta);
    if sigma == 0.0 {
        return Ok((jl, M3::identity()));
    }
    if !dans_carte(&theta) || !w.iter().all(|v| v.is_finite()) || !s.is_finite() {
        return Err("cinématique sigma : tangente hors domaine".into());
    }
    let di = M3::from_columns(&[
        correction(theta, V3::x()),
        correction(theta, V3::y()),
        correction(theta, V3::z()),
    ]);
    let a = M3::identity() - s * derivee_correction(theta, w);
    let b = M3::identity() + sigma * di;
    let k = jl
        * a.lu()
            .solve(&b)
            .ok_or("cinématique sigma : tangente implicite singulière")?;
    let bi = b
        .lu()
        .solve(&M3::identity())
        .ok_or("cinématique sigma : transfert GGL singulier")?;
    if !k.iter().chain(bi.iter()).all(|v| v.is_finite()) {
        return Err("cinématique sigma : dérivées non finies".into());
    }
    Ok((k, bi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{expm, log_so3};

    #[test]
    fn inverse_geometrique_et_limite_au_repos() {
        for t in [0., 1e-10, 1e-5, 0.031, 0.032, 0.7, 2.9] {
            let theta = V3::new(0.3, -0.5, 0.8).normalize() * t;
            let jl = crate::tangent::jac_gauche(&theta);
            for w in [V3::x(), V3::new(-0.7, 0.9, 1.3)] {
                assert!((jl * (w + correction(theta, w)) - w).norm() < 2e-11);
            }
        }
        let w = V3::new(0.3, -0.5, 0.8);
        assert_eq!(rotation(V3::zeros(), w, 0.1).unwrap(), V3::zeros());
        for s in [0., 0.1, 1.0] {
            let b = V3::z() * 0.7;
            assert_eq!(rotation(b, 4.0 * V3::z(), s).unwrap(), b);
        }
    }

    #[test]
    fn coefficient_geometrique_et_signe_spatial() {
        // Vitesse spatiale affine. Magnus donne -h³/12 (w0 x a0).
        // La correction sigma=gamma/(3 beta) doit rendre ce même terme.
        let w0 = V3::new(1.1, -0.4, 0.7);
        let a = V3::new(-0.3, 0.8, 0.2);
        let cible = -w0.cross(&a) / 12.0;
        let mut preced = f64::INFINITY;
        for h in [0.04, 0.02, 0.01, 0.005] {
            let b = h * w0 + 0.5 * h * h * a;
            let th = rotation(b, w0 + h * a, h / 3.0).unwrap();
            let e = ((th - b) / h.powi(3) - cible).norm();
            assert!(e < preced * 0.6, "h={h}, erreur={e}, précédente={preced}");
            preced = e;
        }
        assert!(preced < 3e-4);
    }

    #[test]
    fn tangente_implicite_et_correction_independante() {
        let (h, beta, gamma, sigma) = (0.13, 0.28, 0.56, 0.7);
        let b = V3::new(0.21, -0.37, 0.18);
        let w = V3::new(0.7, 1.1, -0.4);
        let s = sigma * h * beta / gamma;
        let th = rotation(b, w, s).unwrap();
        let r = expm(&th);
        let (k, bi) = tangentes(th, w, s, sigma).unwrap();
        for independant in [false, true] {
            let exact = if independant { k * bi } else { k };
            for i in 0..3 {
                let mut d = V3::zeros();
                d[i] = 1e-6;
                let dw = if independant {
                    V3::zeros()
                } else {
                    gamma / (beta * h) * d
                };
                let p = expm(&rotation(b + d, w + dw, s).unwrap());
                let m = expm(&rotation(b - d, w - dw, s).unwrap());
                let df = (log_so3(&(p * r.transpose())) - log_so3(&(m * r.transpose()))) / 2e-6;
                assert!((df - exact.column(i)).norm() < 5e-9);
            }
        }
        // Supprimer la dérivée du Jacobien inverse ne donne PAS notre tangente.
        assert!((k - M3::identity()).norm() > 1e-3);
    }

    #[test]
    fn covariance_et_refus() {
        let b = V3::new(0.7, -0.9, 1.3);
        let w = V3::new(-0.2, 0.5, 1.1);
        let q = expm(&V3::new(-0.5, 0.4, 0.2));
        let th = rotation(b, w, 0.2).unwrap();
        assert!((rotation(q * b, q * w, 0.2).unwrap() - q * th).norm() < 1e-14);
        let (k, bi) = tangentes(th, w, 0.2, 0.8).unwrap();
        let (kq, biq) = tangentes(q * th, q * w, 0.2, 0.8).unwrap();
        assert!((kq - q * k * q.transpose()).norm() < 2e-14);
        assert!((biq - q * bi * q.transpose()).norm() < 2e-14);
        assert!(rotation(V3::x() * 4.0, w, 0.2).is_err());
        assert!(rotation(b, w, f64::NAN).is_err());
        assert!(rotation(b, V3::repeat(f64::INFINITY), 0.2).is_err());
    }
}
