//! Projection sur SO(3) : itération polaire près du groupe, SVD ailleurs.
//! X⁺ = X(I − E/2), E = XᵀX − I, donne E⁺ = −3E²/4 + E³/4.
//! Les produits sont de taille 3 fixe, sans allocation sur le chemin proche.
use crate::{DMatrix, M3};

fn proche(a: &M3) -> Option<M3> {
    let mut e = a.transpose() * a - M3::identity();
    if !(e.norm_squared() <= 0.125_f64.powi(2)) || !(a.determinant() > 0.) {
        return None;
    }
    let mut x = *a;
    // Avec ||E||_F ≤ 1/8, quatre étapes suffisent en arithmétique exacte
    // pour passer sous epsilon. Le contrôle porte sur le résultat calculé.
    for _ in 0..5 {
        x -= (x * e) * 0.5;
        e = x.transpose() * x - M3::identity();
        if e.norm_squared() <= (8. * f64::EPSILON).powi(2) {
            return Some(x);
        }
    }
    None
}

pub(crate) fn projette(a: &M3) -> Result<M3, String> {
    crate::fini(a.as_slice(), "orientation à projeter")?;
    if let Some(x) = proche(a) {
        return Ok(x);
    }
    let s = crate::svd_sure(DMatrix::from_column_slice(3, 3, a.as_slice()))?;
    let mut u = s.u.unwrap();
    let vt = s.v_t.unwrap();
    // Procrustes propre : corriger la direction de plus petite valeur
    // singulière si le facteur orthogonal serait une réflexion.
    if (&u * &vt).determinant() < 0. {
        u.column_mut(2).neg_mut();
    }
    let x = M3::from_column_slice((u * vt).as_slice());
    if !(x.determinant() > 0.) || (x.transpose() * x - M3::identity()).norm() > 64. * f64::EPSILON {
        return Err("projection d'orientation : facteur hors SO(3)".into());
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{expm, V3};

    #[test]
    fn facteur_connu_et_changements_de_reperes() {
        let q = expm(&V3::new(0.43, -1.11, 2.07));
        let p = expm(&V3::new(-0.83, 0.21, -0.62));
        let gauche = expm(&V3::new(0.51, 0.62, -1.3));
        let droite = expm(&V3::new(1.82, -0.47, 0.15));
        for delta in [0., 1e-15, 1e-9, -0.02, 0.02, 0.2, 0.49] {
            let d = M3::from_diagonal(&V3::new(1. + delta, 1. - 2. * delta, 1. + 0.5 * delta));
            // A=QH avec H symétrique définie positive : le facteur est Q.
            let a = q * p * d * p.transpose();
            let x = projette(&a).unwrap();
            assert!((x - q).norm() < 5e-15, "delta {delta}");
            assert!((x.transpose() * x - M3::identity()).norm() < 4e-15);
            assert!(
                (projette(&(gauche * a * droite)).unwrap() - gauche * x * droite).norm() < 6e-15
            );
        }
        let a = q * (M3::identity() + 1e-9 * p.column(0) * p.column(0).transpose());
        assert!(proche(&a).is_some());
        let plan =
            expm(&V3::new(0., 0., 0.83)) * M3::from_diagonal(&V3::new(1. + 1e-9, 1. - 2e-9, 1.));
        let x = projette(&plan).unwrap();
        // Une symétrie plane exacte ne doit pas créer de couplage parasite.
        assert_eq!([x[(0, 2)], x[(1, 2)], x[(2, 0)], x[(2, 1)]], [0.; 4]);
    }

    #[test]
    fn echelles_reflexions_rang_nul_et_non_finis() {
        let q = expm(&V3::new(0.6, -0.7, 1.8));
        let p = expm(&V3::new(-1.2, 0.4, 0.8));
        for scale in [1e-200, 1., 1e200] {
            for signe in [-1., 1.] {
                let a = q
                    * M3::from_diagonal(&V3::new(2. * scale, scale, signe * 0.5 * scale))
                    * p.transpose();
                let x = projette(&a).unwrap();
                assert!(
                    (x - q * p.transpose()).norm() < 5e-15,
                    "scale={scale}, signe={signe}, erreur={}, sigma={:?}",
                    (x - q * p.transpose()).norm(),
                    crate::svd_sure(DMatrix::from_column_slice(3, 3, a.as_slice()))
                        .unwrap()
                        .singular_values
                        .as_slice()
                );
            }
        }
        let x = projette(&M3::zeros()).unwrap();
        assert!((x.transpose() * x - M3::identity()).norm() < 4e-15);
        assert!((x.determinant() - 1.).abs() < 4e-15);
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(projette(&M3::from_element(v)).is_err());
        }
    }

    #[test]
    fn compositions_reversibles_reste_sur_le_groupe() {
        let initial = expm(&V3::new(-1.2, 0.7, 0.2));
        let pas = expm(&V3::new(0.17, -0.31, 0.23));
        let mut x = initial;
        for _ in 0..20_000 {
            x = projette(&(pas * x)).unwrap();
            x = projette(&(pas.transpose() * x)).unwrap();
            assert!((x.transpose() * x - M3::identity()).norm() < 4e-15);
        }
        assert!((x - initial).norm() < 2e-12);
    }
}
