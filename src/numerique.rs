//! Opérations numériques communes, avec contrôles de finitude et bornes.
use crate::{robustesse::fini, DMatrix};

/// Somme sans perte : s+e représente a+b (hors débordement).
pub(crate) fn deux_sommes(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    (s, (a - (s - v)) + (b - v))
}

/// Accumule un incrément dans une position représentée par deux f64.
pub(crate) fn deplace_compense(
    r: crate::V3,
    bas: crate::V3,
    dr: crate::V3,
) -> (crate::V3, crate::V3) {
    let mut haut = r;
    let mut reste = bas;
    for i in 0..3 {
        let (s, e) = deux_sommes(r[i], dr[i]);
        (haut[i], reste[i]) = deux_sommes(s, bas[i] + e);
    }
    (haut, reste)
}

/// Dernier pas absorbé jusqu'à 1,5 h en pas fixe, sans seuil absolu de temps.
pub(crate) fn fin_pas(t: f64, fin: f64, h: f64, absorbe: bool) -> Result<f64, String> {
    let mut suivant = (t + h).min(fin);
    if absorbe && fin - suivant <= 0.5 * h {
        suivant = fin;
    }
    if !suivant.is_finite() || suivant <= t {
        return Err("le pas ne fait pas progresser le temps".into());
    }
    Ok(suivant)
}

/// Noyau orthonormal de G, sans former GᵀG ni perdre les directions lorsque
/// les lignes ont des échelles différentes. QR avec pivotement de Gᵀ ; les
/// réflecteurs appliqués à I donnent aussi le complément absent du Q mince.
pub(crate) fn noyau(mut g: DMatrix<f64>) -> Result<DMatrix<f64>, String> {
    fini(g.as_slice(), "jacobien des contraintes")?;
    let (m, n) = g.shape();
    if m == 0 || n == 0 {
        return Ok(DMatrix::identity(n, n));
    }
    for mut row in g.row_iter_mut() {
        let s = row.amax();
        if s > 0.0 {
            row /= s;
            let norm = row.norm();
            row /= norm;
        }
    }
    let qr = g.transpose().col_piv_qr();
    let r = qr.r();
    let tol = 64.0 * f64::EPSILON * m.max(n) as f64 * r.amax();
    let rang = (0..m.min(n)).take_while(|&i| r[(i, i)].abs() > tol).count();
    let mut qt = DMatrix::identity(n, n);
    qr.q_tr_mul(&mut qt);
    let z = qt.rows(rang, n - rang).transpose();
    fini(z.as_slice(), "base admissible")?;
    Ok(z)
}

pub(crate) fn propres_sym(
    a: DMatrix<f64>,
) -> Result<nalgebra::linalg::SymmetricEigen<f64, nalgebra::Dyn>, String> {
    fini(a.as_slice(), "matrice spectrale")?;
    let plafond = 1000 * a.nrows().max(1);
    let e = nalgebra::linalg::SymmetricEigen::try_new(a, f64::EPSILON, plafond)
        .ok_or("valeurs propres : plafond d'itérations atteint")?;
    fini(e.eigenvalues.as_slice(), "valeurs propres")?;
    fini(e.eigenvectors.as_slice(), "vecteurs propres")?;
    Ok(e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn position_petits_increments_et_annulation() {
        use crate::V3;
        let initial = V3::new(1., 2_f64.powi(40), -1.);
        let pas = V3::new(2_f64.powi(-60), 2_f64.powi(-20), 2_f64.powi(-60));
        let (mut haut, mut bas) = (initial, V3::zeros());
        for _ in 0..128 {
            (haut, bas) = deplace_compense(haut, bas, pas);
        }
        let (difference, reste) = deplace_compense(haut, bas, -initial);
        assert_eq!(difference, 128. * pas);
        assert_eq!(reste, V3::zeros());
        for _ in 0..128 {
            (haut, bas) = deplace_compense(haut, bas, -pas);
        }
        assert_eq!(haut, initial);
        assert_eq!(bas, V3::zeros());
    }

    #[test]
    fn noyau_rectangulaire_et_echelles() {
        // Deux équations indépendantes, une répétée, un complément de dim 2.
        let g = DMatrix::from_row_slice(
            3,
            4,
            &[1e200, 0., 2e200, 0., 0., 1e-200, 0., 0., 3., 0., 6., 0.],
        );
        let z = noyau(g).unwrap();
        assert_eq!(z.shape(), (4, 2));
        let contraintes = DMatrix::from_row_slice(2, 4, &[1., 0., 2., 0., 0., 1., 0., 0.]);
        assert!((contraintes * &z).norm() < 1e-14);
        assert!((z.transpose() * z - DMatrix::identity(2, 2)).norm() < 1e-14);
        assert_eq!(noyau(DMatrix::identity(5, 3)).unwrap().shape(), (3, 0));
        let libre = noyau(DMatrix::zeros(5, 3)).unwrap();
        assert_eq!(libre.shape(), (3, 3));
        assert!((libre.transpose() * libre - DMatrix::identity(3, 3)).norm() < 1e-14);
        assert_eq!(
            noyau(DMatrix::zeros(0, 3)).unwrap(),
            DMatrix::identity(3, 3)
        );
        assert!(noyau(DMatrix::from_element(1, 1, f64::NAN)).is_err());
    }

    #[test]
    fn spectre_fini_et_borne() {
        assert!(propres_sym(DMatrix::from_element(2, 2, f64::INFINITY)).is_err());
        let a = DMatrix::from_row_slice(2, 2, &[2., 1., 1., 2.]);
        let e = propres_sym(a.clone()).unwrap();
        assert!(
            (&a * &e.eigenvectors - &e.eigenvectors * DMatrix::from_diagonal(&e.eigenvalues))
                .norm()
                < 1e-14
        );
        assert!((e.eigenvalues.min() - 1.).abs() < 1e-14);
        assert!((e.eigenvalues.max() - 3.).abs() < 1e-14);
    }

    #[test]
    fn tables_rust_finies() {
        use crate::Loi;
        for pts in [vec![(f64::NAN, 0.)], vec![(0., f64::INFINITY)]] {
            assert!(Loi::table(pts).is_err());
        }
        let loi = Loi::table(vec![(0., 2.), (1., 4.)]).unwrap();
        assert_eq!(loi.valeur(0.25), 2.5);
        assert_eq!(loi.derivee(0.25), 2.);
    }
}
