//! SVD mince avec faer, sans changer le type de résultat utilisé par le noyau.
//! La SVD précédente perdait des chiffres sur un rang déficient 7×5 :
//! reconstruction et solution fabriquée sont contrôlées indépendamment.
use crate::{DMatrix, DVector};
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::svd::{svd, svd_scratch, ComputeSvdVectors};

pub(crate) fn decompose(
    a: &DMatrix<f64>,
) -> Result<nalgebra::linalg::SVD<f64, nalgebra::Dyn, nalgebra::Dyn>, String> {
    let (m, n) = a.shape();
    let k = m.min(n);
    // Une puissance de deux normalise sans arrondi les coefficients qui
    // restent représentables. Le bidiagonaliseur perdait le spectre d'une
    // matrice dense 3×3 multipliée par 1e-200 ; le cas diagonal le masquait.
    let maximum = a.amax();
    let bits = maximum.to_bits();
    let exposant = bits & 0x7ff0_0000_0000_0000;
    let echelle = if maximum == 0. {
        1.
    } else if exposant != 0 {
        f64::from_bits(exposant)
    } else {
        f64::from_bits(1_u64 << (63 - bits.leading_zeros()))
    };
    let normalisee = (echelle != 1.).then(|| a / echelle);
    let entree = normalisee.as_ref().unwrap_or(a);
    if normalisee
        .as_ref()
        .is_some_and(|scaled| scaled.iter().zip(a.iter()).any(|(&x, &y)| x * echelle != y))
    {
        return Err(
            "SVD : étendue des coefficients incompatible avec la normalisation en f64".into(),
        );
    }
    let mut u = DMatrix::zeros(m, k);
    let mut v = DMatrix::zeros(n, k);
    let mut valeurs = DVector::zeros(k);
    let mut mem = MemBuffer::try_new(svd_scratch::<f64>(
        m,
        n,
        ComputeSvdVectors::Thin,
        ComputeSvdVectors::Thin,
        faer::Par::Seq,
        Default::default(),
    ))
    .map_err(|e| format!("mémoire SVD : {e:?}"))?;
    svd(
        faer::MatRef::from_column_major_slice(entree.as_slice(), m, n),
        faer::diag::DiagMut::from_slice_mut(valeurs.as_mut_slice()),
        Some(faer::MatMut::from_column_major_slice_mut(
            u.as_mut_slice(),
            m,
            k,
        )),
        Some(faer::MatMut::from_column_major_slice_mut(
            v.as_mut_slice(),
            n,
            k,
        )),
        faer::Par::Seq,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .map_err(|e| format!("SVD : {e:?}"))?;
    valeurs *= echelle;
    if !u
        .iter()
        .chain(v.iter())
        .chain(valeurs.iter())
        .all(|x| x.is_finite())
    {
        return Err("SVD : facteurs non finis".into());
    }
    Ok(nalgebra::linalg::SVD {
        u: Some(u),
        v_t: Some(v.transpose()),
        singular_values: valeurs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn echelles_extremes_rang_nul_et_dimensions_vides() {
        for echelle in [1e-200, 1., 1e200] {
            let a = DMatrix::from_diagonal(&DVector::from_vec(vec![
                echelle,
                -3. * echelle,
                2. * echelle,
            ]));
            let s = crate::svd_sure(a).unwrap();
            let expected = DVector::from_vec(vec![3., 2., 1.]);
            assert!((s.singular_values / echelle - expected).amax() < 1e-14);
        }
        for (m, n) in [(0, 5), (5, 0), (3, 3)] {
            let s = crate::svd_sure(DMatrix::zeros(m, n)).unwrap();
            assert_eq!(
                s.solve(&DVector::from_element(m, 1.), 1e-12).unwrap(),
                DVector::zeros(n)
            );
        }
        let min = f64::from_bits(1);
        let s = crate::svd_sure(DMatrix::from_element(1, 1, min)).unwrap();
        assert_eq!(s.singular_values[0], min);
        let etendue = DVector::from_vec(vec![
            2.0_f64.powi(500),
            3.0_f64.next_up() * 2.0_f64.powi(-574),
        ]);
        assert!(crate::svd_sure(DMatrix::from_diagonal(&etendue)).is_err());
        assert!(
            crate::svd_sure(DMatrix::from_diagonal(&DVector::from_vec(vec![
                f64::MAX,
                min
            ])))
            .is_err()
        );
        assert!(crate::svd_sure(DMatrix::from_element(3, 3, 1e308)).is_err());
    }

    #[test]
    fn spectre_dense_independant_de_lechelle() {
        use crate::{expm, V3};
        let gauche = expm(&V3::new(0.6, -0.7, 1.8));
        let droite = expm(&V3::new(-1.2, 0.4, 0.8));
        for scale in [1e-300, 1e-200, 1., 1e200, 1e300] {
            let reference = DVector::from_vec(vec![2., 1., 0.5]);
            let a = gauche
                * crate::M3::from_diagonal(&V3::new(2. * scale, scale, -0.5 * scale))
                * droite.transpose();
            let s = crate::svd_sure(DMatrix::from_column_slice(3, 3, a.as_slice())).unwrap();
            assert!(
                (&s.singular_values / scale - reference).norm() < 3e-15,
                "{scale}"
            );
            let recompose = s.u.unwrap()
                * DMatrix::from_diagonal(&(s.singular_values / scale))
                * s.v_t.unwrap();
            assert!(
                (recompose - DMatrix::from_column_slice(3, 3, (a / scale).as_slice())).norm()
                    < 4e-15
            );
        }
    }

    #[test]
    fn rang_deficient_reconstruction_et_norme_minimale() {
        for (m, n, rang) in [(7, 5, 3), (5, 7, 3), (16, 16, 11), (80, 101, 33)] {
            let gauche = DMatrix::from_fn(m, rang, |i, j| {
                ((i + 1) as f64 * (j + 1) as f64 * 0.731).sin()
            });
            let droite = DMatrix::from_fn(rang, n, |i, j| {
                ((i + 1) as f64 * (j + 1) as f64 * 0.537).cos()
            });
            let a = gauche * droite;
            let s = crate::svd_sure(a.clone()).unwrap();
            let u = s.u.as_ref().unwrap();
            let vt = s.v_t.as_ref().unwrap();
            let k = m.min(n);
            let recomposition = u * DMatrix::from_diagonal(&s.singular_values) * vt;
            assert!((&recomposition - &a).norm() < 2e-14 * a.norm(), "{m}×{n}");
            assert!((u.transpose() * u - DMatrix::identity(k, k)).norm() < 2e-13);
            assert!((vt * vt.transpose() - DMatrix::identity(k, k)).norm() < 2e-13);
            let y = DVector::from_fn(m, |i, _| (0.21 * i as f64).sin());
            let reference = a.transpose() * y;
            let b = &a * &reference;
            let seuil = f64::EPSILON * m.max(n) as f64 * a.norm();
            let x = s.solve(&b, seuil).unwrap();
            assert!(
                (&x - &reference).norm() < 2e-12 * reference.norm(),
                "{m}×{n}"
            );
        }
    }
}
