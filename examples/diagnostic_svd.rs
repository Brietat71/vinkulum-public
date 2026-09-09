//! Cas fabriqué : x=Aᵀy appartient à l'espace des lignes, b=Ax.
//! La solution de norme minimale est donc connue sans décomposition.
//! Sortie JSON pour vérification indépendante avec LAPACK.
use nalgebra::{DMatrix, DVector, Dyn};

fn facteurs(s: &nalgebra::linalg::SVD<f64, Dyn, Dyn>, b: &DVector<f64>, seuil: f64) -> String {
    format!(
        "{{\"sigma\":{:?},\"x\":{:?},\"u_col\":{:?},\"vt_col\":{:?}}}",
        s.singular_values.as_slice(),
        s.solve(b, seuil).unwrap().as_slice(),
        s.u.as_ref().unwrap().as_slice(),
        s.v_t.as_ref().unwrap().as_slice()
    )
}

fn main() {
    let (m, n, rang) = (7, 5, 3);
    let gauche = DMatrix::from_fn(m, rang, |i, j| {
        ((i + 1) as f64 * (j + 1) as f64 * 0.731).sin()
    });
    let droite = DMatrix::from_fn(rang, n, |i, j| {
        ((i + 1) as f64 * (j + 1) as f64 * 0.537).cos()
    });
    let a = gauche * droite;
    let y = DVector::from_fn(m, |i, _| (0.21 * i as f64).sin());
    let reference = a.transpose() * &y;
    let b = &a * &reference;
    let seuil = f64::EPSILON * m.max(n) as f64 * a.norm();
    let ancienne = a
        .clone()
        .try_svd(true, true, f64::EPSILON, 1000 * m.max(n))
        .unwrap();
    let courante = vinkulum::svd_sure(a.clone()).unwrap();
    println!(
        "{{\"m\":{m},\"n\":{n},\"rang\":{rang},\"a_col\":{:?},\"y\":{:?},\"b\":{:?},\"reference\":{:?},\"seuil\":{seuil:?},\"nalgebra\":{},\"faer\":{}}}",
        a.as_slice(),
        y.as_slice(),
        b.as_slice(),
        reference.as_slice(),
        facteurs(&ancienne, &b, seuil),
        facteurs(&courante, &b, seuil)
    );
}
