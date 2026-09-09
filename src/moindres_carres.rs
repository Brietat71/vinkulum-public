//! Moindres carrés de norme minimale par factorisation orthogonale complète.
//! A P = Q R puis R_rᵀ = Z T : x = P Z [T⁻ᵀ (Qᵀb)_r ; 0].
//! Les deux matrices orthogonales restent implicites. Un rang mal séparé
//! ou un contrôle numérique insuffisant conserve le repli SVD antérieur.
use crate::{DMatrix, DVector};
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::linalg::{
    householder as hh,
    qr::{col_pivoting::factor as cpqr, no_pivoting::factor as qr},
};

fn memoire(req: StackReq) -> Result<MemBuffer, String> {
    MemBuffer::try_new(req).map_err(|e| format!("mémoire de factorisation orthogonale : {e:?}"))
}

fn taille_bloc(m: usize, n: usize) -> usize {
    // Petites matrices : appliquer les réflecteurs individuellement évite
    // de préparer des produits matriciels par blocs pour quelques dizaines
    // de coordonnées. Les décisions de rang restent identiques.
    if m.max(n) <= 64 {
        1
    } else {
        qr::recommended_block_size::<f64>(m, n)
    }
}

// Pour un triangle R, les résolutions avec |diag R| et -|hors-diag R|
// majorent les sommes absolues des lignes/colonnes de R⁻¹. Donc
// σ_min(R) >= 1/sqrt(borne_1 * borne_inf), sans inverse ni matrice normale.
// Les opérations positives sont arrondies vers le haut, puis la borne
// finale vers le bas. Cette garantie concerne le triangle calculé.
fn borne_sigma(r: faer::MatRef<'_, f64>, rang: usize) -> f64 {
    if rang == 0 {
        return f64::INFINITY;
    }
    let mut lignes = vec![0.; rang];
    let mut colonnes = vec![0.; rang];
    for i in (0..rang).rev() {
        let mut somme = 1.;
        for j in i + 1..rang {
            somme = (somme + (r[(i, j)].abs() * lignes[j]).next_up()).next_up();
        }
        lignes[i] = (somme / r[(i, i)].abs()).next_up();
    }
    for j in 0..rang {
        let mut somme = 1.;
        for i in 0..j {
            somme = (somme + (r[(i, j)].abs() * colonnes[i]).next_up()).next_up();
        }
        colonnes[j] = (somme / r[(j, j)].abs()).next_up();
    }
    let ni = lignes.into_iter().fold(0., f64::max);
    let n1 = colonnes.into_iter().fold(0., f64::max);
    let majorant = (ni.sqrt().next_up() * n1.sqrt().next_up()).next_up();
    (1. / majorant).next_down().max(0.)
}

fn essaie(a: &DMatrix<f64>, b: &DVector<f64>, seuil: f64) -> Result<Option<DVector<f64>>, String> {
    let (m, n) = a.shape();
    if m == 0 || n == 0 {
        return Ok(Some(DVector::zeros(n)));
    }
    let echelle = a.amax();
    if echelle == 0. {
        return Ok(Some(DVector::zeros(n)));
    }
    let seuil = seuil / echelle;
    if !seuil.is_finite() || seuil <= 0. {
        return Ok(None);
    }
    let mut f = DVector::from_iterator(m, b.iter().map(|v| v / echelle));
    if f.iter()
        .zip(b.iter())
        .any(|(&x, &y)| !x.is_finite() || (x == 0. && y != 0.))
    {
        return Ok(None);
    }
    let f0 = f.clone();
    let mut q = faer::Mat::from_fn(m, n, |i, j| a[(i, j)] / echelle);
    if (0..n).any(|j| (0..m).any(|i| q[(i, j)] == 0. && a[(i, j)] != 0.)) {
        return Ok(None);
    }
    let k = m.min(n);
    let bs = taille_bloc(m, n);
    let mut h = faer::Mat::zeros(bs, k);
    let (mut p, mut pi) = (vec![0usize; n], vec![0usize; n]);
    let mut mem = memoire(cpqr::qr_in_place_scratch::<usize, f64>(
        m,
        n,
        bs,
        faer::Par::Seq,
        Default::default(),
    ))?;
    cpqr::qr_in_place(
        q.as_mut(),
        h.as_mut(),
        &mut p,
        &mut pi,
        faer::Par::Seq,
        MemStack::new(&mut mem),
        Default::default(),
    );
    if (0..k).any(|i| (i..n).any(|j| !q[(i, j)].is_finite())) {
        return Ok(None);
    }
    let rang = (0..k).find(|&i| q[(i, i)].abs() <= seuil).unwrap_or(k);
    let mut reste = 0.0_f64;
    for i in rang..k {
        for j in i..n {
            let v = q[(i, j)];
            if v != 0. {
                reste = (reste + (v * v).next_up()).next_up();
            }
        }
    }
    let reste = reste.sqrt().next_up();
    // Bande de séparation autour du seuil de la SVD. Les bornes portent
    // sur le R calculé ; le contrôle des moindres carrés suit sur A original.
    if !reste.is_finite() || reste > 0.25 * seuil || !(borne_sigma(q.as_ref(), rang) > 4. * seuil) {
        return Ok(None);
    }
    if rang == 0 {
        return Ok(Some(DVector::zeros(n)));
    }
    let mut mem = memoire(
        hh::apply_block_householder_sequence_transpose_on_the_left_in_place_scratch::<f64>(
            m, bs, 1,
        ),
    )?;
    hh::apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj(
        q.as_ref(),
        h.as_ref(),
        faer::Conj::No,
        faer::MatMut::from_column_major_slice_mut(f.as_mut_slice(), m, 1),
        faer::Par::Seq,
        MemStack::new(&mut mem),
    );
    // R_rᵀ, en omettant les réflecteurs situés sous la diagonale de q.
    let mut z = faer::Mat::from_fn(n, rang, |i, j| if i >= j { q[(j, i)] } else { 0. });
    let bz = taille_bloc(n, rang);
    let mut hz = faer::Mat::zeros(bz, rang);
    let mut mem = memoire(qr::qr_in_place_scratch::<f64>(
        n,
        rang,
        bz,
        faer::Par::Seq,
        Default::default(),
    ))?;
    qr::qr_in_place(
        z.as_mut(),
        hz.as_mut(),
        faer::Par::Seq,
        MemStack::new(&mut mem),
        Default::default(),
    );
    let mut y = DVector::zeros(n);
    // Tᵀ y = (Qᵀ b)_r ; les autres composantes sont nulles.
    for i in 0..rang {
        let mut value = f[i];
        for j in 0..i {
            value -= z[(j, i)] * y[j];
        }
        let pivot = z[(i, i)];
        if !pivot.is_finite() || pivot.abs() <= seuil {
            return Ok(None);
        }
        y[i] = value / pivot;
    }
    let mut mem = memoire(
        hh::apply_block_householder_sequence_on_the_left_in_place_scratch::<f64>(n, bz, 1),
    )?;
    hh::apply_block_householder_sequence_on_the_left_in_place_with_conj(
        z.as_ref(),
        hz.as_ref(),
        faer::Conj::No,
        faer::MatMut::from_column_major_slice_mut(y.as_mut_slice(), n, 1),
        faer::Par::Seq,
        MemStack::new(&mut mem),
    );
    let mut x = DVector::zeros(n);
    for j in 0..n {
        x[p[j]] = y[j];
    }
    if !x.iter().all(|v| v.is_finite()) {
        return Ok(None);
    }
    // Stationnarité Aᵀ(Ax-b)=0, sans former AᵀA. La normalisation scalaire
    // commune évite les produits démesurés ; aucun coefficient n'est coupé.
    let mut r = -&f0;
    let mut norme_a = 0.0_f64;
    for j in 0..n {
        for i in 0..m {
            let v = a[(i, j)] / echelle;
            r[i] += v * x[j];
            norme_a = norme_a.hypot(v);
        }
    }
    let mut normal = DVector::<f64>::zeros(n);
    for j in 0..n {
        for i in 0..m {
            normal[j] += (a[(i, j)] / echelle) * r[i];
        }
    }
    let borne = 64. * f64::EPSILON * m.max(n) as f64 * norme_a * (norme_a * x.norm() + f0.norm());
    if !borne.is_finite() || !normal.iter().all(|v| v.is_finite()) || normal.amax() > borne {
        return Ok(None);
    }
    Ok(Some(x))
}

pub(crate) fn resout(
    a: DMatrix<f64>,
    b: &DVector<f64>,
    seuil: f64,
) -> Result<DVector<f64>, String> {
    crate::fini(a.as_slice(), "matrice de moindres carrés")?;
    crate::fini(b.as_slice(), "second membre de moindres carrés")?;
    if let Some(x) = essaie(&a, b, seuil)? {
        return Ok(x);
    }
    crate::svd_sure(a)?
        .solve(b, seuil)
        .map_err(|e| format!("moindres carrés SVD : {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rang_general_et_norme_minimale() {
        // Colonnes et lignes denses, permutations différentes, rangs variés.
        // x=Aᵀy et b=Ax garantissent une référence de norme minimale sans SVD.
        for (m, n, rang) in [(7, 5, 3), (5, 9, 4), (16, 16, 11), (80, 101, 33)] {
            let u = DMatrix::from_fn(m, rang, |i, j| {
                ((i + 1) as f64 * (j + 1) as f64 * 0.731).sin()
            });
            let v = DMatrix::from_fn(rang, n, |i, j| {
                ((i + 1) as f64 * (j + 1) as f64 * 0.537).cos()
            });
            let a = u * v;
            let y = DVector::from_fn(m, |i, _| (0.21 * i as f64).sin());
            let reference = a.transpose() * y;
            let b = &a * &reference;
            let seuil = f64::EPSILON * m.max(n) as f64 * a.norm();
            let x = essaie(&a, &b, seuil).unwrap().expect("rang bien séparé");
            assert!(
                (&x - &reference).norm() < 2e-12 * reference.norm(),
                "{m}×{n}"
            );
            let decomposition = crate::svd_sure(a.clone()).unwrap();
            let svd = decomposition.solve(&b, seuil).unwrap();
            assert!(
                (&x - &svd).norm() < 3e-12 * reference.norm(),
                "{m}×{n}, rang {rang}: COD-réf {}, SVD-réf {}, norme réf {}",
                (&x - &reference).norm(),
                (&svd - &reference).norm(),
                reference.norm()
            );
        }
    }

    #[test]
    fn borne_spectrale_conservative_sur_triangle() {
        for a in [0.0_f64, 0.1, 1., 100., 1e12] {
            let r = faer::Mat::from_fn(2, 2, |i, j| {
                if i == j {
                    1.
                } else if i < j {
                    a
                } else {
                    0.
                }
            });
            let sigma = 2. / ((a * a + 4.).sqrt() + a);
            let borne = borne_sigma(r.as_ref(), 2);
            assert!(
                borne > 0. && borne <= sigma,
                "{a}: borne {borne}, sigma {sigma}"
            );
        }
    }

    #[test]
    fn forces_incompatibles_et_directions_presque_dependantes() {
        for delta in [0., 1e-12, 1e-9, 1e-6, 1.] {
            let a = DMatrix::from_row_slice(3, 2, &[1., 1., 0., delta, 0., 0.]);
            let b = DVector::from_vec(vec![2., 1., 3.]);
            let seuil = f64::EPSILON * 3. * a.norm();
            let expected = if delta == 0. {
                vec![1., 1.]
            } else {
                vec![2. - 1. / delta, 1. / delta]
            };
            let expected = DVector::from_vec(expected);
            let x = resout(a, &b, seuil).unwrap();
            assert!(
                (&x - &expected).norm() <= 2e-14 * expected.norm().max(1.),
                "{delta}: {x:?}"
            );
        }
    }

    #[test]
    fn rang_ambigu_replie_et_valeurs_non_finies_refusees() {
        for delta in [0.5e-12, 1e-12, 2e-12] {
            let a = DMatrix::from_diagonal(&DVector::from_vec(vec![1., delta]));
            let b = DVector::from_vec(vec![1., 1.]);
            assert!(essaie(&a, &b, 1e-12).unwrap().is_none());
            let reference = crate::svd_sure(a.clone())
                .unwrap()
                .solve(&b, 1e-12)
                .unwrap();
            assert_eq!(resout(a, &b, 1e-12).unwrap(), reference);
        }
        for v in [f64::NAN, f64::INFINITY] {
            assert!(resout(DMatrix::from_element(2, 2, v), &DVector::zeros(2), 1e-12).is_err());
            assert!(resout(DMatrix::identity(2, 2), &DVector::from_element(2, v), 1e-12).is_err());
        }
    }
}
