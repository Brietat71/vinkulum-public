//! Résolution d'un système de Newton par raffinement d'une factorisation
//! auxiliaire. La diagonale des multiplicateurs n'est modifiée que dans le
//! préconditionneur ; tous les résidus portent sur le système original.
use crate::{numerique::deux_sommes, DVector, Facto, Triplets};

#[derive(Default)]
pub(crate) struct Cache {
    dimension: usize,
    motif: Vec<(usize, usize)>,
    symbolique: Option<std::sync::Arc<faer::sparse::linalg::lu::SymbolicLu<usize>>>,
}

pub(crate) struct Facteur {
    lu: Facto,
}

fn norme(v: &DVector<f64>) -> f64 {
    v.iter().fold(0.0_f64, |s, &x| s.hypot(x))
}

// Les petits résidus ne doivent pas être jugés par la soustraction de
// grands produits déjà arrondis. TwoSum et le reste du produit par FMA
// accumulent b-Ax en deux composantes, puis le résultat est arrondi en f64.
pub(crate) fn residu(a: &Triplets, b: &DVector<f64>, x: &DVector<f64>) -> DVector<f64> {
    let mut haut = b.clone();
    let mut bas = DVector::zeros(b.len());
    for t in a {
        let p = -t.val * x[t.col];
        let e = (-t.val).mul_add(x[t.col], -p);
        let (s, r) = deux_sommes(haut[t.row], p);
        (haut[t.row], bas[t.row]) = deux_sommes(s, bas[t.row] + r + e);
    }
    haut += bas;
    haut
}

impl Facteur {
    pub fn nouveau(
        a: &Triplets,
        physiques: usize,
        dimension: usize,
        cache: &mut Cache,
    ) -> Result<Self, String> {
        let mut b = a.clone();
        let maximum = a.iter().map(|t| t.val.abs()).fold(0.0_f64, f64::max);
        let decalage = f64::EPSILON.sqrt() * maximum;
        if !decalage.is_finite() || decalage == 0. {
            return Err("raffinement : échelle du préconditionneur invalide".into());
        }
        for i in physiques..dimension {
            b.push(faer::sparse::Triplet {
                row: i,
                col: i,
                val: -decalage,
            });
        }
        crate::analyse::regroupe(&mut b);
        let motif: Vec<_> = b.iter().map(|t| (t.row, t.col)).collect();
        if cache.dimension != dimension || cache.motif != motif {
            cache.dimension = dimension;
            cache.motif = motif;
            cache.symbolique = None;
        }
        let lu = Facto::creux(&b, dimension, &mut cache.symbolique)?;
        Ok(Self { lu })
    }

    pub fn resout(&self, a: &Triplets, b: &DVector<f64>) -> Option<DVector<f64>> {
        let nb = norme(b);
        if nb == 0. {
            return Some(DVector::zeros(b.len()));
        }
        if !nb.is_finite() {
            return None;
        }
        let mut x = self.lu.solve(b);
        let mut precedent = f64::INFINITY;
        for _ in 0..8 {
            if !x.iter().all(|v| v.is_finite()) {
                return None;
            }
            let r = residu(a, b, &x);
            let nr = norme(&r);
            if !nr.is_finite() {
                return None;
            }
            if nr <= 1e-12 * nb {
                return Some(x);
            }
            if nr >= 0.75 * precedent {
                return None;
            }
            precedent = nr;
            x += self.lu.solve(&r);
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DMatrix;

    fn triplets(a: &DMatrix<f64>) -> Triplets {
        let mut t = Triplets::new();
        for j in 0..a.ncols() {
            for i in 0..a.nrows() {
                if a[(i, j)] != 0. {
                    t.push(faer::sparse::Triplet {
                        row: i,
                        col: j,
                        val: a[(i, j)],
                    });
                }
            }
        }
        t
    }

    #[test]
    fn residu_produits_et_sommes_compenses() {
        let a = DMatrix::from_row_slice(1, 3, &[1., 2.0_f64.powi(53), -2.0_f64.powi(53)]);
        let b = DVector::from_element(1, 0.);
        let x = DVector::from_element(3, 1.);
        assert_eq!(residu(&triplets(&a), &b, &x)[0], -1.);
        let a = DMatrix::from_element(1, 1, 1. - 2.0_f64.powi(-27));
        let b = DVector::from_element(1, 1.);
        let x = DVector::from_element(1, 1. + 2.0_f64.powi(-27));
        assert_eq!(residu(&triplets(&a), &b, &x)[0], 2.0_f64.powi(-54));
    }

    #[test]
    fn systeme_singulier_compatible_et_reutilisation() {
        // G a quatre lignes dont deux combinaisons générales. La partie
        // multiplicateur de la référence est G*y, donc de norme minimale.
        let g = DMatrix::from_row_slice(
            4,
            6,
            &[
                1., 2., 0., 0., 0., 0., 0., 0., 1., -1., 0., 0., 1., 2., 1., -1., 0., 0., 2., 4.,
                -1., 1., 0., 0.,
            ],
        );
        let mut a = DMatrix::zeros(10, 10);
        for i in 0..6 {
            a[(i, i)] = 2. + i as f64;
        }
        a[(4, 5)] = 0.7; // K non symétrique.
        a.view_mut((6, 0), (4, 6)).copy_from(&g);
        a.view_mut((0, 6), (6, 4)).copy_from(&g.transpose());
        let t = triplets(&a);
        let mut cache = Cache::default();
        let f = Facteur::nouveau(&t, 6, 10, &mut cache).unwrap();
        for scale in [1e-20, 1., 1e20] {
            let q = DVector::from_iterator(6, (0..6).map(|i| (i as f64 + 1.).sin() * scale));
            let mut attendu = DVector::zeros(10);
            attendu.rows_mut(0, 6).copy_from(&q);
            attendu.rows_mut(6, 4).copy_from(&(&g * &q));
            let b = &a * &attendu;
            let x = f.resout(&t, &b).expect("raffinement compatible");
            assert!((&a * &x - &b).norm() <= 2e-12 * b.norm());
            assert!((x - attendu).norm() < 2e-7 * scale);
        }
    }

    #[test]
    fn contraintes_contradictoires_et_direction_physique_libre() {
        let a = DMatrix::from_row_slice(3, 3, &[2., 1., 1., 1., 0., 0., 1., 0., 0.]);
        let t = triplets(&a);
        let f = Facteur::nouveau(&t, 1, 3, &mut Cache::default()).unwrap();
        let b = DVector::from_vec(vec![0., 1., -1.]);
        assert!(f.resout(&t, &b).is_none());
        let a = DMatrix::from_row_slice(3, 3, &[0., 0., 0., 0., 0., 1., 0., 1., 0.]);
        let t = triplets(&a);
        let b = DVector::from_vec(vec![1., 0., 0.]);
        assert!(Facteur::nouveau(&t, 2, 3, &mut Cache::default())
            .ok()
            .and_then(|f| f.resout(&t, &b))
            .is_none());
    }
}
