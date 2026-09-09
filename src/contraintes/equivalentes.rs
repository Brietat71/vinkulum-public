//! Compression orthonormale des lignes exactement équivalentes (au signe près).
//! G = U G_r, Φ = U Φ_r et UᵀU = I : λ = U μ conserve la norme minimale.
//! Les contraintes originales restent présentes dans le modèle et ses résidus.
use super::Gradient;
use crate::{DVector, Triplets};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

// La clé emprunte une ligne CSR ; aucune allocation de clé par contrainte.
// La table compare aussi les coefficients exacts, même si les hachages coïncident.
struct Cle<'a> {
    phi: u64,
    signe: f64,
    termes: &'a [(usize, f64)],
}

impl PartialEq for Cle<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.phi == other.phi
            && self.termes.len() == other.termes.len()
            && self
                .termes
                .iter()
                .zip(other.termes)
                .all(|(&(j, v), &(k, w))| {
                    j == k && (self.signe * v).to_bits() == (other.signe * w).to_bits()
                })
    }
}

impl Eq for Cle<'_> {}

impl Hash for Cle<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.phi.hash(state);
        self.termes.len().hash(state);
        for &(j, v) in self.termes {
            j.hash(state);
            (self.signe * v).to_bits().hash(state);
        }
    }
}

pub(crate) struct Equivalentes {
    pub gradient: Gradient,
    pub phi: DVector<f64>,
    // Chaque ligne originale a une seule entrée non nulle dans U.
    lignes: Vec<(usize, f64)>,
    representatives: Vec<usize>,
    racines: Vec<f64>,
}

impl Equivalentes {
    pub fn nouvelles(g: &Gradient, phi: &DVector<f64>) -> Option<Self> {
        let mut offsets = vec![0; g.m + 1];
        for t in &g.trip {
            offsets[t.row + 1] += 1;
        }
        for i in 0..g.m {
            offsets[i + 1] += offsets[i];
        }
        let mut curseurs = offsets[..g.m].to_vec();
        let mut termes = vec![(0, 0.); g.trip.len()];
        // Les triplets sont déjà regroupés et triés par (colonne, ligne).
        for t in &g.trip {
            termes[curseurs[t.row]] = (t.col, t.val);
            curseurs[t.row] += 1;
        }
        let mut groupes = HashMap::new();
        let mut representatives = Vec::new();
        let mut tailles = Vec::new();
        let mut injection = Vec::with_capacity(g.m);
        for i in 0..g.m {
            let ligne = &termes[offsets[i]..offsets[i + 1]];
            // Une ligne nulle, notamment contradictoire, reste distincte.
            let signe = ligne
                .first()
                .map_or(1.0, |t| if t.1 < 0.0 { -1.0 } else { 1.0 });
            let cle = Cle {
                phi: if phi[i] == 0.0 {
                    0
                } else {
                    (signe * phi[i]).to_bits()
                },
                signe,
                termes: ligne,
            };
            let (k, relatif) = if let Some(&(k, premier_signe)) = groupes.get(&cle) {
                (k, signe * premier_signe)
            } else {
                let k = representatives.len();
                representatives.push(i);
                tailles.push(0usize);
                if !ligne.is_empty() {
                    groupes.insert(cle, (k, signe));
                }
                (k, 1.0)
            };
            tailles[k] += 1;
            injection.push((k, relatif));
        }
        if representatives.len() == g.m {
            return None;
        }
        let racines: Vec<_> = tailles.iter().map(|&k| (k as f64).sqrt()).collect();
        let mut trip = Triplets::new();
        for (k, &i) in representatives.iter().enumerate() {
            for &(j, v) in &termes[offsets[i]..offsets[i + 1]] {
                let val = racines[k] * v;
                if !val.is_finite() {
                    return None;
                }
                trip.push(faer::sparse::Triplet {
                    row: k,
                    col: j,
                    val,
                });
            }
        }
        crate::analyse::regroupe(&mut trip);
        let phi = DVector::from_iterator(
            representatives.len(),
            representatives
                .iter()
                .enumerate()
                .map(|(k, &i)| racines[k] * phi[i]),
        );
        if !phi.iter().all(|v| v.is_finite()) {
            return None;
        }
        Some(Self {
            gradient: Gradient {
                trip,
                n: g.n,
                m: representatives.len(),
            },
            phi,
            lignes: injection,
            representatives,
            racines,
        })
    }

    pub fn releve(&self, mu: &DVector<f64>) -> DVector<f64> {
        DVector::from_iterator(
            self.lignes.len(),
            self.lignes
                .iter()
                .map(|&(k, signe)| signe * (mu[k] / self.racines[k])),
        )
    }

    // Pendant une recherche de pas, deux fonctions ayant la même valeur et
    // tangente au départ peuvent se séparer. Le mérite naturel ne doit pas
    // ignorer cette partie de Φ hors de l'image de U. Aucun seuil de fusion.
    pub fn reduit_phi(&self, phi: &DVector<f64>) -> Option<DVector<f64>> {
        if self
            .lignes
            .iter()
            .enumerate()
            .any(|(i, &(k, signe))| phi[i] != signe * phi[self.representatives[k]])
        {
            return None;
        }
        Some(DVector::from_iterator(
            self.representatives.len(),
            self.representatives
                .iter()
                .enumerate()
                .map(|(k, &i)| self.racines[k] * phi[i]),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{contraintes::Projection, DMatrix};

    fn gradient(rows: &[&[f64]]) -> Gradient {
        let mut trip = Triplets::new();
        for (i, row) in rows.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                trip.push(faer::sparse::Triplet {
                    row: i,
                    col: j,
                    val,
                });
            }
        }
        crate::analyse::regroupe(&mut trip);
        Gradient {
            trip,
            n: rows[0].len(),
            m: rows.len(),
        }
    }

    #[test]
    fn reactions_minimales_multiplicites_signes_et_autre_dependance() {
        let g = gradient(&[
            &[1., 2., 0.],
            &[0., 1., 0.],
            &[-1., -2., 0.],
            &[1., 3., 0.],
            &[0., -1., 0.],
            &[1., 2., 0.],
            &[0., 0., 0.],
        ]);
        let phi = DVector::from_vec(vec![-0.375, -0.25, 0.375, -0.625, 0.25, -0.375, 0.]);
        let e = Equivalentes::nouvelles(&g, &phi).unwrap();
        assert_eq!(e.gradient.m, 4);
        // λ = G y appartient à l'image de G : c'est la solution de norme
        // minimale de Gᵀλ=f. La ligne r0+r1 reste une dépendance distincte.
        let reference = DVector::from_vec(vec![0.8, 0.3, -0.8, 1.1, -0.3, 0.8, 0.]);
        let mut f = g.transpose_fois(&reference);
        f[2] = 5.;
        let mu = e
            .gradient
            .reactions(&f, &mut Projection::default())
            .unwrap();
        let lambda = e.releve(&mu);
        assert!((&lambda - reference).amax() < 2e-15);
        assert!((lambda.norm() - mu.norm()).abs() < 1e-15);
        let residual = g.residu(&f, &lambda);
        assert!((residual - DVector::from_vec(vec![0., 0., 5.])).amax() < 5e-15);
    }

    #[test]
    fn newton_non_symetrique_equivaut_a_la_svd_complete() {
        let g = gradient(&[
            &[1., 0., 0., 1.],
            &[0., -2., -1., 0.],
            &[-1., 0., 0., -1.],
            &[0., 0., 1., 0.],
            &[1., 0., 0., 1.],
            &[0., 2., 1., 0.],
        ]);
        let phi = DVector::from_vec(vec![0.25, -0.75, -0.25, 0.25, 0.25, 0.75]);
        let e = Equivalentes::nouvelles(&g, &phi).unwrap();
        let (n, m, mr) = (g.n, g.m, e.gradient.m);
        let k = DMatrix::from_row_slice(
            4,
            4,
            &[
                5., 0.1, 0.2, -0.5, -0.3, 3., 0., 0.2, 0., 0.1, -2., 0.4, 0.5, 0., 0.1, 8.,
            ],
        );
        let d = [2., 3., 5., 7.];
        let assemble = |gr: &Gradient| {
            let mut a = DMatrix::zeros(n + gr.m, n + gr.m);
            for i in 0..n {
                for j in 0..n {
                    a[(i, j)] = k[(i, j)] / (d[i] * d[j]);
                }
            }
            for t in &gr.trip {
                a[(t.col, n + t.row)] = t.val / d[t.col];
                a[(n + t.row, t.col)] = t.val / d[t.col];
            }
            a
        };
        let a = assemble(&g);
        let ar = assemble(&e.gradient);
        let mut b = DVector::from_element(n + m, 0.);
        b.rows_mut(0, n)
            .copy_from(&DVector::from_vec(vec![1., -2., 3., -4.]));
        b.rows_mut(n, m).copy_from(&(-phi));
        let mut br = DVector::zeros(n + mr);
        br.rows_mut(0, n).copy_from(&b.rows(0, n));
        br.rows_mut(n, mr).copy_from(&(-&e.phi));
        let xr = ar.lu().solve(&br).unwrap();
        let mut x = DVector::zeros(n + m);
        x.rows_mut(0, n).copy_from(&xr.rows(0, n));
        x.rows_mut(n, m)
            .copy_from(&e.releve(&xr.rows(n, mr).into_owned()));
        let reference = crate::svd_sure(a.clone())
            .unwrap()
            .solve(&b, 1e-12)
            .unwrap();
        assert!((&x - &reference).amax() < 2e-13);
        assert!((&a * &x - b).amax() < 2e-14);
        assert!((x.norm() - xr.norm()).abs() < 3. * f64::EPSILON * x.norm());
    }

    #[test]
    fn ne_fusionne_ni_contradictions_ni_directions_proches_ni_lignes_nulles() {
        for delta in [1e-300, 1e-12, 1e-6] {
            let proche = gradient(&[&[1., 0.], &[1., delta]]);
            assert!(Equivalentes::nouvelles(&proche, &DVector::zeros(2)).is_none());
            let identique = gradient(&[&[1., 0.], &[1., 0.]]);
            let phi = DVector::from_vec(vec![0., delta]);
            assert!(Equivalentes::nouvelles(&identique, &phi).is_none());
        }
        let nul = gradient(&[&[0., 0.], &[0., 0.]]);
        assert!(Equivalentes::nouvelles(&nul, &DVector::zeros(2)).is_none());
        let g = gradient(&[&[1., 0.], &[-1., 0.]]);
        assert!(Equivalentes::nouvelles(&g, &DVector::from_vec(vec![1., 1.])).is_none());
        let e = Equivalentes::nouvelles(&g, &DVector::from_vec(vec![-0., 0.])).unwrap();
        assert!(e.reduit_phi(&DVector::from_vec(vec![0.5, -0.5])).is_some());
        assert!(e
            .reduit_phi(&DVector::from_vec(vec![0.5, -0.5_f64.next_up()]))
            .is_none());
        let immense = gradient(&[&[f64::MAX, 0.], &[-f64::MAX, 0.]]);
        assert!(Equivalentes::nouvelles(&immense, &DVector::zeros(2)).is_none());
    }

    #[test]
    fn reevalue_groupes_valeurs_et_motif() {
        let mut cache = Projection::default();
        for rows in [
            vec![vec![1., 0., 0.], vec![1., 0., 0.], vec![0., 1., 0.]],
            vec![vec![1., 0., 0.], vec![0., -1., 0.], vec![0., 1., 0.]],
            vec![vec![1., 1., 0.], vec![-1., -1., 0.], vec![0., 0., 1.]],
            vec![vec![1., 1., 0.], vec![-1., -1., 0.], vec![1., 1., 0.]],
        ] {
            let g = gradient(&rows.iter().map(Vec::as_slice).collect::<Vec<_>>());
            let e = Equivalentes::nouvelles(&g, &DVector::zeros(3)).unwrap();
            let reference =
                DVector::from_iterator(3, rows.iter().map(|r| r[0] + 2. * r[1] + 3. * r[2]));
            let f = g.transpose_fois(&reference);
            let mu = e.gradient.reactions(&f, &mut cache).unwrap();
            assert!((e.releve(&mu) - reference).amax() < 2e-15);
        }
    }
}
