//! Projection statique sur les contraintes, sans matrice normale globale.
//! Les composantes du graphe (lignes, degrés physiques) sont indépendantes.
//! On résout Gᵀλ ≃ f par QR ; les rangs déficients gardent la solution SVD
//! de norme minimale, sur Gᵀ directement, sans mettre son conditionnement au carré.
use crate::{DMatrix, DVector, Elem, Modele, Triplets};

mod qr_creux;
use qr_creux::{QrSymbolique, Resolution};
mod equivalentes;
pub(crate) use equivalentes::Equivalentes;

// Le QR dense de composantes bandées coûtait plus que le creux dès 64
// colonnes. Ce seuil est propre au projecteur, distinct du seuil LU.
const QR_DENSE_MAX: usize = 32;

pub(crate) struct Gradient {
    pub trip: Triplets,
    n: usize,
    m: usize,
}

#[derive(Default)]
struct CacheBloc {
    motif: Vec<(usize, usize)>,
    dim: (usize, usize),
    qr: Option<QrSymbolique>,
}

#[derive(Default)]
pub(crate) struct Projection {
    blocs: Vec<CacheBloc>,
}

struct Bloc {
    lignes: Vec<usize>,
    colonnes: Vec<usize>,
    // Transposée locale : ligne physique, colonne de multiplicateur.
    a: Triplets,
}

pub(crate) struct Composantes {
    parents: Vec<usize>,
    tailles: Vec<usize>,
}

impl Composantes {
    pub(crate) fn nouvelle(n: usize) -> Self {
        Self {
            parents: (0..n).collect(),
            tailles: vec![1; n],
        }
    }

    pub(crate) fn racine(&mut self, mut i: usize) -> usize {
        while self.parents[i] != i {
            self.parents[i] = self.parents[self.parents[i]];
            i = self.parents[i];
        }
        i
    }

    pub(crate) fn unit(&mut self, i: usize, j: usize) {
        let (mut i, mut j) = (self.racine(i), self.racine(j));
        if i == j {
            return;
        }
        if self.tailles[i] < self.tailles[j] {
            std::mem::swap(&mut i, &mut j);
        }
        self.parents[j] = i;
        self.tailles[i] += self.tailles[j];
    }
}

impl Gradient {
    pub fn est_fini(&self) -> bool {
        self.trip.iter().all(|t| t.val.is_finite())
    }

    pub fn transpose_fois(&self, lambda: &DVector<f64>) -> DVector<f64> {
        let mut y = DVector::zeros(self.n);
        for t in &self.trip {
            y[t.col] += t.val * lambda[t.row];
        }
        y
    }

    pub fn residu(&self, force: &DVector<f64>, lambda: &DVector<f64>) -> DVector<f64> {
        force - self.transpose_fois(lambda)
    }

    fn composantes(&self) -> Vec<Bloc> {
        let mut graphe = Composantes::nouvelle(self.n + self.m);
        let mut presentes = vec![false; self.m];
        for t in &self.trip {
            graphe.unit(t.row, self.m + t.col);
            presentes[t.row] = true;
        }
        let mut numero = vec![usize::MAX; self.n + self.m];
        let mut ligne_locale = vec![0; self.m];
        let mut colonne_locale = vec![0; self.n];
        let mut blocs: Vec<Bloc> = Vec::new();
        for i in 0..self.m {
            if !presentes[i] {
                continue;
            }
            let r = graphe.racine(i);
            if numero[r] == usize::MAX {
                numero[r] = blocs.len();
                blocs.push(Bloc {
                    lignes: Vec::new(),
                    colonnes: Vec::new(),
                    a: Triplets::new(),
                });
            }
            let b = &mut blocs[numero[r]];
            ligne_locale[i] = b.lignes.len();
            b.lignes.push(i);
        }
        for (j, local) in colonne_locale.iter_mut().enumerate() {
            let ib = numero[graphe.racine(self.m + j)];
            if ib == usize::MAX {
                continue;
            }
            *local = blocs[ib].colonnes.len();
            blocs[ib].colonnes.push(j);
        }
        for t in &self.trip {
            let ib = numero[graphe.racine(t.row)];
            blocs[ib].a.push(faer::sparse::Triplet {
                row: colonne_locale[t.col],
                col: ligne_locale[t.row],
                val: t.val,
            });
        }
        for b in &mut blocs {
            b.a.sort_by_key(|v| (v.col, v.row));
        }
        blocs
    }

    pub fn reactions(
        &self,
        f: &DVector<f64>,
        cache: &mut Projection,
    ) -> Result<DVector<f64>, String> {
        crate::fini(f.as_slice(), "forces statiques")?;
        if !self.est_fini() {
            return Err("gradient de contrainte non fini".into());
        }
        let blocs = self.composantes();
        cache.blocs.resize_with(blocs.len(), CacheBloc::default);
        let mut lambda = DVector::zeros(self.m);
        for (b, c) in blocs.iter().zip(&mut cache.blocs) {
            let rhs = DVector::from_iterator(b.colonnes.len(), b.colonnes.iter().map(|&i| f[i]));
            let x = b.resout(&rhs, c)?;
            for (j, &i) in b.lignes.iter().enumerate() {
                lambda[i] = x[j];
            }
        }
        crate::fini(lambda.as_slice(), "réactions statiques")?;
        Ok(lambda)
    }
}

impl Bloc {
    // Des colonnes identiques ne peuvent augmenter le rang. Cette borne
    // exacte évite un QR transposé voué au repli, notamment au premier appel
    // sur des liaisons dupliquées. Elle ne décide jamais qu'un rang est plein.
    fn assez_de_colonnes_distinctes(&self) -> bool {
        let mut distinctes = std::collections::HashSet::new();
        let mut start = 0;
        while start < self.a.len() {
            let col = self.a[start].col;
            let mut end = start + 1;
            while end < self.a.len() && self.a[end].col == col {
                end += 1;
            }
            let key: Vec<_> = self.a[start..end]
                .iter()
                .map(|t| (t.row, t.val.to_bits()))
                .collect();
            distinctes.insert(key);
            if distinctes.len() >= self.colonnes.len() {
                return true;
            }
            start = end;
        }
        false
    }

    fn dense(&self) -> DMatrix<f64> {
        let mut a = DMatrix::zeros(self.colonnes.len(), self.lignes.len());
        for t in &self.a {
            a[(t.row, t.col)] += t.val;
        }
        a
    }

    fn resout(&self, rhs: &DVector<f64>, cache: &mut CacheBloc) -> Result<DVector<f64>, String> {
        let (nr, nc) = (self.colonnes.len(), self.lignes.len());
        if nr == 1 && nc == 1 {
            return Ok(DVector::from_element(1, rhs[0] / self.a[0].val));
        }
        let norme = self.a.iter().fold(0.0_f64, |s, v| s.hypot(v.val));
        let seuil_rang = f64::EPSILON * nr.max(nc) as f64 * norme;
        if nr.min(nc) > QR_DENSE_MAX && (nr >= nc || self.assez_de_colonnes_distinctes()) {
            let motif: Vec<_> = self.a.iter().map(|t| (t.row, t.col)).collect();
            if cache.dim != (nr, nc) || cache.motif != motif {
                cache.qr = None;
                cache.dim = (nr, nc);
                cache.motif = motif;
            }
            let a =
                faer::sparse::SparseColMat::<usize, f64>::try_new_from_triplets(nr, nc, &self.a)
                    .map_err(|e| format!("contraintes creuses : {e:?}"))?;
            let at = a
                .transpose()
                .to_col_major()
                .map_err(|e| format!("transposée creuse : {e:?}"))?;
            // Un bloc large de rang plein se résout en norme minimale
            // par le QR de sa transposée, sans SVD ni équations normales.
            let (b, bt, sens) = if nr >= nc {
                (&a, &at, Resolution::MoindresCarres)
            } else {
                (&at, &a, Resolution::NormeMinimale)
            };
            if cache.qr.is_none() {
                cache.qr = Some(QrSymbolique::nouveau(b, bt)?);
            }
            if let Some(x) = cache
                .qr
                .as_ref()
                .unwrap()
                .resout(bt, rhs, seuil_rang, sens)?
            {
                return Ok(x);
            }
        }
        let a = self.dense();
        if nr >= nc {
            let qr = a.clone().col_piv_qr();
            let r = qr.r();
            if (0..nc).all(|i| r[(i, i)].abs() > seuil_rang) {
                let mut y = rhs.clone();
                qr.q_tr_mul(&mut y);
                if let Some(mut x) = r.solve_upper_triangular(&y.rows(0, nc).into_owned()) {
                    qr.p().inv_permute_rows(&mut x);
                    return Ok(x);
                }
            }
        }
        crate::moindres_carres::resout(a, rhs, seuil_rang)
    }
}

impl Modele {
    /// Même champ que phi_g : toutes les liaisons, lignes K nulles, sans
    /// appliquer le masque d'activité de la dynamique à la statique.
    pub(crate) fn phi_g_locale(&self, t: f64) -> Result<(DVector<f64>, Gradient), String> {
        let mut phi = DVector::zeros(self.m());
        let mut trip = Triplets::new();
        let mut row = 0;
        for e in &self.elems {
            let (p, ge) = match e {
                Elem::L(l) => l.phi_g(&self.corps, t),
                Elem::D(d) => d.phi_g(&self.corps),
                Elem::E(en) => en.phi_g(&self.corps)?,
                Elem::V(vi) => vi.phi_g(&self.corps)?,
                Elem::C(c) => c.phi_g(&self.corps),
                Elem::K(_) => (DVector::zeros(1), DMatrix::zeros(1, 12)),
            };
            for i in 0..e.n() {
                phi[row + i] = p[i];
                for (kb, cb) in e.corps().iter().enumerate() {
                    if let Some(ib) = cb {
                        for j in 0..6 {
                            let val = ge[(i, 6 * kb + j)];
                            if val != 0.0 {
                                trip.push(faer::sparse::Triplet {
                                    row: row + i,
                                    col: 6 * ib + j,
                                    val,
                                });
                            }
                        }
                    }
                }
            }
            row += e.n();
        }
        crate::analyse::regroupe(&mut trip);
        // Comme phi_g, laisser le solveur traiter les valeurs non finies :
        // un essai de recherche linéaire doit pouvoir être raccourci.
        Ok((
            phi,
            Gradient {
                trip,
                n: self.n(),
                m: self.m(),
            },
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Corps, Liaison, M3, V3};

    fn gradient(n: usize, m: usize, entries: &[(usize, usize, f64)]) -> Gradient {
        let mut trip = entries
            .iter()
            .map(|&(row, col, val)| faer::sparse::Triplet { row, col, val })
            .collect();
        crate::analyse::regroupe(&mut trip);
        Gradient { trip, n, m }
    }

    #[test]
    fn rang_proche_et_solution_minimale() {
        for delta in [0.0, 1e-12, 1e-9, 1e-6, 1.0] {
            let g = gradient(3, 2, &[(0, 0, 1.0), (1, 0, 1.0), (1, 1, delta)]);
            let f = DVector::from_vec(vec![2.0, 1.0, 3.0]);
            let lambda = g.reactions(&f, &mut Projection::default()).unwrap();
            let reference = if delta == 0.0 {
                vec![1.0, 1.0]
            } else {
                vec![2.0 - 1.0 / delta, 1.0 / delta]
            };
            let reference = DVector::from_vec(reference);
            assert!(
                (&lambda - &reference).norm() < 1e-14 * reference.norm().max(1.0),
                "{delta}: {lambda:?}"
            );
            let residual = g.residu(&f, &lambda);
            let reference = DVector::from_vec(vec![0.0, if delta == 0.0 { 1.0 } else { 0.0 }, 3.0]);
            assert!((residual - reference).amax() < 1e-12, "delta {delta}");
        }
        let g = gradient(2, 3, &[(0, 0, 1.0), (1, 1, 1.0), (2, 0, 1.0), (2, 1, 1.0)]);
        let f = DVector::from_vec(vec![2.0, 3.0]);
        let l = g.reactions(&f, &mut Projection::default()).unwrap();
        assert!((l - DVector::from_vec(vec![1.0 / 3.0, 4.0 / 3.0, 5.0 / 3.0])).amax() < 2e-15);
    }

    #[test]
    fn qr_creux_permutation_valeurs_et_motif() {
        let nc = 180;
        let mut cache = Projection::default();
        for changement in 0..3 {
            let mut terms = Vec::new();
            let x = DVector::from_fn(nc, |i, _| (i as f64 + 0.7).sin());
            let mut f = DVector::zeros(nc + 3);
            for i in 0..nc {
                let d = 2.0 + 0.01 * i as f64 + 0.1 * changement as f64;
                let j = i + 1 + usize::from(changement == 2);
                terms.extend([(i, i, d), (i, j, -0.25)]);
                f[i] += d * x[i];
                f[j] -= 0.25 * x[i];
            }
            f[nc + 2] = 3.0;
            let g = gradient(nc + 3, nc, &terms);
            let got = g.reactions(&f, &mut cache).unwrap();
            assert!((&got - x).amax() < 2e-14);
            let mut rf = g.residu(&f, &got);
            assert_eq!(rf[nc + 2], 3.0);
            rf[nc + 2] = 0.0;
            assert!(rf.amax() < 2e-14);
            if changement < 2 {
                assert!(cache.blocs[0].qr.is_some());
            }
        }
        let g = gradient(3, 2, &[]);
        assert_eq!(
            g.reactions(&DVector::from_element(3, 1.0), &mut cache)
                .unwrap(),
            DVector::zeros(2)
        );
        assert!(cache.blocs.is_empty());
    }

    #[test]
    fn qr_creux_redondant_refuse_avant_division() {
        let nc = 180;
        let mut terms = Vec::new();
        for j in 0..nc + 1 {
            terms.extend([(0, j, 1.0), (nc - 1, j, 1.0)]);
        }
        for i in 1..nc - 1 {
            terms.push((i, i, 1.0));
        }
        let g = gradient(nc + 1, nc, &terms);
        let reference = DVector::from_fn(nc, |i, _| {
            if i == 0 || i == nc - 1 {
                0.5
            } else {
                (i as f64).sin()
            }
        });
        let f = DVector::from_fn(nc + 1, |i, _| {
            if i > 0 && i < nc - 1 {
                1.0 + (i as f64).sin()
            } else {
                1.0
            }
        });
        let mut cache = Projection::default();
        let got = g.reactions(&f, &mut cache).unwrap();
        assert!(cache.blocs[0].qr.is_some());
        assert!((&got - reference).amax() < 2e-12, "norme minimale perdue");
        assert!(g.residu(&f, &got).amax() < 2e-12);
    }

    #[test]
    fn qr_transpose_norme_minimale_et_permutations() {
        let mut cache = Projection::default();
        // G est haut. Choisir λ dans son image garantit la norme minimale
        // pour f = Gᵀ λ ; les coefficients changent, puis le motif devient dense.
        for changement in 0..3 {
            let dense = changement == 2;
            let (n, m) = (72, 101);
            let mut terms = Vec::new();
            let mut reference = DVector::zeros(m);
            let y = DVector::from_fn(n, |i, _| (0.37 * i as f64).cos());
            for i in 0..m {
                for j in 0..n {
                    let diag = if i == (17 * j) % m { 3.0 } else { 0.0 };
                    let band = if (i + j) % m == 0 || (i + j + 1) % m == 0 {
                        0.25 + 0.1 * changement as f64
                    } else {
                        0.0
                    };
                    let val = diag
                        + band
                        + if dense {
                            0.05 * ((i + 3 * j) as f64).sin()
                        } else {
                            0.0
                        };
                    if val != 0.0 {
                        terms.push((i, j, val));
                        reference[i] += val * y[j];
                    }
                }
            }
            let mut f = DVector::zeros(n);
            for &(i, j, val) in &terms {
                f[j] += val * reference[i];
            }
            let g = gradient(n, m, &terms);
            let got = g.reactions(&f, &mut cache).unwrap();
            assert!(cache.blocs.iter().any(|b| b.qr.is_some()));
            assert!((&got - reference).amax() < 2e-13);
            assert!(g.residu(&f, &got).amax() < 5e-13);
        }
    }

    #[test]
    fn qr_transpose_deficient_conserve_residu_libre() {
        let (n, m) = (80, 121);
        let mut terms = Vec::new();
        for i in 0..m {
            terms.extend([(i, 0, 1.0), (i, n - 1, 1.0)]);
            if i > 0 && i < n - 1 {
                terms.push((i, i, 2.0));
            }
        }
        let reference = DVector::from_fn(m, |i, _| {
            1.0 + if i > 0 && i < n - 1 {
                (i as f64).sin()
            } else {
                0.0
            }
        });
        let mut f = DVector::zeros(n);
        for &(i, j, val) in &terms {
            f[j] += val * reference[i];
        }
        // Ces forces opposées appartiennent exactement au noyau de G.
        f[0] += 1.0;
        f[n - 1] -= 1.0;
        let g = gradient(n, m, &terms);
        let mut cache = Projection::default();
        let got = g.reactions(&f, &mut cache).unwrap();
        // Les colonnes identiques suffisent à prouver le défaut de rang.
        assert!(cache.blocs[0].qr.is_none());
        assert!((&got - reference).amax() < 3e-12);
        let r = g.residu(&f, &got);
        for i in 0..n {
            let expected = if i == 0 {
                1.0
            } else if i == n - 1 {
                -1.0
            } else {
                0.0
            };
            assert!((r[i] - expected).abs() < 3e-12);
        }
    }

    #[test]
    fn gradient_local_alias_et_masque_dynamique() {
        let mut m = Modele::new(V3::zeros());
        for i in 0..3 {
            m.corps.push(Corps {
                nom: i.to_string(),
                m: 1.0,
                j: M3::identity(),
                r: V3::new(i as f64, 0.2, -0.3),
                r_bas: V3::new(1e-17, -1e-17, 2e-17),
                rot: crate::expm(&V3::new(0.1, 0.2 * i as f64, -0.3)),
                v: V3::zeros(),
                w: V3::zeros(),
            });
        }
        for (a, b) in [(Some(2), Some(0)), (None, Some(1)), (Some(0), Some(0))] {
            m.elems.push(Elem::L(Liaison {
                nom: "j".into(),
                a,
                b,
                pa: V3::new(0.1, 0.2, 0.3),
                pb: V3::new(0.2, -0.1, 0.3),
                ra: M3::identity(),
                rb: M3::identity(),
                bt: vec![0, 2],
                br: vec![1],
                cible_t: None,
                cible_r: None,
                nh: false,
            }));
        }
        m.actif = vec![false; m.m()];
        let (p, dense) = m.phi_g(0.0).unwrap();
        let (q, sparse) = m.phi_g_locale(0.0).unwrap();
        assert_eq!(p, q);
        let mut actual = DMatrix::zeros(m.m(), m.n());
        for t in &sparse.trip {
            actual[(t.row, t.col)] += t.val;
        }
        assert_eq!(actual, dense);
    }
}
