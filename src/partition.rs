//! Restriction algébrique des systèmes multi-rythmes. Les états restent
//! globaux ; seuls les triplets factorisés et le second membre sont locaux.
use crate::*;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Indices {
    pub global: Vec<usize>,
    local: Vec<usize>,
}

impl Indices {
    pub fn du_modele(m: &Modele) -> Option<Self> {
        if m.gel.is_empty() {
            return None;
        }
        let n = m.n();
        Some(Self::nouveaux((0..n + m.m()).map(|i| {
            if i < n {
                !m.gele(i / 6)
            } else {
                m.actif.get(i - n).copied().unwrap_or(true)
            }
        })))
    }

    fn nouveaux(garde: impl Iterator<Item = bool>) -> Self {
        let mut global = Vec::new();
        let local = garde
            .enumerate()
            .map(|(i, garde)| {
                if garde {
                    global.push(i);
                    global.len() - 1
                } else {
                    usize::MAX
                }
            })
            .collect();
        Self { global, local }
    }

    pub fn dim(&self) -> usize {
        self.global.len()
    }
}

pub(crate) struct Reduction {
    pub indices: Indices,
    // A_le : ligne locale conservée, colonne globale éliminée.
    couplage: Triplets,
}

impl Reduction {
    pub fn nouvelle(indices: Indices, trip: Triplets) -> (Self, Triplets) {
        let mut reduit = Vec::with_capacity(trip.len());
        let mut couplage = Vec::new();
        for t in trip {
            let row = indices.local[t.row];
            if row == usize::MAX {
                continue;
            }
            let col = indices.local[t.col];
            if col == usize::MAX {
                couplage.push(faer::sparse::Triplet {
                    row,
                    col: t.col,
                    val: t.val,
                });
            } else {
                reduit.push(faer::sparse::Triplet {
                    row,
                    col,
                    val: t.val,
                });
            }
        }
        (Self { indices, couplage }, reduit)
    }

    /// b_l − A_le x_e : x_e est connu, mais n'est pas nécessairement nul.
    pub fn rhs(&self, b: &DVector<f64>, impose: &DVector<f64>) -> DVector<f64> {
        let mut rhs = DVector::from_iterator(
            self.indices.dim(),
            self.indices.global.iter().map(|&i| b[i]),
        );
        for t in &self.couplage {
            rhs[t.row] -= t.val * impose[t.col];
        }
        rhs
    }

    pub fn etend(&self, x: &DVector<f64>, mut impose: DVector<f64>) -> DVector<f64> {
        for (j, &i) in self.indices.global.iter().enumerate() {
            impose[i] = x[j];
        }
        impose
    }

    /// Jacobien global pour le diagnostic optionnel de rotation. Les lignes
    /// éliminées de NEWTON sont des identités (l'acc_init n'utilise pas ceci).
    pub fn triplets_globaux(&self, trip: &Triplets) -> Triplets {
        let mut out: Triplets = trip
            .iter()
            .map(|t| faer::sparse::Triplet {
                row: self.indices.global[t.row],
                col: self.indices.global[t.col],
                val: t.val,
            })
            .collect();
        out.extend(self.couplage.iter().map(|t| faer::sparse::Triplet {
            row: self.indices.global[t.row],
            col: t.col,
            val: t.val,
        }));
        out.extend(
            self.indices
                .local
                .iter()
                .enumerate()
                .filter(|(_, j)| **j == usize::MAX)
                .map(|(i, _)| faer::sparse::Triplet {
                    row: i,
                    col: i,
                    val: 1.0,
                }),
        );
        out
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::cell::RefCell;

    thread_local! { pub static FACTOS: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) }; }

    fn body(i: usize) -> Corps {
        Corps {
            r_bas: crate::V3::zeros(),
            nom: format!("c{i}"),
            m: 1.2,
            j: M3::from_diagonal(&V3::new(1.0, 2.0, 3.0)),
            r: V3::new(i as f64, 0.0, 0.0),
            rot: M3::identity(),
            v: V3::zeros(),
            w: V3::new(0.1, 0.2, 0.3),
        }
    }

    #[test]
    fn elimination_non_nulle_dense_creuse_et_rotation() {
        for nb in [8, 80] {
            let dim = 6 * nb;
            let indices = Indices::nouveaux((0..dim).map(|i| (i / 6) % 2 == 0));
            let mut trip = Vec::new();
            let tp = |row, col, val| faer::sparse::Triplet { row, col, val };
            for i in 0..dim {
                let garde = indices.local[i] != usize::MAX;
                trip.push(tp(
                    i,
                    i,
                    if garde {
                        2.0 + i as f64 / dim as f64
                    } else {
                        1.0
                    },
                ));
                if garde {
                    // Colonnes supprimées à second membre NON NUL + couplage
                    // entre blocs mobiles non contigus, matrice non symétrique.
                    trip.push(tp(i, (i + 6) % dim, 0.7));
                    trip.push(tp(i, (i + 12) % dim, -0.2));
                }
            }
            let mut cs: Vec<_> = (0..nb).map(body).collect();
            let full = Facto::dense_seule(&trip, dim);
            let mut reference = JacGarde::neuf(trip.clone(), Some(full), 0.01, dim, &cs, None);
            let (red, tr) = Reduction::nouvelle(indices, trip);
            let dr = red.indices.dim();
            let f = if dr <= DENSE_MAX {
                Facto::dense_seule(&tr, dr)
            } else {
                Facto::creux(&tr, dr, &mut None).unwrap()
            };
            let mut garde = JacGarde::neuf(tr, Some(f), 0.01, dim, &cs, Some(red));
            for (i, c) in cs.iter_mut().enumerate() {
                c.rot = expm(&(V3::new(0.2, 0.3, -0.1) * (i + 1) as f64));
            }
            let b = DVector::from_fn(dim, |i, _| (i as f64 + 0.4).sin());
            let (xr, sr) = garde.resout(&cs, &[], dim, &b);
            let (xf, sf) = reference.resout(&cs, &[], dim, &b);
            assert!(!sr && !sf);
            assert!((xr - xf).norm() < 1e-12);
            assert_eq!(garde.trip.iter().map(|t| t.row).max(), Some(dr - 1));
        }
    }

    #[test]
    fn acceleration_masses_gelees_et_partition_vide() {
        for gel in [vec![true, false, true], vec![true; 3], vec![false; 3]] {
            let mut m = Modele::new(V3::new(0.0, -9.81, 0.0));
            m.corps = (0..3).map(body).collect();
            m.corps[0].rot = expm(&V3::new(0.1, 0.3, 0.2));
            m.efforts.push((0, V3::x(), V3::new(2.0, -1.0, 3.0)));
            let reference = m.acc_init(0.0).unwrap().0;
            m.gel = gel.clone();
            FACTOS.with(|f| f.borrow_mut().clear());
            let x = m.acc_init(0.0).unwrap().0;
            assert!((x - reference).norm() < 1e-12);
            let dim = 6 * gel.iter().filter(|g| !**g).count();
            FACTOS.with(|f| assert_eq!(*f.borrow(), if dim == 0 { vec![] } else { vec![dim] }));
        }
    }

    #[test]
    fn repli_svd_reduit_et_reaction_redondante() {
        let cs: Vec<_> = (0..8).map(body).collect();
        let indices = Indices::nouveaux((0..50).map(|i| i < 48 && (i / 6) % 2 == 0 || i == 48));
        let trip: Triplets = (0..50)
            .filter(|&i| i != 2)
            .map(|i| faer::sparse::Triplet {
                row: i,
                col: i,
                val: 1.0,
            })
            .collect();
        let b = DVector::from_fn(50, |i, _| if i == 2 { 0.0 } else { 1.0 + i as f64 });
        let mut complet = JacGarde::neuf(
            trip.clone(),
            Some(Facto::dense_seule(&trip, 50)),
            0.01,
            50,
            &cs,
            None,
        );
        let (red, tr) = Reduction::nouvelle(indices, trip);
        let f = Facto::dense_seule(&tr, red.indices.dim());
        let mut reduit = JacGarde::neuf(tr, Some(f), 0.01, 50, &cs, Some(red));
        let (xr, svd) = reduit.resout(&cs, &[], 50, &b);
        let (xf, _) = complet.resout(&cs, &[], 50, &b);
        assert!(svd);
        assert!((xr - xf).norm() < 1e-11);
    }

    #[test]
    fn dimensions_chaine_avec_liaisons_et_cache() {
        for n in [8, 16, 32] {
            let mut m = Modele::new(V3::new(0.0, -9.81, 0.0));
            m.corps = (0..n + 2).map(body).collect();
            for i in 0..n {
                m.corps[i].w = V3::new(0.0, 0.0, 0.3);
                m.elems.push(Elem::L(Liaison {
                    nom: format!("p{i}"),
                    a: None,
                    b: Some(i),
                    pa: m.corps[i].r,
                    pb: V3::zeros(),
                    ra: M3::identity(),
                    rb: M3::identity(),
                    bt: vec![0, 1, 2],
                    br: vec![0, 1],
                    cible_t: None,
                    cible_r: None,
                    nh: false,
                }));
            }
            let rapides: Vec<_> = (0..n + 2).map(|i| i >= n).collect();
            FACTOS.with(|f| f.borrow_mut().clear());
            m.simule_multi(0.003, 0.001, 10, &rapides, 0.9, 1e-12, 20, |_, _| {})
                .unwrap();
            FACTOS.with(|f| {
                let f = f.borrow();
                assert!(f.contains(&12) && f.contains(&(11 * n)), "{f:?}");
                assert!(f.iter().all(|&d| d == 12 || d == 11 * n), "{f:?}");
                assert!(f.len() < 20, "cache perdu entre micro-pas : {f:?}");
            });
            assert!(m.phi_holonome_max(m.t).unwrap() < 1e-12);
            // Les réactions lentes survivent au dernier micro-pas rapide.
            for i in 0..n {
                assert!((m.lam[5 * i + 1] + 1.2 * 9.81).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn dimensions_reellement_factorisees_et_changement_partition() {
        let mut m = Modele::new(V3::new(0.0, -9.81, 0.0));
        m.corps = (0..8).map(body).collect();
        for rapides in [
            vec![true, false, true, false, true, false, false, false],
            vec![true; 8],
            vec![false; 8],
        ] {
            FACTOS.with(|f| f.borrow_mut().clear());
            m.simule_multi(m.t + 0.0025, 0.001, 2, &rapides, 0.9, 1e-12, 20, |_, _| {})
                .unwrap();
            let nr = 6 * rapides.iter().filter(|r| **r).count();
            FACTOS.with(|f| {
                let f = f.borrow();
                assert!(!f.is_empty());
                assert!(
                    f.iter().all(|&d| d > 0 && (d == nr || d == 48 - nr)),
                    "{f:?}"
                );
            });
            assert!(m.gel.is_empty() && m.gel_elem.is_empty());
        }
    }
}
