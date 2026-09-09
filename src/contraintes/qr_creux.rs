//! QR supernodal avec accès à la diagonale de R avant toute résolution.
//! APIs publiques de faer 0.24 : COLAMD, arbre d'élimination et réflecteurs.
//! Le chemin simplicial n'est pas utilisé sur les matrices de rang inconnu.
use crate::DVector;
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::perm::PermRef;
use faer::sparse::linalg::qr::supernodal::{self, SymbolicSupernodalQr};
use faer::sparse::linalg::{colamd, qr};
use faer::sparse::SparseColMat;

pub(super) enum Resolution {
    MoindresCarres,
    NormeMinimale,
}

struct Numerique {
    hrow: Vec<usize>,
    taille: Vec<usize>,
    nrows: Vec<usize>,
    ncols: Vec<usize>,
    r: Vec<f64>,
    hv: Vec<f64>,
    tau: Vec<f64>,
}

fn memoire(req: StackReq) -> Result<MemBuffer, String> {
    MemBuffer::try_new(req).map_err(|e| format!("mémoire QR des contraintes : {e:?}"))
}

pub(super) struct QrSymbolique {
    symbolique: SymbolicSupernodalQr<usize>,
    permutation: Vec<usize>,
    inverse: Vec<usize>,
}

impl QrSymbolique {
    pub fn nouveau(
        a: &SparseColMat<usize, f64>,
        at: &SparseColMat<usize, f64>,
    ) -> Result<Self, String> {
        let (nr, nc) = (a.nrows(), a.ncols());
        let mut permutation = vec![0; nc];
        let mut inverse = vec![0; nc];
        let mut tree = vec![0; nc];
        let mut post = vec![0; nc];
        let mut counts = vec![0; nc];
        let mut min_col = vec![0; nr];
        let req = StackReq::any_of(&[
            colamd::order_scratch::<usize>(nr, nc, a.compute_nnz()),
            qr::col_etree_scratch::<usize>(nr, nc),
            qr::postorder_scratch::<usize>(nc),
            qr::column_counts_aat_scratch::<usize>(nc, nr),
            supernodal::factorize_supernodal_symbolic_qr_scratch::<usize>(nr, nc),
        ]);
        let mut mem = memoire(req)?;
        colamd::order(
            &mut permutation,
            &mut inverse,
            a.symbolic(),
            Default::default(),
            MemStack::new(&mut mem),
        )
        .map_err(|e| format!("ordre QR des contraintes : {e:?}"))?;
        let p = PermRef::new_checked(&permutation, &inverse, nc);
        let etree = qr::col_etree(a.symbolic(), Some(p), &mut tree, MemStack::new(&mut mem));
        qr::postorder(&mut post, etree, MemStack::new(&mut mem));
        qr::column_counts_ata(
            &mut counts,
            &mut min_col,
            at.symbolic(),
            Some(p),
            etree,
            &post,
            MemStack::new(&mut mem),
        );
        let symbolique = supernodal::factorize_supernodal_symbolic_qr(
            a.symbolic(),
            Some(p),
            min_col,
            etree,
            &counts,
            MemStack::new(&mut mem),
            Default::default(),
        )
        .map_err(|e| format!("QR symbolique des contraintes : {e:?}"))?;
        Ok(Self {
            symbolique,
            permutation,
            inverse,
        })
    }

    /// None demande le repli avec révélation du rang et norme minimale.
    /// Le choix se fait sur R, avant une division par un pivot presque nul.
    pub fn resout(
        &self,
        at: &SparseColMat<usize, f64>,
        rhs: &DVector<f64>,
        seuil_rang: f64,
        sens: Resolution,
    ) -> Result<Option<DVector<f64>>, String> {
        let s = &self.symbolique;
        let h = s.householder();
        let mut hrow = vec![0; h.len_householder_row_idx()];
        let nh = h.len_householder_row_idx() + h.n_supernodes();
        let mut taille = vec![0; nh];
        let mut nrows = vec![0; nh];
        let mut ncols = vec![0; nh];
        let mut r = vec![0.0; s.R_adjoint().len_val()];
        let mut hv = vec![0.0; h.len_householder_val()];
        let mut tau = vec![0.0; h.len_tau_val()];
        let par = faer::Par::Seq;
        let mut mem = memoire(supernodal::factorize_supernodal_numeric_qr_scratch::<
            usize,
            f64,
        >(s, par, Default::default()))?;
        let p = PermRef::new_checked(&self.permutation, &self.inverse, self.inverse.len());
        let q = supernodal::factorize_supernodal_numeric_qr(
            &mut hrow,
            &mut taille,
            &mut nrows,
            &mut ncols,
            &mut r,
            &mut hv,
            &mut tau,
            at.as_ref(),
            Some(p),
            s,
            par,
            MemStack::new(&mut mem),
            Default::default(),
        );
        let rt = faer::sparse::linalg::cholesky::supernodal::SupernodalLltRef::new(
            s.R_adjoint(),
            q.R_val(),
        );
        for i in 0..s.R_adjoint().n_supernodes() {
            let block = rt.supernode(i).val();
            for j in 0..block.ncols() {
                let v: f64 = block[(j, j)];
                if !v.is_finite() || v.abs() <= seuil_rang {
                    return Ok(None);
                }
            }
        }
        if let Resolution::NormeMinimale = sens {
            let f = Numerique {
                hrow,
                taille,
                nrows,
                ncols,
                r,
                hv,
                tau,
            };
            let x = f.norme_minimale(self, rhs)?;
            return Ok(x.iter().all(|v| v.is_finite()).then_some(x));
        }
        let mut y = rhs.clone();
        let mut mem = memoire(s.solve_in_place_scratch::<f64>(1, par))?;
        q.solve_in_place_with_conj(
            faer::Conj::No,
            faer::MatMut::from_column_major_slice_mut(y.as_mut_slice(), rhs.len(), 1),
            par,
            MemStack::new(&mut mem),
        );
        let x = DVector::from_iterator(self.inverse.len(), self.inverse.iter().map(|&i| y[i]));
        Ok(x.iter().all(|v| v.is_finite()).then_some(x))
    }
}

impl Numerique {
    /// B P = Q R ; résoudre Bᵀ x = f avec x = Q R⁻ᵀ Pᵀ f.
    /// Une seule résolution triangulaire, puis Q appliqué par réflecteurs.
    fn norme_minimale(
        &self,
        qr: &QrSymbolique,
        rhs: &DVector<f64>,
    ) -> Result<DVector<f64>, String> {
        let s = &qr.symbolique;
        let l = faer::sparse::linalg::cholesky::supernodal::SupernodalLltRef::new(
            s.R_adjoint(),
            &self.r,
        );
        let mut y = DVector::from_iterator(rhs.len(), qr.permutation.iter().map(|&i| rhs[i]));
        // Rᵀ est stocké par colonnes supernodales. La partie hors diagonale
        // met à jour le second membre des panneaux encore à résoudre.
        for panel in 0..s.R_adjoint().n_supernodes() {
            let p = l.supernode(panel);
            let v = p.val();
            let start = p.start();
            for j in 0..v.ncols() {
                y[start + j] /= v[(j, j)];
                let z = y[start + j];
                for i in j + 1..v.ncols() {
                    y[start + i] -= v[(i, j)] * z;
                }
                for (i, &row) in p.pattern().iter().enumerate() {
                    y[row] -= v[(v.ncols() + i, j)] * z;
                }
            }
        }
        self.applique_q(s, &y)
    }

    /// Adjoint du parcours comprimé de Qᵀ : injecter chaque segment de y
    /// aux pivots de son panneau, puis remonter les réflecteurs en ordre
    /// inverse. Aucun Q explicite et aucun produit par la matrice d'origine.
    fn applique_q(
        &self,
        s: &SymbolicSupernodalQr<usize>,
        y: &DVector<f64>,
    ) -> Result<DVector<f64>, String> {
        use faer::linalg::householder::{
            apply_block_householder_sequence_on_the_left_in_place_scratch as scratch,
            apply_block_householder_sequence_on_the_left_in_place_with_conj as applique,
        };
        let h = s.householder();
        let rt = s.R_adjoint();
        // Les métadonnées numériques énumèrent les sous-panneaux en ordre
        // direct ; mémoriser seulement leur découpage pour le parcours inverse.
        let mut blocs = Vec::new();
        let mut bornes = vec![0];
        let mut index = 0;
        let mut max_lignes = 0;
        let mut max_bloc = 0;
        for panel in 0..h.n_supernodes() {
            let nr = h.col_ptr_for_householder_row_idx()[panel + 1]
                - h.col_ptr_for_householder_row_idx()[panel];
            let nv =
                h.col_ptr_for_householder_val()[panel + 1] - h.col_ptr_for_householder_val()[panel];
            let nc = nv / nr;
            max_lignes = max_lignes.max(nr);
            let mut start = 0;
            while start < nr.min(nc) {
                blocs.push((start, index));
                max_bloc = max_bloc.max(self.taille[index]);
                start += self.ncols[index];
                index += 1;
            }
            bornes.push(blocs.len());
        }
        let mut mem = memoire(scratch::<f64>(max_lignes, max_bloc, 1))?;
        let mut work = vec![0.0; max_lignes];
        let mut x = DVector::zeros(h.nrows());
        for panel in (0..h.n_supernodes()).rev() {
            let rows = &self.hrow[h.col_ptr_for_householder_row_idx()[panel]
                ..h.col_ptr_for_householder_row_idx()[panel + 1]];
            let values = &self.hv[h.col_ptr_for_householder_val()[panel]
                ..h.col_ptr_for_householder_val()[panel + 1]];
            let nc = values.len() / rows.len();
            let basis = faer::MatRef::from_column_major_slice(values, rows.len(), nc);
            let tv = &self.tau[h.col_ptr_for_tau_val()[panel]..h.col_ptr_for_tau_val()[panel + 1]];
            let tau_cols = nc.min(rows.len());
            let tau = faer::MatRef::from_column_major_slice(tv, tv.len() / tau_cols, tau_cols);
            for (i, &row) in rows.iter().enumerate() {
                work[i] = x[row];
            }
            let begin = rt.supernode_begin()[panel];
            for i in 0..rt.supernode_end()[panel] - begin {
                work[i] += y[begin + i];
            }
            for &(start, index) in blocs[bornes[panel]..bornes[panel + 1]].iter().rev() {
                let (nr, nc, bs) = (self.nrows[index], self.ncols[index], self.taille[index]);
                applique(
                    basis.submatrix(start, start, nr, nc),
                    tau.submatrix(0, start, bs, nc),
                    faer::Conj::No,
                    faer::MatMut::from_column_major_slice_mut(&mut work[start..start + nr], nr, 1),
                    faer::Par::Seq,
                    MemStack::new(&mut mem),
                );
            }
            for (i, &row) in rows.iter().enumerate() {
                x[row] = work[i];
            }
        }
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::sparse::Triplet;

    #[test]
    fn norme_minimale_directions_proches() {
        // Bᵀ x = f impose x0=2 et x0+δ(x1+x2)=1. La norme minimale
        // répartit exactement le dernier effort entre x1 et x2.
        // Le calcul B(R⁻¹R⁻ᵀf) perdrait x0 par soustraction de grands termes.
        let (m, n) = (74, 72);
        for delta in [1.0, 1e-3, 1e-6, 1e-9, 1e-12, 0.0] {
            let mut terms = vec![
                Triplet {
                    row: 0,
                    col: 0,
                    val: 1.0,
                },
                Triplet {
                    row: 0,
                    col: 1,
                    val: 1.0,
                },
                Triplet {
                    row: 1,
                    col: 1,
                    val: delta,
                },
                Triplet {
                    row: 2,
                    col: 1,
                    val: delta,
                },
            ];
            let mut rhs = DVector::zeros(n);
            rhs[0] = 2.0;
            rhs[1] = 1.0;
            for j in 2..n {
                terms.push(Triplet {
                    row: j + 1,
                    col: j,
                    val: 1.0,
                });
                rhs[j] = (j as f64).sin();
            }
            let b = SparseColMat::try_new_from_triplets(m, n, &terms).unwrap();
            let bt = b.transpose().to_col_major().unwrap();
            let qr = QrSymbolique::nouveau(&b, &bt).unwrap();
            let norm = terms.iter().fold(0.0_f64, |s, v| s.hypot(v.val));
            let x = qr
                .resout(
                    &bt,
                    &rhs,
                    f64::EPSILON * m as f64 * norm,
                    Resolution::NormeMinimale,
                )
                .unwrap();
            if delta == 0.0 {
                assert!(x.is_none());
                continue;
            }
            let x = x.expect("le QR de rang plein doit résoudre directement");
            assert!((x[0] - 2.0).abs() < 1e-14);
            for i in [1, 2] {
                assert!((x[i] * delta + 0.5).abs() < 2e-15);
            }
            for j in 2..n {
                assert!((x[j + 1] - rhs[j]).abs() < 1e-14);
            }
            assert_eq!(x[m - 1], 0.0);
        }
    }
}
