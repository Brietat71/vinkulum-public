//! Accélération de Gauss pour les contraintes redondantes.
//!
//! M = L Lᵀ, W = L⁻¹ Gᵀ. La SVD de W donne directement les espaces
//! admissible et de réaction, sans former G M⁻¹ Gᵀ ni sa condition au carré.
//! Les composantes indépendantes sont traitées séparément. Ce repli dense
//! par composante conserve les lignes actives et minimise la norme de λ.
use crate::{contraintes::Composantes, DMatrix, DVector, Triplets, M3, V3};

mod contact;
pub(crate) use contact::complementarite;

/// Copie des contributions réellement présentées au solveur, pour le
/// certificat facultatif. Aucune reconstruction depuis une trajectoire.
pub(crate) struct TemoinLineaire {
    pub a: Triplets,
    pub b: DVector<f64>,
    pub x: DVector<f64>,
    pub g_complet: Vec<(usize, usize, f64)>,
    pub c_complet: DVector<f64>,
    pub actives: Vec<bool>,
}

pub(crate) fn lineaire(a: &Triplets, n: usize, b: &DVector<f64>) -> Result<DVector<f64>, String> {
    let dim = b.len();
    let inactives = a
        .iter()
        .filter(|t| t.row >= n && t.col == t.row && t.val != 0.)
        .count();
    if dim == 0 {
        return Ok(DVector::zeros(0));
    }
    let direct = if dim - n - inactives > n {
        None // Plus de contraintes que de degrés : LU nécessairement singulier.
    } else if dim <= crate::DENSE_MAX {
        Some(crate::Facto::dense_seule(a, dim))
    } else {
        let mut symbolique = None;
        crate::Facto::creux(a, dim, &mut symbolique).ok()
    };
    if let Some(f) = direct {
        let mut x = f.solve(b);
        // Une accélération presque nulle doit néanmoins respecter G a = c.
        // Le résidu compensé corrige notamment le pendule vertical URDF,
        // sans ajouter un plancher arbitraire à la tolérance des contraintes.
        for essai in 0..4 {
            if accepte(a, n, b, &x) {
                return Ok(x);
            }
            if essai == 3 || x.iter().any(|v| !v.is_finite()) {
                break;
            }
            x += f.solve(&crate::raffinement::residu(a, b, &x));
        }
    }
    resout(a, n, b)
}

#[derive(Default)]
struct Bloc {
    masses: Vec<usize>,
    contraintes: Vec<usize>,
    gt: Triplets,
}

fn descend(l: &M3, b: V3) -> V3 {
    let x = b.x / l[(0, 0)];
    let y = (b.y - l[(1, 0)] * x) / l[(1, 1)];
    let z = (b.z - l[(2, 0)] * x - l[(2, 1)] * y) / l[(2, 2)];
    V3::new(x, y, z)
}

fn remonte(l: &M3, b: V3) -> V3 {
    let z = b.z / l[(2, 2)];
    let y = (b.y - l[(2, 1)] * z) / l[(1, 1)];
    let x = (b.x - l[(1, 0)] * y - l[(2, 0)] * z) / l[(0, 0)];
    V3::new(x, y, z)
}

/// Contrôle en arithmétique compensée, séparément pour la dynamique et
/// pour chaque contrainte. La norme du second membre global ne doit pas
/// permettre à une force importante de masquer une accélération imposée.
pub(crate) fn accepte(a: &Triplets, n: usize, b: &DVector<f64>, x: &DVector<f64>) -> bool {
    if !x.iter().all(|v| v.is_finite()) {
        return false;
    }
    let r = crate::raffinement::residu(a, b, x);
    let mut echelle = b.map(f64::abs);
    let accel = x.rows(0, n).amax();
    for t in a {
        echelle[t.row] += if t.row >= n && t.col < n {
            t.val.abs() * accel
        } else {
            (t.val * x[t.col]).abs()
        };
    }
    let physique = echelle.rows(0, n).amax();
    // Borne d'erreur arrière, pas une borne sur l'erreur de solution.
    let tol = 128. * f64::EPSILON * (b.len() as f64).max(1.);
    r.iter().enumerate().all(|(i, &v)| {
        let s = if i < n { physique } else { echelle[i] };
        v.is_finite() && s.is_finite() && v.abs() <= tol * s
    })
}

pub(crate) fn resout(a: &Triplets, n: usize, b: &DVector<f64>) -> Result<DVector<f64>, String> {
    let m = b.len() - n;
    if n % 3 != 0 {
        return Err("initialisation : blocs de masse incomplets".into());
    }
    let nb = n / 3;
    let mut masses = vec![M3::zeros(); nb];
    let mut diagonales = DVector::<f64>::zeros(m);
    let mut graphe = Composantes::nouvelle(nb + m);
    for t in a {
        if t.row < n && t.col < n {
            if t.row / 3 != t.col / 3 && t.val != 0. {
                return Err("initialisation : masse hors blocs 3×3".into());
            }
            if t.row / 3 == t.col / 3 {
                masses[t.row / 3][(t.row % 3, t.col % 3)] += t.val;
            }
        } else if t.row < n && t.val != 0. {
            graphe.unit(t.row / 3, nb + t.col - n);
        } else if t.row >= n && t.col >= n {
            if t.row != t.col && t.val != 0. {
                return Err("initialisation : couplage inattendu entre multiplicateurs".into());
            }
            diagonales[t.row - n] += t.val;
        }
    }
    let facteurs: Vec<M3> = masses
        .into_iter()
        .map(|mass| {
            mass.cholesky()
                .map(|c| c.l())
                .ok_or_else(|| "initialisation : masse non définie positive".to_string())
        })
        .collect::<Result<_, _>>()?;
    let mut blocs: Vec<Bloc> = Vec::new();
    let mut numero = vec![usize::MAX; nb + m];
    let mut local = vec![0; nb + m];
    for i in 0..nb + m {
        let racine = graphe.racine(i);
        if numero[racine] == usize::MAX {
            numero[racine] = blocs.len();
            blocs.push(Bloc::default());
        }
        let bloc = &mut blocs[numero[racine]];
        let indices = if i < nb {
            &mut bloc.masses
        } else {
            &mut bloc.contraintes
        };
        local[i] = indices.len();
        indices.push(if i < nb { i } else { i - nb });
    }
    for t in a.iter().filter(|t| t.row < n && t.col >= n && t.val != 0.) {
        let ib = numero[graphe.racine(t.row / 3)];
        blocs[ib].gt.push(faer::sparse::Triplet {
            row: 3 * local[t.row / 3] + t.row % 3,
            col: local[nb + t.col - n],
            val: t.val,
        });
    }
    let mut x = DVector::zeros(b.len());
    for bloc in blocs {
        let nr = 3 * bloc.masses.len();
        let nc = bloc.contraintes.len();
        if nr == 0 {
            for &j in &bloc.contraintes {
                x[n + j] = if diagonales[j] != 0. {
                    b[n + j] / diagonales[j]
                } else {
                    0.
                };
                if diagonales[j] == 0. && b[n + j] != 0. {
                    return Err(
                        "initialisation : accélérations de contraintes incompatibles".into(),
                    );
                }
            }
            continue;
        }
        if bloc.contraintes.iter().any(|&j| diagonales[j] != 0.) {
            return Err("initialisation : multiplicateur inactif couplé à une masse".into());
        }
        let mut z0 = DVector::zeros(nr);
        let mut w = DMatrix::zeros(nr, nc);
        for t in &bloc.gt {
            w[(t.row, t.col)] += t.val;
        }
        for (i, &j) in bloc.masses.iter().enumerate() {
            let l = &facteurs[j];
            z0.rows_mut(3 * i, 3)
                .copy_from(&descend(l, b.fixed_rows::<3>(3 * j).into_owned()));
            for c in 0..nc {
                let v = descend(l, w.fixed_view::<3, 1>(3 * i, c).into_owned());
                w.fixed_view_mut::<3, 1>(3 * i, c).copy_from(&v);
            }
        }
        let z = if nc == 0 {
            z0
        } else {
            let norme = w.iter().fold(0.0_f64, |s, &v| s.hypot(v));
            let seuil = f64::EPSILON * nr.max(nc) as f64 * norme;
            let s = crate::svd_sure(w)?;
            let u = s.u.unwrap();
            let vt = s.v_t.unwrap();
            let rang = s.singular_values.iter().take_while(|&&v| v > seuil).count();
            let ur = u.columns(0, rang);
            let vr = vt.rows(0, rang);
            let c = DVector::from_iterator(nc, bloc.contraintes.iter().map(|&j| b[n + j]));
            let p = ur.transpose() * &z0;
            let mut particulier = vr * &c;
            let incompatible = &c - vr.transpose() * &particulier;
            if incompatible.amax() > 128. * f64::EPSILON * nr.max(nc) as f64 * c.amax() {
                return Err("initialisation : accélérations de contraintes incompatibles avec le rang numérique".into());
            }
            for i in 0..rang {
                particulier[i] /= s.singular_values[i];
            }
            // Construire le sous-espace libre séparément. Soustraire Wλ
            // perdrait l'accélération lorsque deux contraintes sont proches.
            let libre = if rang == nr {
                DVector::zeros(nr)
            } else if u.ncols() == nr {
                let nul = u.columns(rang, nr - rang);
                nul * (nul.transpose() * &z0)
            } else {
                let mut libre = &z0 - ur * &p;
                libre -= ur * (ur.transpose() * &libre);
                libre
            };
            let z = libre + ur * &particulier;
            let mut react = p - particulier;
            for i in 0..rang {
                react[i] /= s.singular_values[i];
            }
            let lambda = vr.transpose() * react;
            for (i, &j) in bloc.contraintes.iter().enumerate() {
                x[n + j] = lambda[i];
            }
            z
        };
        for (i, &j) in bloc.masses.iter().enumerate() {
            let acc = remonte(&facteurs[j], z.fixed_rows::<3>(3 * i).into_owned());
            x.rows_mut(3 * j, 3).copy_from(&acc);
        }
    }
    if !accepte(a, n, b, &x) {
        return Err("initialisation : équilibre ou accélérations de contraintes incompatibles avec le rang numérique".into());
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn triplets(m: &DMatrix<f64>) -> Triplets {
        (0..m.ncols())
            .flat_map(|j| {
                (0..m.nrows()).filter_map(move |i| {
                    (m[(i, j)] != 0.).then_some(faer::sparse::Triplet {
                        row: i,
                        col: j,
                        val: m[(i, j)],
                    })
                })
            })
            .collect()
    }

    #[test]
    fn gauss_masse_anisotrope_redondances_et_composantes() {
        // a choisi ; λ = G y est dans im G, donc de norme minimale.
        // Deux blocs indépendants, une rotation d'inertie et un λ inactif.
        for scale in [1e-80, 1., 1e80] {
            let n = 12;
            let g = DMatrix::from_row_slice(
                4,
                6,
                &[
                    1., 2., 0., 0., 0., 0., 0., 0., 1., -1., 0., 0., 1., 2., 1., -1., 0., 0., 2.,
                    4., -1., 1., 0., 0.,
                ],
            );
            let mut a = DMatrix::zeros(21, 21);
            let q = crate::expm(&V3::new(0.3, -0.2, 0.5));
            for i in 0..4 {
                let mass = q * M3::from_diagonal(&V3::new(0.01, 2., 13.)) * q.transpose() * scale;
                a.fixed_view_mut::<3, 3>(3 * i, 3 * i).copy_from(&mass);
            }
            let mut attendu = DVector::zeros(21);
            for i in 0..n {
                attendu[i] = (0.17 * (i + 1) as f64).sin();
            }
            for k in 0..2 {
                a.view_mut((n + 4 * k, 6 * k), (4, 6)).copy_from(&g);
                a.view_mut((6 * k, n + 4 * k), (6, 4))
                    .copy_from(&g.transpose());
                let lambda = &g * attendu.rows(6 * k, 6) * scale;
                attendu.rows_mut(n + 4 * k, 4).copy_from(&lambda);
            }
            a[(20, 20)] = 1.;
            attendu[20] = 0.7 * scale;
            let b = &a * &attendu;
            let x = resout(&triplets(&a), n, &b).unwrap();
            assert!(
                (x.rows(0, n) - attendu.rows(0, n)).amax() < 2e-11,
                "{scale}: {x:?}"
            );
            assert!(
                ((x.rows(n, 9) - attendu.rows(n, 9)) / scale).amax() < 2e-11,
                "{scale}"
            );
            let mut faux = b.clone();
            faux[n + 2] += 0.1;
            assert!(resout(&triplets(&a), n, &faux)
                .unwrap_err()
                .contains("incompatibles"));
        }
    }

    #[test]
    fn gauss_contraintes_proches_acceleration_et_reactions() {
        // Une matrice normale perd δ² ; les deux contraintes fixent pourtant
        // a_x = a_y = 0. Une troisième copie les rend redondantes.
        for delta in [1e-4, 1e-8, 1e-12] {
            let mut a = DMatrix::identity(6, 6);
            a.view_mut((3, 3), (3, 3)).fill(0.);
            let g = DMatrix::from_row_slice(3, 3, &[1., 0., 0., 1., delta, 0., 1., 0., 0.]);
            a.view_mut((3, 0), (3, 3)).copy_from(&g);
            a.view_mut((0, 3), (3, 3)).copy_from(&g.transpose());
            let b = DVector::from_vec(vec![0., 1., 0.3, 0., 0., 0.]);
            let x = resout(&triplets(&a), 3, &b).unwrap();
            assert!(x[0].abs().max(x[1].abs()) < 1e-14, "{delta}: {x:?}");
            assert!((x[2] - 0.3).abs() < 1e-14);
            for (i, v) in [-0.5, 1., -0.5].iter().enumerate() {
                assert!((x[3 + i] * delta - v).abs() < 2e-14, "{delta}: {x:?}");
            }
        }
    }

    #[test]
    fn gauss_masse_indefinie_et_contrainte_nulle() {
        let mut a = DMatrix::identity(4, 4);
        a[(3, 3)] = 0.;
        let mut b = DVector::from_vec(vec![1., 2., 3., 0.]);
        let x = resout(&triplets(&a), 3, &b).unwrap();
        assert_eq!(x, b);
        b[3] = 1e-100;
        assert!(resout(&triplets(&a), 3, &b)
            .unwrap_err()
            .contains("incompatibles"));
        a[(0, 0)] = -1.;
        assert!(resout(&triplets(&a), 3, &b)
            .unwrap_err()
            .contains("positive"));
    }
}
