//! Sélection duale des contacts dans le problème de Gauss convexe.
//! On augmente un multiplicateur jusqu'à fermer sa contrainte OU jusqu'à
//! annuler une réaction active. Ce second événement permet de remplacer
//! une face par une autre même si leurs égalités sont incompatibles.
//! Les directions utilisent la projection de masse, pas une matrice normale.
use crate::{DVector, Triplets};

fn selection(a: &Triplets, b: &DVector<f64>, garde: &[bool]) -> (Triplets, DVector<f64>) {
    let mut trip: Triplets = a
        .iter()
        .filter(|t| garde[t.row] && garde[t.col])
        .copied()
        .collect();
    let mut rhs = b.clone();
    for (i, &oui) in garde.iter().enumerate() {
        if !oui {
            trip.push(faer::sparse::Triplet {
                row: i,
                col: i,
                val: 1.,
            });
            rhs[i] = 0.;
        }
    }
    (trip, rhs)
}

fn energie(a: &Triplets, n: usize, x: &DVector<f64>) -> f64 {
    let (mut haut, mut bas) = (0., 0.);
    for t in a.iter().filter(|t| t.row < n && t.col < n) {
        let (s, r) = crate::numerique::deux_sommes(haut, t.val * x[t.row] * x[t.col]);
        haut = s;
        bas += r;
    }
    (haut + bas).max(0.)
}

pub(crate) fn complementarite(
    a: &Triplets,
    n: usize,
    b: &DVector<f64>,
    contacts: &[(usize, bool)],
) -> Result<DVector<f64>, String> {
    if contacts.is_empty() {
        return super::lineaire(a, n, b);
    }
    let mut garde = vec![true; b.len()];
    for &(i, eligible) in contacts {
        garde[i] = eligible;
    }
    // Si toutes les réactions admissibles sont positives, les conditions
    // KKT du problème convexe sont déjà satisfaites : conserver le LU rapide.
    let (essai, rhs_essai) = selection(a, b, &garde);
    if let Ok(x) = super::lineaire(&essai, n, &rhs_essai) {
        if contacts.iter().all(|&(i, _)| x[i] >= 0.) {
            return Ok(x);
        }
    }
    for &(i, _) in contacts {
        garde[i] = false;
    }
    let (mut trip, rhs) = selection(a, b, &garde);
    let mut x = super::lineaire(&trip, n, &rhs)?;
    let tol = 128. * f64::EPSILON * b.len() as f64;
    let mut pivots = 0;
    let budget = (32 + 8 * contacts.len().saturating_mul(contacts.len())).min(10_000);
    loop {
        let residu = crate::raffinement::residu(a, b, &x);
        let accel = x.rows(0, n).amax();
        let mut echelle = b.map(f64::abs);
        for t in a.iter().filter(|t| t.row >= n && t.col < n) {
            echelle[t.row] += t.val.abs() * accel;
        }
        let candidat = contacts
            .iter()
            .filter(|&&(i, eligible)| eligible && !garde[i])
            .map(|&(i, _)| (i, -residu[i], echelle[i]))
            .filter(|&(_, v, s)| v > tol * s)
            .max_by(|a, b| (a.1 / a.2).total_cmp(&(b.1 / b.2)));
        let Some((p, _, _)) = candidat else {
            if !super::accepte(&trip, n, &selection(a, b, &garde).1, &x)
                || contacts.iter().any(|&(i, _)| x[i] < 0.)
            {
                return Err(
                    "initialisation contact : équilibre ou signe des réactions non vérifié".into(),
                );
            }
            return Ok(x);
        };
        let mut force = DVector::zeros(b.len());
        for t in a.iter().filter(|t| t.row == p && t.col < n) {
            force[t.col] -= t.val;
        }
        let sans_contraintes: Vec<_> = (0..b.len()).map(|i| i < n).collect();
        let (masse, rhs_libre) = selection(a, &force, &sans_contraintes);
        let libre = super::resout(&masse, n, &rhs_libre)?;
        let reference = energie(a, n, &libre);
        loop {
            pivots += 1;
            if pivots > budget {
                return Err("initialisation contact : budget de pivots dépassé".into());
            }
            let direction = super::resout(&trip, n, &force)?;
            let courbure = energie(a, n, &direction);
            let violation = -crate::raffinement::residu(a, b, &x)[p];
            let seuil = (f64::EPSILON * b.len() as f64).powi(2) * reference;
            let primal = if courbure > seuil {
                violation.max(0.) / courbure
            } else {
                f64::INFINITY
            };
            let borne = contacts
                .iter()
                .filter(|&&(i, _)| garde[i] && direction[i] < 0.)
                .map(|&(i, _)| (i, -x[i] / direction[i]))
                .min_by(|a, b| a.1.total_cmp(&b.1));
            let dual = borne.map_or(f64::INFINITY, |(_, t)| t.max(0.));
            let pas = primal.min(dual);
            if !pas.is_finite() {
                return Err(
                    "initialisation contact : accélérations unilatérales incompatibles".into(),
                );
            }
            x += direction * pas;
            x[p] += pas;
            if primal <= dual {
                garde[p] = true;
                let (nouvelle, rhs) = selection(a, b, &garde);
                trip = nouvelle;
                x = super::resout(&trip, n, &rhs)?;
                // Une réaction négative après recomposition traduit une perte
                // de précision ; elle ne doit pas devenir une force attractive.
                for &(i, _) in contacts {
                    if garde[i] && x[i] < 0. {
                        if x[i] < -tol * x.rows(n, b.len() - n).amax() {
                            return Err(
                                "initialisation contact : réaction négative après projection"
                                    .into(),
                            );
                        }
                        x[i] = 0.;
                    }
                }
                break;
            }
            let (sortant, _) = borne.expect("pas dual fini");
            garde[sortant] = false;
            x[sortant] = 0.;
            trip = selection(a, b, &garde).0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DMatrix;

    fn probleme(
        g: &DMatrix<f64>,
        masse: &DMatrix<f64>,
        f: &[f64],
        c: &[f64],
    ) -> (Triplets, DVector<f64>) {
        let mut a = DMatrix::zeros(3 + g.nrows(), 3 + g.nrows());
        a.view_mut((0, 0), (3, 3)).copy_from(masse);
        a.view_mut((3, 0), g.shape()).copy_from(g);
        a.view_mut((0, 3), (3, g.nrows())).copy_from(&g.transpose());
        let mut trip = Triplets::new();
        for j in 0..a.ncols() {
            for i in 0..a.nrows() {
                if a[(i, j)] != 0. {
                    trip.push(faer::sparse::Triplet {
                        row: i,
                        col: j,
                        val: a[(i, j)],
                    });
                }
            }
        }
        (
            trip,
            DVector::from_iterator(a.nrows(), f.iter().chain(c).copied()),
        )
    }

    #[test]
    fn faces_paralleles_remplacement_et_incompatibilite() {
        let g = DMatrix::from_row_slice(3, 3, &[1., 0., 0., 1., 0., 0., -1., 0., 0.]);
        let contacts = [(3, true), (4, true), (5, true)];
        // x <= 1, x <= 0, x >= -1 : deux faces parallèles, une seule active.
        let (a, b) = probleme(
            &g,
            &DMatrix::identity(3, 3),
            &[2., 0.3, -0.2],
            &[1., 0., 1.],
        );
        let x = complementarite(&a, 3, &b, &contacts).unwrap();
        assert!((x - DVector::from_vec(vec![0., 0.3, -0.2, 0., 2., 0.])).amax() < 1e-13);
        let mut faux = b.clone();
        faux[5] = -1.; // x >= 1 et x <= 0.
        assert!(complementarite(&a, 3, &faux, &contacts)
            .unwrap_err()
            .contains("incompatibles"));
        // Contact en décollement : aucune force, même sous une charge entrante.
        let x = complementarite(&a, 3, &b, &[(3, false), (4, false), (5, false)]).unwrap();
        assert_eq!(x, DVector::from_vec(vec![2., 0.3, -0.2, 0., 0., 0.]));
    }

    #[test]
    fn contacts_couples_reference_par_enumeration() {
        // Contre-calcul indépendant : toutes les 2^5 faces, LU nalgebra,
        // puis filtrage des conditions KKT. Pas de reprise de l'algorithme dual.
        for seed in 0..32 {
            let g = DMatrix::from_fn(5, 3, |i, j| {
                let argument = if seed % 2 == 0 {
                    i * 3 + j + 1
                } else {
                    (i + 1) * (j + 1)
                };
                ((argument + seed) as f64 * 0.73).sin()
            });
            let masse =
                DMatrix::from_row_slice(3, 3, &[2., 0.3, -0.1, 0.3, 1., 0.2, -0.1, 0.2, 3.]);
            let force = [(seed as f64 * 0.3).cos(), -1., 0.4];
            let c: Vec<_> = (0..5)
                .map(|i| 0.05 + 0.2 * ((i + seed) as f64).sin().abs())
                .collect();
            let (a, b) = probleme(&g, &masse, &force, &c);
            let mut reference = None;
            for mask in 0_u32..32 {
                let actifs: Vec<_> = (0..5).filter(|i| mask & (1 << i) != 0).collect();
                let mut kkt = DMatrix::zeros(3 + actifs.len(), 3 + actifs.len());
                kkt.view_mut((0, 0), (3, 3)).copy_from(&masse);
                let mut rhs = DVector::zeros(kkt.nrows());
                rhs.rows_mut(0, 3)
                    .copy_from(&DVector::from_column_slice(&force));
                for (j, &i) in actifs.iter().enumerate() {
                    for k in 0..3 {
                        kkt[(3 + j, k)] = g[(i, k)];
                        kkt[(k, 3 + j)] = g[(i, k)];
                    }
                    rhs[3 + j] = c[i];
                }
                let Some(y) = kkt.clone().lu().solve(&rhs) else {
                    continue;
                };
                if (&kkt * &y - &rhs).amax() > 1e-12 * (1. + rhs.amax()) {
                    continue;
                }
                if y.rows(3, actifs.len()).iter().any(|&v| v < -1e-10) {
                    continue;
                }
                let acc = y.rows(0, 3).into_owned();
                if (&g * &acc - DVector::from_vec(c.clone()))
                    .iter()
                    .any(|&v| v > 1e-10)
                {
                    continue;
                }
                reference = Some(acc);
                break;
            }
            let x = complementarite(
                &a,
                3,
                &b,
                &[(3, true), (4, true), (5, true), (6, true), (7, true)],
            )
            .unwrap();
            let attendu = reference.unwrap();
            assert!(
                (x.rows(0, 3) - &attendu).amax() < 2e-11,
                "seed {seed}: {x:?} vs {attendu:?}"
            );
            assert!(x.rows(3, 5).iter().all(|&v| v >= 0.));
        }
    }

    #[test]
    fn contacts_presque_paralleles_direction_faible() {
        for delta in [1e-4, 1e-8, 1e-12] {
            let g = DMatrix::from_row_slice(2, 3, &[1., delta, 0., -1., 0., 0.]);
            let (a, b) = probleme(&g, &DMatrix::identity(3, 3), &[0., 1., 0.3], &[0., 0.]);
            let x = complementarite(&a, 3, &b, &[(3, true), (4, true)]).unwrap();
            assert!(x[0].abs().max(x[1].abs()) < 1e-13, "{delta}: {x:?}");
            assert!((x[2] - 0.3).abs() < 1e-14);
            assert!((x[3] * delta - 1.).abs() < 2e-13);
            assert!((x[4] * delta - 1.).abs() < 2e-13);
        }
    }
}
