//! Équilibrage par congruence, avec une échelle commune par vecteur spatial.
//! Les blocs de trois coordonnées physiques tournent ensemble. Leurs normes
//! de Frobenius ne dépendent pas du choix des axes du repère.
use crate::{DVector, Triplets};

// Exposant d'une puissance de deux normale. Construire directement le
// flottant évite une exponentiation pour chaque coefficient creux.
fn puissance_normale(exposant: i32) -> f64 {
    debug_assert!((-1022..=1023).contains(&exposant));
    f64::from_bits(((exposant + 1023) as u64) << 52)
}

fn change_exposant(mut x: f64, mut exposant: i32) -> f64 {
    // L'échelle totale peut dépasser celle d'un f64 alors que x mis à
    // l'échelle reste représentable. Les étapes vont toutes dans le même
    // sens ; aucun produit intermédiaire des deux échelles n'est formé.
    while exposant > 1023 {
        x *= puissance_normale(1023);
        exposant -= 1023;
    }
    while exposant < -1022 {
        x *= puissance_normale(-1022);
        exposant += 1022;
    }
    x * puissance_normale(exposant)
}

#[derive(Clone, Copy, Default)]
struct Norme {
    maximum: f64,
    inverse: f64,
    somme: f64,
}

impl Norme {
    fn ajoute(&mut self, v: f64) {
        let a = v.abs();
        if a > self.maximum {
            self.somme = 1. + self.somme * (self.maximum / a).powi(2);
            self.maximum = a;
            self.inverse = a.recip();
        } else if a != 0. {
            let ratio = if self.inverse.is_finite() {
                a * self.inverse
            } else {
                a / self.maximum
            };
            self.somme += ratio * ratio;
        }
    }

    // Racine de la norme de Frobenius moyenne par ligne du bloc.
    // Ne former ni les carrés non normalisés, ni une norme qui déborderait
    // alors que sa racine est représentable.
    fn facteur(&self, taille: usize) -> f64 {
        self.maximum.sqrt() * (self.somme / taille as f64).sqrt().sqrt()
    }
}

#[derive(Default)]
pub(crate) struct Cache {
    physiques: usize,
    exposants: Vec<i32>,
}

fn coefficient(v: f64, exposant: i32) -> Option<f64> {
    let x = change_exposant(v, -exposant);
    (x.is_finite() && change_exposant(x, exposant) == v).then_some(x)
}

pub(crate) fn vecteur(mut v: DVector<f64>, d: &DVector<f64>) -> Result<DVector<f64>, String> {
    for i in 0..v.len() {
        let x = v[i] / d[i];
        if !x.is_finite() || x * d[i] != v[i] {
            return Err("équilibrage statique : perte d'une composante du vecteur".into());
        }
        v[i] = x;
    }
    Ok(v)
}

pub(crate) fn statique(
    trip: &mut Triplets,
    physiques: usize,
    dimension: usize,
    cache: &mut Cache,
) -> Result<DVector<f64>, String> {
    debug_assert!(physiques.is_multiple_of(3));
    let np = physiques / 3;
    let nb = np + dimension - physiques;
    let bloc = |i: usize| {
        if i < physiques {
            i / 3
        } else {
            np + i - physiques
        }
    };
    let mut d = DVector::from_element(dimension, 1.);
    if cache.physiques != physiques || cache.exposants.len() != dimension {
        cache.physiques = physiques;
        cache.exposants = vec![0; dimension];
    }
    // Newton change peu les échelles entre deux configurations voisines.
    // Repartir de celles du pas précédent économise les premières passes.
    // Si elles perdraient un coefficient du nouveau système, repartir de
    // l'identité, avant toute modification de la matrice.
    if trip
        .iter()
        .any(|t| coefficient(t.val, cache.exposants[t.row] + cache.exposants[t.col]).is_none())
    {
        cache.exposants.fill(0);
    }
    for t in trip.iter_mut() {
        t.val = change_exposant(t.val, -cache.exposants[t.row] - cache.exposants[t.col]);
    }
    for i in 0..dimension {
        d[i] = change_exposant(1., cache.exposants[i]);
    }
    let mut lignes = vec![Norme::default(); nb];
    let mut colonnes = lignes.clone();
    let mut exposants = vec![0; nb];
    for _ in 0..8 {
        lignes.fill(Norme::default());
        colonnes.fill(Norme::default());
        for t in trip.iter() {
            if !t.val.is_finite() {
                return Err("équilibrage statique : coefficient non fini".into());
            }
            lignes[bloc(t.row)].ajoute(t.val);
            colonnes[bloc(t.col)].ajoute(t.val);
        }
        let mut change = false;
        for i in 0..nb {
            let taille = if i < np { 3 } else { 1 };
            let facteur = lignes[i].facteur(taille).max(colonnes[i].facteur(taille));
            // Puissances de deux : pas d'arrondi de coefficient tant que
            // l'exposant reste dans le domaine représentable. Une grille
            // grossière suffit pour équilibrer, sans rechercher l'égalité.
            exposants[i] = if facteur == 0. {
                0
            } else {
                let bits = facteur.to_bits();
                let exposant = (bits >> 52) as i32 - 1023;
                let mantisse = f64::from_bits((bits & ((1_u64 << 52) - 1)) | (1023_u64 << 52));
                exposant + i32::from(mantisse >= std::f64::consts::SQRT_2)
            };
            change |= exposants[i] != 0;
        }
        if !change {
            break;
        }
        for i in 0..dimension {
            d[i] = change_exposant(d[i], exposants[bloc(i)]);
            if !d[i].is_finite() || d[i] == 0. {
                return Err("équilibrage statique : échelle non représentable".into());
            }
            cache.exposants[i] += exposants[bloc(i)];
        }
        for t in trip.iter_mut() {
            let e = exposants[bloc(t.row)] + exposants[bloc(t.col)];
            t.val = coefficient(t.val, e).ok_or("équilibrage statique : perte d'un coefficient")?;
        }
    }
    Ok(d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DMatrix, M3, V3};

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
    fn diagonale_presque_nulle_et_changements_de_repere() {
        let a = DMatrix::from_fn(9, 9, |i, j| {
            if i >= 6 && j >= 6 {
                0.
            } else if i == 0 && j == 0 {
                1e-24
            } else {
                ((i * j + i + j + 1) as f64).sin() * 3.
            }
        });
        let mut reference = triplets(&a);
        let d = statique(&mut reference, 6, 9, &mut Cache::default()).unwrap();
        assert!(d.min() >= 0.5 && d.max() <= 8.);
        let equilibree = crate::analyse::dense(&reference, 9);
        for j in 0..9 {
            for i in 0..9 {
                assert_eq!((equilibree[(i, j)] * d[i]) * d[j], a[(i, j)]);
            }
        }
        for axe in [V3::new(0.31, -0.47, 0.22), V3::new(-1.1, 0.6, 2.2)] {
            let q: M3 = nalgebra::Rotation3::from_scaled_axis(axe).into_inner();
            let mut s = DMatrix::identity(9, 9);
            s.view_mut((0, 0), (3, 3)).copy_from(&q);
            s.view_mut((3, 3), (3, 3)).copy_from(&q);
            let tournee = &s * &a * s.transpose();
            let mut t = triplets(&tournee);
            let dt = statique(&mut t, 6, 9, &mut Cache::default()).unwrap();
            assert_eq!(dt, d);
            assert!(
                (crate::analyse::dense(&t, 9) - &s * &equilibree * s.transpose()).amax() < 2e-15
            );
        }
    }

    #[test]
    fn zeros_et_etendues_extremes() {
        let mut t = Triplets::new();
        assert_eq!(
            statique(&mut t, 6, 9, &mut Cache::default()).unwrap(),
            DVector::from_element(9, 1.)
        );
        let a = DMatrix::from_diagonal(&DVector::from_vec(vec![f64::MAX, f64::from_bits(1)]));
        let mut t = triplets(&a);
        let mut cache = Cache::default();
        let d = statique(&mut t, 0, 2, &mut cache).unwrap();
        assert!(t
            .iter()
            .all(|v| v.val.is_finite() && v.val >= 0.5 && v.val <= 2.));
        for v in t {
            let e = 2 * d[v.row].log2() as i32;
            assert_eq!(change_exposant(v.val, e), a[(v.row, v.col)]);
        }
        // Le changement extrême rend les anciennes échelles inutilisables.
        // La tentative doit repartir de l'identité sans perdre le petit terme.
        let inverse = DMatrix::from_diagonal(&DVector::from_vec(vec![f64::from_bits(1), f64::MAX]));
        let mut t = triplets(&inverse);
        let d = statique(&mut t, 0, 2, &mut cache).unwrap();
        for v in t {
            let e = 2 * d[v.row].log2() as i32;
            assert_eq!(change_exposant(v.val, e), inverse[(v.row, v.col)]);
        }
        assert_eq!(
            change_exposant(2.0_f64.powi(-500), 1500),
            2.0_f64.powi(1000)
        );
        assert_eq!(
            change_exposant(2.0_f64.powi(1000), -1500),
            2.0_f64.powi(-500)
        );
        let petit = DVector::from_vec(vec![f64::from_bits(3)]);
        assert!(vecteur(petit.clone(), &DVector::from_element(1, 2.)).is_err());
        assert_eq!(
            vecteur(petit.clone(), &DVector::from_element(1, 1.)).unwrap(),
            petit
        );
        assert!(vecteur(
            DVector::from_element(1, f64::MAX),
            &DVector::from_element(1, 0.5)
        )
        .is_err());
    }
}
