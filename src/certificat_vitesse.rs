//! Décision exacte sur |Σ g_j u_j + dt| ≤ tol + 2⁻⁴⁶ (|dt| + Σ |g_j u_j|).
//!
//! Les binary64 finis sont interprétés comme des rationnels dyadiques exacts.
//! Tous les produits sont des entiers en unités de 2⁻²¹⁴⁸. Ni opération
//! flottante ni arrondi dans la décision ; voir docs/CERTIFICATION_NOYAU.md.
use num_bigint::{BigInt, BigUint, Sign};
use pyo3::prelude::*;

struct Dyadique {
    negatif: bool,
    mantisse: u64,
    exposant: i32,
}

fn decode(x: f64) -> Result<Dyadique, String> {
    let bits = x.to_bits();
    let champ = ((bits >> 52) & 0x7ff) as i32;
    if champ == 0x7ff {
        return Err("certificat de vitesse : donnée non finie".into());
    }
    let fraction = bits & ((1_u64 << 52) - 1);
    Ok(Dyadique {
        negatif: bits >> 63 != 0,
        mantisse: fraction | if champ == 0 { 0 } else { 1_u64 << 52 },
        exposant: if champ == 0 { -1074 } else { champ - 1075 },
    })
}

fn entier(mantisse: u128, exposant: i32, negatif: bool) -> BigInt {
    // Pour un produit de deux binary64 finis : -2148 ≤ exposant ≤ 1942.
    // Les mantisses ont au plus 53 bits : leur produit tient dans u128.
    let magnitude = BigUint::from(mantisse) << (exposant + 2148) as usize;
    BigInt::from_biguint(if negatif { Sign::Minus } else { Sign::Plus }, magnitude)
}

pub(crate) struct CertificatLigne {
    residu: BigInt,
    echelle: BigUint,
    tolerance: BigUint,
}

impl CertificatLigne {
    pub(crate) fn accepte(&self) -> bool {
        // 64 * EPSILON = 2⁻⁴⁶. Comparaison entière, y compris à la frontière.
        (self.residu.magnitude() << 46_usize) <= (&self.tolerance << 46_usize) + &self.echelle
    }
}

pub(crate) fn ligne(
    termes: impl IntoIterator<Item = (f64, f64)>,
    dt: f64,
    tol: f64,
) -> Result<CertificatLigne, String> {
    let t = decode(tol)?;
    if t.negatif && t.mantisse != 0 {
        return Err("certificat de vitesse : tolérance négative".into());
    }
    let tolerance = entier(t.mantisse.into(), t.exposant, false)
        .magnitude()
        .clone();
    let d = decode(dt)?;
    let mut residu = entier(d.mantisse.into(), d.exposant, d.negatif);
    let mut echelle = residu.magnitude().clone();
    for (g, u) in termes {
        let (a, b) = (decode(g)?, decode(u)?);
        if a.mantisse == 0 || b.mantisse == 0 {
            continue;
        }
        let produit = entier(
            u128::from(a.mantisse) * u128::from(b.mantisse),
            a.exposant + b.exposant,
            a.negatif != b.negatif,
        );
        echelle += produit.magnitude();
        residu += produit;
    }
    Ok(CertificatLigne {
        residu,
        echelle,
        tolerance,
    })
}

/// Entrée privée de contre-épreuve : même calcul que le garde du noyau.
/// Entiers hexadécimaux signés, en unités de 2⁻²¹⁴⁸ ; pas de f64 de sortie.
#[pyfunction]
fn _certificat_ligne_vitesse(
    g: Vec<f64>,
    u: Vec<f64>,
    dt: f64,
    tol: f64,
) -> PyResult<(bool, String, String, String)> {
    if g.len() != u.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "certificat de vitesse : dimensions incompatibles",
        ));
    }
    let c =
        ligne(g.into_iter().zip(u), dt, tol).map_err(pyo3::exceptions::PyValueError::new_err)?;
    Ok((
        c.accepte(),
        c.residu.to_str_radix(16),
        c.echelle.to_str_radix(16),
        c.tolerance.to_str_radix(16),
    ))
}

pub(crate) fn enregistrer(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_certificat_ligne_vitesse, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sous_debordement_et_annulation_exacte() {
        let petit = f64::from_bits(1);
        let c = ligne([(petit, petit)], 0., 0.).unwrap();
        assert_eq!(c.residu, BigInt::from(1));
        assert!(!c.accepte());
        assert!(ligne([(f64::MAX, f64::MAX), (-f64::MAX, f64::MAX)], 0., 0.)
            .unwrap()
            .accepte());
    }

    #[test]
    fn frontiere_non_arrondie() {
        let u = f64::from_bits(1_f64.to_bits() - 1);
        let tol = u * (1. - 64. * f64::EPSILON);
        assert!(u <= tol + 64. * f64::EPSILON * u);
        assert!(!ligne([(1., u)], 0., tol).unwrap().accepte());
    }
}
