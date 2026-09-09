//! Assemblage des positions holonomes et projection des vitesses.
//! Les corrections portent sur G directement ; former GGᵀ perdait des
//! contraintes indépendantes avec les grands bras de levier.
use crate::{expm, fini, reglages, rotation, DMatrix, DVector, Elem, Modele, Sauvegarde, V3};

fn correction(g: &DMatrix<f64>, r: &DVector<f64>) -> Result<DVector<f64>, String> {
    fini(g.as_slice(), "jacobien d'assemblage")?;
    fini(r.as_slice(), "résidu d'assemblage")?;
    let (mut a, mut b) = (g.clone(), r.clone());
    // Normaliser chaque équation préserve son ensemble de solutions et la
    // correction de norme minimale quand le système est compatible.
    // Le maximum évite le débordement du carré de la norme d'une ligne.
    for i in 0..a.nrows() {
        let maximum = a.row(i).amax();
        if maximum > 0. {
            a.row_mut(i).unscale_mut(maximum);
            b[i] /= maximum;
            let norme = a.row(i).norm();
            a.row_mut(i).unscale_mut(norme);
            b[i] /= norme;
        }
    }
    let x = crate::moindres_carres::resout(a, &b, 1e-12)?;
    fini(x.as_slice(), "correction d'assemblage")?;
    Ok(x)
}

/// Résidu sur les équations originales, avec accumulation compensée.
/// La seconde sortie mesure l'échelle des produits de chaque ligne, pour
/// distinguer la tolérance demandée des arrondis des données cinématiques.
fn residu(
    g: &DMatrix<f64>,
    u: &DVector<f64>,
    dt: &DVector<f64>,
) -> Result<(DVector<f64>, DVector<f64>), String> {
    let (mut r, mut echelle) = (dt.clone(), dt.map(f64::abs));
    for i in 0..g.nrows() {
        let mut bas = 0.;
        for j in 0..g.ncols() {
            let p = g[(i, j)] * u[j];
            let e = g[(i, j)].mul_add(u[j], -p);
            let (s, reste) = crate::numerique::deux_sommes(r[i], p);
            (r[i], bas) = crate::numerique::deux_sommes(s, bas + reste + e);
            echelle[i] += p.abs();
        }
        r[i] += bas;
    }
    fini(r.as_slice(), "vitesse des contraintes")?;
    fini(echelle.as_slice(), "échelle des vitesses des contraintes")?;
    Ok((r, echelle))
}

impl Modele {
    /// Masque les éléments gelés et, en position, les lignes non holonomes.
    /// Les redondances restent présentes : une ligne écartée de la
    /// factorisation dynamique doit néanmoins être physiquement satisfaite.
    fn masque_assemblage(&self, r: &mut DVector<f64>, g: &mut DMatrix<f64>, position: bool) {
        let mut row = 0;
        for (ie, e) in self.elems.iter().enumerate() {
            let ne = e.n();
            if self.gel_elem.get(ie).copied().unwrap_or(false)
                || (position && matches!(e, Elem::L(l) if l.nh))
            {
                r.rows_mut(row, ne).fill(0.);
                g.rows_mut(row, ne).fill(0.);
            }
            row += ne;
        }
    }

    fn gele_correction(&self, g: &mut DMatrix<f64>) {
        for i in 0..self.corps.len() {
            if self.gele(i) {
                g.columns_mut(6 * i, 6).fill(0.);
            }
        }
    }

    /// Résout les contraintes holonomes depuis une pose approchée.
    /// Newton amorti sur SO(3), correction de norme minimale par QR/SVD
    /// sur le jacobien équilibré. Rend (résidu maximal, corrections).
    /// En cas d'échec, restitue l'état antérieur à l'appel.
    pub fn assemble(&mut self, t: f64, tol: f64, iters: usize) -> Result<(f64, usize), String> {
        reglages(t, t, 1.0, tol, iters)?;
        let sauve = Sauvegarde::new(self);
        let r = self.assemble_interne(t, tol, iters);
        if r.is_err() {
            sauve.restaure(self);
        }
        r
    }

    pub(super) fn assemble_interne(
        &mut self,
        t: f64,
        tol: f64,
        iters: usize,
    ) -> Result<(f64, usize), String> {
        reglages(t, t, 1.0, tol, iters)?;
        self.verifie_etat_fini()?;
        // it compte les corrections déjà effectuées : la dernière mérite
        // le même contrôle de convergence que toutes les précédentes.
        for it in 0..=iters {
            let (mut phi, mut g) = self.phi_g(t)?;
            self.masque_assemblage(&mut phi, &mut g, true);
            fini(phi.as_slice(), "contraintes assemblage")?;
            fini(g.as_slice(), "jacobien assemblage")?;
            let norme = phi.amax();
            if norme <= tol {
                self.verifie_bassin(t)?;
                return Ok((norme, it));
            }
            if it == iters {
                return Err(format!(
                    "assemblage : pas convergé en {iters} itérations, ‖Φ‖ = {norme:.3e} — \
                     le mécanisme est-il assemblable depuis cette pose ?"
                ));
            }
            self.gele_correction(&mut g);
            let dq = -correction(&g, &phi)?;
            // Borne angulaire du Newton amorti, identique à l'assemblage
            // historique. La projection de R empêche la dérive hors SO(3).
            let ang = (0..self.corps.len())
                .map(|i| dq[6 * i + 3].hypot(dq[6 * i + 4]).hypot(dq[6 * i + 5]))
                .fold(0.0f64, f64::max);
            let sc = if ang > 0.5 { 0.5 / ang } else { 1.0 };
            for i in 0..self.corps.len() {
                if self.gele(i) {
                    continue;
                }
                let c = &mut self.corps[i];
                c.deplace(sc * V3::new(dq[6 * i], dq[6 * i + 1], dq[6 * i + 2]));
                c.rot = expm(&(sc * V3::new(dq[6 * i + 3], dq[6 * i + 4], dq[6 * i + 5]))) * c.rot;
                c.rot = rotation::projette(&c.rot)?;
            }
            self.verifie_etat_fini()?;
        }
        unreachable!("le dernier contrôle retourne un résultat")
    }

    fn phi_dt_assemblage(&self, t: f64) -> Result<DVector<f64>, String> {
        let mut dt = DVector::zeros(self.m());
        let mut row = 0;
        for e in &self.elems {
            if let Elem::L(l) = e {
                dt.rows_mut(row, e.n()).copy_from(&l.phi_dt(&self.corps, t));
            }
            row += e.n();
        }
        fini(dt.as_slice(), "dérivée temporelle des contraintes")?;
        Ok(dt)
    }

    /// Projette au plus près sur G u + Φ_t = 0. Aucune vitesse n'est écrite
    /// avant le contrôle des équations originales, y compris redondantes.
    pub(super) fn projette_vitesse(&mut self, t: f64, tol: f64) -> Result<(), String> {
        reglages(t, t, 1., tol, 1)?;
        self.verifie_etat_fini()?;
        let (_, mut g) = self.phi_g(t)?;
        let mut dt = self.phi_dt_assemblage(t)?;
        self.masque_assemblage(&mut dt, &mut g, false);
        fini(g.as_slice(), "jacobien de projection des vitesses")?;
        let mut gc = g.clone();
        self.gele_correction(&mut gc);
        let mut un = self.u();
        for it in 0..=3 {
            let (r, echelle) = residu(&g, &un, &dt)?;
            if r.iter()
                .zip(echelle.iter())
                .all(|(r, e)| r.abs() <= tol + 64. * f64::EPSILON * e)
            {
                // Le résidu compensé propose ; l'arithmétique entière
                // décide sur chaque équation originale avant toute écriture.
                let mut certifie = true;
                for i in 0..g.nrows() {
                    let c = crate::certificat_vitesse::ligne(
                        (0..g.ncols()).map(|j| (g[(i, j)], un[j])),
                        dt[i],
                        tol,
                    )?;
                    certifie &= c.accepte();
                }
                if certifie {
                    self.set_u(&un);
                    return Ok(());
                }
            }
            if it == 3 {
                return Err(format!(
                    "assemblage : vitesses incompatibles ou précision insuffisante, \
                     ‖G·u + ∂Φ/∂t‖ = {:.3e} (tolérance {tol:.3e})",
                    r.amax()
                ));
            }
            // Le contrôle suivant permet aussi de raffiner une première
            // correction affectée par l'annulation u - δu.
            un -= correction(&gc, &r)?;
            fini(un.as_slice(), "vitesses projetées")?;
        }
        unreachable!("le dernier contrôle retourne un résultat")
    }

    /// Diagnostic à état fixé, avec la dérivée de la loi utilisée par GGL
    /// et les liaisons non holonomes, sans différence temporelle auxiliaire.
    pub(super) fn phi_dot_a(&self, t: f64) -> Result<DVector<f64>, String> {
        fini(&[t], "temps du diagnostic")?;
        self.verifie_etat_fini()?;
        let (_, g) = self.phi_g(t)?;
        fini(g.as_slice(), "jacobien du diagnostic")?;
        Ok(residu(&g, &self.u(), &self.phi_dt_assemblage(t)?)?.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Corps, Liaison, M3};

    fn modele() -> Modele {
        let mut m = Modele::new(V3::zeros());
        for i in 0..2 {
            m.corps.push(Corps {
                nom: format!("c{i}"),
                m: 1.,
                j: M3::identity(),
                r: V3::new(1. - 0.75 * i as f64, 0., 0.),
                r_bas: V3::zeros(),
                rot: M3::identity(),
                v: V3::new(3. - i as f64, 0., 0.),
                w: V3::zeros(),
            });
        }
        m.gel = vec![true, false];
        m
    }

    fn liaison(a: Option<usize>, b: usize) -> Elem {
        Elem::L(Liaison {
            nom: format!("l{b}"),
            a,
            b: Some(b),
            pa: V3::zeros(),
            pb: V3::zeros(),
            ra: M3::identity(),
            rb: M3::identity(),
            bt: vec![0],
            br: vec![],
            cible_t: None,
            cible_r: None,
            nh: false,
        })
    }

    #[test]
    fn correction_ne_deplace_pas_un_corps_gele() {
        let mut m = modele();
        m.elems.push(liaison(Some(0), 1));
        m.actif = vec![false]; // masque de factorisation, pas de physique
        m.assemble(0., 1e-12, 1).unwrap();
        m.projette_vitesse(0., 1e-12).unwrap();
        assert_eq!(m.corps[0].r.x, 1.);
        assert_eq!(m.corps[0].v.x, 3.);
        assert!((m.corps[1].r.x - 1.).abs() < 1e-12);
        assert!((m.corps[1].v.x - 3.).abs() < 1e-12);
    }

    #[test]
    fn assemblage_ignore_les_elements_geles() {
        let mut m = modele();
        m.elems = vec![liaison(None, 0), liaison(None, 1)];
        m.gel_elem = vec![true, false];
        // Une autre partition peut être hors bassin à la date courante
        // d'une commande ; son contrôle appartient à son propre sous-pas.
        if let Elem::L(l) = &mut m.elems[0] {
            l.br = vec![0, 1, 2];
        }
        m.corps[0].rot = expm(&V3::new(0., 0., std::f64::consts::PI));
        let rotation_gelee = m.corps[0].rot;
        m.assemble(0., 1e-12, 1).unwrap();
        m.projette_vitesse(0., 1e-12).unwrap();
        assert_eq!(m.corps[0].r.x, 1.);
        assert_eq!(m.corps[0].v.x, 3.);
        assert_eq!(m.corps[0].rot, rotation_gelee);
        assert!(m.corps[1].r.x.abs() < 1e-12);
        assert!(m.corps[1].v.x.abs() < 1e-12);
    }
}
