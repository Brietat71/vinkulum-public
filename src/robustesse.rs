//! Validation commune et sauvegarde des états évolutifs, sans copier les
//! maillages ni les matrices de raideur à chaque pas.
use crate::*;

pub(crate) fn fini(v: &[f64], nom: &str) -> Result<(), String> {
    if v.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(format!("{nom} : valeur non finie"))
    }
}
pub(crate) fn reglages(t0: f64, t1: f64, h: f64, tol: f64, iters: usize) -> Result<(), String> {
    fini(&[t0, t1, h, tol], "temps, pas ou tolérance")?;
    if h <= 0.0 || tol <= 0.0 || iters == 0 || t1 < t0 {
        return Err(
            "temps final ≥ temps courant, pas et tolérance > 0, itérations ≥ 1 requis".into(),
        );
    }
    Ok(())
}

#[derive(Clone)]
pub(crate) struct Sauvegarde {
    corps: Vec<(V3, M3, V3, V3, V3)>,
    elems: Vec<(usize, (f64, f64))>,
    couples: Vec<f64>,
    pales: Vec<Pale>,
    inflows: Vec<Inflow>,
    contacts: Vec<Contact>,
    t: f64,
    h: f64,
    lam: DVector<f64>,
    actif: Vec<bool>,
    gel: Vec<bool>,
    gel_elem: Vec<bool>,
    k_pas: Vec<(bool, f64)>,
    moy: Option<(f64, V3, V3, f64)>,
    hist_len: Option<usize>,
    schema_init: Option<(DVector<f64>, DVector<f64>, DVector<f64>)>,
    schema_fin: Option<(DVector<f64>, DVector<f64>, DVector<f64>)>,
    n_paires: usize,
    s_fb: f64,
    l_ref: Option<V3>,
    e_ref: Option<f64>,
}
impl Sauvegarde {
    pub(crate) fn new(m: &Modele) -> Self {
        Self {
            corps: m
                .corps
                .iter()
                .map(|c| (c.r, c.rot, c.v, c.w, c.r_bas))
                .collect(),
            elems: m
                .elems
                .iter()
                .enumerate()
                .filter_map(|(i, e)| match e {
                    Elem::E(e) => Some((i, e.prev)),
                    Elem::V(v) => Some((i, (v.prev, 0.0))),
                    _ => None,
                })
                .collect(),
            couples: m.couples.iter().map(|c| c.prev).collect(),
            pales: m.pales.clone(),
            inflows: m.inflows.clone(),
            contacts: m.contacts.clone(),
            t: m.t,
            h: m.h_dernier,
            lam: m.lam.clone(),
            actif: m.actif.clone(),
            gel: m.gel.clone(),
            gel_elem: m.gel_elem.clone(),
            k_pas: m.k_pas.clone(),
            moy: m.moy_aero,
            hist_len: m.hist_schema.as_ref().map(Vec::len),
            schema_init: m.schema_init.clone(),
            schema_fin: m.schema_fin.clone(),
            n_paires: m.n_paires,
            s_fb: m.s_fb,
            l_ref: m.l_ref,
            e_ref: m.e_ref,
        }
    }
    pub(crate) fn actualise(&mut self, m: &Modele) {
        for (q, c) in self.corps.iter_mut().zip(&m.corps) {
            *q = (c.r, c.rot, c.v, c.w, c.r_bas);
        }
        for (i, p) in &mut self.elems {
            match &m.elems[*i] {
                Elem::E(e) => *p = e.prev,
                Elem::V(v) => p.0 = v.prev,
                _ => {}
            }
        }
        for (p, c) in self.couples.iter_mut().zip(&m.couples) {
            *p = c.prev;
        }
        self.pales.clone_from(&m.pales);
        self.inflows.clone_from(&m.inflows);
        self.contacts.clone_from(&m.contacts);
        self.t = m.t;
        self.h = m.h_dernier;
        self.lam.clone_from(&m.lam);
        self.actif.clone_from(&m.actif);
        self.gel.clone_from(&m.gel);
        self.gel_elem.clone_from(&m.gel_elem);
        self.k_pas.clone_from(&m.k_pas);
        self.moy = m.moy_aero;
        self.hist_len = m.hist_schema.as_ref().map(Vec::len);
        self.schema_init.clone_from(&m.schema_init);
        self.schema_fin.clone_from(&m.schema_fin);
        self.n_paires = m.n_paires;
        self.s_fb = m.s_fb;
        self.l_ref = m.l_ref;
        self.e_ref = m.e_ref;
    }
    pub(crate) fn restaure(self, m: &mut Modele) {
        for (c, (r, rot, v, w, bas)) in m.corps.iter_mut().zip(self.corps) {
            c.r = r;
            c.r_bas = bas;
            c.rot = rot;
            c.v = v;
            c.w = w;
        }
        for (i, p) in self.elems {
            match &mut m.elems[i] {
                Elem::E(e) => e.prev = p,
                Elem::V(v) => v.prev = p.0,
                _ => {}
            }
        }
        for (c, p) in m.couples.iter_mut().zip(self.couples) {
            c.prev = p;
        }
        m.pales = self.pales;
        m.inflows = self.inflows;
        m.contacts = self.contacts;
        m.t = self.t;
        m.h_dernier = self.h;
        m.lam = self.lam;
        m.actif = self.actif;
        m.gel = self.gel;
        m.gel_elem = self.gel_elem;
        m.k_pas = self.k_pas;
        m.moy_aero = self.moy;
        if let (Some(h), Some(n)) = (&mut m.hist_schema, self.hist_len) {
            h.truncate(n);
        }
        m.schema_init = self.schema_init;
        m.schema_fin = self.schema_fin;
        m.n_paires = self.n_paires;
        m.s_fb = self.s_fb;
        m.l_ref = self.l_ref;
        m.e_ref = self.e_ref;
        m.symb = None;
        m.symb_motif = 0;
    }
}

impl Modele {
    pub(crate) fn verifie_etat_fini(&self) -> Result<(), String> {
        fini(self.g.as_slice(), "gravité")?;
        fini(self.lam.as_slice(), "multiplicateurs")?;
        for c in &self.corps {
            for x in [
                c.r.as_slice(),
                c.r_bas.as_slice(),
                c.rot.as_slice(),
                c.v.as_slice(),
                c.w.as_slice(),
                c.j.as_slice(),
            ] {
                fini(x, &c.nom)?;
            }
            if !c.m.is_finite() || c.m <= 0.0 {
                return Err(format!("{} : masse invalide", c.nom));
            }
        }
        for (_, f, m) in &self.efforts {
            fini(f.as_slice(), "effort")?;
            fini(m.as_slice(), "moment")?;
        }
        for charge in &self.efforts_temporels {
            charge.verifie(&self.corps)?;
        }
        for i in &self.inflows {
            fini(
                &[i.v_i, i.v1s, i.v1c, i.poussee, i.m_roul, i.m_tang],
                "état inflow",
            )?;
        }
        for p in &self.pales {
            for &(a, b) in &p.z {
                fini(&[a, b], "état Wagner")?;
            }
            for e in &p.etats_lb {
                fini(
                    &[e.dp, e.cn_prev, e.ff, e.tau_v, e.cv_prev, e.cnv, e.al_prev],
                    "état LB",
                )?;
            }
        }
        for se in &self.supers {
            fini(se.k.as_slice(), "raideur superélément")?;
            se.verifie_branche(&self.corps)?;
        }
        for p in &self.poutres {
            let rel = (self.corps[p.a].rot * p.ra).transpose() * self.corps[p.b].rot * p.rb;
            if (rel.trace() + 1.0).abs() < 1e-12 {
                return Err(format!(
                    "poutre « {} » : rotation relative à π, tangente non différentiable",
                    p.nom
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn verifie_options(&self) -> Result<(), String> {
        if !(0.0..=1.0).contains(&self.sigma_lie) {
            return Err("sigma_lie : valeur finie dans [0, 1] requise".into());
        }
        if self.sigma_lie != 0.0 && self.adapt.is_some() {
            return Err("sigma_lie : estimateur adaptatif non validé, utiliser un pas fixe".into());
        }
        if let Some(a) = self.adapt {
            fini(&[a, self.adapt_bornes.0, self.adapt_bornes.1], "adaptatif")?;
            if a <= 0.0
                || self.adapt_bornes.0 <= 0.0
                || self.adapt_bornes.0 > 1.0
                || self.adapt_bornes.1 < 1.0
            {
                return Err(
                    "adaptatif : tolérance > 0 et 0 < borne min ≤ 1 ≤ borne max requis".into(),
                );
            }
        }
        if let Some((l, c, m)) = self.pas_contact {
            fini(&[l, c, m], "pas de contact")?;
            if c <= 0.0 || l < c || m < 0.0 {
                return Err("contact : 0 < pas contact ≤ pas libre, marge ≥ 0 requis".into());
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn body(i: usize) -> Corps {
        Corps {
            r_bas: crate::V3::zeros(),
            nom: format!("c{i}"),
            m: 1.0,
            j: M3::identity(),
            r: V3::new(i as f64, 0.0, 0.0),
            rot: M3::identity(),
            v: V3::zeros(),
            w: V3::zeros(),
        }
    }

    #[test]
    fn log_pi_et_signe() {
        for axis in [V3::x(), V3::new(-0.8, 0.2, 0.1).normalize()] {
            for delta in [0.0, 1e-8, 1e-6] {
                let r = expm(&((std::f64::consts::PI - delta) * axis));
                assert!((expm(&log_so3(&r)) - r).norm() < 5e-8);
                let th = tangent::logm::<f64>(&ad::m_from(&r));
                assert!((expm(&V3::from(th)) - r).norm() < 1e-10);
            }
        }
        let r = M3::from_diagonal(&V3::new(1.0, -1.0, -1.0));
        assert!(
            (V3::from(tangent::logm::<f64>(&ad::m_from(&r))).norm() - std::f64::consts::PI).abs()
                < 1e-14
        );
    }

    #[test]
    fn svd_revient_au_repere_courant() {
        let mut cs: Vec<_> = (0..8).map(body).collect();
        let trip: Triplets = (0..48)
            .filter(|&i| i != 1)
            .map(|i| faer::sparse::Triplet {
                row: i,
                col: i,
                val: 1.0,
            })
            .collect();
        let mut garde = JacGarde::neuf(
            trip.clone(),
            Some(Facto::dense_seule(&trip, 48)),
            0.01,
            48,
            &cs,
            None,
        );
        let rot = expm(&(V3::z() * 0.7));
        for c in &mut cs {
            c.rot = rot;
        }
        let mut b = DVector::zeros(48);
        b.rows_mut(0, 3).copy_from(&(rot * V3::x()));
        let (x, svd) = garde.resout(&cs, &[], 48, &b);
        assert!(svd);
        assert!((x - b).norm() < 1e-12);
    }

    #[test]
    fn erreur_restaure_dernier_pas_et_historique() {
        let mut m = Modele::new(V3::new(0.0, 0.0, -9.81));
        m.corps.push(body(0));
        m.hist_schema = Some(vec![]);
        m.corps[0].r_bas.x = 2_f64.powi(-60);
        let before = m.corps.clone();
        // Overflow from finite inputs must not turn into a converged NaN state.
        m.efforts.push((0, V3::new(1e308, 1e308, 0.0), V3::zeros()));
        assert!(m
            .simule(0.1, 0.01, 0.9, 1e-12, 3, |_, _| panic!("pas accepté"))
            .is_err());
        assert_eq!(m.t, 0.0);
        assert_eq!(m.corps[0].r, before[0].r);
        assert_eq!(m.corps[0].r_bas, before[0].r_bas);
        assert_eq!(m.hist_schema.unwrap().len(), 0);
    }
}
