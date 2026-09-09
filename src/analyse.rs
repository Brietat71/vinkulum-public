//! Tangentes d'analyse assemblées à partir des dépendances des éléments.
//! K = -∂f/∂q + ∂(Gᵀλ)/∂q, C = -∂f/∂u, à états internes figés.
//! Les rotations sont perturbées à gauche, en coordonnées monde.
use crate::{DMatrix, DVector, Elem, Modele, Ref, Triplets};

/// Dimensions, pointeurs de colonnes, indices de lignes, coefficients.
/// Les tableaux possèdent leurs données ; aucune vue sur un état mutable.
pub(crate) type CscPython = (usize, usize, Vec<usize>, Vec<usize>, Vec<f64>);

/// Extraction partielle : D = dz/dq des seules poutres, sans tangente globale.
/// Les six lignes de chaque poutre sont contiguës ; les colonnes gardent les
/// six coordonnées physiques de tous les corps, même ceux encastrés.
pub(crate) struct FacteursMaterielsPoutres {
    pub d: CscPython,
    pub deformations: Vec<f64>,
    pub masse: CscPython,
    pub contraintes: CscPython,
    pub phi: Vec<f64>,
    pub poutres: Vec<(String, usize, usize)>,
    pub infos_domaine: Vec<String>,
}

pub(crate) fn direction_poutre(quoi: &str) -> Result<([f64; 3], [f64; 3]), String> {
    match quoi {
        "ea" => Ok(([1.0, 0.0, 0.0], [0.0; 3])),
        "ga" => Ok(([0.0, 1.0, 1.0], [0.0; 3])),
        "gj" => Ok(([0.0; 3], [1.0, 0.0, 0.0])),
        // Directions historiques : ga et ei varient leurs deux composantes.
        "ei" => Ok(([0.0; 3], [0.0, 1.0, 1.0])),
        "ei_e2" => Ok(([0.0; 3], [0.0, 1.0, 0.0])),
        "ei_e3" => Ok(([0.0; 3], [0.0, 0.0, 1.0])),
        _ => Err("quoi : \"ea\", \"ga\", \"gj\", \"ei\", \"ei_e2\" ou \"ei_e3\"".into()),
    }
}

fn csc_python(mut trip: Triplets, rows: usize, cols: usize) -> Result<CscPython, String> {
    regroupe(&mut trip);
    if trip.iter().any(|v| !v.val.is_finite()) {
        return Err("linéarisation creuse : coefficient non fini après assemblage".into());
    }
    let mut pointers = vec![0; cols + 1];
    let mut indices = Vec::with_capacity(trip.len());
    let mut values = Vec::with_capacity(trip.len());
    for v in trip {
        pointers[v.col + 1] += 1;
        indices.push(v.row);
        values.push(v.val);
    }
    for i in 0..cols {
        pointers[i + 1] += pointers[i];
    }
    Ok((rows, cols, pointers, indices, values))
}

fn ajoute(t: &mut Triplets, row: usize, col: usize, val: f64) {
    if val != 0.0 {
        t.push(faer::sparse::Triplet { row, col, val });
    }
}

fn disperse(t: &mut Triplets, cs: &[Ref], a: &DMatrix<f64>, signe: f64) {
    for (il, i) in cs.iter().enumerate() {
        let Some(i) = i else { continue };
        for (jl, j) in cs.iter().enumerate() {
            let Some(j) = j else { continue };
            for p in 0..6 {
                for q in 0..6 {
                    ajoute(t, 6 * i + p, 6 * j + q, signe * a[(6 * il + p, 6 * jl + q)]);
                }
            }
        }
    }
}

pub(crate) fn dense(trip: &Triplets, dim: usize) -> DMatrix<f64> {
    let mut a = DMatrix::zeros(dim, dim);
    for v in trip {
        a[(v.row, v.col)] += v.val;
    }
    a
}

/// Regrouper les contributions avant équilibrage, dans l'ordre de leur
/// dispersion. Les petits coefficients restent présents ; seul zéro est omis.
pub(crate) fn regroupe(trip: &mut Triplets) {
    trip.sort_by_key(|v| (v.col, v.row)); // tri stable : ordre des sommes
    let mut out = 0;
    for read in 0..trip.len() {
        let v = trip[read];
        if out > 0 && (trip[out - 1].row, trip[out - 1].col) == (v.row, v.col) {
            trip[out - 1].val += v.val;
        } else {
            trip[out] = v;
            out += 1;
        }
    }
    trip.truncate(out);
    trip.retain(|v| v.val != 0.0);
}

impl Modele {
    /// Masse spatiale seule, sans forces, réactions, tangentes ni projection.
    fn masse_spatiale_creuse(&self, symetrise: bool) -> Result<CscPython, String> {
        let mut m = Triplets::with_capacity(12 * self.corps.len());
        for (i, body) in self.corps.iter().enumerate() {
            if !body.m.is_finite() || body.m <= 0.0 {
                return Err(format!("{} : masse invalide", body.nom));
            }
            let mut js = body.rot * body.j * body.rot.transpose();
            if symetrise {
                // La nouvelle extraction définit une masse flottante exactement
                // symétrique. L'analyse historique garde ses propres arrondis.
                for a in 0..3 {
                    for b in a + 1..3 {
                        let x = 0.5 * js[(a, b)] + 0.5 * js[(b, a)];
                        js[(a, b)] = x;
                        js[(b, a)] = x;
                    }
                }
            }
            for a in 0..3 {
                ajoute(&mut m, 6 * i + a, 6 * i + a, body.m);
                for b in 0..3 {
                    ajoute(&mut m, 6 * i + 3 + a, 6 * i + 3 + b, js[(a, b)]);
                }
            }
        }
        csc_python(m, self.n(), self.n())
    }

    /// Facteurs matériels natifs au point courant. La précontrainte des
    /// poutres reste dans z ; aucun terme géométrique ni d'autre élément
    /// n'est ajouté à D. G suit phi_g_locale, y compris ses lignes K nulles.
    /// Ne pas appeler la linéarisation : elle initialise éventuellement λ
    /// et calcule un autre objet (la tangente du champ complet).
    pub(crate) fn facteurs_materiels_poutres(
        &self,
        t: f64,
    ) -> Result<FacteursMaterielsPoutres, String> {
        crate::fini(&[t], "temps des facteurs matériels")?;
        crate::fini(self.g.as_slice(), "gravité")?;
        crate::fini(self.lam.as_slice(), "multiplicateurs")?;
        for c in &self.corps {
            for x in [
                c.r.as_slice(),
                c.r_bas.as_slice(),
                c.rot.as_slice(),
                c.v.as_slice(),
                c.w.as_slice(),
                c.j.as_slice(),
            ] {
                crate::fini(x, &c.nom)?;
            }
        }
        let mut d = Triplets::with_capacity(72 * self.poutres.len());
        let mut deformations = Vec::with_capacity(6 * self.poutres.len());
        let mut poutres = Vec::with_capacity(self.poutres.len());
        for (b, poutre) in self.poutres.iter().enumerate() {
            let local = poutre.facteur_materiel(&self.corps)?;
            deformations.extend_from_slice(local.deformation_ponderee.as_slice());
            poutres.push((poutre.nom.clone(), local.noeuds[0], local.noeuds[1]));
            for r in 0..6 {
                for c in 0..12 {
                    ajoute(
                        &mut d,
                        6 * b + r,
                        6 * local.noeuds[c / 6] + c % 6,
                        local.d[(r, c)],
                    );
                }
            }
        }
        let (phi, g) = self.phi_g_locale(t)?;
        crate::fini(phi.as_slice(), "contraintes des facteurs matériels")?;
        Ok(FacteursMaterielsPoutres {
            d: csc_python(d, 6 * self.poutres.len(), self.n())?,
            deformations,
            masse: self.masse_spatiale_creuse(true)?,
            contraintes: csc_python(g.trip, self.m(), self.n())?,
            phi: phi.as_slice().to_vec(),
            poutres,
            infos_domaine: self.infos_domaine_facteurs(),
        })
    }

    /// Diagnostics de périmètre, pas une validation de linéarité ni de SPD.
    fn infos_domaine_facteurs(&self) -> Vec<String> {
        let mut out = Vec::new();
        if self
            .corps
            .iter()
            .any(|c| c.v.iter().chain(c.w.iter()).any(|x| *x != 0.0))
        {
            out.push("vitesses non nulles : gyroscopie et termes de mouvement absents de D".into());
        }
        if self.g.iter().any(|x| *x != 0.0) {
            out.push("pesanteur présente : charge hors des facteurs de poutres".into());
        }
        if !self.efforts.is_empty() {
            out.push(format!(
                "{} efforts constants au CdM hors des facteurs de poutres",
                self.efforts.len()
            ));
        }
        for (type_, noms) in [
            (
                "couples à loi",
                self.couples
                    .iter()
                    .map(|e| e.nom.as_str())
                    .collect::<Vec<_>>(),
            ),
            (
                "superéléments",
                self.supers.iter().map(|e| e.nom.as_str()).collect(),
            ),
            ("pales", self.pales.iter().map(|e| e.nom.as_str()).collect()),
            (
                "contacts",
                self.contacts.iter().map(|e| e.nom.as_str()).collect(),
            ),
        ] {
            if !noms.is_empty() {
                out.push(format!(
                    "{type_} hors des facteurs de poutres : {}",
                    noms.join(", ")
                ));
            }
        }
        if !self.inflows.is_empty() {
            out.push("états de fluide présents : absents de D et de la masse mécanique".into());
        }
        if !self.spheres.is_empty() || self.loi_contact.is_some() {
            out.push(
                "appariement de contacts configuré : changements de paires absents de D et G"
                    .into(),
            );
        }
        if self.lam.iter().any(|x| *x != 0.0) {
            out.push("multiplicateurs non nuls : tangente des réactions absente de D".into());
        }
        if self.lam.len() != self.m() {
            out.push(
                "multiplicateurs non initialisés : extraction sans calcul de réactions".into(),
            );
        }
        for e in &self.elems {
            match e {
                Elem::K(_) => out.push(format!("{} : ligne de contact complémentaire nulle dans G et phi", e.nom())),
                Elem::L(l) => {
                    if l.nh {
                        out.push(format!("{} : contrainte non holonome, G agit sur les vitesses", l.nom));
                    }
                    if l.cible_t.is_some() || l.cible_r.is_some() {
                        out.push(format!("{} : contrainte pilotée, extraction à temps figé", l.nom));
                    }
                }
                _ => out.push(format!("{} : contrainte géométrique non linéaire, seule sa différentielle G est fournie", e.nom())),
            }
        }
        if self.actif.iter().any(|a| !a) {
            out.push("masque d'activité dynamique non appliqué à G ni phi".into());
        }
        out
    }

    /// Même champ que K,C,M,G d'analyse, sans construire les matrices
    /// denses ni une base du noyau des contraintes. Une évaluation locale
    /// commune fournit les dérivées en pose et en vitesse.
    pub(crate) fn linearisation_creuse(
        &self,
        t: f64,
    ) -> Result<(CscPython, CscPython, CscPython, CscPython), String> {
        let lam = self.multiplicateurs_analyse(t)?;
        let (mut k, c) = self.tangentes_forces_locales(t, true, true)?;
        self.precontrainte_locale(t, &lam, &mut k)?;
        let (_, g) = self.phi_g_locale(t)?;
        Ok((
            csc_python(k, self.n(), self.n())?,
            csc_python(c, self.n(), self.n())?,
            self.masse_spatiale_creuse(false)?,
            csc_python(g.trip, self.m(), self.n())?,
        ))
    }

    /// Contributions du champ `forces_a`, sans les réactions. L'ordre de
    /// dispersion suit celui du résidu ; les indices partagés s'additionnent.
    pub(crate) fn tangentes_forces_locales(
        &self,
        t: f64,
        pose: bool,
        vitesse: bool,
    ) -> Result<(Triplets, Triplets), String> {
        use crate::ad::Dual;
        crate::fini(&[t], "temps d'analyse")?;
        self.verifie_etat_fini()?;
        let (mut k, mut c) = (Triplets::new(), Triplets::new());
        for (i, cb) in self.corps.iter().enumerate() {
            type D = Dual<f64, 6>;
            // -f_rot = ω × J_s ω ; aucune accélération M(q)·a dans
            // cette définition de K, et aucun second terme gyroscopique.
            let d = crate::tangent::residu_rot(
                cb,
                std::array::from_fn(|a| D::var(0.0, a)),
                std::array::from_fn(|a| D::var(cb.w[a], 3 + a)),
                [0.0; 3],
            );
            for p in 0..3 {
                for q in 0..3 {
                    if pose {
                        ajoute(&mut k, 6 * i + 3 + p, 6 * i + 3 + q, d[p].d[q]);
                    }
                    if vitesse {
                        ajoute(&mut c, 6 * i + 3 + p, 6 * i + 3 + q, d[p].d[3 + q]);
                    }
                }
            }
        }
        // Pesanteur et efforts constants au CdM : tangentes nulles.
        for cp in &self.couples {
            let d = cp.tangentes(&self.corps, t);
            for (target, sg) in [(cp.a, 1.0), (cp.b, -1.0)] {
                let Some(i) = target else { continue };
                for (jlocal, source) in [cp.a, cp.b].iter().enumerate() {
                    let Some(j) = source else { continue };
                    for p in 0..3 {
                        for q in 0..3 {
                            if pose {
                                ajoute(
                                    &mut k,
                                    6 * i + 3 + p,
                                    6 * j + 3 + q,
                                    sg * d[(p, 3 * jlocal + q)],
                                );
                            }
                            if vitesse {
                                ajoute(
                                    &mut c,
                                    6 * i + 3 + p,
                                    6 * j + 3 + q,
                                    sg * d[(p, 6 + 3 * jlocal + q)],
                                );
                            }
                        }
                    }
                }
            }
        }
        for pale in &self.pales {
            let d = pale.tangentes(&self.corps, &self.inflows, self.vent);
            let base = 6 * pale.corps;
            for p in 0..6 {
                for q in 0..6 {
                    if pose {
                        ajoute(&mut k, base + p, base + q, -d[(p, q)]);
                    }
                    if vitesse {
                        ajoute(&mut c, base + p, base + q, -d[(p, 6 + q)]);
                    }
                }
            }
        }
        if pose {
            for poutre in &self.poutres {
                disperse(
                    &mut k,
                    &[Some(poutre.a), Some(poutre.b)],
                    &poutre.tangente(&self.corps),
                    -1.0,
                );
            }
        }
        for se in &self.supers {
            let (ke, ce) = se.tangentes(&self.corps)?;
            let cs: Vec<_> = se.noeuds.iter().copied().map(Some).collect();
            if pose {
                disperse(&mut k, &cs, &ke, -1.0);
            }
            if vitesse {
                disperse(&mut c, &cs, &ce, -1.0);
            }
        }
        for ct in &self.contacts {
            // Même λ figé que forces_a : la normale non lisse est hors
            // champ, le frottement dépend du dernier multiplicateur.
            let mut ct = ct.clone();
            if ct.nonlisse {
                ct.fn_impose = Some(self.lam.get(ct.row).copied().unwrap_or(0.0).max(0.0));
            }
            let d = ct.tangentes(&self.corps);
            let cs = [Some(ct.corps), ct.b];
            for (il, i) in cs.iter().enumerate() {
                let Some(i) = i else { continue };
                for (jl, j) in cs.iter().enumerate() {
                    let Some(j) = j else { continue };
                    for p in 0..6 {
                        for q in 0..6 {
                            if pose {
                                ajoute(&mut k, 6 * i + p, 6 * j + q, -d[(6 * il + p, 12 * jl + q)]);
                            }
                            if vitesse {
                                ajoute(
                                    &mut c,
                                    6 * i + p,
                                    6 * j + q,
                                    -d[(6 * il + p, 12 * jl + 6 + q)],
                                );
                            }
                        }
                    }
                }
            }
        }
        if k.iter().chain(&c).any(|v| !v.val.is_finite()) {
            return Err("tangentes locales : valeur non finie".into());
        }
        Ok((k, c))
    }

    fn multiplicateurs_analyse(&self, t: f64) -> Result<DVector<f64>, String> {
        crate::fini(&[t], "temps d'analyse")?;
        self.verifie_etat_fini()?;
        // Les réactions initialisées voient le modèle entier. Une copie
        // n'est nécessaire que lorsque λ n'a pas encore été calculé.
        if self.lam.len() == self.m() {
            Ok(self.lam.clone())
        } else {
            Ok(self.clone().acc_init(t)?.1)
        }
    }

    pub(crate) fn raideur_locale(&self, t: f64) -> Result<Triplets, String> {
        let lam = self.multiplicateurs_analyse(t)?;
        let (mut k, _) = self.tangentes_forces_locales(t, true, false)?;
        self.precontrainte_locale(t, &lam, &mut k)?;
        Ok(k)
    }

    fn precontrainte_locale(
        &self,
        t: f64,
        lam: &DVector<f64>,
        k: &mut Triplets,
    ) -> Result<(), String> {
        let mut row = 0;
        for e in &self.elems {
            let ne = e.n();
            // phi_g omet les lignes de contact complémentaire dans
            // l'analyse statique ; conserver ce même champ de réactions.
            if !matches!(e, Elem::K(_)) {
                let (_, _, ke) = e.eval_ad(&self.corps, t, &lam.as_slice()[row..row + ne])?;
                disperse(k, &e.corps(), &ke, 1.0);
            }
            row += ne;
        }
        if k.iter().any(|v| !v.val.is_finite()) {
            return Err("précontrainte locale : valeur non finie".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        expm, Cardan, Corps, Couple, Distance, Engrenage, Liaison, Loi, LoiCouple, Modele, Vis, M3,
        V3,
    };

    fn modele_facteurs(integree: bool) -> Modele {
        let mut mo = Modele::new(V3::zeros());
        // Ordre des corps et des éléments volontairement non géométrique.
        for (i, x) in [1.4, 0.0, 0.7].into_iter().enumerate() {
            mo.corps.push(Corps {
                nom: format!("noeud{i}"),
                m: i as f64 + 1.0,
                j: M3::from_diagonal(&V3::new(0.2, 0.3, 0.4)),
                r: V3::new(x, 0.0, 0.0),
                r_bas: V3::zeros(),
                rot: M3::identity(),
                v: V3::zeros(),
                w: V3::zeros(),
            });
        }
        for (a, b) in [(2, 0), (1, 2)] {
            mo.poutres.push(crate::Poutre {
                nom: format!("poutre{a}{b}"),
                integree,
                a,
                b,
                l: 0.7,
                cn: [4e5, 8e4, 5e4],
                cm: [120., 350., 900.],
                ra: M3::identity(),
                rb: M3::identity(),
            });
        }
        mo
    }

    fn dense_csc(c: &CscPython) -> DMatrix<f64> {
        let mut a = DMatrix::zeros(c.0, c.1);
        assert_eq!(c.2.len(), c.1 + 1);
        assert_eq!(c.3.len(), c.4.len());
        for j in 0..c.1 {
            let rows = &c.3[c.2[j]..c.2[j + 1]];
            assert!(rows.windows(2).all(|w| w[0] < w[1]));
            for p in c.2[j]..c.2[j + 1] {
                a[(c.3[p], j)] = c.4[p];
            }
        }
        a
    }

    fn encastrement_facteurs() -> Elem {
        Elem::L(Liaison {
            nom: "racine".into(),
            a: None,
            b: Some(1),
            pa: V3::zeros(),
            pb: V3::zeros(),
            ra: M3::identity(),
            rb: M3::identity(),
            bt: vec![0, 1, 2],
            br: vec![0, 1, 2],
            nh: false,
            cible_t: None,
            cible_r: None,
        })
    }

    #[test]
    fn facteurs_poutres_dispersion_et_raideur_au_repos() {
        for integree in [false, true] {
            let mo = modele_facteurs(integree);
            let f = mo.facteurs_materiels_poutres(0.0).unwrap();
            assert_eq!((f.d.0, f.d.1), (12, 18));
            assert_eq!(
                f.poutres,
                vec![("poutre20".into(), 2, 0), ("poutre12".into(), 1, 2)]
            );
            let d = dense_csc(&f.d);
            for (b, p) in mo.poutres.iter().enumerate() {
                let local = p.facteur_materiel(&mo.corps).unwrap();
                for r in 0..6 {
                    for c in 0..12 {
                        assert_eq!(
                            d[(6 * b + r, 6 * local.noeuds[c / 6] + c % 6)],
                            local.d[(r, c)]
                        );
                    }
                }
            }
            assert!(f.deformations.iter().all(|z| *z == 0.0));
            let (k, _, _, _) = mo.linearisation_creuse(0.0).unwrap();
            let k = dense_csc(&k);
            assert!((&d.transpose() * &d - &k).norm() < 3e-15 * k.norm());
            assert_eq!(dense_csc(&f.masse), mo.masse_spatiale());
            assert_eq!((f.contraintes.0, f.contraintes.1), (0, 18));
        }
    }

    #[test]
    fn facteurs_poutres_masse_g_sans_initialiser_ni_factoriser() {
        let mut mo = modele_facteurs(true);
        mo.elems.push(encastrement_facteurs());
        if let Elem::L(l) = &mut mo.elems[0] {
            l.nh = true;
            l.cible_t = Some((V3::x(), Loi::Lineaire { a0: 0.0, taux: 0.2 }));
        }
        mo.corps[0].rot = expm(&V3::new(0.1, 0.2, -0.3));
        mo.actif = vec![false; mo.m()];
        mo.t = 0.2;
        let avant = format!("{:?}", mo.corps);
        crate::partition::tests::FACTOS.with(|v| v.borrow_mut().clear());
        let f = mo.facteurs_materiels_poutres(0.8).unwrap();
        crate::partition::tests::FACTOS.with(|v| assert!(v.borrow().is_empty()));
        assert!(mo.lam.is_empty());
        assert_eq!(mo.t, 0.2);
        assert_eq!(mo.actif, vec![false; mo.m()]);
        assert_eq!(format!("{:?}", mo.corps), avant);
        let (phi, g) = mo.phi_g(0.8).unwrap();
        assert_eq!(f.phi, phi.as_slice());
        assert_eq!(dense_csc(&f.contraintes), g);
        let m = dense_csc(&f.masse);
        assert_eq!(m, m.transpose());
        assert!((&m - mo.masse_spatiale()).norm() < 3e-16 * m.norm());
        for motif in ["non holonome", "pilotée", "non initialisés", "masque"] {
            assert!(f.infos_domaine.iter().any(|s| s.contains(motif)), "{motif}");
        }
    }

    #[test]
    fn facteurs_poutres_precontraints_restent_partiels_et_conjugues() {
        let mut mo = modele_facteurs(true);
        mo.corps[0].r += V3::new(-0.02, 0.03, 0.01);
        mo.corps[0].rot = expm(&V3::new(0.04, -0.02, 0.13));
        mo.corps[2].v = V3::y();
        mo.corps[2].w = V3::new(0.1, 0.2, -0.1);
        mo.g = V3::new(0.0, 0.0, -9.81);
        mo.efforts.push((0, V3::x(), V3::y()));
        let f = mo.facteurs_materiels_poutres(0.0).unwrap();
        let d = dense_csc(&f.d);
        let z = DVector::from_vec(f.deformations);
        assert!(z.norm() > 1.0);
        let mut forces = DVector::zeros(mo.n());
        let mut k = Triplets::new();
        for p in &mo.poutres {
            let local = p.forces(&mo.corps);
            for (c, i) in [p.a, p.b].into_iter().enumerate() {
                for r in 0..6 {
                    forces[6 * i + r] += local[6 * c + r];
                }
            }
            disperse(
                &mut k,
                &[Some(p.a), Some(p.b)],
                &p.tangente(&mo.corps),
                -1.0,
            );
        }
        assert!((&d.transpose() * z + &forces).norm() < 3e-13 * forces.norm());
        let k = dense(&k, mo.n());
        assert!((d.transpose() * d - &k).norm() > 1e-3 * k.norm());
        for motif in ["vitesses", "pesanteur", "efforts constants"] {
            assert!(f.infos_domaine.iter().any(|s| s.contains(motif)), "{motif}");
        }
    }

    #[test]
    fn facteurs_poutres_contact_signale_ligne_nulle_et_conserve_etat() {
        let mut mo = modele_facteurs(true);
        mo.contacts.push(crate::Contact {
            nom: "sol".into(),
            corps: 0,
            p0: V3::zeros(),
            p1: V3::zeros(),
            rayon: 0.1,
            b: None,
            pb: V3::zeros(),
            p1b: V3::zeros(),
            rayon_b: 0.0,
            normale: V3::z(),
            origine: V3::zeros(),
            k: 1e5,
            expo: 1.5,
            d_hat: 1e-3,
            c: 0.0,
            mu: 0.2,
            v_eps: 1e-3,
            demi: None,
            cylindre: None,
            maille: None,
            q_maille: V3::new(0.1, 0.2, 0.3),
            s_axe: 0.3,
            cable: false,
            sortie: (0.1, 2.0),
            nonlisse: true,
            row: 0,
            fn_impose: Some(2.0),
            restitution: 0.0,
            auto: false,
        });
        mo.elems.push(Elem::K(crate::Unilateral {
            nom: "sol".into(),
            ct: 0,
            a: 0,
            b: None,
        }));
        mo.lam = DVector::from_element(1, 2.0);
        mo.actif = vec![true];
        mo.spheres.push((0, V3::zeros(), 0.1));
        let avant = format!("{:?}", mo.contacts);
        let f = mo.facteurs_materiels_poutres(0.0).unwrap();
        assert_eq!(f.phi, vec![0.0]);
        assert_eq!((f.contraintes.0, f.contraintes.1), (1, 18));
        assert!(f.contraintes.4.is_empty());
        assert_eq!(mo.lam[0], 2.0);
        assert_eq!(format!("{:?}", mo.contacts), avant);
        for motif in [
            "contacts",
            "complémentaire nulle",
            "appariement",
            "multiplicateurs non nuls",
        ] {
            assert!(f.infos_domaine.iter().any(|s| s.contains(motif)), "{motif}");
        }
    }

    #[test]
    fn facteurs_poutres_vides_et_entrees_invalides() {
        let vide = Modele::new(V3::zeros())
            .facteurs_materiels_poutres(0.0)
            .unwrap();
        assert_eq!(vide.d, (0, 0, vec![0], vec![], vec![]));
        assert!(vide.poutres.is_empty() && vide.deformations.is_empty());
        let mut mo = modele_facteurs(true);
        mo.poutres.clear();
        let f = mo.facteurs_materiels_poutres(0.0).unwrap();
        assert_eq!(f.d, (0, 18, vec![0; 19], vec![], vec![]));
        assert_eq!(dense_csc(&f.masse), mo.masse_spatiale());
        assert!(mo.facteurs_materiels_poutres(f64::NAN).is_err());
        let mut mo = modele_facteurs(true);
        mo.poutres[0].a = 99;
        assert!(mo.facteurs_materiels_poutres(0.0).is_err());
        let mut mo = modele_facteurs(true);
        mo.corps[0].v.x = f64::INFINITY;
        assert!(mo.facteurs_materiels_poutres(0.0).is_err());
        let mut mo = modele_facteurs(true);
        mo.corps[0].rot = expm(&V3::new(std::f64::consts::PI, 0.0, 0.0));
        assert!(mo.facteurs_materiels_poutres(0.0).is_err());
    }

    #[test]
    fn csc_indices_sommes_et_finitude() {
        let trip = vec![
            faer::sparse::Triplet {
                row: 1,
                col: 2,
                val: 1e16,
            },
            faer::sparse::Triplet {
                row: 0,
                col: 0,
                val: 2.0f64.powi(-70),
            },
            faer::sparse::Triplet {
                row: 1,
                col: 2,
                val: -1e16,
            },
            faer::sparse::Triplet {
                row: 1,
                col: 2,
                val: 1.0,
            },
        ];
        let c = csc_python(trip, 2, 4).unwrap();
        assert_eq!(c.2, vec![0, 1, 1, 2, 2]);
        assert_eq!(c.3, vec![0, 1]);
        assert_eq!(c.4, vec![2.0f64.powi(-70), 1.0]);
        assert_eq!(csc_python(vec![], 0, 4).unwrap().2, vec![0; 5]);
        let overflow = vec![
            faer::sparse::Triplet {
                row: 0,
                col: 0,
                val: f64::MAX
            };
            2
        ];
        assert!(csc_python(overflow, 1, 1).is_err());
    }

    #[test]
    fn gouverneur_sans_couple_ne_cree_pas_amortissement() {
        let mut mo = Modele::new(V3::zeros());
        mo.corps.push(Corps {
            nom: "arbre".into(),
            m: 1.0,
            j: M3::identity(),
            r: V3::zeros(),
            r_bas: V3::zeros(),
            rot: M3::identity(),
            v: V3::zeros(),
            w: V3::zeros(),
        });
        mo.couples.push(Couple {
            nom: "gouverneur nul".into(),
            a: None,
            b: Some(0),
            axe: V3::z(),
            ref_rel: M3::identity(),
            prev: 0.0,
            loi: LoiCouple::Gouverneur {
                kg: 2.0,
                q_max: 0.0,
                cible: Loi::Constante(0.0),
            },
        });
        for w in [-1.0, 0.0, 1.0] {
            mo.corps[0].w = V3::new(0.0, 0.0, w);
            assert_eq!(mo.forces_a(0.0).amax(), 0.0);
            assert_eq!(mo.raideur(0.0).unwrap().amax(), 0.0);
            assert_eq!(mo.amortissement(0.0).unwrap().amax(), 0.0);
        }
    }

    #[test]
    fn regroupement_preserve_petits_termes_et_ordre() {
        let mut t = Triplets::new();
        for (i, j, x) in [
            (0, 1, 1e16),
            (1, 0, 2.0f64.powi(-70)),
            (0, 1, -1e16),
            (1, 1, 3.0),
            (0, 1, 1.0),
            (1, 1, -3.0),
        ] {
            ajoute(&mut t, i, j, x);
        }
        let avant = dense(&t, 2);
        regroupe(&mut t);
        assert_eq!(dense(&t, 2), avant);
        assert_eq!(t.len(), 2);
        assert_eq!(t[0].val, 2.0f64.powi(-70));
        assert_eq!(t[1].val, 1.0);
    }

    #[test]
    fn assemblage_precontraint_contre_champ_global() {
        let mut mo = Modele::new(V3::new(0.0, 0.0, -9.81));
        for i in 0..4 {
            let a = i as f64 + 1.0;
            mo.corps.push(Corps {
                nom: i.to_string(),
                m: a,
                j: M3::from_diagonal(&V3::new(0.2, 0.3, 0.4)),
                r: V3::new(0.3 * a, -0.2 * a, 0.1 * a),
                r_bas: V3::new(1e-17, -2e-17, 3e-17),
                rot: expm(&V3::new(0.03 * a, -0.02 * a, 0.04 * a)),
                v: V3::new(0.1, 0.2, -0.1),
                w: V3::new(0.2 * a, -0.3, 0.4),
            });
        }
        // Indices dispersés, porteur mobile, bâti, corps commun aux deux
        // extrémités et liaison non holonome : aucune symétrisation de K.
        for (a, b, nh) in [
            (Some(3), Some(1), false),
            (None, Some(2), false),
            (Some(0), Some(2), true),
            (Some(2), Some(2), false),
        ] {
            mo.elems.push(Elem::L(Liaison {
                nom: "liaison".into(),
                a,
                b,
                pa: V3::new(0.1, -0.02, 0.03),
                pb: V3::new(0.02, 0.04, -0.01),
                ra: expm(&V3::new(0.1, -0.2, 0.3)),
                rb: M3::identity(),
                bt: vec![0, 2],
                br: vec![0, 1],
                nh,
                cible_t: Some((V3::x(), Loi::Constante(0.03))),
                cible_r: Some((V3::z(), Loi::Lineaire { a0: 0.1, taux: 0.2 })),
            }));
        }
        mo.elems.push(Elem::D(Distance {
            nom: "distance".into(),
            a: Some(2),
            b: Some(0),
            pa: V3::new(0.02, 0.03, 0.01),
            pb: V3::zeros(),
            l: 0.7,
        }));
        mo.elems.push(Elem::E(Engrenage {
            nom: "engrenage".into(),
            a: Some(3),
            b: Some(0),
            c: Some(1),
            na: V3::z(),
            nb: V3::z(),
            rapport: -2.0,
            ref_a: M3::identity(),
            ref_b: M3::identity(),
            prev: (0.0, 0.0),
        }));
        mo.elems.push(Elem::V(Vis {
            nom: "vis".into(),
            a: Some(1),
            b: Some(3),
            n: V3::x(),
            pas: 0.03,
            ref_rel: M3::identity(),
            c0: 0.0,
            prev: 0.0,
        }));
        mo.elems.push(Elem::C(Cardan {
            nom: "cardan".into(),
            a: Some(0),
            b: Some(3),
            na: V3::x(),
            nb: V3::y(),
        }));
        mo.lam = DVector::from_iterator(mo.m(), (0..mo.m()).map(|i| (i as f64 + 0.3).sin()));
        let t = 0.7;
        let k = mo.raideur(t).unwrap();
        let champ = |m: &Modele| -m.forces_a(t) + m.phi_g(t).unwrap().1.transpose() * &mo.lam;
        let mut fd = DMatrix::zeros(mo.n(), mo.n());
        let h = 1e-6;
        for j in 0..mo.n() {
            let mut values = Vec::new();
            for s in [1.0, -1.0] {
                let mut m = mo.clone();
                let c = &mut m.corps[j / 6];
                let mut d = V3::zeros();
                d[j % 3] = s * h;
                if j % 6 < 3 {
                    c.deplace(d);
                } else {
                    c.rot = expm(&d) * c.rot;
                }
                values.push(champ(&m));
            }
            fd.set_column(j, &((&values[0] - &values[1]) / (2.0 * h)));
        }
        assert!((&k - &fd).amax() < 2e-8, "écart {}", (&k - &fd).amax());
        assert!(
            (&k - k.transpose()).norm() > 0.01,
            "hors équilibre, K n'est pas symétrique"
        );
    }
}
