//! tangent — le jacobien du résidu de Newton, par différentiation
//! automatique élément par élément (nombres duaux de `ad`), assemblé
//! analytiquement avec la tangente de l'exponentielle sur SO(3).
//!
//! Résidu (index 3, α-généralisé de Lie) en x = (u̇, λ) :
//!     r_dyn = M(q)·u̇ + G(q)ᵀλ − f(q, u, t)
//!     r_c   = Φ(q) / (β h²)
//! À sigma=0 : q = q_n ⊕ h·(u_n + (½−β)h a_n + βh a₁), u = u_n + (1−γ)h a_n + γh a₁,
//! a₁ = c·u̇ + cte, c = (1−α_f)/(1−α_m). Donc, par corps :
//!     ∂q/∂u̇ = βh²c · diag(I₃, J_l(θ))     (J_l : jacobien gauche de exp, R = exp(θ)R₀)
//!     ∂u/∂u̇ = γhc · I₆
//! et le jacobien :
//!     [ M + (K_M + K_c − K_f)·βh²c·J_l − C_f·γhc ,  Gᵀ ]
//!     [ G·J_l·c                                ,  0  ]
//! où K_M = ∂(M u̇)/∂q, K_c = ∂(Gᵀλ)/∂q (la Hessienne de λ·Φ — la « raideur
//! géométrique » des liaisons), K_f = ∂f/∂q, C_f = ∂f/∂u. K_c et G sortent des
//! duaux EMBOÎTÉS sur chaque élément ; K_M, K_f, C_f des duaux sur chaque corps.
//! Les dérivées locales sont automatiques ou analytiques. Les blocs de poutre
//! entre rotations restent approchés ; les cassures et les états internes
//! figés limitent également les garanties de convergence quadratique.
//! Pour sigma!=0, la configuration et sa tangente sont celles de schema_lie :
//! le facteur J_l est remplacé par J_l A^-1 B dans les colonnes u̇,
//! et par J_l A^-1 dans les colonnes de configuration indépendantes GGL.

use crate::ad::{self, Dual, Scalar, M, V};
use crate::{Cardan, Corps, Distance, Elem, Engrenage, Liaison, Modele, Pas, Vis};
use nalgebra::{DMatrix, DVector, SMatrix, SVector};

const PI: f64 = std::f64::consts::PI;
/// en dessous de ce nombre d'éléments, tout est séquentiel (mesuré)
pub const PARALLELE_MIN: usize = 48;

/// pose perturbée : r = r₀ + δr, R = exp([δθ]×)·R₀ — δ porte les tangentes.
fn pose_pert<T: Scalar>(r0: &crate::V3, rot0: &crate::M3, dr: V<T>, dth: V<T>) -> (V<T>, M<T>) {
    let r = [
        T::from_f64(r0.x) + dr[0],
        T::from_f64(r0.y) + dr[1],
        T::from_f64(r0.z) + dr[2],
    ];
    // Duaux du premier ordre à valeur nulle — le cas de toutes les tangentes :
    // exp([δθ]×)·R₀ = R₀ + [δθ]×·R₀ exactement (cf. `ad::expm`), et [δθ]× a une
    // diagonale nulle : 18 produits par constante au lieu d'un expm et de 27
    // produits de duaux. Mesuré 7 sept. : expm<Dual> + m_mul pesaient 24 % du
    // résidu de Princeton.
    if T::ORDRE == 1 && dth.iter().all(|c| c.val() == 0.0) {
        let k = ad::skew(dth);
        let mut rot = ad::m_from::<T>(rot0);
        for i in 0..3 {
            for j in 0..3 {
                for l in 0..3 {
                    if l != i {
                        rot[i][j] = rot[i][j] + k[i][l].scale(rot0[(l, j)]);
                    }
                }
            }
        }
        return (r, rot);
    }
    let e = ad::expm(dth);
    let rot = ad::m_mul(&e, &ad::m_from::<T>(rot0));
    (r, rot)
}

type PoseElement<T> = (V<T>, M<T>, V<T>);

/// Différence haute et basse de deux positions. L'erreur d'arrondi de la
/// soustraction des parties hautes n'est pas une déformation physique.
pub(crate) fn difference_positions<T: Scalar>(
    a: V<T>,
    ab: V<T>,
    b: V<T>,
    bb: V<T>,
) -> (V<T>, V<T>) {
    let h = ad::v_sub(b, a);
    let bas = std::array::from_fn(|i| {
        let e = crate::numerique::deux_sommes(b[i].val(), -a[i].val()).1;
        (bb[i] - ab[i]) + T::from_f64(e)
    });
    (h, bas)
}

/// ||haut+bas||-L : ajouter la petite différence de normes après -L.
pub(crate) fn distance_ecart<T: Scalar>(haut: V<T>, bas: V<T>, longueur: f64) -> T {
    let l = ad::v_norm(haut);
    let complet = ad::v_norm(ad::v_add(haut, bas));
    if l.val() == 0.0 {
        return complet - T::from_f64(longueur);
    }
    (l - T::from_f64(longueur))
        + (ad::v_dot(haut, bas).scale(2.0) + ad::v_dot(bas, bas)) / (complet + l)
}

fn pose_de<T: Scalar>(
    corps: &[Corps],
    c: crate::Ref,
    base: usize,
    seed: &dyn Fn(usize) -> T,
) -> PoseElement<T> {
    match c {
        Some(i) => {
            let p = pose_pert(
                &corps[i].r,
                &corps[i].rot,
                [seed(base), seed(base + 1), seed(base + 2)],
                [seed(base + 3), seed(base + 4), seed(base + 5)],
            );
            (p.0, p.1, ad::v_from(corps[i].r_bas.into()))
        }
        None => ([T::zero(); 3], ad::m_ident(), [T::zero(); 3]),
    }
}

// ── Φ générique par élément ──────────────────────────────────────────────────
impl Liaison {
    pub fn phi_t<T: Scalar>(&self, pa: &PoseElement<T>, pb: &PoseElement<T>, t: f64) -> Vec<T> {
        let cr = match &self.cible_r {
            Some((axe, loi)) => crate::expm(&(axe * loi.valeur(t))),
            None => crate::M3::identity(),
        };
        let ct = self
            .cible_t
            .as_ref()
            .map_or(crate::V3::zeros(), |(dir, loi)| dir * loi.valeur(t));
        self.phi_cibles(pa, pb, &ad::m_from::<T>(&cr), ad::v_from(ct.into()))
    }

    fn phi_cibles<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        cr: &M<T>,
        ct: V<T>,
    ) -> Vec<T> {
        let ua = ad::m_mul(&pa.1, &ad::m_from::<T>(&self.ra));
        let ub = ad::m_mul(&pb.1, &ad::m_from::<T>(&self.rb));
        let qa = ad::m_v(&pa.1, ad::v_from::<T>([self.pa.x, self.pa.y, self.pa.z]));
        // NON HOLONOME : le bras vers B est le point de contact GÉOMÉTRIQUE,
        // comme dans `phi_g` — sans quoi `d_gu` dérivait une autre fonction
        // que celle du résidu (audit du 5 sept. : 9e-4 sur la ligne nh)
        let qb = if self.nh {
            ad::v_add(ad::v_sub(ad::v_add(pa.0, qa), pb.0), ad::v_sub(pa.2, pb.2))
        } else {
            ad::m_v(&pb.1, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]))
        };
        let dw = ad::v_sub(ad::v_add(pb.0, qb), ad::v_add(pa.0, qa));
        let db = ad::m_tv(&ua, ad::v_sub(pb.2, pa.2));
        let d = ad::m_tv(&ua, dw);
        let e = ad::m_mul(&ad::m_mul(&ad::m_t(&ua), &ub), &ad::m_t(cr));
        let half = T::from_f64(0.5);
        let et = ad::m_t(&e);
        let mut sk = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sk[i][j] = half * (e[i][j] - et[i][j]);
            }
        }
        let phi_r = ad::vee(&sk);
        let mut out = Vec::with_capacity(self.bt.len() + self.br.len());
        for &i in &self.bt {
            out.push((d[i] - ct[i]) + db[i]);
        }
        for &i in &self.br {
            out.push(phi_r[i]);
        }
        out
    }
}

impl Liaison {
    /// G(q)·u d'une liaison NON HOLONOME, générique en T : la vitesse relative
    /// du point de contact GÉOMÉTRIQUE dans le repère de la liaison, et la
    /// vitesse angulaire relative projetée — exactement ce que `phi_g` assemble
    /// ligne à ligne. Une ligne nh n'a PAS de Φ (dw ≡ 0 par construction) :
    /// c'est cette fonction qu'il faut dériver en pose, pas Φ.
    pub fn gu_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        t: f64,
        u: [[f64; 3]; 4],
    ) -> Vec<T> {
        let cr = match &self.cible_r {
            Some((axe, loi)) => crate::expm(&(axe * loi.valeur(t))),
            None => crate::M3::identity(),
        };
        self.gu_cible(pa, pb, &ad::m_from::<T>(&cr), u)
    }

    fn gu_cible<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        cr: &M<T>,
        u: [[f64; 3]; 4],
    ) -> Vec<T> {
        let ua = ad::m_mul(&pa.1, &ad::m_from::<T>(&self.ra));
        let ub = ad::m_mul(&pb.1, &ad::m_from::<T>(&self.rb));
        let qa = ad::m_v(&pa.1, ad::v_from::<T>([self.pa.x, self.pa.y, self.pa.z]));
        let qb = ad::v_add(ad::v_sub(ad::v_add(pa.0, qa), pb.0), ad::v_sub(pa.2, pb.2));
        let (va, wa, vb, wb) = (
            ad::v_from::<T>(u[0]),
            ad::v_from::<T>(u[1]),
            ad::v_from::<T>(u[2]),
            ad::v_from::<T>(u[3]),
        );
        let vpa = ad::v_add(va, ad::v_cross(wa, qa));
        let vpb = ad::v_add(vb, ad::v_cross(wb, qb));
        let d = ad::m_tv(&ua, ad::v_sub(vpb, vpa));
        let e = ad::m_mul(&ad::m_mul(&ad::m_t(&ua), &ub), &ad::m_t(cr));
        let half = T::from_f64(0.5);
        let wrel = ad::m_tv(&ua, ad::v_sub(wb, wa));
        // dwm : colonne k = vee(skew([e_k]× E)) ; gr = dwm · uaᵀ(ω_b − ω_a)
        let mut gr = [T::zero(); 3];
        for k in 0..3 {
            let mut ek = [T::zero(); 3];
            ek[k] = T::one();
            let sk = ad::m_mul(&ad::skew(ek), &e);
            let skt = ad::m_t(&sk);
            let mut a = [[T::zero(); 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] = half * (sk[i][j] - skt[i][j]);
                }
            }
            let col = ad::vee(&a);
            for i in 0..3 {
                gr[i] = gr[i] + col[i] * wrel[k];
            }
        }
        let mut out = Vec::with_capacity(self.bt.len() + self.br.len());
        for &i in &self.bt {
            out.push(d[i]);
        }
        for &i in &self.br {
            out.push(gr[i]);
        }
        out
    }
}

impl Liaison {
    /// ∂Φ/∂t d'une cible RHÉONOME, générique en T — le miroir de `phi_dt`
    /// (lib.rs) : la partie rotation dépend de la pose par E.
    pub fn phit_t<T: Scalar>(&self, pa: &PoseElement<T>, pb: &PoseElement<T>, t: f64) -> Vec<T> {
        let cr = match &self.cible_r {
            Some((axe, loi)) => crate::expm(&(axe * loi.valeur(t))),
            None => crate::M3::identity(),
        };
        self.phit_cible(pa, pb, &ad::m_from::<T>(&cr), t)
    }

    fn phit_cible<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        cr: &M<T>,
        t: f64,
    ) -> Vec<T> {
        let n = self.bt.len() + self.br.len();
        let mut out = vec![T::zero(); n];
        if let Some((dir, loi)) = &self.cible_t {
            for (row, &i) in self.bt.iter().enumerate() {
                out[row] = T::from_f64(-dir[i] * loi.derivee(t));
            }
        }
        if let Some((axe, loi)) = &self.cible_r {
            let ua = ad::m_mul(&pa.1, &ad::m_from::<T>(&self.ra));
            let ub = ad::m_mul(&pb.1, &ad::m_from::<T>(&self.rb));
            let e = ad::m_mul(&ad::m_mul(&ad::m_t(&ua), &ub), &ad::m_t(cr));
            let sa = ad::m_from::<T>(&(crate::skew(axe) * loi.derivee(t)));
            let de0 = ad::m_mul(&e, &sa);
            let mut de = [[T::zero(); 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    de[i][j] = T::zero() - de0[i][j];
                }
            }
            let det = ad::m_t(&de);
            let half = T::from_f64(0.5);
            let mut a = [[T::zero(); 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] = half * (de[i][j] - det[i][j]);
                }
            }
            let v = ad::vee(&a);
            for (k, &i) in self.br.iter().enumerate() {
                out[self.bt.len() + k] = v[i];
            }
        }
        out
    }

    // Jet temporel dans le segment courant de la loi (constante, affine
    // ou table affine par morceaux). Une cassure de pente n'a pas de
    // dérivée seconde classique ; on conserve la branche de Loi::derivee.
    fn cibles_jet<T: Scalar>(&self, t: f64, dt: T) -> (M<T>, V<T>) {
        let loi = |l: &crate::Loi| T::from_f64(l.valeur(t)) + dt.scale(l.derivee(t));
        let cr = self.cible_r.as_ref().map_or(ad::m_ident(), |(axe, l)| {
            ad::expm(std::array::from_fn(|i| loi(l).scale(axe[i])))
        });
        let ct = self.cible_t.as_ref().map_or([T::zero(); 3], |(dir, l)| {
            std::array::from_fn(|i| loi(l).scale(dir[i]))
        });
        (cr, ct)
    }
}

impl Modele {
    /// Terme convectif de l'accélération des contraintes, à u fixé :
    /// d²/ds² Φ(q ⊕ su, t+s), ou d/ds (G(q ⊕ su,t+s)u + Φ_t)
    /// pour les lignes non holonomes. Un seul axe dual, sans Hessienne
    /// globale ni perturbation des poses du modèle.
    pub(crate) fn biais_acceleration(&self, t: f64) -> Result<DVector<f64>, String> {
        crate::fini(&[t], "temps d'initialisation")?;
        self.verifie_etat_fini()?;
        type D = Dual<f64, 1>;
        type J = Dual<D, 1>;
        let s = J {
            v: D::var(0., 0),
            d: [D::cst(1.)],
        };
        let mut out = DVector::zeros(self.m());
        let mut row = 0;
        for e in &self.elems {
            let refs = e.corps();
            let immobile = refs.iter().flatten().all(|&i| {
                self.corps[i].v == crate::V3::zeros() && self.corps[i].w == crate::V3::zeros()
            });
            let commande = matches!(e, Elem::L(l) if l.cible_r.iter().chain(l.cible_t.iter())
                .any(|(_, loi)| loi.derivee(t) != 0.));
            if immobile && !commande
                || matches!(e, Elem::K(k) if !matches!(self.k_pas.get(k.ct), Some((true, _))))
            {
                row += e.n();
                continue;
            }
            let non_holonome = match e {
                Elem::L(l) if l.nh => Some(l),
                _ => None,
            };
            let valeurs = if let Some(l) = non_holonome {
                let dt = D::var(0., 0);
                let pose = |c: crate::Ref| {
                    pose_de(&self.corps, c, 0, &|i| {
                        dt.scale(c.map_or(0., |k| self.corps[k].u6(i)))
                    })
                };
                let (pa, pb) = (pose(l.a), pose(l.b));
                let vitesse = |c: crate::Ref, start| {
                    std::array::from_fn(|i| c.map_or(0., |k| self.corps[k].u6(start + i)))
                };
                let (cr, _) = l.cibles_jet(t, dt);
                let gu = l.gu_cible(
                    &pa,
                    &pb,
                    &cr,
                    [
                        vitesse(l.a, 0),
                        vitesse(l.a, 3),
                        vitesse(l.b, 0),
                        vitesse(l.b, 3),
                    ],
                );
                let pt = l.phit_cible(&pa, &pb, &cr, t);
                gu.into_iter()
                    .zip(pt)
                    .map(|(g, p)| (g + p).d[0])
                    .collect::<Vec<_>>()
            } else if let Elem::K(k) = e {
                let ct = &self.contacts[k.ct];
                let ca = &self.corps[ct.corps];
                let cb = ct.b.map(|i| &self.corps[i]);
                let vitesse = |c: Option<&Corps>, start| {
                    std::array::from_fn(|i| J::from_f64(c.map_or(0., |b| b.u6(start + i))))
                };
                let (va, wa, vb, wb) = (
                    vitesse(Some(ca), 0),
                    vitesse(Some(ca), 3),
                    vitesse(cb, 0),
                    vitesse(cb, 3),
                );
                let (_, _, _, _, phi, _) = ct.eval_t(
                    ca,
                    cb,
                    ad::v_scale(s, va),
                    ad::v_scale(s, wa),
                    va,
                    wa,
                    ad::v_scale(s, vb),
                    ad::v_scale(s, wb),
                    vb,
                    wb,
                );
                vec![phi.d[0].d[0]]
            } else {
                let poses: Vec<_> = refs
                    .iter()
                    .map(|&c| {
                        pose_de(&self.corps, c, 0, &|i| {
                            s.scale(c.map_or(0., |k| self.corps[k].u6(i)))
                        })
                    })
                    .collect();
                let phi = if let Elem::L(l) = e {
                    let (cr, ct) = l.cibles_jet(t, s);
                    l.phi_cibles(&poses[0], &poses[1], &cr, ct)
                } else {
                    e.phi_poses(&poses, t)
                };
                phi.into_iter().map(|p| p.d[0].d[0]).collect()
            };
            out.rows_mut(row, e.n())
                .copy_from(&DVector::from_vec(valeurs));
            row += e.n();
        }
        crate::fini(out.as_slice(), "courbure des contraintes initiales")?;
        Ok(out)
    }
}

impl Distance {
    pub fn phi_t<T: Scalar>(&self, pa: &PoseElement<T>, pb: &PoseElement<T>) -> Vec<T> {
        let qa = ad::m_v(&pa.1, ad::v_from::<T>([self.pa.x, self.pa.y, self.pa.z]));
        let qb = ad::m_v(&pb.1, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]));
        let dw = ad::v_sub(ad::v_add(pb.0, qb), ad::v_add(pa.0, qa));
        vec![distance_ecart(dw, ad::v_sub(pb.2, pa.2), self.l)]
    }
}

impl Engrenage {
    /// angle déroulé autour de n : la branche d'atan2 la plus proche de `prev`
    pub(crate) fn angle_t<T: Scalar>(e: &M<T>, n: [f64; 3], prev: f64) -> T {
        let half = T::from_f64(0.5);
        let et = ad::m_t(e);
        let mut sk = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sk[i][j] = half * (e[i][j] - et[i][j]);
            }
        }
        let s = ad::v_dot(ad::vee(&sk), ad::v_from::<T>(n));
        let c = half * (ad::m_trace(e) - T::one());
        let w = s.atan2(c);
        let tours = ((prev - w.val()) / (2.0 * PI)).round();
        w + T::from_f64(2.0 * PI * tours)
    }
    pub fn phi_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        pc: &PoseElement<T>,
    ) -> Vec<T> {
        let rct = ad::m_t(&pc.1);
        let ea = ad::m_mul(
            &ad::m_mul(&rct, &pa.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_a)),
        );
        let eb = ad::m_mul(
            &ad::m_mul(&rct, &pb.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_b)),
        );
        let tha = Self::angle_t(&ea, [self.na.x, self.na.y, self.na.z], self.prev.0);
        let thb = Self::angle_t(&eb, [self.nb.x, self.nb.y, self.nb.z], self.prev.1);
        vec![tha - T::from_f64(self.rapport) * thb]
    }
}

impl Engrenage {
    /// GRADIENT EXACT de l'angle déroulé autour de n, pour E → exp(ε̂)·E :
    ///
    /// ```text
    ///   dθ/dε = [½c·(tr E·n − Eᵀn) + s·a] / (s² + c²),
    ///   a = vee(skew E), s = a·n, c = (tr E − 1)/2
    /// ```
    ///
    /// par l'identité vee(ε̂M + Mᵀε̂) = (tr M·I − M)ε. Pour une rotation PURE
    /// autour de n il vaut n (la formule à la main) ; incliné, il porte des
    /// composantes transverses que l'audit du 5 sept. a mesurées (7e-4). En
    /// forme fermée, il coûte trois produits 3×3 — les duaux emboîtés 18×18
    /// essayés d'abord multipliaient le pas par 15 (mesuré : 0,022 → 0,338 ms).
    pub(crate) fn dtheta_t<T: Scalar>(e: &M<T>, n: [f64; 3]) -> V<T> {
        let half = T::from_f64(0.5);
        let et = ad::m_t(e);
        let mut sk = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sk[i][j] = half * (e[i][j] - et[i][j]);
            }
        }
        let a = ad::vee(&sk);
        let nn = ad::v_from::<T>(n);
        let s = ad::v_dot(a, nn);
        let tr = ad::m_trace(e);
        let c = half * (tr - T::one());
        let etn = ad::m_v(&et, nn);
        let w = ad::v_add(
            ad::v_scale(half * c, ad::v_sub(ad::v_scale(tr, nn), etn)),
            ad::v_scale(s, a),
        );
        ad::v_scale(T::one() / (s * s + c * c), w)
    }

    /// G exact générique (1 × 18 : a 0..6, b 6..12, porteur 12..18)
    pub fn g_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        pc: &PoseElement<T>,
    ) -> [T; 18] {
        let rct = ad::m_t(&pc.1);
        let ea = ad::m_mul(
            &ad::m_mul(&rct, &pa.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_a)),
        );
        let eb = ad::m_mul(
            &ad::m_mul(&rct, &pb.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_b)),
        );
        // ε = R_cᵀ δ  (E → exp(ε̂)E) ⇒ gradient en δ : R_c · dθ/dε
        let wa = ad::m_v(
            &pc.1,
            Self::dtheta_t(&ea, [self.na.x, self.na.y, self.na.z]),
        );
        let wb = ad::m_v(
            &pc.1,
            Self::dtheta_t(&eb, [self.nb.x, self.nb.y, self.nb.z]),
        );
        let r = T::from_f64(self.rapport);
        let mut g = [T::zero(); 18];
        for c in 0..3 {
            g[3 + c] = wa[c];
            g[9 + c] = T::zero() - r * wb[c];
            g[15 + c] = r * wb[c] - wa[c];
        }
        g
    }

    pub fn g_exact(&self, corps: &[Corps]) -> DMatrix<f64> {
        let z = |_: usize| 0.0f64;
        let pa = pose_de::<f64>(corps, self.a, 0, &z);
        let pb = pose_de::<f64>(corps, self.b, 6, &z);
        let pc = pose_de::<f64>(corps, self.c, 12, &z);
        let gt = self.g_t(&pa, &pb, &pc);
        let mut g = DMatrix::zeros(1, 18);
        for k in 0..18 {
            g[(0, k)] = gt[k];
        }
        g
    }
}

impl Vis {
    /// G exact générique (1 × 12) : translation ±n_w, rotation de a
    /// n_w × d − R_a·dθ, rotation de b + R_a·dθ, moins pas·(…) sur l'angle —
    /// cf. `Engrenage::dtheta_t`.
    pub fn g_t<T: Scalar>(&self, pa: &PoseElement<T>, pb: &PoseElement<T>) -> [T; 12] {
        let e = ad::m_mul(
            &ad::m_mul(&ad::m_t(&pa.1), &pb.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_rel)),
        );
        let n = [self.n.x, self.n.y, self.n.z];
        // E → exp(ε̂)E avec ε = R_aᵀ δ_b (et −R_aᵀ δ_a) ⇒ gradient R_a·dθ/dε
        let w = ad::m_v(&pa.1, Engrenage::dtheta_t(&e, n));
        let nw = ad::m_v(&pa.1, ad::v_from::<T>(n));
        let d = ad::v_add(ad::v_sub(pb.0, pa.0), ad::v_sub(pb.2, pa.2));
        let croix = ad::v_cross(nw, d);
        let pas = T::from_f64(self.pas);
        let mut g = [T::zero(); 12];
        for c in 0..3 {
            g[c] = T::zero() - nw[c];
            g[3 + c] = croix[c] + pas * w[c];
            g[6 + c] = nw[c];
            g[9 + c] = T::zero() - pas * w[c];
        }
        g
    }

    pub fn g_exact(&self, corps: &[Corps]) -> DMatrix<f64> {
        let z = |_: usize| 0.0f64;
        let pa = pose_de::<f64>(corps, self.a, 0, &z);
        let pb = pose_de::<f64>(corps, self.b, 6, &z);
        let gt = self.g_t(&pa, &pb);
        let mut g = DMatrix::zeros(1, 12);
        for k in 0..12 {
            g[(0, k)] = gt[k];
        }
        g
    }

    /// Gᵀλ en duaux du premier ordre sur la pose (12) — K par `eval_ad`
    pub fn gtl_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        lam: &[f64],
    ) -> [T; 12] {
        let g = self.g_t(pa, pb);
        let l = T::from_f64(lam[0]);
        let mut out = [T::zero(); 12];
        for k in 0..12 {
            out[k] = l * g[k];
        }
        out
    }

    /// Φ = n·(r_b − r_a) − pas·θ − c₀, en duaux : même déroulement d'angle que
    /// l'engrenage, plus le terme de translation.
    pub fn phi_t<T: Scalar>(&self, pa: &PoseElement<T>, pb: &PoseElement<T>) -> Vec<T> {
        let e = ad::m_mul(
            &ad::m_mul(&ad::m_t(&pa.1), &pb.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_rel)),
        );
        let th = Engrenage::angle_t(&e, [self.n.x, self.n.y, self.n.z], self.prev);
        let nw = ad::m_v(&pa.1, ad::v_from::<T>([self.n.x, self.n.y, self.n.z]));
        let (d, db) = difference_positions(pa.0, pa.2, pb.0, pb.2);
        vec![
            (ad::v_dot(nw, d) - T::from_f64(self.c0)) - T::from_f64(self.pas) * th
                + ad::v_dot(nw, db),
        ]
    }
}

impl Cardan {
    /// Φ = (R_a·n_a)·(R_b·n_b) — deux rotations, un produit scalaire.
    pub fn phi_t<T: Scalar>(&self, pa: &PoseElement<T>, pb: &PoseElement<T>) -> Vec<T> {
        let ua = ad::m_v(&pa.1, ad::v_from::<T>([self.na.x, self.na.y, self.na.z]));
        let ub = ad::m_v(&pb.1, ad::v_from::<T>([self.nb.x, self.nb.y, self.nb.z]));
        vec![ad::v_dot(ua, ub)]
    }
}

fn m_add<T: Scalar>(a: &M<T>, b: &M<T>) -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = a[i][j] + b[i][j];
        }
    }
    o
}
fn m_neg<T: Scalar>(a: &M<T>) -> M<T> {
    let mut o = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            o[i][j] = -a[i][j];
        }
    }
    o
}
/// ligne i d'une matrice 3×3
fn ligne<T: Scalar>(a: &M<T>, i: usize) -> V<T> {
    a[i]
}

impl Liaison {
    /// Gᵀλ (12 composantes : δr_a, δθ_a, δr_b, δθ_b), le G ANALYTIQUE de
    /// `phi_g` écrit sur tout Scalar — dérivé une fois par duaux du premier
    /// ordre, il donne K = ∂(Gᵀλ)/∂δ exact pour 12 tangentes au lieu de 144.
    pub fn gtl_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        t: f64,
        lam: &[f64],
    ) -> [T; 12] {
        let ua = ad::m_mul(&pa.1, &ad::m_from::<T>(&self.ra));
        let ub = ad::m_mul(&pb.1, &ad::m_from::<T>(&self.rb));
        let qa = ad::m_v(&pa.1, ad::v_from::<T>([self.pa.x, self.pa.y, self.pa.z]));
        // Même point géométrique mobile que phi_g pour le roulement :
        // Gᵀλ dépend de ce bras, et sa dérivée doit également le suivre.
        let qb = if self.nh {
            ad::v_add(ad::v_sub(ad::v_add(pa.0, qa), pb.0), ad::v_sub(pa.2, pb.2))
        } else {
            ad::m_v(&pb.1, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]))
        };
        let dw = ad::v_add(
            ad::v_sub(ad::v_add(pb.0, qb), ad::v_add(pa.0, qa)),
            ad::v_sub(pb.2, pa.2),
        );
        let cr = match &self.cible_r {
            Some((axe, loi)) => crate::expm(&(axe * loi.valeur(t))),
            None => crate::M3::identity(),
        };
        let uat = ad::m_t(&ua);
        let gt_ra = m_neg(&uat);
        let gt_ta = ad::m_mul(&uat, &m_add(&ad::skew(dw), &ad::skew(qa)));
        let gt_rb = uat;
        let gt_tb = m_neg(&ad::m_mul(&uat, &ad::skew(qb)));
        let e = ad::m_mul(&ad::m_mul(&uat, &ub), &ad::m_t(&ad::m_from::<T>(&cr)));
        let half = T::from_f64(0.5);
        let mut dwm = [[T::zero(); 3]; 3];
        for k in 0..3 {
            let mut ek = [T::zero(); 3];
            ek[k] = T::one();
            let sk = ad::m_mul(&ad::skew(ek), &e);
            let skt = ad::m_t(&sk);
            let mut a = [[T::zero(); 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] = half * (sk[i][j] - skt[i][j]);
                }
            }
            let col = ad::vee(&a);
            for i in 0..3 {
                dwm[i][k] = col[i];
            }
        }
        let gr_ta = m_neg(&ad::m_mul(&dwm, &uat));
        let gr_tb = ad::m_mul(&dwm, &uat);
        let mut out = [T::zero(); 12];
        let mut row = 0;
        for &i in &self.bt {
            let l = T::from_f64(lam[row]);
            let (r0, r1, r2, r3) = (
                ligne(&gt_ra, i),
                ligne(&gt_ta, i),
                ligne(&gt_rb, i),
                ligne(&gt_tb, i),
            );
            for c in 0..3 {
                out[c] = out[c] + l * r0[c];
                out[3 + c] = out[3 + c] + l * r1[c];
                out[6 + c] = out[6 + c] + l * r2[c];
                out[9 + c] = out[9 + c] + l * r3[c];
            }
            row += 1;
        }
        for &i in &self.br {
            let l = T::from_f64(lam[row]);
            let (r1, r3) = (ligne(&gr_ta, i), ligne(&gr_tb, i));
            for c in 0..3 {
                out[3 + c] = out[3 + c] + l * r1[c];
                out[9 + c] = out[9 + c] + l * r3[c];
            }
            row += 1;
        }
        out
    }
}

impl Distance {
    pub fn gtl_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        lam: &[f64],
    ) -> [T; 12] {
        let qa = ad::m_v(&pa.1, ad::v_from::<T>([self.pa.x, self.pa.y, self.pa.z]));
        let qb = ad::m_v(&pb.1, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]));
        let dw = ad::v_add(
            ad::v_sub(ad::v_add(pb.0, qb), ad::v_add(pa.0, qa)),
            ad::v_sub(pb.2, pa.2),
        );
        let l = T::from_f64(lam[0]);
        let nrm = ad::v_scale(T::one() / ad::v_norm(dw), dw);
        let ta = ad::v_cross(nrm, qa); // ([qa]×)ᵀ n = n × qa
        let tb = ad::v_cross(qb, nrm); // −([qb]×)ᵀ n = qb × n
        let mut out = [T::zero(); 12];
        for c in 0..3 {
            out[c] = -l * nrm[c];
            out[3 + c] = l * ta[c];
            out[6 + c] = l * nrm[c];
            out[9 + c] = l * tb[c];
        }
        out
    }
}

impl Engrenage {
    pub fn gtl_t<T: Scalar>(
        &self,
        pa: &PoseElement<T>,
        pb: &PoseElement<T>,
        pc: &PoseElement<T>,
        lam: &[f64],
    ) -> [T; 18] {
        let g = self.g_t(pa, pb, pc);
        let l = T::from_f64(lam[0]);
        let mut out = [T::zero(); 18];
        for k in 0..18 {
            out[k] = l * g[k];
        }
        out
    }
}

// ── couples à loi ─────────────────────────────────────────────────────────────
/// Projection sur un intervalle : à une cassure, choisir la moyenne des
/// deux pentes (élément du jacobien généralisé). La valeur est inchangée.
/// Cela conserve la convention centrée des tangentes d'analyse historiques,
/// notamment au repos lorsque le plafond de vitesse d'un servo s'active.
fn borne_couple<T: Scalar>(x: T, lo: f64, hi: f64) -> T {
    if lo == hi {
        // Un intervalle réduit à un point est constant, même à la borne.
        // Le gouverneur autorise q_max = 0 : sa dérivée doit rester nulle.
        T::from_f64(lo)
    } else if x.val() == lo || x.val() == hi {
        let valeur = T::from_f64(x.val());
        valeur + (x - valeur).scale(0.5)
    } else {
        x.max(T::from_f64(lo)).min(T::from_f64(hi))
    }
}

impl crate::Couple {
    /// (θ déroulé, ω relative, τ) génériques ; θ mesuré autour de l'axe dans le
    /// repère de a depuis la pose initiale, branche d'atan2 la plus proche de `prev`.
    pub fn eval_t<T: Scalar>(
        &self,
        pa: &(V<T>, M<T>),
        pb: &(V<T>, M<T>),
        wa: V<T>,
        wb: V<T>,
        t: f64,
    ) -> (T, T, T, V<T>) {
        let e = ad::m_mul(
            &ad::m_mul(&ad::m_t(&pa.1), &pb.1),
            &ad::m_t(&ad::m_from::<T>(&self.ref_rel)),
        );
        let axe = [self.axe.x, self.axe.y, self.axe.z];
        let th = Engrenage::angle_t(&e, axe, self.prev);
        let nw = ad::m_v(&pa.1, ad::v_from::<T>(axe));
        let w = ad::v_dot(nw, ad::v_sub(wb, wa));
        let tau = match &self.loi {
            crate::LoiCouple::Constant(c) => T::from_f64(*c),
            crate::LoiCouple::Ressort { k, c, theta0 } => {
                -T::from_f64(*k) * (th - T::from_f64(*theta0)) - T::from_f64(*c) * w
            }
            crate::LoiCouple::Butee {
                k,
                c,
                theta_min,
                theta_max,
            } => {
                // dépassement SIGNÉ : positif au-delà de θ_max, négatif sous
                // θ_min, nul entre — et l'amortissement suit |d|, pas ω seul.
                let d_hi = borne_couple(th - T::from_f64(*theta_max), 0.0, f64::INFINITY);
                let d_lo = borne_couple(th - T::from_f64(*theta_min), f64::NEG_INFINITY, 0.0);
                let d = d_hi + d_lo;
                -T::from_f64(*k) * d - T::from_f64(*c) * (d_hi - d_lo) * w
            }
            crate::LoiCouple::PD {
                kp,
                kd,
                q_max,
                w_nl,
                cible,
            } => {
                let q = T::from_f64(*q_max);
                let pd =
                    T::from_f64(*kp) * (T::from_f64(cible.valeur(t)) - th) - T::from_f64(*kd) * w;
                let u = (pd / q).tanh();
                let un = T::one();
                let qp = q * borne_couple(un - w / T::from_f64(*w_nl), 0.0, 1.0);
                let qn = q * borne_couple(un + w / T::from_f64(*w_nl), 0.0, 1.0);
                borne_couple(u, 0.0, f64::INFINITY) * qp
                    + borne_couple(u, f64::NEG_INFINITY, 0.0) * qn
            }
            crate::LoiCouple::Gouverneur { kg, q_max, cible } => borne_couple(
                T::from_f64(*kg) * (T::from_f64(cible.valeur(t)) - w),
                -*q_max,
                *q_max,
            ),
        };
        (th, w, tau, nw)
    }

    fn poses<T: Scalar>(
        &self,
        corps: &[Corps],
        seed: &dyn Fn(usize) -> T,
    ) -> ((V<T>, M<T>), (V<T>, M<T>), V<T>, V<T>) {
        // colonnes : 0..3 δθ_a, 3..6 δθ_b, 6..9 ω_a, 9..12 ω_b (pas de translation)
        let z = || [T::zero(); 3];
        let pa = match self.a {
            Some(i) => pose_pert(&corps[i].r, &corps[i].rot, z(), [seed(0), seed(1), seed(2)]),
            None => (z(), ad::m_ident()),
        };
        let pb = match self.b {
            Some(i) => pose_pert(&corps[i].r, &corps[i].rot, z(), [seed(3), seed(4), seed(5)]),
            None => (z(), ad::m_ident()),
        };
        let wv = |c: crate::Ref, base: usize| -> V<T> {
            match c {
                Some(i) => {
                    let w = corps[i].w;
                    [
                        seed(base) + T::from_f64(w.x),
                        seed(base + 1) + T::from_f64(w.y),
                        seed(base + 2) + T::from_f64(w.z),
                    ]
                }
                None => z(),
            }
        };
        (pa, pb, wv(self.a, 6), wv(self.b, 9))
    }

    /// couple sur b, en monde (f64)
    pub fn valeur(&self, corps: &[Corps], t: f64) -> crate::V3 {
        let seed = |_k: usize| 0.0_f64;
        let (pa, pb, wa, wb) = self.poses::<f64>(corps, &seed);
        let (_, _, tau, nw) = self.eval_t(&pa, &pb, wa, wb, t);
        crate::V3::new(tau * nw[0], tau * nw[1], tau * nw[2])
    }
    pub fn angle(&self, corps: &[Corps]) -> f64 {
        let seed = |_k: usize| 0.0_f64;
        let (pa, pb, wa, wb) = self.poses::<f64>(corps, &seed);
        self.eval_t(&pa, &pb, wa, wb, 0.0).0
    }
    pub fn etat(&self, corps: &[Corps], t: f64) -> (f64, f64, f64) {
        let seed = |_k: usize| 0.0_f64;
        let (pa, pb, wa, wb) = self.poses::<f64>(corps, &seed);
        let (th, w, tau, _) = self.eval_t(&pa, &pb, wa, wb, t);
        (th, w, tau)
    }
    /// ∂(τ·n)/∂(δθ_a, δθ_b, ω_a, ω_b) : 3 × 12, par duaux
    pub fn tangentes(&self, corps: &[Corps], t: f64) -> DMatrix<f64> {
        type D = Dual<f64, 12>;
        let seed = |k: usize| D::var(0.0, k);
        let (pa, pb, wa, wb) = self.poses::<D>(corps, &seed);
        let (_, _, tau, nw) = self.eval_t(&pa, &pb, wa, wb, t);
        let mut out = DMatrix::zeros(3, 12);
        for i in 0..3 {
            let ti = tau * nw[i];
            for k in 0..12 {
                out[(i, k)] = ti.d[k];
            }
        }
        out
    }
}

// ── pales : théorie des tranches quasi-stationnaire ─────────────────────────
pub(crate) fn interp<T: Scalar>(xs: &[f64], ys: &[f64], x: T) -> T {
    // linéaire par morceaux ; hors table, on prolonge par la valeur de bord
    let xv = x.val();
    if xv <= xs[0] {
        return T::from_f64(ys[0]);
    }
    let n = xs.len();
    if xv >= xs[n - 1] {
        return T::from_f64(ys[n - 1]);
    }
    let lo = xs.partition_point(|&v| v <= xv) - 1;
    let (x0, x1, y0, y1) = (xs[lo], xs[lo + 1], ys[lo], ys[lo + 1]);
    T::from_f64(y0) + (x - T::from_f64(x0)) * T::from_f64((y1 - y0) / (x1 - x0))
}

const GAUSS5: [(f64, f64); 5] = [
    (-0.906179845938664, 0.236926885056189),
    (-0.538469310105683, 0.478628670499366),
    (0.0, 0.568888888888889),
    (0.538469310105683, 0.478628670499366),
    (0.906179845938664, 0.236926885056189),
];
const GAUSS3: [(f64, f64); 3] = [
    (-0.774596669241483, 0.555555555555556),
    (0.0, 0.888888888888889),
    (0.774596669241483, 0.555555555555556),
];

impl crate::Pale {
    /// (F, M au CdM, en monde) génériques en (δr, δθ, v, ω) ; inflow figé.
    #[allow(clippy::too_many_arguments)]
    pub fn eval_t<T: Scalar>(
        &self,
        c: &Corps,
        dr: V<T>,
        dth: V<T>,
        v: V<T>,
        w: V<T>,
        inflow: Option<&crate::Inflow>,
        vent: crate::V3,
    ) -> (V<T>, V<T>, T, T) {
        let (position, rot) = pose_pert(&c.r, &c.rot, dr, dth);
        let es = ad::m_v(&rot, ad::v_from::<T>([self.es.x, self.es.y, self.es.z]));
        let ec = ad::m_v(&rot, ad::v_from::<T>([self.ec.x, self.ec.y, self.ec.z]));
        let en = ad::m_v(&rot, ad::v_from::<T>([self.en.x, self.en.y, self.en.z]));
        let p0 = ad::m_v(&rot, ad::v_from::<T>([self.p0.x, self.p0.y, self.p0.z]));
        // vent induit : air vers −axe à v_i ⇒ la pale voit +v_i·axe en vitesse propre
        let (ax, vi): (V<T>, T) = match inflow {
            Some(i) => (
                ad::v_from::<T>([i.axe.x, i.axe.y, i.axe.z]),
                T::from_f64(i.v_i),
            ),
            None => ([T::zero(); 3], T::zero()),
        };
        // Pitt–Peters : l'inflow n'est plus une constante sur le disque. Il
        // faut donc situer CHAQUE station en (r̄, ψ), ce qui demande le centre
        // du disque et une référence d'azimut — d'où ces trois cotes sur
        // `Inflow`. Le repère (e1, e2) est direct autour de l'axe.
        let pitt = inflow.is_some_and(|i| i.pitt);
        let (ctr, e1v, e2v, r_max) = match inflow {
            Some(i) if i.pitt => {
                let e2 = i.axe.cross(&i.e1);
                (
                    ad::v_from::<T>([i.centre.x, i.centre.y, i.centre.z]),
                    ad::v_from::<T>([i.e1.x, i.e1.y, i.e1.z]),
                    ad::v_from::<T>([e2.x, e2.y, e2.z]),
                    (i.aire / std::f64::consts::PI).sqrt().max(1e-9),
                )
            }
            _ => ([T::zero(); 3], [T::zero(); 3], [T::zero(); 3], 1.0),
        };
        let (v1s, v1c) = match inflow {
            Some(i) => (T::from_f64(i.v1s), T::from_f64(i.v1c)),
            None => (T::zero(), T::zero()),
        };
        let mut f_tot = [T::zero(); 3];
        let mut m_tot = [T::zero(); 3];
        let mut pousse = T::zero();
        let mut couple = T::zero();
        let pts: Vec<(f64, f64)> = if self.gauss <= 3 {
            GAUSS3.to_vec()
        } else {
            GAUSS5.to_vec()
        };
        let half_l = T::from_f64(0.5 * self.longueur);
        for (i_st, (xi, wg)) in pts.into_iter().enumerate() {
            let sk = half_l * (T::one() + T::from_f64(xi));
            let r_k = ad::v_add(p0, ad::v_scale(sk, es)); // position de la station, depuis le CdM
            let vk = ad::v_add(v, ad::v_cross(w, r_k)); // vitesse de la station (monde)
                                                        // vent induit : uniforme, ou la distribution 1/rev de Pitt–Peters
            let vi_k = if pitt {
                let pw = ad::v_add(position, r_k);
                let d = ad::v_sub(pw, ctr);
                let (x1, x2) = (ad::v_dot(d, e1v), ad::v_dot(d, e2v));
                // λ(r̄,ψ) = λ₀ + r̄(λ₁c cos ψ + λ₁s sin ψ), et r̄ cos ψ = x1/R
                // exactement : inutile de passer par r̄ et ψ séparément, ce qui
                // évite un atan2 et sa singularité au centre du disque.
                let vh = vi + (x1 * v1c + x2 * v1s) / T::from_f64(r_max);
                // + PROFIL RADIAL imposé (sillage libre), interpolé en r̄
                let vh = match inflow {
                    Some(i) if !i.profil.is_empty() => {
                        let rb =
                            (x1 * x1 + x2 * x2 + T::from_f64(1e-24)).sqrt() / T::from_f64(r_max);
                        let (rs, vs): (Vec<f64>, Vec<f64>) = i.profil.iter().cloned().unzip();
                        vh + interp(&rs, &vs, rb)
                    }
                    _ => vh,
                };
                // + CARTE w(r̄, ψ) imposée (Peters–He, sillage en avancement) :
                // bilinéaire, périodique en ψ — ψ par atan2, gardé de l'origine
                match inflow {
                    Some(i) if !i.carte.2.is_empty() => {
                        let (rs, ps, ws) = (&i.carte.0, &i.carte.1, &i.carte.2);
                        let rb =
                            (x1 * x1 + x2 * x2 + T::from_f64(1e-24)).sqrt() / T::from_f64(r_max);
                        let psi = x2.atan2(x1 + T::from_f64(1e-300));
                        let two_pi = T::from_f64(2.0 * std::f64::consts::PI);
                        let psi = if psi.val() < 0.0 { psi + two_pi } else { psi };
                        // rangée r̄ inférieure et supérieure : interpolation en ψ
                        // (périodique) sur chacune, puis en r̄
                        let np_ = ps.len();
                        let ligne = |ir: usize, psi: T| -> T {
                            let mut pp: Vec<f64> = ps.clone();
                            let mut ww: Vec<f64> = ws[ir * np_..(ir + 1) * np_].to_vec();
                            pp.push(ps[0] + 2.0 * std::f64::consts::PI);
                            ww.push(ws[ir * np_]);
                            let psi = if psi.val() < ps[0] { psi + two_pi } else { psi };
                            interp(&pp, &ww, psi)
                        };
                        let rbv = rb.val().clamp(rs[0], *rs.last().unwrap());
                        let k = rs
                            .iter()
                            .rposition(|&r| r <= rbv)
                            .unwrap_or(0)
                            .min(rs.len() - 2);
                        let (r0, r1) = (rs[k], rs[k + 1]);
                        let t = ((rb - T::from_f64(r0)) / T::from_f64(r1 - r0))
                            .max(T::zero())
                            .min(T::one());
                        let (w0, w1) = (ligne(k, psi), ligne(k + 1, psi));
                        vh + w0 + t * (w1 - w0)
                    }
                    _ => vh,
                }
            } else {
                vi
            };
            let vk = ad::v_add(vk, ad::v_scale(vi_k, ax)); // + vent induit vu comme vitesse propre
                                                           // vent AMBIANT : l'air se déplace à `vent` (monde), la section voit −vent
            let vk = ad::v_sub(vk, ad::v_from::<T>([vent.x, vent.y, vent.z]));
            let ut = ad::v_dot(vk, ec); // avance (bord d'attaque devant)
            let up = ad::v_dot(vk, en); // montée, sur la ligne de référence (c/4)
                                        // 3/4 de corde : seule la part de ROTATION change, ω × (−c/2·ec)
            let up_c = if self.a34 {
                let r34 = ad::v_scale(T::from_f64(-0.5 * self.corde), ec);
                up + ad::v_dot(ad::v_cross(w, r34), en)
            } else {
                up
            };
            let alpha_geo = -up_c.atan2(ut); // rad
            let alpha = alpha_geo;
            // INSTATIONNAIRE : ce qui alimente la table n'est pas l'incidence
            // géométrique mais l'incidence EFFECTIVE de Wagner — l'incidence
            // que verrait une section stationnaire portant autant. Les deux
            // états z sont figés sur le pas, comme l'inflow : le sillage est un
            // état de fluide, pas un degré de liberté mécanique.
            let alpha = if self.instat {
                let (z1, z2) = self.z[i_st];
                alpha * T::from_f64(1.0 - crate::W_A1 - crate::W_A2)
                    + T::from_f64(crate::W_B1 * crate::W_A1 * z1 + crate::W_B2 * crate::W_A2 * z2)
            } else {
                alpha
            };
            let alpha_deg = alpha * T::from_f64(180.0 / std::f64::consts::PI);
            let cl = if self.lb {
                // Leishman–Beddoes : f'' et C_N^v figés sur le pas, forme d'Øye
                let e = &self.etats_lb[i_st];
                let f = T::from_f64(e.ff);
                let att = T::from_f64(self.a0) * (alpha_deg - T::from_f64(self.alpha0));
                f * att
                    + (T::one() - f) * interp(&self.polaire.alpha, &self.cl_fs, alpha_deg)
                    + T::from_f64(e.cnv)
            } else if self.oye {
                // Øye : f figé sur le pas (état de fluide, comme z et l'inflow)
                let f = T::from_f64(self.f[i_st]);
                let att = T::from_f64(self.a0) * (alpha_deg - T::from_f64(self.alpha0));
                f * att + (T::one() - f) * interp(&self.polaire.alpha, &self.cl_fs, alpha_deg)
            } else if self.instat && self.retard_lineaire {
                // retard sur la partie linéaire seule ; la non-linéarité statique
                // Δ(α) = c81(α) − a₀(α − α₀) à l'incidence INSTANTANÉE
                let ag = alpha_geo * T::from_f64(180.0 / std::f64::consts::PI);
                let a0 = T::from_f64(self.a0);
                let al0 = T::from_f64(self.alpha0);
                a0 * (alpha_deg - al0) + interp(&self.polaire.alpha, &self.polaire.cl, ag)
                    - a0 * (ag - al0)
            } else {
                interp(&self.polaire.alpha, &self.polaire.cl, alpha_deg)
            };
            let cd = interp(&self.polaire.alpha, &self.polaire.cd, alpha_deg);
            let cm = interp(&self.polaire.alpha, &self.polaire.cm, alpha_deg);
            let v2 = ut * ut + up * up;
            // pale à l'arrêt : |V| = 0 rend lx, dx indéterminés et la dérivée de
            // √ infinie — plancher à 1 mm/s, sans effet dès que ça tourne
            let vn = (v2 + T::from_f64(1e-6)).sqrt();
            let q = T::from_f64(0.5 * self.rho * self.corde) * v2;
            // portance ⟂ au vent relatif (−ut, −up) côté +en ; traînée le long du vent
            let (lx, ln) = (-up / vn, ut / vn);
            let (dx, dn) = (-ut / vn, -up / vn);
            let fc = q * (cl * lx + cd * dx);
            let fn_ = q * (cl * ln + cd * dn);
            // MASSE AJOUTÉE : L_nc = πρb²·ẇ, w = vitesse normale de l'AIR vue au
            // milieu de corde (c/4 derrière la référence), ẇ en différence arrière
            // sur le pas ; force le long de la normale, comme chez Theodorsen
            let fn_ = if self.nc && self.nc_h > 0.0 {
                let rm = ad::v_scale(T::from_f64(-0.25 * self.corde), ec);
                let w_mid = -(up + ad::v_dot(ad::v_cross(w, rm), en));
                let b = 0.5 * self.corde;
                fn_ + T::from_f64(std::f64::consts::PI * self.rho * b * b / self.nc_h)
                    * (w_mid - T::from_f64(self.nc_w[i_st]))
            } else {
                fn_
            };
            let fk = ad::v_add(ad::v_scale(fc, ec), ad::v_scale(fn_, en));
            let mk = ad::v_add(
                ad::v_cross(r_k, fk),
                ad::v_scale(q * cm * T::from_f64(self.corde), es),
            );
            let wl = T::from_f64(wg) * half_l;
            f_tot = ad::v_add(f_tot, ad::v_scale(wl, fk));
            m_tot = ad::v_add(m_tot, ad::v_scale(wl, mk));
            if inflow.is_some() {
                pousse = pousse + wl * ad::v_dot(fk, ax);
                couple = couple + wl * ad::v_dot(ad::v_cross(r_k, fk), ax);
            }
        }
        (f_tot, m_tot, pousse, couple)
    }

    /// (α géométrique — au 3/4 de corde si `a34` —, |V|, w au milieu de corde)
    /// par station — ce qu'il faut pour faire avancer les états de Wagner, de
    /// LB et de la masse ajoutée. On les relit plutôt que de les faire remonter
    /// par `eval_t` : cette fonction est appelée par l'AD sur des duaux, et un
    /// effet de bord y serait une faute.
    pub fn stations(
        &self,
        corps: &[Corps],
        inflows: &[crate::Inflow],
        vent: crate::V3,
    ) -> Vec<(f64, f64, f64)> {
        let c = &corps[self.corps];
        let inf = self.inflow.map(|i| &inflows[i]);
        let (rot, p0) = (c.rot, c.rot * self.p0);
        let (es, ec, en) = (rot * self.es, rot * self.ec, rot * self.en);
        let (ax, vi) = match inf {
            Some(i) => (i.axe, i.v_i),
            None => (crate::V3::zeros(), 0.0),
        };
        let pts: Vec<(f64, f64)> = if self.gauss <= 3 {
            GAUSS3.to_vec()
        } else {
            GAUSS5.to_vec()
        };
        pts.iter()
            .map(|(xi, _)| {
                let sk = 0.5 * self.longueur * (1.0 + xi);
                let r_k = p0 + es * sk;
                let vk = c.v + c.w.cross(&r_k) + ax * vi - vent;
                let (ut, up) = (vk.dot(&ec), vk.dot(&en));
                let up_c = if self.a34 {
                    up + c.w.cross(&(ec * (-0.5 * self.corde))).dot(&en)
                } else {
                    up
                };
                let w_mid = -(up + c.w.cross(&(ec * (-0.25 * self.corde))).dot(&en));
                (-up_c.atan2(ut), (ut * ut + up * up).sqrt(), w_mid)
            })
            .collect()
    }

    /// (F, M) f64 et (poussée, couple) autour de l'axe d'inflow
    pub fn valeur(
        &self,
        corps: &[Corps],
        inflows: &[crate::Inflow],
        vent: crate::V3,
    ) -> (crate::V3, crate::V3, (f64, f64)) {
        let c = &corps[self.corps];
        let inf = self.inflow.map(|i| &inflows[i]);
        let (f, m, p, q) = self.eval_t::<f64>(
            c,
            [0.0; 3],
            [0.0; 3],
            [c.v.x, c.v.y, c.v.z],
            [c.w.x, c.w.y, c.w.z],
            inf,
            vent,
        );
        (
            crate::V3::new(f[0], f[1], f[2]),
            crate::V3::new(m[0], m[1], m[2]),
            (p, q),
        )
    }

    /// ∂(F, M)/∂(δr, δθ, v, ω) : 6 × 12, par duaux.
    /// Les translations changent l'inflow spatial de Pitt–Peters, les
    /// profils radiaux et les cartes imposées, même à état d'inflow figé.
    pub fn tangentes(
        &self,
        corps: &[Corps],
        inflows: &[crate::Inflow],
        vent: crate::V3,
    ) -> DMatrix<f64> {
        let inf = self.inflow.map(|i| &inflows[i]);
        let spatial = inf.is_some_and(|i| {
            i.pitt
                && (i.v1c != 0.0 || i.v1s != 0.0 || !i.profil.is_empty() || !i.carte.2.is_empty())
        });
        if spatial {
            self.tangentes_n::<12>(&corps[self.corps], inf, vent)
        } else {
            // Trois colonnes exactement nulles : ne pas les transporter
            // dans chaque opération de l'inflow uniforme.
            self.tangentes_n::<9>(&corps[self.corps], inf, vent)
        }
    }

    fn tangentes_n<const N: usize>(
        &self,
        c: &Corps,
        inf: Option<&crate::Inflow>,
        vent: crate::V3,
    ) -> DMatrix<f64> {
        let offset = N - 9;
        let var = |x, k| Dual::<f64, N>::var(x, k);
        let dr = std::array::from_fn(|k| if N == 12 { var(0.0, k) } else { Dual::zero() });
        let (f, m, _, _) = self.eval_t(
            c,
            dr,
            std::array::from_fn(|k| var(0.0, offset + k)),
            std::array::from_fn(|k| var(c.v[k], offset + 3 + k)),
            std::array::from_fn(|k| var(c.w[k], offset + 6 + k)),
            inf,
            vent,
        );
        let mut out = DMatrix::zeros(6, 12);
        for i in 0..3 {
            for k in 0..N {
                out[(i, k + 12 - N)] = f[i].d[k];
                out[(3 + i, k + 12 - N)] = m[i].d[k];
            }
        }
        out
    }
}

impl crate::Contact {
    /// (F_a, M_a, F_b, M_b, δ, F_n) génériques en (δθ, v, ω) des DEUX corps.
    /// Contre un demi-espace fixe, seuls a et δθ_a/v_a/ω_a comptent.
    /// Tout est C¹ en δ, y compris à δ = 0 : c'est ce qui permet à Newton de
    /// traverser l'instant du choc.
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    pub fn eval_t<T: Scalar>(
        &self,
        ca: &Corps,
        cb: Option<&Corps>,
        dra: V<T>,
        dta: V<T>,
        va: V<T>,
        wa: V<T>,
        drb: V<T>,
        dtb: V<T>,
        vb: V<T>,
        wb: V<T>,
    ) -> (V<T>, V<T>, V<T>, V<T>, T, T) {
        // Géométrie dans un repère translaté : utiliser aussi la partie basse
        // avant les opérations locales, plutôt que les seules coordonnées monde.
        let (d, bas) = match cb {
            Some(c2) => ca.difference_precise(c2),
            None => crate::numerique::deplace_compense(ca.r, ca.r_bas, -self.origine),
        };
        let (ra_pos, rota) = pose_pert(&(d + bas), &ca.rot, dra, dta);
        let rb_pos = match cb {
            Some(c2) => pose_pert(&crate::V3::zeros(), &c2.rot, drb, dtb).0,
            None => [T::zero(); 3],
        };
        // CAPSULE contre boîte, cylindre ou maillage : le point de l'axe figé
        // sur le pas (`s_axe`, préparé avec `q_maille`) tient lieu de centre ;
        // s_axe = 0 partout ailleurs, donc p0 au bit près
        let p_ax = self.p0 + (self.p1 - self.p0) * self.s_axe;
        let ra = ad::m_v(&rota, ad::v_from::<T>([p_ax.x, p_ax.y, p_ax.z]));
        let pa = ad::v_add(ra_pos, ra);
        let zero = [T::zero(); 3];
        // n = normale du contact, dirigée du second vers le premier ; d =
        // enfoncement ; rb_ = bras de levier sur le second corps
        let (n, d, rb_, vb_pt) = match cb {
            None => {
                let nn = ad::v_from::<T>([self.normale.x, self.normale.y, self.normale.z]);
                let o = [T::zero(); 3];
                (
                    nn,
                    T::from_f64(self.rayon) - ad::v_dot(ad::v_sub(pa, o), nn),
                    zero,
                    zero,
                )
            }
            // MAILLAGE : géométrie quelconque. Le point le plus proche sort du
            // BVH — donc d'une recherche NON différentiable — puis on traite
            // ce point comme fixe dans le repère du corps porteur. C'est le
            // procédé standard : la distance est différentiable À TRIANGLE
            // FIXE, et le triangle ne change qu'aux frontières de région.
            Some(c2) if self.maille.is_some() => {
                let (_, rotb) = pose_pert(&c2.r, &c2.rot, [T::zero(); 3], dtb);
                let ql = self.q_maille; // point le plus proche, repère b
                let rq = ad::m_v(&rotb, ad::v_from::<T>([ql.x, ql.y, ql.z]));
                let qw = ad::v_add(rb_pos, rq);
                let dv = ad::v_sub(pa, qw);
                let dist = (ad::v_dot(dv, dv) + T::from_f64(1e-24)).sqrt();
                let nn = ad::v_scale(T::one() / dist, dv);
                let dd = T::from_f64(self.rayon) - dist;
                (nn, dd, rq, ad::v_add(vb, ad::v_cross(wb, rq)))
            }
            // CYLINDRE FINI porté par b : dans son repère, on sépare l'axial z
            // (clampé à ±L) du radial ρ (ramené à R). Quatre régimes, un seul
            // code : flanc (|z| ≤ L, ρ ≥ R), fond (|z| > L, ρ < R), arête
            // circulaire (|z| > L, ρ ≥ R) et, non couvert, l'intérieur.
            Some(c2) if self.cylindre.is_some() => {
                let (axe, rc, lc) = self.cylindre.unwrap();
                let (_, rotb) = pose_pert(&c2.r, &c2.rot, [T::zero(); 3], dtb);
                let rb = ad::m_v(&rotb, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]));
                let ctr = ad::v_add(rb_pos, rb);
                let rel = ad::v_sub(pa, ctr);
                let rt = ad::m_t(&rotb);
                let loc = ad::m_v(&rt, rel);
                let e = ad::v_from::<T>([axe.x, axe.y, axe.z]);
                let z = ad::v_dot(loc, e);
                let rv = ad::v_sub(loc, ad::v_scale(z, e));
                let rho = (ad::v_dot(rv, rv) + T::from_f64(1e-24)).sqrt();
                let zc = z.max(T::from_f64(-lc)).min(T::from_f64(lc));
                let q_loc = if rho.val() >= rc {
                    // flanc ou arête : radial ramené à R, axial clampé
                    ad::v_add(ad::v_scale(zc, e), ad::v_scale(T::from_f64(rc) / rho, rv))
                } else if z.val().abs() > lc {
                    // fond : le disque, radial inchangé
                    ad::v_add(ad::v_scale(zc, e), rv)
                } else {
                    // intérieur : non couvert — distance nulle, pas de signe faux
                    loc
                };
                let q = ad::m_v(&rotb, q_loc);
                let qw = ad::v_add(ctr, q);
                let dv = ad::v_sub(pa, qw);
                let dist = (ad::v_dot(dv, dv) + T::from_f64(1e-24)).sqrt();
                let nn = ad::v_scale(T::one() / dist, dv);
                let dd = T::from_f64(self.rayon) - dist;
                let r_cb = ad::v_add(rb, q);
                (nn, dd, r_cb, ad::v_add(vb, ad::v_cross(wb, r_cb)))
            }
            // BOÎTE (OBB) portée par b : le point le plus proche est un CLAMP
            // dans son repère. Trois régimes selon le nombre d'axes saturés :
            // face (1), arête (2), coin (3) — et c'est le même code.
            Some(c2) if self.demi.is_some() => {
                let hd = self.demi.unwrap();
                let (_, rotb) = pose_pert(&c2.r, &c2.rot, [T::zero(); 3], dtb);
                let rb = ad::m_v(&rotb, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]));
                let ctr = ad::v_add(rb_pos, rb);
                let rel = ad::v_sub(pa, ctr);
                let rt = ad::m_t(&rotb);
                let loc = ad::m_v(&rt, rel); // dans le repère boîte
                let h = [T::from_f64(hd.x), T::from_f64(hd.y), T::from_f64(hd.z)];
                let mut cl = [T::zero(); 3];
                for k in 0..3 {
                    cl[k] = loc[k].max(-h[k]).min(h[k]);
                }
                let q = ad::m_v(&rotb, cl); // point le plus proche, monde
                let qw = ad::v_add(ctr, q);
                let dv = ad::v_sub(pa, qw);
                let dist = (ad::v_dot(dv, dv) + T::from_f64(1e-24)).sqrt();
                let nn = ad::v_scale(T::one() / dist, dv);
                let dd = T::from_f64(self.rayon) - dist;
                let r_cb = ad::v_add(rb, q);
                (nn, dd, r_cb, ad::v_add(vb, ad::v_cross(wb, r_cb)))
            }
            Some(c2) if self.rayon_b > 0.0 => {
                let (_, rotb) = pose_pert(&c2.r, &c2.rot, [T::zero(); 3], dtb);
                let rb = ad::m_v(&rotb, ad::v_from::<T>([self.pb.x, self.pb.y, self.pb.z]));
                let pb = ad::v_add(rb_pos, rb);
                // CAPSULE : les deux primitives sont des SEGMENTS dilatés. Le
                // point le plus proche se résout en forme fermée (Ericson,
                // *Real-Time Collision Detection*) ; la sphère est le cas
                // dégénéré p1 = p0, où toutes les branches rendent s = 0.
                //
                // ⚠ Les `clamp` rendent la distance C⁰ mais pas C¹ AUX
                // TRANSITIONS de région (extrémité ↔ intérieur du segment).
                // L'AD y prend la dérivée de la branche active, ce qui est
                // correct presque partout ; le cas dégénéré des segments
                // PARALLÈLES est traité par le repli `denom ≈ 0 → s = 0`.
                let (pa, pb, _ra, rb) = if self.p1 != self.p0 || self.p1b != self.pb {
                    let ra1 = ad::m_v(&rota, ad::v_from::<T>([self.p1.x, self.p1.y, self.p1.z]));
                    let pa1 = ad::v_add(ra_pos, ra1);
                    let rb1 = ad::m_v(&rotb, ad::v_from::<T>([self.p1b.x, self.p1b.y, self.p1b.z]));
                    let pb1 = ad::v_add(rb_pos, rb1);
                    let (d1, d2) = (ad::v_sub(pa1, pa), ad::v_sub(pb1, pb));
                    let r0 = ad::v_sub(pa, pb);
                    let (aa, ee) = (ad::v_dot(d1, d1), ad::v_dot(d2, d2));
                    let (ff, cc, bb) = (ad::v_dot(d2, r0), ad::v_dot(d1, r0), ad::v_dot(d1, d2));
                    let un = T::one();
                    let z = T::zero();
                    let cl = |x: T| x.max(z).min(un);
                    let den = aa * ee - bb * bb;
                    let sp = if den.val().abs() > 1e-14 {
                        cl((bb * ff - cc * ee) / den)
                    } else {
                        z
                    };
                    let tp = if ee.val() > 1e-14 {
                        (bb * sp + ff) / ee
                    } else {
                        z
                    };
                    let (tp, sp) = if tp.val() < 0.0 {
                        (z, if aa.val() > 1e-14 { cl(-cc / aa) } else { z })
                    } else if tp.val() > 1.0 {
                        (
                            un,
                            if aa.val() > 1e-14 {
                                cl((bb - cc) / aa)
                            } else {
                                z
                            },
                        )
                    } else {
                        (tp, sp)
                    };
                    let qa = ad::v_add(pa, ad::v_scale(sp, d1));
                    let qb = ad::v_add(pb, ad::v_scale(tp, d2));
                    (
                        qa,
                        qb,
                        ad::v_add(ra, ad::v_scale(sp, d1)),
                        ad::v_add(rb, ad::v_scale(tp, d2)),
                    )
                } else {
                    (pa, pb, ra, rb)
                };
                let dv = ad::v_sub(pa, pb);
                let dist = (ad::v_dot(dv, dv) + T::from_f64(1e-24)).sqrt();
                let nn = ad::v_scale(T::one() / dist, dv);
                // CÂBLE : d change de signe — il « pénètre » quand la distance
                // DÉPASSE la longueur au repos, et la normale s'inverse avec,
                // donc la force TIRE au lieu de pousser. Le même élément donne
                // ainsi le câble, le ressort unilatéral et la butée de
                // translation ; le « rayon » y devient la longueur au repos.
                let dd = if self.cable {
                    dist - T::from_f64(self.rayon + self.rayon_b)
                } else {
                    T::from_f64(self.rayon + self.rayon_b) - dist
                };
                let nn = if self.cable {
                    ad::v_scale(-T::one(), nn)
                } else {
                    nn
                };
                // point de contact vu du corps b : centre b + R_b·n
                let r_cb = ad::v_add(rb, ad::v_scale(T::from_f64(self.rayon_b), nn));
                (nn, dd, r_cb, ad::v_add(vb, ad::v_cross(wb, r_cb)))
            }
            // PLAN PORTÉ PAR UN CORPS (rayon_b = 0) : la normale et l'origine
            // sont dans le repère de b et le suivent. C'est ce qui rend un
            // contact utile à une butée, une glissière ou une came — un plan
            // figé dans le monde ne sert qu'au sol.
            Some(c2) => {
                let (_, rotb) = pose_pert(&c2.r, &c2.rot, [T::zero(); 3], dtb);
                let nn = ad::m_v(
                    &rotb,
                    ad::v_from::<T>([self.normale.x, self.normale.y, self.normale.z]),
                );
                let ro = ad::m_v(
                    &rotb,
                    ad::v_from::<T>([self.origine.x, self.origine.y, self.origine.z]),
                );
                let o = ad::v_add(rb_pos, ro);
                let dd = T::from_f64(self.rayon) - ad::v_dot(ad::v_sub(pa, o), nn);
                // point de contact ramené au corps b : projection du centre a
                let proj = ad::v_sub(pa, ad::v_scale(ad::v_dot(ad::v_sub(pa, o), nn), nn));
                let r_cb = ad::v_sub(proj, rb_pos);
                (nn, dd, r_cb, ad::v_add(vb, ad::v_cross(wb, r_cb)))
            }
        };
        // BARRIÈRE IPC (Li & al. 2020) : au lieu d'une force qui croît AVEC
        // l'enfoncement, un potentiel qui DIVERGE quand la distance tend vers
        // zéro. Trois propriétés que la pénalité de Hertz n'a pas :
        //   · l'interpénétration devient impossible, pas seulement petite ;
        //   · le support est COMPACT (rien au-delà de d̂), donc le
        //     conditionnement n'est pas dégradé loin du contact ;
        //   · B est C², donc la force est C¹ — l'AD passe au travers.
        //     B(d) = −(d − d̂)² ln(d/d̂)    B'(d) = −2(d−d̂)ln(d/d̂) − (d−d̂)²/d
        if let Some(l) = self.fn_impose {
            // NON LISSE : la normale est Gᵀλ (hors d'ici) ; ne reste que Coulomb sur λ
            let r_ca = ad::v_sub(ra, ad::v_scale(T::from_f64(self.rayon), n));
            let v_rel = ad::v_sub(ad::v_add(va, ad::v_cross(wa, r_ca)), vb_pt);
            let vt = ad::v_sub(v_rel, ad::v_scale(ad::v_dot(v_rel, n), n));
            let vn = (ad::v_dot(vt, vt) + T::from_f64(1e-24)).sqrt();
            let ech = (vn / T::from_f64(self.v_eps)).tanh() * T::from_f64(self.mu * l) / vn;
            let f = ad::v_scale(-ech, vt);
            let m = ad::v_cross(r_ca, f);
            let (fb, mb) = match cb {
                None => (zero, zero),
                Some(_) => {
                    let fneg = ad::v_scale(-T::one(), f);
                    (fneg, ad::v_cross(rb_, fneg))
                }
            };
            return (f, m, fb, mb, d, T::from_f64(l));
        }
        if self.expo < 0.0 {
            let dist = -d; // d > 0 = séparation
            let dh = T::from_f64(self.d_hat);
            if dist.val() >= self.d_hat {
                return (zero, zero, zero, zero, d, T::zero());
            }
            // plancher : sous 1e-4·d̂ la barrière est déjà colossale, et le
            // logarithme cesse d'être représentable
            let dd = dist.max(T::from_f64(1e-4 * self.d_hat));
            let ecart = dd - dh;
            let bp = -T::from_f64(2.0) * ecart * (dd / dh).ln() - ecart * ecart / dd;
            let r_ca = ad::v_sub(ra, ad::v_scale(T::from_f64(self.rayon), n));
            let v_ca = ad::v_add(va, ad::v_cross(wa, r_ca));
            let v_rel = ad::v_sub(v_ca, vb_pt);
            // amortissement PROPORTIONNEL à la force de barrière, comme
            // Hunt–Crossley l'est à δ^e : il s'éteint avec elle, donc il ne
            // colle pas à l'entrée et ne tire pas à la sortie
            let dd_pt = -ad::v_dot(v_rel, n);
            let fn_ = (-T::from_f64(self.k) * bp * (T::one() + T::from_f64(self.c) * dd_pt))
                .max(T::zero());
            let vt = ad::v_sub(v_rel, ad::v_scale(ad::v_dot(v_rel, n), n));
            let vn = (ad::v_dot(vt, vt) + T::from_f64(1e-24)).sqrt();
            let ech = (vn / T::from_f64(self.v_eps)).tanh() * T::from_f64(self.mu) * fn_ / vn;
            let f = ad::v_sub(ad::v_scale(fn_, n), ad::v_scale(ech, vt));
            let m = ad::v_cross(r_ca, f);
            let (fb, mb) = match cb {
                None => (zero, zero),
                Some(_) => {
                    let fneg = ad::v_scale(-T::one(), f);
                    (fneg, ad::v_cross(rb_, fneg))
                }
            };
            return (f, m, fb, mb, d, fn_);
        }
        if d.val() <= 0.0 {
            // force nulle, mais on rend le vrai δ (NÉGATIF) : c'est la
            // DISTANCE à la collision, et le pas piloté par le contact s'en
            // sert pour resserrer AVANT le choc. La renvoyer à zéro rendait le
            // critère « proche » vrai partout, donc le pas fin partout —
            // mesuré : aucun gain, et rien ne criait.
            return (zero, zero, zero, zero, d, T::zero());
        }
        // point de contact vu du corps a : centre a − R_a·n
        let r_ca = ad::v_sub(ra, ad::v_scale(T::from_f64(self.rayon), n));
        let v_ca = ad::v_add(va, ad::v_cross(wa, r_ca));
        let v_rel = ad::v_sub(v_ca, vb_pt);
        let dd_pt = -ad::v_dot(v_rel, n); // δ̇ = −(v_rel · n)
                                          // HUNT–CROSSLEY : l'amortissement est proportionnel à δ^e, donc la
                                          // force part de zéro ET y revient continûment. Un amortisseur linéaire
                                          // c·δ̇ colle à l'impact (force non nulle à δ = 0) et TIRE au
                                          // décollement — c'est la discontinuité que toute la littérature du
                                          // contact régularisé cherche à éviter.
        let dp = d.powf(self.expo);
        let fn_ =
            (T::from_f64(self.k) * dp * (T::one() + T::from_f64(self.c) * dd_pt)).max(T::zero());
        // frottement de Coulomb LISSÉ : tanh(‖v_t‖/v_ε) au lieu du signe, qui
        // est non dérivable en zéro et fait osciller Newton au collage
        let vt = ad::v_sub(v_rel, ad::v_scale(ad::v_dot(v_rel, n), n));
        let vn = (ad::v_dot(vt, vt) + T::from_f64(1e-24)).sqrt();
        let ech = (vn / T::from_f64(self.v_eps)).tanh() * T::from_f64(self.mu) * fn_ / vn;
        let f = ad::v_sub(ad::v_scale(fn_, n), ad::v_scale(ech, vt));
        let m = ad::v_cross(r_ca, f);
        let (fb, mb) = match cb {
            None => (zero, zero),
            Some(_) => {
                let fneg = ad::v_scale(-T::one(), f);
                (fneg, ad::v_cross(rb_, fneg))
            }
        };
        (f, m, fb, mb, d, fn_)
    }

    #[allow(clippy::type_complexity)]
    pub fn valeur(
        &self,
        corps: &[Corps],
    ) -> (crate::V3, crate::V3, crate::V3, crate::V3, (f64, f64)) {
        let ca = &corps[self.corps];
        let cb = self.b.map(|i| &corps[i]);
        let (zb, zv, zw) = match cb {
            None => ([0.0; 3], [0.0; 3], [0.0; 3]),
            Some(c) => ([0.0; 3], [c.v.x, c.v.y, c.v.z], [c.w.x, c.w.y, c.w.z]),
        };
        let (f, m, fb, mb, d, fn_) = self.eval_t::<f64>(
            ca,
            cb,
            [0.0; 3],
            [0.0; 3],
            [ca.v.x, ca.v.y, ca.v.z],
            [ca.w.x, ca.w.y, ca.w.z],
            [0.0; 3],
            zb,
            zv,
            zw,
        );
        (
            crate::V3::new(f[0], f[1], f[2]),
            crate::V3::new(m[0], m[1], m[2]),
            crate::V3::new(fb[0], fb[1], fb[2]),
            crate::V3::new(mb[0], mb[1], mb[2]),
            (d, fn_),
        )
    }

    /// Φ = d (enfoncement, > 0 en pénétration) et G = ∂d/∂(δr_a, δθ_a, δr_b, δθ_b),
    /// 1 × 12 — la ligne du contact NON LISSE. Gᵀλ avec λ ≥ 0 pousse alors vers
    /// l'extérieur : −λ·∂d/∂q est la force normale.
    pub fn phi_g(&self, corps: &[Corps]) -> (DVector<f64>, DMatrix<f64>) {
        type D = Dual<f64, 12>;
        let ca = &corps[self.corps];
        let cb = self.b.map(|i| &corps[i]);
        let vv = |i: usize| D::var(0.0, i);
        let cst = |v: crate::V3| [D::from_f64(v.x), D::from_f64(v.y), D::from_f64(v.z)];
        let (zv, zw) = match cb {
            None => ([D::zero(); 3], [D::zero(); 3]),
            Some(c) => (cst(c.v), cst(c.w)),
        };
        let (_, _, _, _, d, _) = self.eval_t::<D>(
            ca,
            cb,
            [vv(0), vv(1), vv(2)],
            [vv(3), vv(4), vv(5)],
            cst(ca.v),
            cst(ca.w),
            [vv(6), vv(7), vv(8)],
            [vv(9), vv(10), vv(11)],
            zv,
            zw,
        );
        let mut g = DMatrix::zeros(1, 12);
        for k in 0..12 {
            g[(0, k)] = d.d[k];
        }
        (DVector::from_element(1, d.val()), g)
    }

    /// ∂(F_a, M_a, F_b, M_b)/∂(δr_a, δθ_a, v_a, ω_a, δr_b, δθ_b, v_b, ω_b) : 12 × 24, duaux.
    /// La TRANSLATION y est depuis le 5 sept. : sans elle la raideur de
    /// position de la pénalité (k·e·δ^{e−1}, ~3e3 N/m sur le banc) manquait au
    /// jacobien — quasi-Newton déclaré, mesuré à 7e-4 relatif à h = 1e-3, et
    /// qui croît en k·h².
    pub fn tangentes(&self, corps: &[Corps]) -> DMatrix<f64> {
        type D = Dual<f64, 24>;
        let ca = &corps[self.corps];
        let cb = self.b.map(|i| &corps[i]);
        let vv = |x: f64, i: usize| D::var(x, i);
        let (zr, zb, zv, zw) = match cb {
            None => (
                [D::zero(); 3],
                [D::zero(); 3],
                [D::zero(); 3],
                [D::zero(); 3],
            ),
            Some(c) => (
                [vv(0.0, 12), vv(0.0, 13), vv(0.0, 14)],
                [vv(0.0, 15), vv(0.0, 16), vv(0.0, 17)],
                [vv(c.v.x, 18), vv(c.v.y, 19), vv(c.v.z, 20)],
                [vv(c.w.x, 21), vv(c.w.y, 22), vv(c.w.z, 23)],
            ),
        };
        let (f, m, fb, mb, _, _) = self.eval_t::<D>(
            ca,
            cb,
            [vv(0.0, 0), vv(0.0, 1), vv(0.0, 2)],
            [vv(0.0, 3), vv(0.0, 4), vv(0.0, 5)],
            [vv(ca.v.x, 6), vv(ca.v.y, 7), vv(ca.v.z, 8)],
            [vv(ca.w.x, 9), vv(ca.w.y, 10), vv(ca.w.z, 11)],
            zr,
            zb,
            zv,
            zw,
        );
        let mut out = DMatrix::zeros(12, 24);
        for i in 0..3 {
            for k in 0..24 {
                out[(i, k)] = f[i].d[k];
                out[(3 + i, k)] = m[i].d[k];
                out[(6 + i, k)] = fb[i].d[k];
                out[(9 + i, k)] = mb[i].d[k];
            }
        }
        out
    }
}

// ── poutre géométriquement exacte ─────────────────────────────────────────────
/// Six déformations pondérées et leur dérivée première en coordonnées physiques.
///
/// U = ||z||²/2 et f = -Dᵀz. DᵀD est la contribution MATÉRIELLE : à z=0,
/// elle coïncide avec la raideur élastique. Hors repos, elle ne contient pas
/// les termes géométriques ni le transport des efforts sous rotations gauches.
/// L'API facteurs_materiels_poutres assemble ces blocs sans former leur Gram.
#[derive(Clone, Debug)]
pub(crate) struct FacteurMaterielPoutre {
    pub noeuds: [usize; 2],
    /// (sqrt(L C_eff) γ, sqrt(L C_M) κ), dans le repère matériel.
    pub deformation_ponderee: SVector<f64, 6>,
    /// Colonnes (dr_A,dtheta_A,dr_B,dtheta_B), perturbations monde à gauche.
    pub d: SMatrix<f64, 6, 12>,
}

/// log sur SO(3) : le vecteur rotation θ tel que R = exp([θ]×). θ = atan2(|w|, (tr R − 1)/2)
/// avec w = vee(R − Rᵀ)/2 = sinθ·n ; série de θ/sinθ près de zéro pour rester dérivable.
pub(crate) fn logm<T: Scalar>(r: &M<T>) -> V<T> {
    let rt = ad::m_t(r);
    let mut a = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = T::from_f64(0.5) * (r[i][j] - rt[i][j]);
        }
    }
    let w = ad::vee(&a); // sinθ·n
    let s2 = ad::v_dot(w, w);
    let c = T::from_f64(0.5) * (ad::m_trace(r) - T::one()); // cosθ
    if c.val() < -0.999999 {
        // Axe depuis la partie symétrique de R+I, robuste même à π exact.
        // Le signe suit sin(θ)n hors de la coupure. À π, la valeur a une
        // branche déterministe ; ses dérivées ne sont pas définies.
        let j = (0..3)
            .max_by(|&a, &b| r[a][a].val().total_cmp(&r[b][b].val()))
            .unwrap();
        let mut axis = [T::zero(); 3];
        axis[j] = ((r[j][j] - c) / (T::one() - c)).sqrt();
        for a in 0..3 {
            if a != j {
                axis[a] = (r[a][j] + r[j][a]) / ((T::one() - c).scale(2.0) * axis[j]);
            }
        }
        if ad::v_dot(axis, w).val() < 0.0 {
            axis = ad::v_scale(T::from_f64(-1.0), axis);
        }
        if T::ORDRE > 0 && s2.val() < 1e-24 {
            return [T::from_f64(f64::NAN); 3];
        }
        let angle = (s2 + T::from_f64(1e-32)).sqrt().atan2(c);
        return ad::v_scale(angle, axis);
    }
    let s = (s2 + T::from_f64(1e-32)).sqrt();
    let th = s.atan2(c);
    // θ/sinθ : 1 + θ²/6 + 7θ⁴/360 pour θ < 1e-3, sinon le quotient
    let fac = if th.val().abs() < 1e-3 {
        T::one() + th * th / T::from_f64(6.0) + th * th * th * th * T::from_f64(7.0 / 360.0)
    } else {
        th / s
    };
    ad::v_scale(fac, w)
}

impl crate::Poutre {
    /// Extraire le facteur matériel depuis γ,κ, sans former ni factoriser K.
    ///
    /// Les formulations milieu/intégrée partagent la cinématique et les poids
    /// de leurs efforts natifs. Les petites déformations pondérées ne sont
    /// jamais annulées par un seuil : le vecteur z reste disponible pour
    /// constater une précontrainte. Ni masse, ni réaction de liaison, ni
    /// contribution d'un autre type d'élément ne sont incluses.
    pub(crate) fn facteur_materiel(
        &self,
        corps: &[Corps],
    ) -> Result<FacteurMaterielPoutre, String> {
        if self.a == self.b || self.a >= corps.len() || self.b >= corps.len() {
            return Err(format!(
                "poutre « {} » : nœuds invalides pour le facteur matériel",
                self.nom
            ));
        }
        if !self.l.is_finite()
            || self.l <= 0.0
            || self
                .cn
                .iter()
                .chain(&self.cm)
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(format!(
                "poutre « {} » : longueur et rigidités positives finies requises",
                self.nom
            ));
        }
        let cn = self.coefficients();
        let poids: [f64; 6] = std::array::from_fn(|i| {
            self.l.sqrt() * if i < 3 { cn[i] } else { self.cm[i - 3] }.sqrt()
        });
        if poids.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(format!(
                "poutre « {} » : poids du facteur matériel hors domaine flottant",
                self.nom
            ));
        }
        type D = Dual<f64, 12>;
        let (pa, pb, bas) = self.poses_t(corps, &|i| D::var(0.0, i));
        let (gam, kap) = self.deform_t(&pa, &pb, bas);
        let angle2 = kap.iter().map(|x| (x.val() * self.l).powi(2)).sum::<f64>();
        if !angle2.is_finite() || angle2 >= (PI - 1e-12).powi(2) {
            return Err(format!(
                "poutre « {} » : branche de rotation non différentiable pour le facteur matériel",
                self.nom
            ));
        }
        let z: [D; 6] =
            std::array::from_fn(|i| (if i < 3 { gam[i] } else { kap[i - 3] }).scale(poids[i]));
        let deformation_ponderee = SVector::<f64, 6>::from_fn(|i, _| z[i].val());
        let d = SMatrix::<f64, 6, 12>::from_fn(|i, j| z[i].d[j]);
        if deformation_ponderee
            .iter()
            .chain(d.iter())
            .any(|v| !v.is_finite())
        {
            return Err(format!(
                "poutre « {} » : facteur matériel non fini",
                self.nom
            ));
        }
        Ok(FacteurMaterielPoutre {
            noeuds: [self.a, self.b],
            deformation_ponderee,
            d,
        })
    }

    /// Origine commune locale, et partie basse de la corde gardée séparée.
    /// Les perturbations restent celles des deux translations monde.
    fn poses_t<T: Scalar>(
        &self,
        corps: &[Corps],
        seed: &dyn Fn(usize) -> T,
    ) -> ((V<T>, M<T>), (V<T>, M<T>), V<T>) {
        let (a, b) = (&corps[self.a], &corps[self.b]);
        let (d, bas) = b.difference_precise(a);
        let pa = pose_pert(
            &crate::V3::zeros(),
            &a.rot,
            [seed(0), seed(1), seed(2)],
            [seed(3), seed(4), seed(5)],
        );
        let pb = pose_pert(
            &d,
            &b.rot,
            [seed(6), seed(7), seed(8)],
            [seed(9), seed(10), seed(11)],
        );
        (pa, pb, ad::v_from(bas.into()))
    }
    /// (γ, κ) génériques aux poses (r_A, R_A), (r_B, R_B)
    fn deform_t<T: Scalar>(&self, pa: &(V<T>, M<T>), pb: &(V<T>, M<T>), bas: V<T>) -> (V<T>, V<T>) {
        let ma = ad::m_mul_cst(&pa.1, &self.ra); // repère matériel porté par A
        let mb = ad::m_mul_cst(&pb.1, &self.rb); // … et par B
        let rel = ad::m_mul(&ad::m_t(&ma), &mb); // R_Aᵀ R_B (matériel)
        let th = logm(&rel);
        let half = ad::v_scale(T::from_f64(0.5), th);
        let rmoy = ad::m_mul(&ma, &ad::expm(half)); // rotation moyenne
        let d = ad::v_sub(pb.0, pa.0);
        let inv_l = T::from_f64(1.0 / self.l);
        // Intégrer R(s) = R_m exp((s/L - 1/2)[θ]×) : la corde vaut
        // L R_m J_sym(θ)(e1+γ). J_sym agit comme sinc(|θ|/2) dans le
        // plan normal à θ et comme l'identité sur son axe. L'inverser
        // évite de confondre longueur de corde et longueur d'arc.
        let chord = ad::v_scale(inv_l, ad::m_tv(&rmoy, d));
        let chord_bas = ad::v_scale(inv_l, ad::m_tv(&rmoy, bas));
        if !self.integree {
            let mut gam = chord;
            gam[0] = gam[0] - T::one();
            gam = ad::v_add(gam, chord_bas);
            return (gam, ad::v_scale(inv_l, th));
        }
        let mut gam = Self::inverse_arc(th, chord);
        gam[0] = gam[0] - T::one();
        gam = ad::v_add(gam, Self::inverse_arc(th, chord_bas));
        (gam, ad::v_scale(inv_l, th))
    }
    /// a(q), a'(q) pour J_sym^-1 = I + a(q) [theta]^2, q=theta.theta.
    /// La série évite deux annulations successives dans a et sa dérivée.
    fn arc_coefficients<T: Scalar>(q: T) -> (T, T) {
        if q.val() < 0.25 {
            const A: [f64; 8] = [
                -1. / 24.,
                -7. / 5760.,
                -31. / 967680.,
                -127. / 154828800.,
                -73. / 3503554560.,
                -1414477. / 2678117105664000.,
                -8191. / 612141052723200.,
                -16931177. / 49950709902213120000.,
            ];
            let mut value = T::from_f64(A[7]);
            let mut derivative = T::zero();
            for c in A[..7].iter().rev() {
                derivative = derivative * q + value;
                value = value * q + T::from_f64(*c);
            }
            (value, derivative)
        } else {
            let h = q.sqrt().scale(0.5);
            let f = h / h.sin();
            let a = (T::one() - f) / q;
            let b = (T::one() - h * h.cos() / h.sin()) / q;
            (a, -(f * b).scale(0.5) / q - a / q)
        }
    }
    /// Inverse de l'intégrale symétrique de rotation, auto-adjoint.
    fn inverse_arc<T: Scalar>(th: V<T>, value: V<T>) -> V<T> {
        let coeff = Self::arc_coefficients(ad::v_dot(th, th)).0;
        ad::v_add(
            value,
            ad::v_scale(coeff, ad::v_cross(th, ad::v_cross(th, value))),
        )
    }

    /// Gradient à gauche de l'énergie : moments par travail virtuel,
    /// sans dérivation automatique du scalaire d'énergie.
    fn forces_analytiques_t<T: Scalar>(
        &self,
        pa: &(V<T>, M<T>),
        pb: &(V<T>, M<T>),
        bas: V<T>,
        cn: [f64; 3],
        cm: [f64; 3],
    ) -> [T; 12] {
        let (gam, kap) = self.deform_t(pa, pb, bas);
        let th = ad::v_scale(T::from_f64(self.l), kap);
        let q = ad::v_dot(th, th);
        // Le gradient de la branche principale n'existe pas à pi.
        // Même voisinage exclu que logm sous dérivation (|sin(theta)| < 1e-12).
        if q.val() >= (PI - 1e-12).powi(2) {
            return [q.scale(f64::NAN); 12];
        }
        let ma = ad::m_mul_cst(&pa.1, &self.ra);
        let rmoy = ad::m_mul(&ma, &ad::expm(ad::v_scale(T::from_f64(0.5), th)));
        let chord = ad::v_sub(pb.0, pa.0);
        let c = ad::v_scale(T::from_f64(1. / self.l), ad::m_tv(&rmoy, chord));
        let c_bas = ad::v_scale(T::from_f64(1. / self.l), ad::m_tv(&rmoy, bas));
        let n: V<T> = std::array::from_fn(|i| gam[i].scale(cn[i]));
        let mut m: V<T> = std::array::from_fn(|i| kap[i].scale(cm[i]));
        let local_force = if self.integree {
            Self::inverse_arc(th, n)
        } else {
            n
        };
        let f = ad::m_v(&rmoy, local_force);
        if self.integree {
            let (a, ap) = Self::arc_coefficients(q);
            let tn = ad::v_dot(th, n);
            // La contribution de la corde est linéaire : conserver sa
            // partie basse jusque dans les produits du travail virtuel.
            for chord_local in [c, c_bas] {
                let tc = ad::v_dot(th, chord_local);
                let nc = ad::v_dot(n, chord_local);
                let contraction = tn * tc - q * nc;
                for i in 0..3 {
                    let gradient = a * (chord_local[i] * tn + n[i] * tc - (th[i] * nc).scale(2.))
                        + (ap * th[i] * contraction).scale(2.);
                    m[i] = m[i] + gradient.scale(self.l);
                }
            }
        }
        let polynomial = |coefficients: &[f64]| {
            coefficients
                .iter()
                .rev()
                .fold(T::zero(), |v, c| v * q + T::from_f64(*c))
        };
        let (b, c_moyen) = if q.val() < 0.25 {
            (
                polynomial(&[
                    1. / 12.,
                    1. / 720.,
                    1. / 30240.,
                    1. / 1209600.,
                    1. / 47900160.,
                    691. / 1307674368000.,
                    1. / 74724249600.,
                    3617. / 10670622842880000.,
                ]),
                polynomial(&[
                    1. / 8.,
                    1. / 384.,
                    1. / 15360.,
                    17. / 10321920.,
                    31. / 743178240.,
                    691. / 653996851200.,
                    5461. / 204047017574400.,
                    929569. / 1371195958099968000.,
                ]),
            )
        } else {
            let h = q.sqrt().scale(0.5);
            let quarter = h.scale(0.5);
            (
                (T::one() - h * h.cos() / h.sin()) / q,
                quarter.sin() / (quarter.cos() * h).scale(4.),
            )
        };
        // J_l(theta)^-T m = m + theta×m/2 + b theta×(theta×m).
        let cross = ad::v_cross(th, m);
        let transported = ad::v_add(
            m,
            ad::v_add(
                ad::v_scale(T::from_f64(0.5), cross),
                ad::v_scale(b, ad::v_cross(th, cross)),
            ),
        );
        let total_moment = ad::v_add(ad::v_cross(chord, f), ad::v_cross(bas, f));
        let th_world = ad::m_v(&ma, th);
        let difference = ad::v_sub(
            ad::m_v(&ma, transported),
            ad::v_scale(c_moyen, ad::v_cross(th_world, total_moment)),
        );
        std::array::from_fn(|i| match i / 3 {
            0 => f[i],
            1 => total_moment[i - 3].scale(0.5) + difference[i - 3],
            2 => -f[i - 6],
            _ => total_moment[i - 9].scale(0.5) - difference[i - 9],
        })
    }
    /// Flexibilité nodale exacte de Timoshenko dans la limite linéaire.
    /// La correction L²/(12 EI) représente la flexion interne condensée,
    /// et non une modification du module physique GA fourni par l'utilisateur.
    fn coefficients(&self) -> [f64; 3] {
        if !self.integree {
            return self.cn;
        }
        let mut cn = self.cn;
        for (i, c) in cn.iter_mut().enumerate().skip(1) {
            *c = 1.0 / (1.0 / self.cn[i] + self.l * self.l / (12.0 * self.cm[3 - i]));
        }
        cn
    }
    pub(crate) fn coefficients_derives(&self, mut cn: [f64; 3], cm: [f64; 3]) -> [f64; 3] {
        if !self.integree {
            return cn;
        }
        let effectif = self.coefficients();
        for i in 1..3 {
            cn[i] = cn[i] * (effectif[i] / self.cn[i]).powi(2)
                + cm[3 - i] * self.l * self.l / 12.0 * (effectif[i] / self.cm[3 - i]).powi(2);
        }
        cn
    }
    #[cfg(test)]
    fn energie_coefficients_t<T: Scalar>(
        &self,
        pa: &(V<T>, M<T>),
        pb: &(V<T>, M<T>),
        bas: V<T>,
        cn: [f64; 3],
        cm: [f64; 3],
    ) -> T {
        let (g, k) = self.deform_t(pa, pb, bas);
        let mut u = T::zero();
        for i in 0..3 {
            u = u + T::from_f64(cn[i]) * g[i] * g[i] + T::from_f64(cm[i]) * k[i] * k[i];
        }
        u * T::from_f64(0.5 * self.l)
    }
    pub fn deformations(&self, corps: &[Corps]) -> ([f64; 3], [f64; 3]) {
        let z = |_k: usize| 0.0_f64;
        let (pa, pb, bas) = self.poses_t(corps, &z);
        self.deform_t(&pa, &pb, bas)
    }
    /// forces nodales (12 : F_A, M_A, F_B, M_B, monde) = −∇U sous perturbation gauche
    pub fn forces(&self, corps: &[Corps]) -> [f64; 12] {
        self.forces_coefficients(corps, self.coefficients(), self.cm)
    }
    pub(crate) fn forces_coefficients(
        &self,
        corps: &[Corps],
        cn: [f64; 3],
        cm: [f64; 3],
    ) -> [f64; 12] {
        let (pa, pb, bas) = self.poses_t(corps, &|_| 0.);
        self.forces_analytiques_t(&pa, &pb, bas, cn, cm)
    }
    /// Contre-calcul indépendant : dérivée du scalaire d'énergie.
    #[cfg(test)]
    fn forces_energie_coefficients(
        &self,
        corps: &[Corps],
        cn: [f64; 3],
        cm: [f64; 3],
    ) -> [f64; 12] {
        type D = Dual<f64, 12>;
        let seed = |k: usize| D::var(0.0, k);
        let (pa, pb, bas) = self.poses_t(corps, &seed);
        let u = self.energie_coefficients_t(&pa, &pb, bas, cn, cm);
        let mut f = [0.0; 12];
        for k in 0..12 {
            f[k] = -u.d[k];
        }
        f
    }
    /// Contre-calcul par dérivation emboîtée de l'énergie, avec rotations
    /// successives exp(epsilon) exp(delta) R : dérivée du champ à gauche,
    /// et non Hessienne symétrique d'une seule carte exponentielle.
    #[cfg(test)]
    pub fn tangente_exacte(&self, corps: &[Corps]) -> DMatrix<f64> {
        self.tangente_exacte_coefficients(corps, self.coefficients(), self.cm)
    }
    #[cfg(test)]
    fn tangente_exacte_coefficients(
        &self,
        corps: &[Corps],
        cn: [f64; 3],
        cm: [f64; 3],
    ) -> DMatrix<f64> {
        type D1 = Dual<f64, 12>;
        type D2 = Dual<Dual<f64, 12>, 12>;
        let noeuds = [self.a, self.b];
        let (corde, bas) = corps[self.b].difference_precise(&corps[self.a]);
        let pose = |n: usize| -> (V<D2>, M<D2>) {
            let c = &corps[noeuds[n]];
            let b = 6 * n;
            let ext = |k: usize| -> D2 {
                let mut d = D2::cst(D1::cst(0.0));
                d.d[b + k] = D1::cst(1.0);
                d
            };
            let int = |k: usize| -> D2 { D2::cst(D1::var(0.0, b + k)) };
            let r0 = if n == 0 { crate::V3::zeros() } else { corde };
            let r = [
                D2::from_f64(r0.x) + ext(0) + int(0),
                D2::from_f64(r0.y) + ext(1) + int(1),
                D2::from_f64(r0.z) + ext(2) + int(2),
            ];
            let e_ext = ad::expm([ext(3), ext(4), ext(5)]);
            let e_int = ad::expm([int(3), int(4), int(5)]);
            (
                r,
                ad::m_mul(&e_int, &ad::m_mul(&e_ext, &ad::m_from::<D2>(&c.rot))),
            )
        };
        let (pa, pb) = (pose(0), pose(1));
        let u = self.energie_coefficients_t(&pa, &pb, ad::v_from(bas.into()), cn, cm);
        let mut k = DMatrix::zeros(12, 12);
        // f_i = −∂U/∂ε_i ; K[i][j] = ∂f_i/∂δ_j = −∂²U/∂ε_i∂δ_j
        for i in 0..12 {
            for j in 0..12 {
                k[(i, j)] = -u.d[j].d[i];
            }
        }
        k
    }

    /// Tangente exacte du champ de forces : un passage de duaux d'ordre un
    /// sur le gradient analytique. Le bloc rotation/rotation n'est pas
    /// symétrisé : les perturbations gauches ne commutent pas hors équilibre.
    pub fn tangente(&self, corps: &[Corps]) -> DMatrix<f64> {
        self.tangente_coefficients(corps, self.coefficients(), self.cm)
    }
    pub(crate) fn tangente_coefficients(
        &self,
        corps: &[Corps],
        cn: [f64; 3],
        cm: [f64; 3],
    ) -> DMatrix<f64> {
        type D = Dual<f64, 12>;
        let (pa, pb, bas) = self.poses_t(corps, &|i| D::var(0., i));
        let f = self.forces_analytiques_t(&pa, &pb, bas, cn, cm);
        DMatrix::from_fn(12, 12, |i, j| f[i].d[j])
    }
}

impl Elem {
    /// Φ sur des poses DONNÉES (a, b, porteur), génériques en T
    fn phi_poses<T: Scalar>(&self, p: &[PoseElement<T>], t: f64) -> Vec<T> {
        match self {
            Elem::L(l) => l.phi_t(&p[0], &p[1], t),
            Elem::D(d) => d.phi_t(&p[0], &p[1]),
            Elem::E(e) => e.phi_t(&p[0], &p[1], &p[2]),
            Elem::V(v) => v.phi_t(&p[0], &p[1]),
            Elem::C(c) => c.phi_t(&p[0], &p[1]),
            Elem::K(_) => unreachable!("contact non lisse : Φ vit dans Contact::phi_g"),
        }
    }

    /// ∂(G(q)·u)/∂q — la dérivée en POSE de la ligne en VITESSE (G·u), ne × 6k.
    ///
    /// C'est le terme que les lignes GGL et non holonomes laissaient de côté
    /// (« O(h·|u|/L) relatif », disait le commentaire) : l'audit jacobien du
    /// 5 sept. l'a mesuré à 5e-2 relatif sur un cardan sous GGL, et Newton y
    /// faisait 4,6 itérations par pas au lieu de 2. Duaux EMBOÎTÉS : à
    /// l'intérieur une direction s le long de u, à l'extérieur les 6k
    /// perturbations de pose δ — et l'ORDRE de composition compte :
    /// G(q)·u = d/ds Φ(exp(sû)·q), donc la pose est exp(sû)·exp(δ̂)·R₀, pas
    /// exp((sû + δ̂)) qui diffère par ½s[u, δ] au second ordre mixte.
    pub fn d_gu(&self, corps: &[Corps], t: f64, u: &[f64]) -> DMatrix<f64> {
        let refs = self.corps();
        if let Elem::L(l) = self {
            if l.nh {
                // non holonome : G·u est une fonction de la pose à part entière
                type D = Dual<f64, 12>;
                let seed = |k: usize| D::var(0.0, k);
                let pa = pose_de::<D>(corps, l.a, 0, &seed);
                let pb = pose_de::<D>(corps, l.b, 6, &seed);
                let uu = |c: crate::Ref, o: usize| -> [f64; 3] {
                    match c {
                        Some(i) => [u[6 * i + o], u[6 * i + o + 1], u[6 * i + o + 2]],
                        None => [0.0; 3],
                    }
                };
                let gu = l.gu_t(
                    &pa,
                    &pb,
                    t,
                    [uu(l.a, 0), uu(l.a, 3), uu(l.b, 0), uu(l.b, 3)],
                );
                let mut d = DMatrix::zeros(gu.len(), 12);
                for (i, g) in gu.iter().enumerate() {
                    for j in 0..12 {
                        d[(i, j)] = g.d[j];
                    }
                }
                return d;
            }
        }
        let k = refs.len();
        if k == 3 {
            self.d_gu_n::<18>(corps, t, u, &refs)
        } else {
            self.d_gu_n::<12>(corps, t, u, &refs)
        }
    }

    /// ∂(∂Φ/∂t)/∂q d'une liaison à cible rhéonome (ne × 12), zéro sinon — le
    /// terme de pose de la ligne GGL que `d_gu` ne porte pas.
    pub fn d_phit(&self, corps: &[Corps], t: f64) -> Option<DMatrix<f64>> {
        let Elem::L(l) = self else { return None };
        l.cible_r.as_ref()?;
        type D = Dual<f64, 12>;
        let seed = |k: usize| D::var(0.0, k);
        let pa = pose_de::<D>(corps, l.a, 0, &seed);
        let pb = pose_de::<D>(corps, l.b, 6, &seed);
        let ph = l.phit_t(&pa, &pb, t);
        let mut d = DMatrix::zeros(ph.len(), 12);
        for (i, p) in ph.iter().enumerate() {
            for j in 0..12 {
                d[(i, j)] = p.d[j];
            }
        }
        Some(d)
    }

    fn d_gu_n<const N: usize>(
        &self,
        corps: &[Corps],
        t: f64,
        u: &[f64],
        refs: &[crate::Ref],
    ) -> DMatrix<f64> {
        type In = Dual<f64, 1>;
        let s: In = In::var(0.0, 0);
        let poses: Vec<PoseElement<Dual<In, N>>> = refs
            .iter()
            .enumerate()
            .map(|(kb, c)| match c {
                Some(i) => {
                    let base = 6 * kb;
                    let cb = &corps[*i];
                    let lift = |x: f64| Dual::<In, N>::cst(In::cst(x));
                    let sd = Dual::<In, N>::cst(s);
                    let ui = |c: usize| lift(u[6 * i + c]);
                    let r = [
                        lift(cb.r.x) + Dual::<In, N>::var(In::cst(0.0), base) + sd * ui(0),
                        lift(cb.r.y) + Dual::<In, N>::var(In::cst(0.0), base + 1) + sd * ui(1),
                        lift(cb.r.z) + Dual::<In, N>::var(In::cst(0.0), base + 2) + sd * ui(2),
                    ];
                    let dth = [
                        Dual::<In, N>::var(In::cst(0.0), base + 3),
                        Dual::<In, N>::var(In::cst(0.0), base + 4),
                        Dual::<In, N>::var(In::cst(0.0), base + 5),
                    ];
                    let su = [sd * ui(3), sd * ui(4), sd * ui(5)];
                    let rot = ad::m_mul(
                        &ad::expm(su),
                        &ad::m_mul(&ad::expm(dth), &ad::m_from::<Dual<In, N>>(&cb.rot)),
                    );
                    (r, rot, ad::v_from(cb.r_bas.into()))
                }
                None => (
                    [Dual::<In, N>::cst(In::cst(0.0)); 3],
                    ad::m_ident(),
                    [Dual::<In, N>::cst(In::cst(0.0)); 3],
                ),
            })
            .collect();
        let phi = self.phi_poses(&poses, t);
        let mut d = DMatrix::zeros(phi.len(), N);
        for (i, ph) in phi.iter().enumerate() {
            for j in 0..N {
                d[(i, j)] = ph.d[j].d[0];
            }
        }
        d
    }

    /// (Φ, G, K) : Φ et G analytiques (f64, `phi_g`), K = ∂(Gᵀλ)/∂δ (6k × 6k)
    /// par duaux du PREMIER ordre sur le gradient analytique — 12 (ou 18)
    /// tangentes, au lieu des 144 des duaux emboîtés (mesuré ×2,7 trop lent).
    ///
    /// `Err` si l'élément n'est pas défini à cette pose — angle d'engrenage
    /// ou de vis sauté d'un quart de tour sur l'itéré. C'est le cas où
    /// `residu_pose` rend 1e30 pour faire reculer Newton ; le jacobien y
    /// PANIQUAIT (`expect`, jusqu'au 6 sept.), au travers de PyO3.
    pub fn eval_ad(
        &self,
        corps: &[Corps],
        t: f64,
        lam: &[f64],
    ) -> Result<(Vec<f64>, DMatrix<f64>, DMatrix<f64>), String> {
        match self {
            Elem::L(l) => {
                let (phi, g) = l.phi_g(corps, t);
                type D = Dual<f64, 12>;
                let seed = |k: usize| D::var(0.0, k);
                let pa = pose_de::<D>(corps, l.a, 0, &seed);
                let pb = pose_de::<D>(corps, l.b, 6, &seed);
                let gl = l.gtl_t::<D>(&pa, &pb, t, lam);
                let mut k = DMatrix::zeros(12, 12);
                for a in 0..12 {
                    for b in 0..12 {
                        k[(a, b)] = gl[a].d[b];
                    }
                }
                Ok((phi.iter().copied().collect(), g, k))
            }
            Elem::D(d) => {
                let (phi, g) = d.phi_g(corps);
                type D = Dual<f64, 12>;
                let seed = |k: usize| D::var(0.0, k);
                let pa = pose_de::<D>(corps, d.a, 0, &seed);
                let pb = pose_de::<D>(corps, d.b, 6, &seed);
                let gl = d.gtl_t::<D>(&pa, &pb, lam);
                let mut k = DMatrix::zeros(12, 12);
                for a in 0..12 {
                    for b in 0..12 {
                        k[(a, b)] = gl[a].d[b];
                    }
                }
                Ok((phi.iter().copied().collect(), g, k))
            }
            // ENGRENAGE et VIS : G EXACT en forme fermée (`g_t`), Gᵀλ en duaux
            // du premier ordre sur la pose — la composition exp(δ₁)·exp(δ₂)
            // est exacte par construction (G est le gradient à gauche au point
            // perturbé). Les duaux emboîtés 18×18 essayés le 5 sept. coûtaient
            // ×15 sur le pas.
            Elem::E(e) => {
                let (phi, g) = e.phi_g(corps)?;
                type D = Dual<f64, 18>;
                let seed = |k: usize| D::var(0.0, k);
                let pa = pose_de::<D>(corps, e.a, 0, &seed);
                let pb = pose_de::<D>(corps, e.b, 6, &seed);
                let pc = pose_de::<D>(corps, e.c, 12, &seed);
                let gl = e.gtl_t::<D>(&pa, &pb, &pc, lam);
                let mut k = DMatrix::zeros(18, 18);
                for a in 0..18 {
                    for b in 0..18 {
                        k[(a, b)] = gl[a].d[b];
                    }
                }
                Ok((phi.iter().copied().collect(), g, k))
            }
            Elem::V(v) => {
                let (phi, g) = v.phi_g(corps)?;
                type D = Dual<f64, 12>;
                let seed = |k: usize| D::var(0.0, k);
                let pa = pose_de::<D>(corps, v.a, 0, &seed);
                let pb = pose_de::<D>(corps, v.b, 6, &seed);
                let gl = v.gtl_t::<D>(&pa, &pb, lam);
                let mut k = DMatrix::zeros(12, 12);
                for a in 0..12 {
                    for b in 0..12 {
                        k[(a, b)] = gl[a].d[b];
                    }
                }
                Ok((phi.iter().copied().collect(), g, k))
            }
            // CARDAN : duaux emboîtés (composés), une ligne, rare
            Elem::C(_) => Ok(self.eval_k::<12>(corps, t, lam)),
            Elem::K(_) => unreachable!("contact non lisse : jugé dans jacobien_ad"),
        }
    }

    /// La version par duaux EMBOÎTÉS (Hessienne de λ·Φ) — gardée comme
    /// CONTRE-CALCUL de `eval_ad` : les deux K doivent coïncider (test).
    pub fn eval_ad2(
        &self,
        corps: &[Corps],
        t: f64,
        lam: &[f64],
    ) -> (Vec<f64>, DMatrix<f64>, DMatrix<f64>) {
        match self {
            Elem::E(_) => self.eval_k::<18>(corps, t, lam),
            _ => self.eval_k::<12>(corps, t, lam),
        }
    }

    fn eval_k<const K: usize>(
        &self,
        corps: &[Corps],
        t: f64,
        lam: &[f64],
    ) -> (Vec<f64>, DMatrix<f64>, DMatrix<f64>) {
        type D1<const K: usize> = Dual<f64, K>;
        type D2<const K: usize> = Dual<Dual<f64, K>, K>;
        // LA COMPOSITION COMPTE. K = ∂(Gᵀλ)/∂δ₂ où G est le gradient à GAUCHE
        // au point exp(δ₂)q : c'est la dérivée mixte de Φ(exp(δ₁)·exp(δ₂)·q),
        // δ₁ (inner) appliqué à gauche. Une seule exponentielle exp(δ₁ + δ₂)
        // — ce que cette fonction faisait jusqu'au 5 sept. — rend la Hessienne
        // SYMÉTRIQUE, qui diffère de ½[Gᵀλ]× (BCH : ½[δ₂, δ₁]). Mesuré par le
        // test `g_dual_egal_g_analytique` sur l'engrenage : K[4,3] 0 contre
        // −0,40 en DF, exactement ½·λ·G_z. La vis et le cardan la portaient
        // depuis leur naissance (audit : 8e-7 à h = 1e-3, en k·h² au-delà).
        let refs = self.corps();
        let poses: Vec<PoseElement<D2<K>>> = refs
            .iter()
            .enumerate()
            .map(|(kb, c)| match c {
                Some(i) => {
                    let base = 6 * kb;
                    let cb = &corps[*i];
                    let lift = |x: f64| D2::<K>::cst(D1::<K>::cst(x));
                    let d1 = |k: usize| D2::<K>::cst(D1::<K>::var(0.0, k));
                    let d2 = |k: usize| D2::<K>::var(D1::<K>::cst(0.0), k);
                    let r = [
                        lift(cb.r.x) + d1(base) + d2(base),
                        lift(cb.r.y) + d1(base + 1) + d2(base + 1),
                        lift(cb.r.z) + d1(base + 2) + d2(base + 2),
                    ];
                    let e1 = ad::expm([d1(base + 3), d1(base + 4), d1(base + 5)]);
                    let e2 = ad::expm([d2(base + 3), d2(base + 4), d2(base + 5)]);
                    let rot = ad::m_mul(&e1, &ad::m_mul(&e2, &ad::m_from::<D2<K>>(&cb.rot)));
                    (r, rot, ad::v_from(cb.r_bas.into()))
                }
                None => (
                    [D2::<K>::cst(D1::<K>::cst(0.0)); 3],
                    ad::m_ident(),
                    [D2::<K>::cst(D1::<K>::cst(0.0)); 3],
                ),
            })
            .collect();
        let phi = self.phi_poses(&poses, t);
        let n = phi.len();
        let mut val = Vec::with_capacity(n);
        let mut g = DMatrix::zeros(n, K);
        let mut kmat = DMatrix::zeros(K, K);
        for (i, p) in phi.iter().enumerate() {
            val.push(p.v.v);
            for a in 0..K {
                g[(i, a)] = p.v.d[a];
                for b in 0..K {
                    // K[a][b] = ∂(Gᵀλ)_a/∂δ₂_b = λ · ∂²Φ/∂δ₂_b ∂δ₁_a
                    kmat[(a, b)] += lam[i] * p.d[b].d[a];
                }
            }
        }
        (val, g, kmat)
    }
}

// ── corps : résidu inertiel en rotation et ses tangentes ────────────────────
/// R J Rᵀ ẇ + ω × (R J Rᵀ ω) avec R = exp(δθ)·R₀ ; tangentes en (δθ, ω).
pub(crate) fn residu_rot<T: Scalar>(c: &Corps, dth: V<T>, w: V<T>, wdot: [f64; 3]) -> V<T> {
    let (_, rot) = pose_pert(&c.r, &c.rot, [T::zero(); 3], dth);
    let js = ad::m_mul(&ad::m_mul(&rot, &ad::m_from::<T>(&c.j)), &ad::m_t(&rot));
    let a = ad::m_v(&js, ad::v_from::<T>(wdot));
    let gyr = ad::v_cross(w, ad::m_v(&js, w));
    ad::v_add(a, gyr)
}

/// jacobien gauche de exp sur SO(3) : exp(θ + ε) = exp(J_l(θ)ε) exp(θ)
pub(crate) fn jac_gauche(th: &crate::V3) -> crate::M3 {
    let t = th.norm();
    let k = crate::skew(th);
    if t < 1e-6 {
        return crate::M3::identity() + 0.5 * k + k * k / 6.0;
    }
    crate::M3::identity() + ((1.0 - t.cos()) / (t * t)) * k + ((t - t.sin()) / (t * t * t)) * k * k
}

impl Modele {
    /// Le jacobien exact du résidu de `Pas::residu` au point x, le modèle
    /// étant déjà posé à l'état (q₁, u₁) de x (appeler `pas.etat` avant).
    /// Rendu en TRIPLETS ; les blocs d'éléments se calculent EN PARALLÈLE
    /// (rayon), la dispersion est séquentielle et ordonnée (déterminisme).
    pub(crate) fn jacobien_ad(
        &self,
        pas: &Pas,
        x: &DVector<f64>,
        dth: &[crate::V3],
    ) -> Result<(crate::Triplets, f64), String> {
        use rayon::prelude::*;
        let (n, m) = (self.n(), self.m());
        // u₁ courant, lu UNE fois : les lignes en vitesse (nh, GGL) le
        // relisaient par élément — une allocation de n par élément
        let u1 = self.u();
        let c = (1.0 - pas.af) / (1.0 - pas.am);
        let kq = pas.bet * pas.hh * pas.hh * c; // ∂q/∂u̇ avant le facteur de pose
        let ku = pas.gam * pas.hh * c; // ∂u/∂u̇

        // CHRONO DE LA TANGENTE DE exp — c'est le seul chiffre qui tranche le
        // refus de Cayley (note d'architecture de sept. 2026, §1.3 : « le coût
        // par itération de Newton est réduit par rapport à la formulation
        // exponentielle »). Il suppose que les dérivées de exp sont DANS la
        // boucle ; ici `jac_gauche` est appelée une fois par CORPS et par
        // jacobien, pas par élément ni par itération de Newton.
        // Ce qu'il compte : la construction des J_l. Ce qu'il NE compte PAS,
        // et c'est déclaré : les produits `J_l · kq` dispersés dans les blocs
        // AD (une poignée de 3×3 par élément, inséparables du reste). C'est
        // donc la BORNE BASSE de ce que Cayley remplacerait, pour sigma=0.
        // Sigma!=0 inclut la résolution locale et les tangentes implicites.
        let te0 = std::time::Instant::now();
        // Sigma=0 conserve le chemin historique. Pour sigma!=0, le facteur
        // de pose est J_l A^-1 B ; GGL doit aussi conserver B^-1 pour sa
        // variation indépendante de la vitesse (voir schema_lie).
        let (jls, transfert_ggl): (Vec<crate::M3>, Vec<crate::M3>) = if pas.sigma == 0.0 {
            (dth.iter().map(jac_gauche).collect(), Vec::new())
        } else {
            let (_, dq) = pas.a1_dq(x)?;
            let s = pas.sigma * pas.hh * pas.bet / pas.gam;
            let cartes = (0..self.corps.len())
                .map(|i| {
                    let k = 6 * i + 3;
                    let th = pas.hh * crate::V3::new(dq[k], dq[k + 1], dq[k + 2]);
                    let w = crate::V3::new(u1[k], u1[k + 1], u1[k + 2]);
                    crate::schema_lie::tangentes(th, w, s, pas.sigma)
                })
                .collect::<Result<Vec<_>, String>>()?;
            cartes.into_iter().unzip()
        };
        let t_exp = te0.elapsed().as_secs_f64();
        let mut trip: crate::Triplets =
            Vec::with_capacity(36 * self.corps.len() + 400 * self.elems.len());
        let t = |r: usize, cc: usize, v: f64| faer::sparse::Triplet {
            row: r,
            col: cc,
            val: v,
        };
        // ∂r/∂(δq) PAR UNITÉ D'INCRÉMENT DE POSE (δθ = incrément du vecteur
        // rotation, facteur de pose déjà appliqué), gardé à part : c'est la matrice dont
        // les colonnes ζ de GGL sont une COMBINAISON LINÉAIRE — ζ n'agit sur le
        // résidu que par la pose, qu'il déplace de βh²·G₀ᵀζ. Colonnes EXACTES,
        // là où les différences finies (jusqu'au 5 sept.) n'avaient que
        // ~4 chiffres à h = 1e-3 et cassaient au point mort de `srscm`
        // (cond(GGᵀ) ×4 800 : 1e-8 de bruit relatif × 2,8e6 = un incrément faux
        // au pour-cent, Newton stagne). Ne garde rien quand GGL est absent.
        let mut kqt: Vec<(usize, usize, f64)> = Vec::new();
        let garde_kq = pas.ggl.is_some();
        // ── corps : M, K_M, C (gyroscopique), avec la tangente de exp
        // rayon ne paie qu'au-delà d'une taille : sur la bielle (4 éléments) le
        // dispatch coûtait ×4 (0,044 → 0,19 ms/pas). Seuil mesuré, pas deviné.
        // ⚠ LES POUTRES COMPTENT AUSSI. Le critère ne regardait que `elems` :
        // sur un modèle de POUTRES pures il n'y a qu'un élément (l'encastrement),
        // donc les tangentes des N poutres se calculaient en SÉRIE. Trouvé le
        // 3 sept. en balayant les topologies, à côté de la copie O(N²).
        let seuil_par = std::env::var("VINKULUM_PAR_MIN")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(PARALLELE_MIN);
        let par = self.elems.len() + self.poutres.len() >= seuil_par;
        if std::env::var("VINKULUM_CHRONO").is_ok() {
            static UNE: std::sync::Once = std::sync::Once::new();
            UNE.call_once(|| {
                eprintln!(
                    "éléments {} + poutres {} · corps {} · parallèle {}",
                    self.elems.len(),
                    self.poutres.len(),
                    self.corps.len(),
                    par
                )
            });
        }
        let gel = |i: usize| self.gel.get(i).copied().unwrap_or(false);
        let bloc_corps = |(i, cb): (usize, &Corps)| -> [[(f64, f64); 3]; 3] {
            if gel(i) {
                return [[(0.0, 0.0); 3]; 3]; // multi-rythme : ligne identité, posée plus bas
            }
            let k = 6 * i;
            let wdot = [x[k + 3], x[k + 4], x[k + 5]];
            type D6 = Dual<f64, 6>;
            let d = residu_rot::<D6>(
                cb,
                [D6::var(0.0, 0), D6::var(0.0, 1), D6::var(0.0, 2)],
                [D6::var(cb.w.x, 3), D6::var(cb.w.y, 4), D6::var(cb.w.z, 5)],
                wdot,
            );
            let jl = &jls[i];
            let js = cb.rot * cb.j * cb.rot.transpose();
            let mut b = [[(0.0, 0.0); 3]; 3];
            for a in 0..3 {
                for bb in 0..3 {
                    let mut dth_ab = 0.0;
                    for cc in 0..3 {
                        dth_ab += d[a].d[cc] * jl[(cc, bb)];
                    }
                    b[a][bb] = (js[(a, bb)] + dth_ab * kq + d[a].d[3 + bb] * ku, dth_ab);
                }
            }
            b
        };
        let blocs_corps: Vec<[[(f64, f64); 3]; 3]> = if par {
            self.corps.par_iter().enumerate().map(bloc_corps).collect()
        } else {
            self.corps.iter().enumerate().map(bloc_corps).collect()
        };
        for (i, cb) in self.corps.iter().enumerate() {
            let k = 6 * i;
            for a in 0..3 {
                trip.push(t(k + a, k + a, cb.m));
            }
            for a in 0..3 {
                for bb in 0..3 {
                    trip.push(t(k + 3 + a, k + 3 + bb, blocs_corps[i][a][bb].0));
                    if garde_kq && blocs_corps[i][a][bb].1 != 0.0 {
                        kqt.push((k + 3 + a, k + 3 + bb, blocs_corps[i][a][bb].1));
                    }
                }
            }
        }
        // ── éléments : (G, K_c) en parallèle, puis dispersion
        let mut rows = Vec::with_capacity(self.elems.len());
        let mut row = n;
        for e in &self.elems {
            rows.push(row);
            row += e.n();
        }
        let sc = pas.bet * pas.hh * pas.hh;
        let bloc_elem =
            |(ie, e): (usize, &Elem)| -> Result<(DMatrix<f64>, DMatrix<f64>, (f64, f64)), String> {
                let ne = e.n();
                let r0 = rows[ie];
                if self.gel_elem.get(ie).copied().unwrap_or(false) {
                    // multi-rythme : élément entièrement gelé, ses lignes sont masquées
                    return Ok((DMatrix::zeros(ne, 12), DMatrix::zeros(12, 12), (1.0, 0.0)));
                }
                let lam: Vec<f64> = (0..ne)
                    .map(|i| {
                        if self.actif[r0 - n + i] {
                            x[r0 + i]
                        } else {
                            0.0
                        }
                    })
                    .collect();
                if let Elem::K(k) = e {
                    // Fischer–Burmeister : ∂φ/∂a sur la ligne (a = −d/sc, donc facteur
                    // −∂φ/∂a devant G·c), ∂φ/∂b·s/sc sur la diagonale λ. La raideur
                    // de contact λ·∂G/∂q est laissée de côté (ponytail: quasi-Newton,
                    // comme la pénalité qui ignore déjà sa raideur en position).
                    let (_, g) = self.contacts[k.ct].phi_g(&self.corps);
                    let (actif_k, rhs) = self.k_pas.get(k.ct).copied().unwrap_or((false, 0.0));
                    if !actif_k {
                        return Ok((g, DMatrix::zeros(12, 12), (0.0, 1.0)));
                    }
                    let mut gu = 0.0;
                    for (kb, cb) in e.corps().iter().enumerate() {
                        if let Some(ib) = cb {
                            for cc in 0..6 {
                                gu -= g[(0, 6 * kb + cc)] * self.corps[*ib].u6(cc);
                            }
                        }
                    }
                    let scv = pas.gam * pas.hh;
                    let (a, b) = ((gu + rhs) / scv, lam[0] * self.s_fb / sc);
                    let nrm = (a * a + b * b).sqrt();
                    let (fa, fb) = if nrm > 1e-300 {
                        (1.0 - a / nrm, 1.0 - b / nrm)
                    } else {
                        (1.0 - 0.5f64.sqrt(), 1.0 - 0.5f64.sqrt())
                    };
                    return Ok((g, DMatrix::zeros(12, 12), (-fa, fb * self.s_fb / sc)));
                }
                let (_, g, kc) = e.eval_ad(&self.corps, pas.t1, &lam)?;
                Ok((g, kc, (1.0, 0.0)))
            };
        let blocs: Vec<(DMatrix<f64>, DMatrix<f64>, (f64, f64))> = if par {
            self.elems
                .par_iter()
                .enumerate()
                .map(bloc_elem)
                .collect::<Result<Vec<_>, String>>()?
        } else {
            self.elems
                .iter()
                .enumerate()
                .map(bloc_elem)
                .collect::<Result<Vec<_>, String>>()?
        };
        for (ie, e) in self.elems.iter().enumerate() {
            let ne = e.n();
            let row = rows[ie];
            let (g, kc, (fa, fb)) = &blocs[ie];
            for i in 0..ne {
                if !self.actif[row - n + i] {
                    trip.push(t(row + i, row + i, 1.0));
                } else if *fb != 0.0 {
                    trip.push(t(row + i, row + i, *fb));
                }
            }
            let cs = e.corps();
            // la ligne dépend-elle de la POSE (Φ/sc) ou de la VITESSE (G·u) ?
            // Les lignes en vitesse (non holonomes, GGL) portent leur terme
            // ∂(G·u)/∂q EXACT depuis le 5 sept. (`Elem::d_gu`, ci-dessous) ;
            // seules les lignes de contact (K) le laissent de côté.
            let en_pose = !(matches!(e, Elem::L(l) if l.nh) || matches!(e, Elem::K(_)));
            let nh_l = matches!(e, Elem::L(l) if l.nh);
            if nh_l {
                // ∂(G·u₁)/∂q₁ · ∂q₁/∂u̇ / (γh), sur les u̇ des corps de l'élément
                let dg = e.d_gu(&self.corps, pas.t1, u1.as_slice());
                let scv = pas.gam * pas.hh;
                for (ka, ca) in cs.iter().enumerate() {
                    let Some(ia) = ca else { continue };
                    let jl_a = &jls[*ia];
                    for i in 0..ne {
                        if !self.actif[row - n + i] {
                            continue;
                        }
                        for b in 0..3 {
                            let vt = dg[(i, 6 * ka + b)] / scv;
                            if vt != 0.0 {
                                trip.push(t(row + i, 6 * ia + b, vt * kq));
                                if garde_kq {
                                    kqt.push((row + i, 6 * ia + b, vt));
                                }
                            }
                            let mut vr = 0.0;
                            for cc in 0..3 {
                                vr += dg[(i, 6 * ka + 3 + cc)] * jl_a[(cc, b)];
                            }
                            let vr = vr / scv;
                            if vr != 0.0 {
                                trip.push(t(row + i, 6 * ia + 3 + b, vr * kq));
                                if garde_kq {
                                    kqt.push((row + i, 6 * ia + 3 + b, vr));
                                }
                            }
                        }
                    }
                }
            }
            for (ka, ca) in cs.iter().enumerate() {
                let Some(ia) = ca else { continue };
                let jl_a = &jls[*ia];
                for i in 0..ne {
                    if !self.actif[row - n + i] {
                        continue;
                    }
                    for b in 0..3 {
                        let gt = g[(i, 6 * ka + b)];
                        trip.push(t(6 * ia + b, row + i, gt)); // Gᵀ
                        trip.push(t(row + i, 6 * ia + b, gt * c * fa)); // G·∂q/∂u̇ /(βh²)
                        let mut gr = 0.0;
                        for cc in 0..3 {
                            gr += g[(i, 6 * ka + 3 + cc)] * jl_a[(cc, b)];
                        }
                        trip.push(t(6 * ia + 3 + b, row + i, g[(i, 6 * ka + 3 + b)]));
                        // ligne en POSE : ∂Φ/∂u̇ passe par q₁, donc par J_l ;
                        // ligne en VITESSE (nh, contact) : ∂(G·u₁)/∂u̇ passe par
                        // u₁ = … + γh·a₁, SANS J_l — appliqué à tort jusqu'au
                        // 5 sept. (audit : 4e-4 sur la ligne nh, colonnes θ)
                        let gv = if en_pose { gr } else { g[(i, 6 * ka + 3 + b)] };
                        trip.push(t(row + i, 6 * ia + 3 + b, gv * c * fa));
                        if garde_kq && en_pose {
                            kqt.push((row + i, 6 * ia + b, gt / sc));
                            kqt.push((row + i, 6 * ia + 3 + b, gr / sc));
                        }
                    }
                }
                for (kb, cbb) in cs.iter().enumerate() {
                    let Some(ib) = cbb else { continue };
                    let jl_b = &jls[*ib];
                    for p in 0..6 {
                        for q in 0..6 {
                            let mut v = 0.0;
                            if q < 3 {
                                v = kc[(6 * ka + p, 6 * kb + q)];
                            } else {
                                for cc in 0..3 {
                                    v += kc[(6 * ka + p, 6 * kb + 3 + cc)] * jl_b[(cc, q - 3)];
                                }
                            }
                            if v != 0.0 {
                                trip.push(t(6 * ia + p, 6 * ib + q, v * kq));
                                if garde_kq {
                                    kqt.push((6 * ia + p, 6 * ib + q, v));
                                }
                            }
                        }
                    }
                }
            }
        }
        // ── couples à loi : r_b −= τ, r_a += τ  →  ∂r_b = −∂τ, ∂r_a = +∂τ
        for charge in &self.efforts_temporels {
            let i = charge.corps;
            if self.gele(i) {
                continue;
            }
            let d = -charge.tangente_moment(&self.corps, pas.t1) * jls[i];
            for p in 0..3 {
                for q in 0..3 {
                    trip.push(t(6 * i + 3 + p, 6 * i + 3 + q, d[(p, q)] * kq));
                    if garde_kq {
                        kqt.push((6 * i + 3 + p, 6 * i + 3 + q, d[(p, q)]));
                    }
                }
            }
        }
        for cp in &self.couples {
            if [cp.a, cp.b].iter().flatten().all(|&i| gel(i)) {
                continue;
            }
            let d = cp.tangentes(&self.corps, pas.t1);
            let cibles: [(crate::Ref, f64); 2] = [(cp.a, 1.0), (cp.b, -1.0)];
            for (ic, sg) in cibles {
                let Some(i) = ic else { continue };
                for p in 0..3 {
                    // colonnes θ de a (0..3) et b (3..6) via J_l et kq ; ω de a (6..9), b (9..12) via ku
                    for (src, base, coef) in [(cp.a, 0usize, kq), (cp.b, 3, kq)] {
                        let Some(j) = src else { continue };
                        let jl = &jls[j];
                        for q in 0..3 {
                            let mut v = 0.0;
                            for cc in 0..3 {
                                v += d[(p, base + cc)] * jl[(cc, q)];
                            }
                            if v != 0.0 {
                                trip.push(t(6 * i + 3 + p, 6 * j + 3 + q, sg * v * coef));
                                if garde_kq {
                                    kqt.push((6 * i + 3 + p, 6 * j + 3 + q, sg * v));
                                }
                            }
                        }
                    }
                    for (src, base) in [(cp.a, 6usize), (cp.b, 9)] {
                        let Some(j) = src else { continue };
                        for q in 0..3 {
                            let v = d[(p, base + q)];
                            if v != 0.0 {
                                trip.push(t(6 * i + 3 + p, 6 * j + 3 + q, sg * v * ku));
                            }
                        }
                    }
                }
            }
        }
        // ── pales : dérivées en position, orientation et vitesse
        for p in &self.pales {
            let d = p.tangentes(&self.corps, &self.inflows, self.vent);
            let i = p.corps;
            let jl = &jls[i];
            for row_l in 0..6 {
                for q in 0..3 {
                    let vr = d[(row_l, q)];
                    if vr != 0.0 {
                        trip.push(t(6 * i + row_l, 6 * i + q, -vr * kq));
                        if garde_kq {
                            kqt.push((6 * i + row_l, 6 * i + q, -vr));
                        }
                    }
                    let mut v = 0.0;
                    for cc in 0..3 {
                        v += d[(row_l, 3 + cc)] * jl[(cc, q)];
                    }
                    if v != 0.0 {
                        trip.push(t(6 * i + row_l, 6 * i + 3 + q, -v * kq));
                        if garde_kq {
                            kqt.push((6 * i + row_l, 6 * i + 3 + q, -v));
                        }
                    }
                    let vv = d[(row_l, 6 + q)];
                    if vv != 0.0 {
                        trip.push(t(6 * i + row_l, 6 * i + q, -vv * ku));
                    }
                    let vw = d[(row_l, 9 + q)];
                    if vw != 0.0 {
                        trip.push(t(6 * i + row_l, 6 * i + 3 + q, -vw * ku));
                    }
                }
            }
        }
        // ── contacts : deux corps possibles, donc blocs croisés
        for ct in &self.contacts {
            if [Some(ct.corps), ct.b].iter().flatten().all(|&i| gel(i)) {
                continue;
            }
            let ct_nl;
            let ct = if ct.nonlisse {
                ct_nl = crate::Contact {
                    fn_impose: Some(x[n + ct.row].max(0.0)),
                    ..ct.clone()
                };
                &ct_nl
            } else {
                ct
            };
            let d = ct.tangentes(&self.corps);
            let noeuds = [Some(ct.corps), ct.b];
            for (ia, na) in noeuds.iter().enumerate() {
                let Some(ni) = *na else { continue };
                for (ib, nb) in noeuds.iter().enumerate() {
                    let Some(nj) = *nb else { continue };
                    let jl = &jls[nj];
                    for row_l in 0..6 {
                        let rl = 6 * ia + row_l;
                        for q in 0..3 {
                            // translation : raideur de position de la pénalité
                            let vr = d[(rl, 12 * ib + q)];
                            if vr != 0.0 {
                                trip.push(t(6 * ni + row_l, 6 * nj + q, -vr * kq));
                                if garde_kq {
                                    kqt.push((6 * ni + row_l, 6 * nj + q, -vr));
                                }
                            }
                            let mut v = 0.0;
                            for cc in 0..3 {
                                v += d[(rl, 12 * ib + 3 + cc)] * jl[(cc, q)];
                            }
                            if v != 0.0 {
                                trip.push(t(6 * ni + row_l, 6 * nj + 3 + q, -v * kq));
                                if garde_kq {
                                    kqt.push((6 * ni + row_l, 6 * nj + 3 + q, -v));
                                }
                            }
                            let vv = d[(rl, 12 * ib + 6 + q)];
                            if vv != 0.0 {
                                trip.push(t(6 * ni + row_l, 6 * nj + q, -vv * ku));
                            }
                            let vw = d[(rl, 12 * ib + 9 + q)];
                            if vw != 0.0 {
                                trip.push(t(6 * ni + row_l, 6 * nj + 3 + q, -vw * ku));
                            }
                        }
                    }
                }
            }
            // NON LISSE : le frottement est LINÉAIRE en λ (−μ·λ·tanh·v̂_t), donc
            // sa dérivée en λ est la force par unité de λ — absente jusqu'au
            // 5 sept. (mesurée à 0,22 relatif par l'audit jacobien)
            if ct.nonlisse && x[n + ct.row] > 0.0 {
                let ct1 = crate::Contact {
                    fn_impose: Some(1.0),
                    ..ct.clone()
                };
                let (f1, m1, fb1, mb1, _) = ct1.valeur(&self.corps);
                for k in 0..3 {
                    trip.push(t(6 * ct.corps + k, n + ct.row, -f1[k]));
                    trip.push(t(6 * ct.corps + 3 + k, n + ct.row, -m1[k]));
                    if let Some(ib) = ct.b {
                        trip.push(t(6 * ib + k, n + ct.row, -fb1[k]));
                        trip.push(t(6 * ib + 3 + k, n + ct.row, -mb1[k]));
                    }
                }
            }
        }
        // ── poutres : r −= f  →  ∂r = −K·(J_l·kq sur les colonnes θ, kq sur les colonnes r)
        let tang_p = |b: &crate::Poutre| -> DMatrix<f64> {
            if gel(b.a) && gel(b.b) {
                DMatrix::zeros(12, 12)
            } else {
                b.tangente(&self.corps)
            }
        };
        let blocs_p: Vec<DMatrix<f64>> = if par {
            self.poutres.par_iter().map(tang_p).collect()
        } else {
            self.poutres.iter().map(tang_p).collect()
        };
        for (b, k) in self.poutres.iter().zip(blocs_p.iter()) {
            let noeuds = [b.a, b.b];
            for (ia, &ni) in noeuds.iter().enumerate() {
                for (ib, &nj) in noeuds.iter().enumerate() {
                    let jl = &jls[nj];
                    for p in 0..6 {
                        for q in 0..6 {
                            let mut v = 0.0;
                            if q < 3 {
                                v = k[(6 * ia + p, 6 * ib + q)];
                            } else {
                                for cc in 0..3 {
                                    v += k[(6 * ia + p, 6 * ib + 3 + cc)] * jl[(cc, q - 3)];
                                }
                            }
                            if v != 0.0 {
                                trip.push(t(6 * ni + p, 6 * nj + q, -v * kq));
                                if garde_kq {
                                    kqt.push((6 * ni + p, 6 * nj + q, -v));
                                }
                            }
                        }
                    }
                }
            }
        }
        for se in &self.supers {
            let (kf, cf) = se.tangentes(&self.corps)?;
            for (ia, &ni) in se.noeuds.iter().enumerate() {
                for (ib, &nj) in se.noeuds.iter().enumerate() {
                    for p in 0..6 {
                        for q in 0..6 {
                            let row = 6 * ia + p;
                            let v = if q < 3 {
                                kf[(row, 6 * ib + q)]
                            } else {
                                (0..3)
                                    .map(|a| kf[(row, 6 * ib + 3 + a)] * jls[nj][(a, q - 3)])
                                    .sum()
                            };
                            trip.push(t(
                                6 * ni + p,
                                6 * nj + q,
                                -v * kq - cf[(row, 6 * ib + q)] * ku,
                            ));
                            if garde_kq {
                                kqt.push((6 * ni + p, 6 * nj + q, -v));
                            }
                        }
                    }
                }
            }
        }
        if let Some((g0, _)) = &pas.ggl {
            // lignes GGL G(q₁)·u₁ + Φ_t, mises à l'échelle 1/(γh) : ∂/∂u̇ = G·ku/(γh) = G·c,
            // sans J_l (u₁ est une vitesse, pas une pose).
            // ponytail: le terme (∂G/∂q·u₁)·kq est laissé de côté — O(h·|u|/L)
            // relatif, Newton reste convergent (linéaire à ce taux).
            for (ie, e) in self.elems.iter().enumerate() {
                let ne = e.n();
                let r0 = rows[ie] - n;
                let nh = matches!(e, Elem::L(l) if l.nh) || matches!(e, Elem::K(_));
                let g = &blocs[ie].0;
                for i in 0..ne {
                    if !self.actif[r0 + i] || nh {
                        // ligne sans objet : r = ζ, diagonale 1
                        trip.push(t(n + m + r0 + i, n + m + r0 + i, 1.0));
                        continue;
                    }
                    for (ka, ca) in e.corps().iter().enumerate() {
                        let Some(ia) = ca else { continue };
                        for b in 0..6 {
                            let v = g[(i, 6 * ka + b)] * c;
                            if v != 0.0 {
                                trip.push(t(n + m + r0 + i, 6 * ia + b, v));
                            }
                        }
                    }
                }
                // ∂(G(q₁)·u₁)/∂q₁ sur la ligne GGL — exact, par `d_gu` ; les
                // colonnes ζ suivent par `kqt` (même dépendance en pose)
                if ne > 0 && (0..ne).any(|i| self.actif[r0 + i]) {
                    let mut dg = e.d_gu(&self.corps, pas.t1, u1.as_slice());
                    if let Some(dp) = e.d_phit(&self.corps, pas.t1) {
                        dg += dp; // cible rhéonome : Φ_t dépend de la pose aussi
                    }
                    let scv = pas.gam * pas.hh;
                    for (ka, ca) in e.corps().iter().enumerate() {
                        let Some(ia) = ca else { continue };
                        let jl_a = &jls[*ia];
                        for i in 0..ne {
                            if !self.actif[r0 + i] {
                                continue;
                            }
                            for b in 0..3 {
                                let vt = dg[(i, 6 * ka + b)] / scv;
                                if vt != 0.0 {
                                    trip.push(t(n + m + r0 + i, 6 * ia + b, vt * kq));
                                    kqt.push((n + m + r0 + i, 6 * ia + b, vt));
                                }
                                let mut vr = 0.0;
                                for cc in 0..3 {
                                    vr += dg[(i, 6 * ka + 3 + cc)] * jl_a[(cc, b)];
                                }
                                let vr = vr / scv;
                                if vr != 0.0 {
                                    trip.push(t(n + m + r0 + i, 6 * ia + 3 + b, vr * kq));
                                    kqt.push((n + m + r0 + i, 6 * ia + 3 + b, vr));
                                }
                            }
                        }
                    }
                }
            }
            // COLONNES ζ, EXACTES : ∂r/∂ζ_l = (∂r/∂δq) · βh² · G₀ᵀ[:, l]. Un
            // produit creux-creux, indexé par colonne de G₀ — aucune copie du
            // modèle, aucune différence finie (elles coûtaient O(m·(n + nnz))
            // par jacobien ET n'avaient pas la précision du point mort).
            let mut g0_par_col: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
            for &(l, j, v) in g0 {
                if v != 0.0 {
                    if pas.sigma == 0.0 || j % 6 < 3 {
                        g0_par_col[j].push((l, v));
                    } else {
                        // kqt porte J_l A^-1 B. La variation GGL ne modifie
                        // pas omega1 : appliquer B^-1 à G0^T avant ce produit.
                        let corps = j / 6;
                        for c in 0..3 {
                            let gv = transfert_ggl[corps][(c, j % 6 - 3)] * v;
                            if gv != 0.0 {
                                g0_par_col[6 * corps + 3 + c].push((l, gv));
                            }
                        }
                    }
                }
            }
            let bh2 = pas.bet * pas.hh * pas.hh;
            for &(i, j, v) in &kqt {
                for &(l, gv) in &g0_par_col[j] {
                    trip.push(t(i, n + m + l, v * bh2 * gv));
                }
            }
        }
        let _ = m;
        // multi-rythme : les lignes des corps GELÉS ne sont pas des inconnues —
        // identité, et rien d'autre (leur résidu rend u̇, cf. `residu_pose`)
        if self.gel.iter().any(|g| *g) {
            let gele = |r: usize| r < n && self.gel[r / 6];
            trip.retain(|tr| !gele(tr.row));
            for i in 0..self.corps.len() {
                if self.gel[i] {
                    for c in 0..6 {
                        trip.push(t(6 * i + c, 6 * i + c, 1.0));
                    }
                }
            }
        }
        if !trip.iter().all(|t| t.val.is_finite()) {
            return Err(
                "jacobien non fini (débordement ou rotation sur la coupure du logarithme)".into(),
            );
        }
        Ok((trip, t_exp))
    }
}

#[cfg(test)]
mod tests {
    use super::{SMatrix, SVector, PI};
    use crate::*;

    #[test]
    fn courbure_distance_echelles_et_reperes() {
        for q in [M3::identity(), expm(&V3::new(0.31, -0.47, 0.22))] {
            for (longueur, vitesse) in [(1., 1.), (1e-6, 10.), (1., 1e8), (1e6, 1e8)] {
                let mut mo = Modele::new(V3::zeros());
                mo.corps.push(Corps {
                    nom: "masse".into(),
                    m: 1.,
                    j: M3::identity(),
                    r: q * V3::new(longueur, 0., 0.),
                    r_bas: V3::zeros(),
                    rot: q,
                    v: q * V3::new(0., vitesse, 0.),
                    w: V3::zeros(),
                });
                mo.elems.push(Elem::D(Distance {
                    nom: "fil".into(),
                    a: None,
                    b: Some(0),
                    pa: V3::zeros(),
                    pb: V3::zeros(),
                    l: longueur,
                }));
                let reference = vitesse * vitesse / longueur;
                assert!((mo.biais_acceleration(0.).unwrap()[0] / reference - 1.).abs() < 3e-15);
                let before = mo.poses();
                let (a, _) = mo.acc_init(0.).unwrap();
                assert!(
                    (q.transpose() * a.rows(0, 3) + V3::new(reference, 0., 0.)).norm()
                        < 4e-15 * reference
                );
                assert_eq!(before, mo.poses());
            }
        }
    }

    #[test]
    fn courbure_cible_rotation_holonome_et_vitesse() {
        let (theta, angle_cible, omega, taux) = (0.3_f64, 0.1_f64, 2.0_f64, 0.7_f64);
        for nh in [false, true] {
            let mut mo = Modele::new(V3::zeros());
            mo.corps.push(Corps {
                nom: "arbre".into(),
                m: 1.,
                j: M3::identity(),
                r: V3::zeros(),
                r_bas: V3::zeros(),
                rot: expm(&(theta * V3::z())),
                v: V3::zeros(),
                w: omega * V3::z(),
            });
            mo.elems.push(Elem::L(Liaison {
                nom: "commande".into(),
                a: None,
                b: Some(0),
                pa: V3::zeros(),
                pb: V3::zeros(),
                ra: M3::identity(),
                rb: M3::identity(),
                bt: vec![],
                br: vec![2],
                cible_t: None,
                cible_r: Some((
                    V3::z(),
                    Loi::Lineaire {
                        a0: angle_cible,
                        taux,
                    },
                )),
                nh,
            }));
            let expected = -(omega - taux).powi(2) * (theta - angle_cible).sin();
            assert!((mo.biais_acceleration(0.).unwrap()[0] - expected).abs() < 1e-15);
        }
    }

    #[test]
    fn courbure_contact_excentre_et_acceleration_admissible() {
        let mut mo = Modele::new(V3::zeros());
        let p = V3::new(0.4, 0., 0.2);
        let w = V3::new(0., 2., 0.);
        mo.corps.push(Corps {
            nom: "porteur".into(),
            m: 1.,
            j: M3::identity(),
            r: V3::new(0., 0., 0.3),
            r_bas: V3::zeros(),
            rot: M3::identity(),
            v: -w.cross(&p),
            w,
        });
        mo.contacts.push(Contact {
            nom: "appui".into(),
            corps: 0,
            p0: p,
            p1: p,
            rayon: 0.5,
            b: None,
            pb: V3::zeros(),
            p1b: V3::zeros(),
            rayon_b: 0.,
            normale: V3::z(),
            origine: V3::zeros(),
            k: 1e5,
            expo: 1.5,
            d_hat: 1e-3,
            c: 0.,
            mu: 0.,
            v_eps: 1e-3,
            demi: None,
            cylindre: None,
            maille: None,
            q_maille: V3::zeros(),
            s_axe: 0.,
            cable: false,
            sortie: (0., 0.),
            nonlisse: true,
            row: 0,
            fn_impose: None,
            restitution: 0.,
            auto: false,
        });
        mo.elems.push(Elem::K(Unilateral {
            nom: "appui".into(),
            ct: 0,
            a: 0,
            b: None,
        }));
        mo.k_pas = vec![(true, 0.)];
        // Le point de contact est instantanément immobile mais sa courbure
        // ferme le gap. λ = γ/(G M⁻¹ Gᵀ), ici un scalaire analytique.
        let gamma = -w.cross(&w.cross(&p)).z;
        assert!((mo.biais_acceleration(0.).unwrap()[0] - gamma).abs() < 1e-15);
        let lambda = gamma / (1. + p.x * p.x);
        let (a, l) = mo.acc_init(0.).unwrap();
        assert!((l[0] - lambda).abs() < 1e-14);
        let reference = DVector::from_vec(vec![0., 0., lambda, 0., -p.x * lambda, 0.]);
        assert!((a - reference).amax() < 1e-14);
    }

    fn modele_facteur_materiel(integree: bool) -> (Poutre, Vec<Corps>) {
        let b = Poutre {
            nom: "facteur matériel".into(),
            integree,
            a: 0,
            b: 1,
            l: 0.7,
            cn: [4e5, 8e4, 5e4],
            cm: [120., 350., 900.],
            ra: M3::identity(),
            rb: M3::identity(),
        };
        let mut cs = vec![
            Corps {
                nom: "a".into(),
                m: 1.,
                j: M3::identity(),
                r: V3::zeros(),
                r_bas: V3::zeros(),
                rot: M3::identity(),
                v: V3::zeros(),
                w: V3::zeros(),
            };
            2
        ];
        cs[1].nom = "b".into();
        cs[1].r.x = b.l;
        (b, cs)
    }

    #[test]
    fn facteur_materiel_repos_deformations_et_tangente() {
        for integree in [false, true] {
            let (b, cs) = modele_facteur_materiel(integree);
            let facteur = b.facteur_materiel(&cs).unwrap();
            assert_eq!(facteur.noeuds, [0, 1]);
            assert!(facteur.deformation_ponderee.norm() < 1e-12);
            let mut attendu = SMatrix::<f64, 6, 12>::zeros();
            // Référence indépendante : différences et déformations de la
            // poutre droite, dans les deux plans anisotropes de Timoshenko.
            for a in 0..3 {
                let c = if integree && a > 0 {
                    1. / (1. / b.cn[a] + b.l.powi(2) / (12. * b.cm[3 - a]))
                } else {
                    b.cn[a]
                };
                let translation = (c / b.l).sqrt();
                attendu[(a, a)] = -translation;
                attendu[(a, 6 + a)] = translation;
                if a > 0 {
                    let signe = if a == 1 { -1. } else { 1. };
                    let rotation = signe * (b.l * c).sqrt() / 2.;
                    attendu[(a, 6 - a)] = rotation;
                    attendu[(a, 12 - a)] = rotation;
                }
                let courbure = (b.cm[a] / b.l).sqrt();
                attendu[(3 + a, 3 + a)] = -courbure;
                attendu[(3 + a, 9 + a)] = courbure;
            }
            assert!((facteur.d - attendu).norm() < 2e-14 * attendu.norm());
            let materiel = facteur.d.transpose() * facteur.d;
            let tangente = -b.tangente(&cs);
            let ecart = DMatrix::from_fn(12, 12, |i, j| materiel[(i, j)] - tangente[(i, j)]);
            assert!(ecart.norm() < 2e-13 * tangente.norm());
        }
    }

    #[test]
    fn facteur_materiel_objectivite_et_travail_virtuel() {
        let q = expm(&V3::new(-0.7, 0.4, 0.9));
        let translation = V3::new(2., -0.7, 0.4);
        let mut changement = SMatrix::<f64, 12, 12>::zeros();
        for bloc in 0..4 {
            changement
                .fixed_view_mut::<3, 3>(3 * bloc, 3 * bloc)
                .copy_from(&q);
        }
        for integree in [false, true] {
            let (mut b, mut cs) = modele_facteur_materiel(integree);
            // Repères des corps distincts du repère matériel de la poutre.
            cs[0].rot = expm(&V3::new(0.2, -0.1, 0.3));
            cs[1].rot = expm(&V3::new(-0.1, 0.3, 0.15));
            b.ra = cs[0].rot.transpose();
            b.rb = cs[1].rot.transpose();
            cs[1].r += V3::new(0.01, 0.02, -0.008);
            cs[1].r_bas = V3::new(1e-18, -2e-18, 3e-18);
            cs[1].rot = expm(&V3::new(0.04, -0.02, 0.1)) * cs[1].rot;
            let facteur = b.facteur_materiel(&cs).unwrap();
            let force = -facteur.d.transpose() * facteur.deformation_ponderee;
            let reference = SVector::<f64, 12>::from_row_slice(&b.forces(&cs));
            assert!((force - reference).norm() < 2e-12 * reference.norm());
            let (pa, pb, bas) = b.poses_t(&cs, &|_| 0.);
            let energie = b.energie_coefficients_t(&pa, &pb, bas, b.coefficients(), b.cm);
            assert!(
                (0.5 * facteur.deformation_ponderee.norm_squared() - energie).abs()
                    < 2e-14 * energie
            );
            let mut tournes = cs.clone();
            for c in &mut tournes {
                c.r = q * c.r + translation;
                c.r_bas = q * c.r_bas;
                c.rot = q * c.rot;
            }
            let apres = b.facteur_materiel(&tournes).unwrap();
            assert!(
                (apres.deformation_ponderee - facteur.deformation_ponderee).norm()
                    < 2e-12 * facteur.deformation_ponderee.norm()
            );
            assert!((apres.d * changement - facteur.d).norm() < 2e-13 * facteur.d.norm());
            for axe in 0..3 {
                let mut mouvement = SVector::<f64, 12>::zeros();
                mouvement[axe] = 1.;
                mouvement[6 + axe] = 1.;
                assert!((facteur.d * mouvement).norm() < 2e-14 * facteur.d.norm());
                let mut vitesse = V3::zeros();
                vitesse[axe] = 1.;
                let corde = cs[1].r + cs[1].r_bas - cs[0].r - cs[0].r_bas;
                mouvement.fill(0.);
                mouvement
                    .fixed_rows_mut::<3>(6)
                    .copy_from(&vitesse.cross(&corde));
                mouvement[3 + axe] = 1.;
                mouvement[9 + axe] = 1.;
                assert!((facteur.d * mouvement).norm() < 2e-13 * facteur.d.norm());
            }
        }
    }

    #[test]
    fn facteur_materiel_ne_remplace_pas_la_tangente_precontrainte() {
        for integree in [false, true] {
            let (b, mut cs) = modele_facteur_materiel(integree);
            cs[1].r.x *= 0.97;
            cs[1].r.y = 0.02;
            cs[1].rot = expm(&V3::new(0.04, -0.02, 0.13));
            let facteur = b.facteur_materiel(&cs).unwrap();
            assert!(facteur.deformation_ponderee.norm() > 1.);
            let materiel = facteur.d.transpose() * facteur.d;
            let tangente = -b.tangente(&cs);
            let ecart = DMatrix::from_fn(12, 12, |i, j| tangente[(i, j)] - materiel[(i, j)]);
            assert!(ecart.norm() > 1e-3 * tangente.norm());
            assert!((&tangente - tangente.transpose()).norm() > 1e-7 * tangente.norm());
            // Le facteur garde pourtant exactement le gradient d'énergie.
            let force = -facteur.d.transpose() * facteur.deformation_ponderee;
            let reference = SVector::<f64, 12>::from_row_slice(&b.forces(&cs));
            assert!((force - reference).norm() < 2e-12 * reference.norm());
        }
    }

    #[test]
    fn facteur_materiel_refuse_parametres_et_branche_invalides() {
        let (b, mut cs) = modele_facteur_materiel(true);
        for invalide in [
            Poutre { a: 2, ..b.clone() },
            Poutre { a: 1, ..b.clone() },
            Poutre { l: 0., ..b.clone() },
            Poutre {
                cn: [1., -1., 1.],
                ..b.clone()
            },
        ] {
            assert!(invalide.facteur_materiel(&cs).is_err());
        }
        cs[1].rot = expm(&V3::new(PI, 0., 0.));
        assert!(b.facteur_materiel(&cs).is_err());
        cs[1].rot = M3::identity();
        cs[1].r.x = f64::NAN;
        assert!(b.facteur_materiel(&cs).is_err());
    }

    /// La tangente de poutre doit annuler les MODES RIGIDES (une translation
    /// d'ensemble ne change aucune force) et être quasi symétrique (elle dérive
    /// d'une énergie). Deux propriétés que les différences finies ne peuvent pas
    /// se donner à elles-mêmes.
    #[test]
    fn poutre_tangente_modes_rigides() {
        let mut mo = Modele::new(V3::new(0.0, 0.0, 0.0));
        let rot = expm(&V3::new(0.2, -0.1, 0.3));
        let rb = expm(&V3::new(0.1, 0.25, -0.2));
        mo.corps.push(Corps {
            r_bas: crate::V3::zeros(),
            nom: "a".into(),
            m: 0.1,
            j: M3::identity() * 1e-6,
            r: V3::new(0.0, 0.0, 0.0),
            rot,
            v: V3::zeros(),
            w: V3::zeros(),
        });
        mo.corps.push(Corps {
            r_bas: crate::V3::zeros(),
            nom: "b".into(),
            m: 0.1,
            j: M3::identity() * 1e-6,
            r: V3::new(0.09, 0.02, -0.01),
            rot: rb,
            v: V3::zeros(),
            w: V3::zeros(),
        });
        let b = Poutre {
            nom: "p".into(),
            integree: true,
            a: 0,
            b: 1,
            l: 0.1,
            cn: [7e6, 2.7e6, 2.7e6],
            cm: [8.0, 14.6, 233.0],
            ra: rot.transpose(),
            rb: rb.transpose(),
        };
        let k = b.tangente_exacte(&mo.corps);
        let ech = (0..12)
            .map(|i| (0..12).map(|j| k[(i, j)].abs()).fold(0.0, f64::max))
            .fold(0.0, f64::max);
        for tr in 0..3 {
            for row in 0..12 {
                let s: f64 = k[(row, tr)] + k[(row, 6 + tr)];
                assert!(
                    s.abs() < 1e-6 * ech,
                    "mode rigide {tr} : ligne {row} rend {s} (échelle {ech})"
                );
            }
        }
        // la tangente EXACTE (duaux emboîtés) doit coïncider avec les DF du champ
        let kd = b.tangente(&mo.corps);
        // seuil : relatif à la composante, PLUS le bruit de la DF elle-même
        // (plancher d'arrondi des forces ≈ 1e-9 N, divisé par 2ε)
        let bruit = 1e-9 / (2.0 * 1e-5);
        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (k[(i, j)] - kd[(i, j)]).abs() < 1e-4 * kd[(i, j)].abs() + 10.0 * bruit,
                    "K[{i},{j}] : AD {} vs DF {}",
                    k[(i, j)],
                    kd[(i, j)]
                );
            }
        }
        // PAS de contrôle de symétrie : la tangente d'une poutre géométriquement
        // exacte sous perturbation gauche est NON symétrique hors équilibre —
        // c'est un fait établi (Simo–Vu Quoc), pas un défaut. Mesuré ici :
        // K[3,4] = 23 475 contre K[4,3] = −7 825. Ce qui doit tenir est le mode
        // rigide, et il tient.
        let sous_bloc = k[(0, 0)].abs() + k[(6, 6)].abs();
        assert!(sous_bloc > 0.0 && ech.is_finite());
    }

    #[test]
    fn poutre_blocs_mixtes_exacts() {
        let mut corps = vec![
            Corps {
                r_bas: crate::V3::zeros(),
                nom: "a".into(),
                m: 1.0,
                j: M3::identity(),
                r: V3::new(-0.3, 0.2, 0.4),
                rot: expm(&V3::new(0.2, -0.1, 0.3)),
                v: V3::zeros(),
                w: V3::zeros(),
            },
            Corps {
                r_bas: crate::V3::zeros(),
                nom: "b".into(),
                m: 1.0,
                j: M3::identity(),
                r: V3::new(0.6, 0.5, -0.1),
                rot: expm(&V3::new(0.1, 0.25, -0.2)),
                v: V3::zeros(),
                w: V3::zeros(),
            },
        ];
        for integree in [false, true] {
            let b = Poutre {
                nom: "p".into(),
                integree,
                a: 0,
                b: 1,
                l: 1.0,
                cn: [1e5, 2000.0, 3000.0],
                cm: [30.0, 50.0, 20.0],
                ra: M3::identity(),
                rb: M3::identity(),
            };
            // Le grand déplacement amplifie l'annulation qui polluait les
            // colonnes de translation calculées par différences finies.
            for amplitude in [1.0, 1e6] {
                corps[1].r = V3::new(0.6, 0.5, -0.1) * amplitude;
                corps[1].r_bas = V3::new(1e-18, -2e-18, 3e-18) * amplitude;
                let reference = b.tangente_exacte(&corps);
                let k = b.tangente(&corps);
                for i in 0..12 {
                    for j in 0..12 {
                        if i % 6 < 3 || j % 6 < 3 {
                            assert!(
                                (k[(i, j)] - reference[(i, j)]).abs()
                                    < 1e-10 * (1.0 + reference[(i, j)].abs()),
                                "intégrée={integree} amplitude={amplitude} K[{i},{j}] : {} vs {}",
                                k[(i, j)],
                                reference[(i, j)]
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn poutre_moments_analytiques_et_travail_virtuel() {
        use crate::ad::Dual;
        use std::f64::consts::PI;

        let frame = expm(&V3::new(0.4, -0.2, 0.1));
        let axis = V3::new(1., -2., 3.).normalize();
        let mut corps = vec![
            Corps {
                nom: "a".into(),
                m: 1.,
                j: M3::identity(),
                r: V3::zeros(),
                r_bas: V3::zeros(),
                rot: frame,
                v: V3::zeros(),
                w: V3::zeros(),
            };
            2
        ];
        corps[1].r = frame * V3::new(1.01, 0.03, -0.02);
        corps[1].r_bas = V3::new(1e-19, -2e-19, 3e-19);
        for integree in [false, true] {
            let beam = Poutre {
                nom: "p".into(),
                integree,
                a: 0,
                b: 1,
                l: 1.,
                cn: [1e5, 2000., 3000.],
                cm: [30., 50., 20.],
                ra: M3::identity(),
                rb: M3::identity(),
            };
            for angle in [
                0.,
                1e-9,
                1e-5,
                0.1,
                0.4999,
                0.5001,
                1.,
                2.5,
                PI - 0.01,
                PI - 1e-4,
            ] {
                corps[1].rot = frame * expm(&(angle * axis));
                let (pa, pb, bas) = beam.poses_t(&corps, &|_| 0.);
                let f = beam.forces_analytiques_t(&pa, &pb, bas, beam.coefficients(), beam.cm);
                let reference =
                    beam.forces_energie_coefficients(&corps, beam.coefficients(), beam.cm);
                for i in 0..12 {
                    assert!(
                        (f[i] - reference[i]).abs() <= 1e-9 * (1. + reference[i].abs()),
                        "integree={integree} angle={angle} f[{i}]={} energie={}",
                        f[i],
                        reference[i]
                    );
                }
                let force = V3::new(f[0], f[1], f[2]);
                let moment = V3::new(f[3] + f[9], f[4] + f[10], f[5] + f[11]);
                assert!(
                    (moment - (corps[1].r + corps[1].r_bas).cross(&force)).norm()
                        <= 1e-10 * (1. + force.norm())
                );
                type D = Dual<f64, 12>;
                let (pa, pb, bas) = beam.poses_t(&corps, &|i| D::var(0., i));
                let fd = beam.forces_analytiques_t(&pa, &pb, bas, beam.coefficients(), beam.cm);
                let exact = beam.tangente_exacte(&corps);
                for i in 0..12 {
                    for j in 0..12 {
                        assert!(
                            (fd[i].d[j] - exact[(i, j)]).abs() <= 1e-8 * (1. + exact[(i, j)].abs()),
                            "integree={integree} angle={angle} K[{i},{j}]={} energie={}",
                            fd[i].d[j],
                            exact[(i, j)]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn poutre_gradient_reperes_unites_et_coefficients_signes() {
        let q = expm(&V3::new(-0.7, 0.2, 0.5));
        for integree in [false, true] {
            for longueur in [1e-6, 1., 1e6] {
                let beam = Poutre {
                    nom: "p".into(),
                    integree,
                    a: 0,
                    b: 1,
                    l: longueur,
                    cn: [1e5, 2000., 3000.],
                    cm: [30., 50., 20.].map(|c| c * longueur * longueur),
                    ra: expm(&V3::new(0.3, -0.1, 0.2)),
                    rb: expm(&V3::new(-0.1, 0.4, -0.2)),
                };
                let mut corps = vec![
                    Corps {
                        nom: "a".into(),
                        m: 1.,
                        j: M3::identity(),
                        r: V3::zeros(),
                        r_bas: V3::zeros(),
                        rot: beam.ra.transpose(),
                        v: V3::zeros(),
                        w: V3::zeros(),
                    };
                    2
                ];
                corps[1].r = longueur * V3::new(0.9, 0.2, -0.3);
                corps[1].r_bas = longueur * V3::new(1e-12, -2e-12, 3e-12);
                corps[1].rot = expm(&V3::new(0.8, -0.2, 0.5)) * beam.rb.transpose();
                // La linéarité en rigidités sert aux sensibilités : leurs
                // dérivées peuvent être nulles ou négatives.
                for (cn, cm) in [
                    (beam.coefficients(), beam.cm),
                    (
                        [0., -1200., 700.],
                        [12., -15., 0.].map(|c| c * longueur * longueur),
                    ),
                ] {
                    let f = beam.forces_coefficients(&corps, cn, cm);
                    let reference = beam.forces_energie_coefficients(&corps, cn, cm);
                    let k = beam.tangente_coefficients(&corps, cn, cm);
                    let exact = beam.tangente_exacte_coefficients(&corps, cn, cm);
                    let unite = |i: usize| if i % 6 < 3 { 1. } else { longueur };
                    let fnorm = (0..12)
                        .map(|i| (reference[i] / unite(i)).abs())
                        .fold(1., f64::max);
                    let mut knorm = 1_f64;
                    for i in 0..12 {
                        for j in 0..12 {
                            knorm =
                                knorm.max((exact[(i, j)] * longueur / (unite(i) * unite(j))).abs());
                        }
                    }
                    for i in 0..12 {
                        assert!((f[i] - reference[i]).abs() / unite(i) < 2e-11 * fnorm);
                        for j in 0..12 {
                            assert!(
                                (k[(i, j)] - exact[(i, j)]).abs() * longueur
                                    / (unite(i) * unite(j))
                                    < 2e-11 * knorm
                            );
                        }
                    }
                    let mut tournes = corps.clone();
                    for c in &mut tournes {
                        c.r = q * c.r;
                        c.r_bas = q * c.r_bas;
                        c.rot = q * c.rot;
                        // Origine très éloignée, avec addition compensée.
                        c.deplace(longueur * V3::new(1e12, -2e12, 3e12));
                    }
                    let fq = beam.forces_coefficients(&tournes, cn, cm);
                    let kq = beam.tangente_coefficients(&tournes, cn, cm);
                    for block in 0..4 {
                        let refq = q * V3::new(f[3 * block], f[3 * block + 1], f[3 * block + 2]);
                        for i in 0..3 {
                            assert!(
                                (fq[3 * block + i] - refq[i]).abs() / unite(3 * block + i)
                                    < 2e-11 * fnorm
                            );
                        }
                        for other in 0..4 {
                            let km = M3::from_fn(|i, j| k[(3 * block + i, 3 * other + j)]);
                            let refq = q * km * q.transpose();
                            for i in 0..3 {
                                for j in 0..3 {
                                    let row = 3 * block + i;
                                    let col = 3 * other + j;
                                    assert!(
                                        (kq[(row, col)] - refq[(i, j)]).abs() * longueur
                                            / (unite(row) * unite(col))
                                            < 2e-11 * knorm
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn poutre_gradient_coupure_pi_explicite() {
        let mut corps = vec![
            Corps {
                nom: "a".into(),
                m: 1.,
                j: M3::identity(),
                r: V3::zeros(),
                r_bas: V3::zeros(),
                rot: M3::identity(),
                v: V3::zeros(),
                w: V3::zeros(),
            };
            2
        ];
        corps[1].r.x = 1.;
        corps[1].rot = expm(&V3::new(std::f64::consts::PI, 0., 0.));
        for integree in [false, true] {
            let beam = Poutre {
                nom: "p".into(),
                integree,
                a: 0,
                b: 1,
                l: 1.,
                cn: [1.; 3],
                cm: [1.; 3],
                ra: M3::identity(),
                rb: M3::identity(),
            };
            assert!(beam.forces(&corps).iter().all(|v| v.is_nan()));
            assert!(beam.tangente(&corps).iter().all(|v| v.is_nan()));
        }
    }

    #[test]
    fn poutre_allongement_inferieur_ulp() {
        let mut corps = vec![
            Corps {
                r_bas: V3::zeros(),
                nom: "a".into(),
                m: 1.,
                j: M3::identity(),
                r: V3::zeros(),
                rot: M3::identity(),
                v: V3::zeros(),
                w: V3::zeros(),
            };
            2
        ];
        corps[1].r.x = 1.;
        corps[1].r_bas.x = 2_f64.powi(-60);
        for integree in [false, true] {
            let b = Poutre {
                nom: "p".into(),
                integree,
                a: 0,
                b: 1,
                l: 1.,
                cn: [2_f64.powi(30), 1., 1.],
                cm: [1.; 3],
                ra: M3::identity(),
                rb: M3::identity(),
            };
            let (g, _) = b.deformations(&corps);
            assert_eq!(g[0], 2_f64.powi(-60));
            let f = b.forces(&corps);
            assert_eq!(f[0], 2_f64.powi(-30));
            assert_eq!(f[6], -f[0]);
            let k = b.tangente(&corps);
            let exact = b.tangente_exacte(&corps);
            assert_eq!(k[(0, 0)], -2_f64.powi(30));
            for i in 0..12 {
                for j in 0..12 {
                    if i % 6 < 3 || j % 6 < 3 {
                        assert!(
                            (k[(i, j)] - exact[(i, j)]).abs() < 1e-12 * (1. + exact[(i, j)].abs())
                        );
                    }
                }
            }
        }
    }

    /// Un angle d'engrenage ou de vis qui saute d'un quart de tour rend Φ
    /// indéfini : `residu_pose` y met 1e30 et Newton recule. `eval_ad` doit
    /// rendre Err, pas paniquer au travers de PyO3 (6 sept.).
    #[test]
    fn eval_ad_angle_saute_rend_err() {
        let mut mo = Modele::new(V3::zeros());
        for (nom, x) in [("a", 0.0), ("b", 0.1)] {
            mo.corps.push(Corps {
                r_bas: crate::V3::zeros(),
                nom: nom.into(),
                m: 1.0,
                j: M3::identity() * 1e-3,
                r: V3::new(x, 0.0, 0.0),
                rot: M3::identity(),
                v: V3::zeros(),
                w: V3::zeros(),
            });
        }
        let e = Elem::E(Engrenage {
            nom: "e".into(),
            a: Some(0),
            b: Some(1),
            c: None,
            na: V3::z(),
            nb: V3::z(),
            rapport: -2.0,
            ref_a: M3::identity(),
            ref_b: M3::identity(),
            prev: (3.0, 0.0),
        });
        let v = Elem::V(Vis {
            nom: "v".into(),
            a: Some(0),
            b: Some(1),
            n: V3::z(),
            pas: 1e-3,
            ref_rel: M3::identity(),
            c0: 0.0,
            prev: 3.0,
        });
        for el in [&e, &v] {
            let r = el.eval_ad(&mo.corps, 0.0, &[1.0]);
            assert!(r.is_err(), "{} : Err attendu", el.nom());
            assert!(r.unwrap_err().contains("quart de tour"));
        }
    }

    /// G par duaux = G analytique de `phi_g`, sur une liaison générale et une bielle.
    #[test]
    fn g_dual_egal_g_analytique() {
        let mut mo = Modele::new(V3::new(0.0, 0.0, -9.81));
        let rot = expm(&V3::new(0.3, -0.2, 0.5));
        mo.corps.push(Corps {
            r_bas: crate::V3::zeros(),
            nom: "a".into(),
            m: 1.0,
            j: M3::identity() * 1e-3,
            r: V3::new(0.1, 0.2, 0.3),
            rot,
            v: V3::zeros(),
            w: V3::zeros(),
        });
        mo.corps.push(Corps {
            r_bas: crate::V3::zeros(),
            nom: "b".into(),
            m: 2.0,
            j: M3::identity() * 2e-3,
            r: V3::new(-0.2, 0.4, 0.1),
            rot: expm(&V3::new(-0.4, 0.1, 0.2)),
            v: V3::zeros(),
            w: V3::zeros(),
        });
        // Petites parties amplifiées pour rendre toute omission visible
        // dans les gradients, indépendamment de la normalisation du stockage.
        mo.corps[0].r_bas = V3::new(1e-3, -2e-3, 3e-3);
        mo.corps[1].r_bas = V3::new(-2e-3, 1e-3, 2e-3);
        let l = Liaison {
            nom: "l".into(),
            a: Some(0),
            b: Some(1),
            pa: V3::new(0.05, -0.02, 0.03),
            ra: expm(&V3::new(0.1, 0.2, -0.3)),
            pb: V3::new(-0.01, 0.04, 0.02),
            rb: expm(&V3::new(0.2, -0.1, 0.4)),
            bt: vec![0, 1, 2],
            br: vec![0, 1, 2],
            cible_t: None,
            cible_r: Some((V3::new(0.0, 0.0, 1.0), Loi::Lineaire { a0: 0.2, taux: 1.0 })),
            nh: false,
        };
        let (phi, g) = l.phi_g(&mo.corps, 0.7);
        let el = Elem::L(l.clone());
        let (phi2, g2, _) = el.eval_ad2(&mo.corps, 0.7, &[0.0; 6]);
        for i in 0..6 {
            assert!((phi[i] - phi2[i]).abs() < 1e-14);
            for c in 0..12 {
                assert!(
                    (g[(i, c)] - g2[(i, c)]).abs() < 1e-12,
                    "G[{i},{c}] {} vs {}",
                    g[(i, c)],
                    g2[(i, c)]
                );
            }
        }
        let d = Distance {
            nom: "d".into(),
            a: Some(0),
            b: Some(1),
            pa: V3::new(0.05, 0.0, 0.1),
            pb: V3::new(0.0, 0.06, -0.02),
            l: 0.3,
        };
        let (_, gd) = d.phi_g(&mo.corps);
        let (_, gd2, _) = Elem::D(d.clone()).eval_ad2(&mo.corps, 0.0, &[0.0]);
        for c in 0..12 {
            assert!((gd[(0, c)] - gd2[(0, c)]).abs() < 1e-12);
        }
        // K = ∂(Gᵀλ)/∂δ : la dérivée du CHAMP de gradient sous perturbation
        // gauche (ce que Newton exige), contrôlée par différences finies
        // centrées du gradient analytique. NB : la Hessienne de λ·Φ en
        // coordonnées exponentielles (duaux emboîtés, `eval_ad2`) n'est PAS la
        // même matrice — elle porte en plus la courbure de exp (½[δθ]×²) ;
        // les deux ont été confondues une heure, la différence est un fait de
        // géométrie, pas un bug.
        let lam = [0.3, -1.2, 0.7, 2.0, -0.4, 0.9];
        let (_, _, k1) = el.eval_ad(&mo.corps, 0.7, &lam).unwrap();
        let gl_a = |mo2: &Modele| -> Vec<f64> {
            let (_, g) = l.phi_g(&mo2.corps, 0.7);
            (0..12)
                .map(|c| (0..6).map(|i| lam[i] * g[(i, c)]).sum())
                .collect()
        };
        let eps = 1e-6;
        for b in 0..12 {
            let (mut mp, mut mm) = (mo.clone(), mo.clone());
            let (ib, comp) = (b / 6, b % 6);
            let mut dp = [0.0; 6];
            dp[comp] = eps;
            let pert = |m: &mut Modele, s: f64| {
                let c = &mut m.corps[ib];
                c.r += s * V3::new(dp[0], dp[1], dp[2]);
                c.rot = expm(&(s * V3::new(dp[3], dp[4], dp[5]))) * c.rot;
            };
            pert(&mut mp, 1.0);
            pert(&mut mm, -1.0);
            let (gp, gm) = (gl_a(&mp), gl_a(&mm));
            for a in 0..12 {
                let fd = (gp[a] - gm[a]) / (2.0 * eps);
                assert!(
                    (k1[(a, b)] - fd).abs() < 1e-7,
                    "K[{a},{b}] {} vs DF {}",
                    k1[(a, b)],
                    fd
                );
            }
        }
        // idem bielle et engrenage
        let d2 = Elem::D(d.clone());
        let (_, _, kd1) = d2.eval_ad(&mo.corps, 0.0, &[1.7]).unwrap();
        let e = Engrenage {
            nom: "e".into(),
            a: Some(0),
            b: Some(1),
            c: None,
            na: V3::new(0.0, 0.0, 1.0),
            nb: V3::new(0.0, 0.0, 1.0),
            rapport: -4.0,
            ref_a: mo.corps[0].rot,
            ref_b: mo.corps[1].rot,
            prev: (0.0, 0.0),
        };
        let e2 = Elem::E(e.clone());
        let (_, _, ke1) = e2.eval_ad(&mo.corps, 0.0, &[0.8]).unwrap();
        for (el2, kk, lamv, nk) in [(&d2, &kd1, 1.7, 12usize), (&e2, &ke1, 0.8, 18usize)] {
            for b in 0..nk.min(12) {
                let (mut mp, mut mm) = (mo.clone(), mo.clone());
                let (ib, comp) = (b / 6, b % 6);
                let mut dp = [0.0; 6];
                dp[comp] = eps;
                for (m, sg) in [(&mut mp, 1.0), (&mut mm, -1.0)] {
                    let c = &mut m.corps[ib];
                    c.r += sg * V3::new(dp[0], dp[1], dp[2]);
                    c.rot = expm(&(sg * V3::new(dp[3], dp[4], dp[5]))) * c.rot;
                }
                let g_of = |m: &Modele| -> Vec<f64> {
                    let (_, g) = match el2 {
                        Elem::D(dd) => dd.phi_g(&m.corps),
                        Elem::E(ee) => ee.phi_g(&m.corps).unwrap(),
                        _ => unreachable!(),
                    };
                    (0..g.ncols()).map(|c| lamv * g[(0, c)]).collect()
                };
                let (gp, gm) = (g_of(&mp), g_of(&mm));
                for a in 0..gp.len() {
                    let fd = (gp[a] - gm[a]) / (2.0 * eps);
                    assert!(
                        (kk[(a, b)] - fd).abs() < 1e-7,
                        "K[{a},{b}] {} vs DF {}",
                        kk[(a, b)],
                        fd
                    );
                }
            }
        }
    }
}
