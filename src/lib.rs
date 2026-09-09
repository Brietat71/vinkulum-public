//! Vinkulum — noyau multicorps : intégrateur, liaisons, contact, poutres, pales.
//!
//! ## Lints autorisés au niveau du crate, et pourquoi
//!
//! La CI passe `clippy -D warnings` : ce qui suit est donc une DÉCISION, pas
//! un oubli. Chaque lint ci-dessous est stylistique ou contraire à une
//! intention du code ; aucun ne signale un défaut ici.
//!
//! · `type_complexity` — les types du solveur sont des tuples de matrices et
//!   de vecteurs nommés par leur rôle dans la formulation ; un alias les
//!   rendrait plus courts et moins lisibles ;
//! · `too_many_arguments` — une fonction de contact ou de pale porte la
//!   géométrie ET la loi ; les grouper en struct déplacerait le bruit ;
//! · `neg_cmp_op_on_partial_ord` — `!(x > 0.0)` est VOULU : il attrape aussi
//!   NaN, ce que `x <= 0.0` ne fait pas. C'est exactement ce qu'on veut d'une
//!   garde de domaine, et c'est le contraire d'une maladresse ;
//! · `manual_is_multiple_of` — `n % 3 != 0` se lit mieux que la méthode dans
//!   une validation de tableau à plat.
#![allow(
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::manual_is_multiple_of,
    // · `needless_range_loop` — dans une boucle qui indexe PLUSIEURS tableaux
    //   au même rang (résidu, jacobien, blocs de corps), l'itérateur ne rend
    //   pas le code plus clair : il oblige à un `zip` de trois termes là où
    //   `for k in 0..n` dit exactement ce que fait la formule.
    clippy::needless_range_loop
)]
//!
//! Projet DISTINCT de FRELON (décision du 3 sept. 2026) : FRELON s'appuie sur
//! lui, il ne lui appartient pas. Nom : latin *vinculum*, « lien » — ce que
//! le noyau résout, des corps tenus par des contraintes.
//!
//! CLEAN-ROOM (décision projet du 3 sept. 2026) : MBDyn et OCCT se lisent,
//! ils ne se copient pas. La formulation est celle de la maquette Python
//! `python/vinkulum/maquette.py`, qui reste le contre-solveur (R1 : deux implémentations
//! indépendantes de la même formulation doivent coïncider).
//!
//! · corps rigides en coordonnées absolues : r ∈ R³, R ∈ SO(3) (matrice),
//!   vitesses v, ω SPATIALES ; mise à jour de R sur le groupe (Rodrigues) ;
//! · liaisons par contraintes Φ(q) = 0 + multiplicateurs λ : une liaison
//!   générique à deux repères et six degrés relatifs bloqués à la carte,
//!   cible imposée possible (loi de commande) ; bielle `Distance` ;
//! · index 3 résolu directement par α-généralisé (Chung–Hulbert) sur groupe
//!   de Lie (Brüls–Cardona), index 3 selon Arnold–Brüls : Φ tenu à chaque
//!   pas, un seul paramètre ρ∞ ;
//! · Newton sur (u̇, λ), jacobien EXACT du résidu par différentiation
//!   automatique élément par élément (`tangent`, nombres duaux de `ad`),
//!   assemblé en triplets creux — blocs d'éléments EN PARALLÈLE (rayon) à
//!   partir de `tangent::PARALLELE_MIN` ; LU dense sous `DENSE_MAX` inconnues,
//!   creux (faer) au-delà. Contraintes redondantes neutralisées par QR pivoté
//!   au départ, repli SVD sinon : elles ne cassent pas.

pub mod ad;
mod analyse;
mod assemblage;
mod certificat;
mod certificat_vitesse;
mod contraintes;
mod equilibrage;
mod initialisation;
mod moindres_carres;
mod numerique;
mod partition;
mod raffinement;
use partition::{Indices, Reduction};
mod robustesse;
mod rotation;
mod schema_lie;
mod superelement;
mod svd;
pub mod tangent;
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use pyo3::prelude::*;
use robustesse::{fini, reglages, Sauvegarde};

pub type M3 = Matrix3<f64>;
pub type V3 = Vector3<f64>;

// ── SO(3) ────────────────────────────────────────────────────────────────────
pub fn skew(w: &V3) -> M3 {
    M3::new(0.0, -w.z, w.y, w.z, 0.0, -w.x, -w.y, w.x, 0.0)
}

pub fn vee(s: &M3) -> V3 {
    V3::new(s[(2, 1)], s[(0, 2)], s[(1, 0)])
}

/// exp([w]×), Rodrigues — `ad::expm` à T = f64 : une seule définition.
/// log sur SO(3) en f64 : le vecteur rotation θ tel que R = exp([θ]×)
pub fn log_so3(r: &M3) -> V3 {
    V3::from(tangent::logm::<f64>(&ad::m_from(r)))
}

pub fn expm(w: &V3) -> M3 {
    let o = ad::expm::<f64>([w.x, w.y, w.z]);
    M3::from_fn(|i, j| o[i][j])
}

// ── lois de commande ─────────────────────────────────────────────────────────
/// Valeur imposée dans le temps : constante, rampe, ou table interpolée
/// linéairement (le `multilinear` de MBDyn — abscisses STRICTEMENT croissantes,
/// leçon du 2 sept. : deux points à la même abscisse font viser le second).
#[derive(Clone, Debug)]
pub enum Loi {
    Constante(f64),
    Lineaire { a0: f64, taux: f64 },
    Table(Vec<(f64, f64)>),
}

impl Loi {
    /// dérivée en t (0 hors table ; pente du segment dedans)
    pub fn derivee(&self, t: f64) -> f64 {
        match self {
            Loi::Constante(_) => 0.0,
            Loi::Lineaire { taux, .. } => *taux,
            Loi::Table(pts) => pts
                .windows(2)
                .find(|w| t <= w[1].0 && t >= w[0].0)
                .map_or(0.0, |w| (w[1].1 - w[0].1) / (w[1].0 - w[0].0)),
        }
    }
    pub fn valeur(&self, t: f64) -> f64 {
        match self {
            Loi::Constante(c) => *c,
            Loi::Lineaire { a0, taux } => a0 + taux * t,
            Loi::Table(pts) => {
                if t <= pts[0].0 {
                    return pts[0].1;
                }
                for w in pts.windows(2) {
                    let ((t0, v0), (t1, v1)) = (w[0], w[1]);
                    if t <= t1 {
                        return v0 + (v1 - v0) * (t - t0) / (t1 - t0);
                    }
                }
                pts[pts.len() - 1].1
            }
        }
    }

    pub fn table(pts: Vec<(f64, f64)>) -> Result<Loi, String> {
        if pts.is_empty() {
            return Err("table vide".into());
        }
        if pts.iter().any(|(t, v)| !t.is_finite() || !v.is_finite()) {
            return Err("table : abscisses et valeurs finies requises".into());
        }
        if pts.windows(2).any(|w| w[1].0 <= w[0].0) {
            return Err("abscisses non strictement croissantes".into());
        }
        Ok(Loi::Table(pts))
    }
}

// ── corps ────────────────────────────────────────────────────────────────────
#[derive(Clone, Debug)]
pub struct Corps {
    pub nom: String,
    pub m: f64,
    pub j: M3, // au CdM, axes corps
    pub r: V3,
    /// Partie basse des translations, conservée lors des incréments et sauvegardes.
    pub r_bas: V3,
    pub rot: M3,
    pub v: V3,
    pub w: V3,
}

/// `None` = le bâti (pose fixe, aucune inconnue).
pub type Ref = Option<usize>;

// ── liaisons ─────────────────────────────────────────────────────────────────
#[derive(Clone, Debug)]
pub struct Liaison {
    pub nom: String,
    pub a: Ref,
    pub b: Ref,
    pub pa: V3,
    pub ra: M3,
    pub pb: V3,
    pub rb: M3,
    pub bt: Vec<usize>,
    pub br: Vec<usize>,
    /// translation imposée : direction (repère liaison) × loi(t)
    pub cible_t: Option<(V3, Loi)>,
    /// rotation imposée : axe (repère liaison), angle loi(t)
    pub cible_r: Option<(V3, Loi)>,
    /// CONTRAINTE NON HOLONOME : les lignes de cette liaison sont imposées au
    /// niveau VITESSE (`G·u = 0`) et non position (`Φ = 0`).
    ///
    /// C'est le roulement sans glissement, et c'est ce qui manquait à ce noyau
    /// pour qu'un modèle de véhicule ou de robot mobile puisse y entrer : la
    /// vitesse du point de contact est nulle, mais aucune fonction de la
    /// POSITION ne l'exprime — la roue peut atteindre n'importe quelle pose en
    /// roulant. Le jacobien est le même que celui de la liaison holonome
    /// correspondante ; seul le niveau auquel on l'impose change.
    pub nh: bool,
}

#[derive(Clone, Debug)]
pub struct Distance {
    pub nom: String,
    pub a: Ref,
    pub b: Ref,
    pub pa: V3,
    pub pb: V3,
    pub l: f64,
}

/// Engrenage : θ_a − rapport·θ_b = 0, θ mesurés autour des axes `na`/`nb`
/// (repère du PORTEUR `c`, bâti si None) depuis la pose initiale. Rapport
/// NÉGATIF pour deux roues extérieures (elles tournent en sens inverse).
/// C'est la liaison que MBDyn n'a pas : la réaction de denture devient une
/// SORTIE (le multiplicateur λ, couple sur a ; −rapport·λ sur b).
///
/// Les angles sont DÉROULÉS : atan2 rend θ modulo 2π, on choisit la branche
/// la plus proche de l'angle continu du dernier pas accepté (`prev`). Ça
/// suppose moins d'un demi-tour par pas — vérifié, sinon on lève.
#[derive(Clone, Debug)]
pub struct Engrenage {
    pub nom: String,
    pub a: Ref,
    pub b: Ref,
    pub c: Ref,
    pub na: V3,
    pub nb: V3,
    pub rapport: f64,
    pub ref_a: M3, // (R_c₀ᵀ R_a₀) : rotation relative initiale, l'origine des angles
    pub ref_b: M3,
    pub prev: (f64, f64),
}

fn angle_autour(e: &M3, n: &V3) -> f64 {
    // E rotation d'angle θ autour de n : vee(skew E) = sinθ·n, tr E = 1 + 2cosθ
    let s = vee(&(0.5 * (e - e.transpose()))).dot(n);
    let c = 0.5 * (e.trace() - 1.0);
    s.atan2(c)
}

impl Engrenage {
    fn angles(&self, corps: &[Corps]) -> Result<(f64, f64), String> {
        let (_, rc) = pose_de(corps, self.c);
        let (_, ra) = pose_de(corps, self.a);
        let (_, rb) = pose_de(corps, self.b);
        let ea = rc.transpose() * ra * self.ref_a.transpose();
        let eb = rc.transpose() * rb * self.ref_b.transpose();
        let mut out = [0.0; 2];
        for (k, (e, n, p)) in [(ea, self.na, self.prev.0), (eb, self.nb, self.prev.1)]
            .iter()
            .enumerate()
        {
            let w = angle_autour(e, n);
            let tours = ((p - w) / (2.0 * std::f64::consts::PI)).round();
            let th = w + 2.0 * std::f64::consts::PI * tours;
            if (th - p).abs() > std::f64::consts::FRAC_PI_2 {
                return Err(format!("{} : plus d'un quart de tour en un pas ({:.3} rad) — pas de temps trop grand pour dérouler l'angle", self.nom, th - p));
            }
            out[k] = th;
        }
        Ok((out[0], out[1]))
    }

    pub fn phi_g(&self, corps: &[Corps]) -> Result<(DVector<f64>, DMatrix<f64>), String> {
        let (tha, thb) = self.angles(corps)?;
        // G EXACT par duaux (colonnes : 0..6 a, 6..12 b, 12..18 porteur). La
        // formule à la main — n_w sur les rotations — n'est juste que pour un
        // arbre exactement sur son axe ; cf. `Engrenage::g_exact`.
        let g = self.g_exact(corps);
        Ok((DVector::from_element(1, tha - self.rapport * thb), g))
    }
}

/// VIS–ÉCROU (glissière hélicoïdale) — la liaison qui lie une TRANSLATION à
/// une ROTATION autour du même axe, et qui manquait au noyau. Un vérin à vis,
/// une table de machine-outil, un actionneur linéaire électrique : tous en
/// sont faits, et aucun ne se modélise avec une glissière plus un pivot.
///
/// ```text
///     Φ = n·(r_b − r_a) − pas·θ − c₀ ,   pas en m/rad
/// ```
///
/// θ est l'angle relatif DÉROULÉ de b autour de n (même machinerie que
/// l'engrenage : sans déroulement, la contrainte saute d'un tour). `c₀` est
/// la valeur à la création, donc Φ = 0 sur la pose de départ — comme toutes
/// les autres liaisons du noyau.
#[derive(Clone, Debug)]
pub struct Vis {
    pub nom: String,
    pub a: Ref,
    pub b: Ref,
    pub n: V3,       // axe, repère de a
    pub pas: f64,    // m par radian
    pub ref_rel: M3, // rotation relative initiale : l'origine des angles
    pub c0: f64,
    pub prev: f64,
}

impl Vis {
    fn angle(&self, corps: &[Corps]) -> Result<f64, String> {
        let (_, ra) = pose_de(corps, self.a);
        let (_, rb) = pose_de(corps, self.b);
        let e = ra.transpose() * rb * self.ref_rel.transpose();
        let w = angle_autour(&e, &self.n);
        let tours = ((self.prev - w) / (2.0 * std::f64::consts::PI)).round();
        let th = w + 2.0 * std::f64::consts::PI * tours;
        if (th - self.prev).abs() > std::f64::consts::FRAC_PI_2 {
            return Err(format!(
                "{} : plus d'un quart de tour en un pas ({:.3} rad) — \
                                pas de temps trop grand pour dérouler l'angle",
                self.nom,
                th - self.prev
            ));
        }
        Ok(th)
    }

    pub fn phi_g(&self, corps: &[Corps]) -> Result<(DVector<f64>, DMatrix<f64>), String> {
        let th = self.angle(corps)?;
        let (pa_, ra) = pose_de(corps, self.a);
        let (pb_, _) = pose_de(corps, self.b);
        let nw = ra * self.n;
        let (d, db) = tangent::difference_positions::<f64>(
            pa_.into(),
            partie_basse(corps, self.a).into(),
            pb_.into(),
            partie_basse(corps, self.b).into(),
        );
        let (d, db) = (V3::from(d), V3::from(db));
        // G EXACT par duaux (colonnes : 0..6 a, 6..12 b), cf. `Vis::g_exact`
        let g = self.g_exact(corps);
        Ok((
            DVector::from_element(1, (nw.dot(&d) - self.c0) - self.pas * th + nw.dot(&db)),
            g,
        ))
    }
}

/// Joint de **Cardan** (croix de Hooke) : les deux axes de la croix restent
/// perpendiculaires.
///
/// ```text
/// Φ = (R_a·n_a) · (R_b·n_b) = 0
/// ```
///
/// UNE contrainte scalaire — à composer avec une rotule (`liaison` bloquée en
/// translation seule) pour obtenir le joint complet à deux rotations libres.
///
/// Ce n'est PAS ce que rend `liaison(bloque_r = [k])` : celle-là annule la
/// composante k de la partie antisymétrique de RᵀR', qui vaut
/// ½(uₐ₃·u_b₂ − uₐ₂·u_b₃) — le second terme est parasite, et les deux
/// liaisons ne coïncident qu'au premier ordre.
#[derive(Clone, Debug)]
pub struct Cardan {
    pub nom: String,
    pub a: Ref,
    pub b: Ref,
    pub na: V3, // axe de la croix côté a, dans le repère de a
    pub nb: V3, // axe de la croix côté b, dans le repère de b
}

impl Cardan {
    pub fn phi_g(&self, corps: &[Corps]) -> (DVector<f64>, DMatrix<f64>) {
        let (_, ra) = pose_de(corps, self.a);
        let (_, rb) = pose_de(corps, self.b);
        let ua = ra * self.na;
        let ub = rb * self.nb;
        // δΦ = δθ_a·(u_a × u_b) − δθ_b·(u_a × u_b)
        let c = ua.cross(&ub);
        let mut g = DMatrix::zeros(1, 12);
        for k in 0..3 {
            g[(0, 3 + k)] = c[k];
            g[(0, 9 + k)] = -c[k];
        }
        (DVector::from_element(1, ua.dot(&ub)), g)
    }
}

/// Loi d'un couple autour d'un axe, fonction de l'angle relatif θ (déroulé),
/// de la vitesse relative ω et du temps.
#[derive(Clone, Debug)]
pub enum LoiCouple {
    Constant(f64),
    /// −k·(θ − θ₀) − c·ω
    Ressort {
        k: f64,
        c: f64,
        theta0: f64,
    },
    /// servo PD saturé en couple ET en vitesse (le patron du n0 FRELON) :
    /// u = tanh((Kp·(cible(t) − θ) − Kd·ω)/Q), τ = max(u,0)·Q·clamp(1 − ω/ω_nl) + min(u,0)·Q·clamp(1 + ω/ω_nl)
    PD {
        kp: f64,
        kd: f64,
        q_max: f64,
        w_nl: f64,
        cible: Loi,
    },
    /// gouverneur : τ = sat(Kg·(Ω_cible(t) − ω), ±Q)
    Gouverneur {
        kg: f64,
        q_max: f64,
        cible: Loi,
    },
    /// BUTÉE UNILATÉRALE d'articulation — la limite mécanique que tout
    /// mécanisme réel a et qu'aucun code ne peut ignorer : un vérin en fin de
    /// course, un palonnier contre son arrêt, un genou. Elle n'agit QUE
    /// au-delà de la plage [θ_min, θ_max] :
    ///
    /// ```text
    ///     τ = −k·(θ − θ_max) − c·ω  si θ > θ_max ,  symétrique en θ_min ,  0 entre
    /// ```
    ///
    /// L'amortissement suit la PÉNÉTRATION (patron de Hunt–Crossley, comme le
    /// contact) : sans ça un amortisseur linéaire TIRE au décollement et
    /// arrache la pièce de sa butée. C'est la même faute que celle payée sur
    /// le contact, et c'est pourquoi elle est écrite ici.
    Butee {
        k: f64,
        c: f64,
        theta_min: f64,
        theta_max: f64,
    },
}

/// Couple entre deux corps autour d'un axe fixé dans le repère de `a` (monde
/// si a = bâti) : +τ·n sur b, −τ·n sur a. L'angle relatif est mesuré depuis
/// la pose initiale et DÉROULÉ (comme l'engrenage).
#[derive(Clone, Debug)]
pub struct Couple {
    pub nom: String,
    pub a: Ref,
    pub b: Ref,
    pub axe: V3,
    pub ref_rel: M3,
    pub prev: f64,
    pub loi: LoiCouple,
}

/// Polaire de section : Cl, Cd, Cm sur une grille d'incidence (deg, −180…180),
/// interpolation linéaire — la table c81 mono-Mach de FRELON, passée en tableaux.
#[derive(Clone, Debug)]
pub struct Polaire {
    pub alpha: Vec<f64>,
    pub cl: Vec<f64>,
    pub cd: Vec<f64>,
    pub cm: Vec<f64>,
}

/// Inflow uniforme de rotor (théorie du disque). En AVANCEMENT c'est la
/// relation de Glauert qui vaut, pas celle du stationnaire :
///     v_i = T / (2 ρ A · U),   U = ‖V_disque/air + v_i·axe‖
/// résolue par le point fixe amorti à ½ — qui EST l'itération babylonienne en
/// stationnaire (U = v_i ⇒ v ← ½(v + v_h²/v), convergence quadratique vers
/// √(T/2ρA)) et une contraction franche dès qu'il y a de la vitesse d'avance.
/// Le résultat est suivi avec un retard du premier ordre τ — un ÉTAT
/// explicite, mis à jour à chaque pas accepté depuis la poussée du pas (pas de
/// couplage dans le jacobien).
#[derive(Clone, Debug)]
pub struct Inflow {
    pub axe: V3,
    pub aire: f64,
    pub tau: f64,
    pub rho: f64,
    pub v_i: f64,
    pub poussee: f64,
    /// INFLOW DYNAMIQUE DE PITT–PETERS (1981), trois états — le standard des
    /// codes de rotor (CAMRAD, HOST, RCAS). `false` = le modèle uniforme
    /// quasi-statique d'origine, qui reste le défaut.
    ///
    /// La distribution devient  λ(r̄, ψ) = λ₀ + r̄·(λ₁c cos ψ + λ₁s sin ψ),
    /// et les trois états obéissent à
    ///
    /// ```text
    ///     [M] dλ/dτ + [L]⁻¹ λ = (C_T, −C_L, C_M)ᵀ ,   τ = Ω t
    ///     M = diag(8/3π, 16/45π, 16/45π)
    /// ```
    ///
    /// Ce n'est pas un raffinement cosmétique : c'est le seul modèle qui donne
    /// au rotor sa VRAIE constante de temps de réponse (et donc sa dynamique
    /// de commande), et le gradient de sillage en avancement au lieu d'une
    /// valeur moyenne.
    pub pitt: bool,
    /// INFLOW IMPOSÉ DE L'EXTÉRIEUR (sillage libre en Python, table) : ni le
    /// point fixe uniforme ni Pitt–Peters ne le font évoluer — `actualise` le
    /// laisse tel quel. Les harmoniques v1c/v1s sont lues avec `pitt` = true.
    pub impose: bool,
    /// centre du disque (monde, supposé fixe), référence d'azimut dans le plan,
    /// et vitesse de bout ΩR — nécessaires pour situer une station en (r̄, ψ)
    pub centre: V3,
    pub e1: V3,
    pub vtip: f64,
    /// harmoniques 1/rev de l'inflow, en m/s comme `v_i`
    pub v1s: f64,
    pub v1c: f64,
    /// PROFIL RADIAL d'inflow imposé : (r̄, v) croissant en r̄, interpolé
    /// linéairement, AJOUTÉ à v_i + harmoniques à la station. C'est ce qu'un
    /// sillage libre en stationnaire rend (`vinkulum.sillage`) : la perte de
    /// bout et la distribution que le disque uniforme n'a pas. Vide = rien.
    pub profil: Vec<(f64, f64)>,
    /// CARTE d'inflow imposée w(r̄, ψ) (m/s, > 0 vers le bas), grille r̄ × ψ
    /// (ψ en rad depuis e1, périodique), interpolation bilinéaire, AJOUTÉE à
    /// la station comme le profil. C'est ce qu'un modèle à états finis
    /// (Peters–He, `vinkulum.peters_he`) ou un sillage en avancement rend :
    /// toutes les harmoniques, pas seulement la 1/rev. Vide = rien.
    pub carte: (Vec<f64>, Vec<f64>, Vec<f64>),
    /// moments de rouli/tangage aérodynamiques du dernier pas, autour du centre
    pub m_roul: f64,
    pub m_tang: f64,
}

/// LEISHMAN–BEDDOES (1989), constantes de temps en SEMI-CORDES : retard de
/// pression T_p (le point de séparation répond à une incidence EFFECTIVE
/// retardée), retard de séparation T_f, décroissance du tourbillon T_v et
/// temps de convection du tourbillon sur la corde T_vl. Valeurs de la
/// littérature (profil NACA 0012 en subsonique) — jamais calibrées ici, et
/// c'est déclaré : elles se mesurent en soufflerie sur LA section.
pub const LB_TP: f64 = 1.7;
pub const LB_TF: f64 = 3.0;
pub const LB_TV: f64 = 6.0;
pub const LB_TVL: f64 = 11.0;

/// États de Leishman–Beddoes d'une station : déficit de pression D_p, dernier
/// C_N potentiel, point de séparation retardé f'', temps de tourbillon τ_v,
/// dernier C_v et portance tourbillonnaire C_N^v.
#[derive(Clone, Copy, Debug, Default)]
pub struct EtatLb {
    pub dp: f64,
    pub cn_prev: f64,
    pub ff: f64,
    pub tau_v: f64,
    pub cv_prev: f64,
    pub cnv: f64,
    /// incidence du pas précédent, pour détecter le DÉBUT d'une montée
    pub al_prev: f64,
    pub descend: bool,
}

/// Constantes de R. T. Jones pour la fonction de Wagner
/// φ(s) = 1 − A₁e^{−b₁s} − A₂e^{−b₂s}. A₁ + A₂ = 0,5 EXACTEMENT : c'est ce qui
/// donne φ(0) = 0,5, la moitié de portance à l'instant d'un échelon.
pub const W_A1: f64 = 0.165;
pub const W_A2: f64 = 0.335;
pub const W_B1: f64 = 0.0455;
pub const W_B2: f64 = 0.300;

/// UN PAS D'INFLOW DYNAMIQUE DE PITT–PETERS (1981).
///
/// Trois états — uniforme et les deux harmoniques 1/rev — au lieu d'une
/// moyenne. C'est le modèle standard des codes de rotor, et ce qu'il apporte
/// n'est pas un raffinement : c'est la CONSTANTE DE TEMPS du rotor. Le modèle
/// uniforme d'origine la posait à la main (`tau`) ; ici elle SORT de la
/// masse apparente du fluide, qui est une donnée de la théorie.
///
/// ```text
///     [M] dλ/dτ + [L]⁻¹ λ = (C_T, −C_L, C_M)ᵀ ,   τ = Ω t
///     M = diag(8/3π, 16/45π, 16/45π)
/// ```
///
/// avec, χ étant l'inclinaison du sillage (0 au stationnaire, π/2 en vol
/// rapide) et V_T, V_m les vitesses de masse :
///
/// ```text
///     L = [ 1/(2V_T)              0                (15π/64)tan(χ/2)/V_T ]
///         [ 0                    −4/(V_m(1+cos χ))  0                   ]
///         [ (15π/64)tan(χ/2)/V_T  0                −4cos χ/(V_m(1+cos χ))]
/// ```
///
/// TROIS CONTRÔLES ANALYTIQUES le gardent (`bancs.pitt_peters`) :
///  1. au stationnaire l'équilibre doit retomber EXACTEMENT sur Froude,
///     λ₀ = √(C_T/2) — c'est L₁₁ = 1/(2λ₀) qui le fait, sans qu'on l'écrive ;
///  2. la constante de temps linéarisée vaut τ_λ = (8/3π)/(4λ₀) en temps
///     adimensionné, soit 0,2122/λ₀ — une valeur qu'aucune ligne ne pose ;
///  3. en vol rapide λ₀ → C_T/(2μ), la formule de Glauert.
///
/// L'intégration est EXPLICITE et séparée du multicorps, comme dans les codes
/// du domaine : l'inflow est un état de fluide, pas un degré de liberté
/// mécanique, et le coupler implicitement demanderait de dériver la table c81.
fn pitt_peters(inf: &mut Inflow, vent: V3, h: f64) {
    if !(h > 0.0) || inf.vtip <= 0.0 {
        return;
    }
    let r = (inf.aire / std::f64::consts::PI).sqrt();
    let om = inf.vtip / r; // Ω
    let q = inf.rho * inf.aire * inf.vtip * inf.vtip;
    // coefficients adimensionnés ; le moment est rapporté à R de plus
    let ct = inf.poussee / q;
    let cl = inf.m_roul / (q * r);
    let cm = inf.m_tang / (q * r);
    // λ (axial, positif vers le bas du disque) et μ (avancement)
    let vz = vent.dot(&inf.axe);
    let mu = (vent - inf.axe * vz).norm() / inf.vtip;
    let l0 = inf.v_i / inf.vtip;
    let lz = -vz / inf.vtip; // composante axiale du vent, même signe que λ₀
    let lam = lz + l0; // inflow total vu par le disque
                       // vitesses de masse (Peters–HaQuang) : V_T totale, V_m « moyenne »
    let vt = (mu * mu + lam * lam).sqrt().max(1e-6);
    let vm = ((mu * mu + lam * (lam + l0)) / vt).max(1e-6);
    // χ : 0 au stationnaire (sillage sous le disque), π/2 en vol rasant
    let chi = mu.atan2(lam.abs().max(1e-9));
    let (cc, tc) = (chi.cos(), (0.5 * chi).tan());
    let k = 15.0 * std::f64::consts::PI / 64.0 * tc;
    // L, puis λ̇ = M⁻¹(C − L⁻¹λ). Le 2×2 (0,2) s'inverse à la main.
    // SIGNES : forme de Peters & HaQuang (1988). L₂₂ et L₃₃ sont POSITIFS,
    // seul le couplage L₃₁ porte le signe — c'est ce qui rend −L⁻¹λ
    // stabilisant. Essayé avec L₂₂, L₃₃ négatifs : det(L) < 0, donc
    // (L⁻¹)₃₃ < 0, donc λ₁c croît exponentiellement — mesuré, λ₀ à 1e48 en
    // vol d'avancement. Contrôle qui le tranche : à χ = 0 (stationnaire) le
    // couplage s'annule et les trois termes diagonaux doivent être > 0.
    let (l11, l13, l31, l33) = (0.5 / vt, k / vt, -k / vt, 4.0 * cc / (vm * (1.0 + cc)));
    let l22 = 4.0 / (vm * (1.0 + cc));
    let det = l11 * l33 - l13 * l31;
    let (i11, i13, i31, i33) = if det.abs() > 1e-12 {
        (l33 / det, -l13 / det, -l31 / det, l11 / det)
    } else {
        (1.0 / l11, 0.0, 0.0, 0.0)
    };
    let i22 = if l22.abs() > 1e-12 { 1.0 / l22 } else { 0.0 };
    let (m0, m1) = (
        8.0 / (3.0 * std::f64::consts::PI),
        16.0 / (45.0 * std::f64::consts::PI),
    );
    let (a0, a1s, a1c) = (l0, inf.v1s / inf.vtip, inf.v1c / inf.vtip);
    let d0 = (ct - (i11 * a0 + i13 * a1c)) / m0;
    let d1s = (-cl - i22 * a1s) / m1;
    let d1c = (cm - (i31 * a0 + i33 * a1c)) / m1;
    let dtau = om * h; // τ = Ω t
    inf.v_i = (a0 + dtau * d0).max(0.0) * inf.vtip;
    inf.v1s = (a1s + dtau * d1s) * inf.vtip;
    inf.v1c = (a1c + dtau * d1c) * inf.vtip;
}

/// Pale en théorie des tranches quasi-stationnaire : stations de Gauss le long
/// de l'envergure (axe `es` du repère de section, dans le corps), corde `c`,
/// polaire c81, inflow optionnel. Repère de section (repère corps) : `ec` =
/// direction d'AVANCE de la section (bord d'attaque devant), `en` = normale
/// (portance positive), `es` = envergure = ec × en.
#[derive(Clone, Debug)]
pub struct Pale {
    pub nom: String,
    pub corps: usize,
    pub p0: V3, // pied d'envergure, repère corps (au CdM)
    pub es: V3,
    pub ec: V3,
    pub en: V3,
    pub longueur: f64,
    pub corde: f64,
    pub polaire: std::sync::Arc<Polaire>,
    pub inflow: Option<usize>,
    pub rho: f64,
    pub gauss: usize,
    /// dernières sorties : (poussée le long de l'axe d'inflow, couple autour de cet axe)
    pub sortie: (f64, f64),
    /// AÉRODYNAMIQUE INSTATIONNAIRE DE SECTION (Wagner / Theodorsen, forme
    /// indicielle de R. T. Jones). `false` = quasi-stationnaire, le défaut.
    ///
    /// Une section de pale de rotor voit son incidence varier à 1/rev en
    /// avancement. En quasi-stationnaire elle répond instantanément ; en
    /// réalité le sillage déposé derrière elle retarde et atténue la portance.
    /// La fonction de Wagner dit de combien :
    ///
    /// ```text
    ///     φ(s) = 1 − A₁e^{−b₁s} − A₂e^{−b₂s} ,  s = 2Ut/c  (semi-cordes)
    ///     A₁ = 0,165  b₁ = 0,0455   A₂ = 0,335  b₂ = 0,300
    /// ```
    ///
    /// φ(0) = 1 − A₁ − A₂ = **0,5** : à l'instant d'un échelon d'incidence, une
    /// section ne porte que la MOITIÉ de sa valeur stationnaire. C'est le
    /// résultat classique, et c'est le contrôle qui garde le modèle.
    ///
    /// Mise en état (deux états par station) : ż = −b z + α en temps
    /// semi-corde, puis α_eff = α(1−A₁−A₂) + b₁A₁z₁ + b₂A₂z₂. C'est α_eff qui
    /// alimente la table c81 — donc le retard porte sur la partie linéaire et
    /// la table garde le décrochage, ce qui est l'usage du domaine.
    pub instat: bool,
    /// deux états de Wagner par station de Gauss
    pub z: Vec<(f64, f64)>,
    /// DÉCROCHAGE DYNAMIQUE d'Øye (un état par station) : la portance est
    /// f·CL_att + (1 − f)·CL_fs, CL_att la droite attachée a₀(α − α₀), CL_fs
    /// la table « entièrement décollée » déduite de la polaire statique par
    /// f_st = (2√(CL/CL_att) − 1)², et f relaxe vers f_st(α) avec τ = 4c/|V|.
    /// Le c81 statique reste la limite t → ∞ par construction. Pas de
    /// Leishman–Beddoes : ni retard de pression, ni tourbillon de bord
    /// d'attaque — déclaré.
    pub oye: bool,
    pub a0: f64,
    pub alpha0: f64,
    pub cl_fs: Vec<f64>,
    pub fst: Vec<f64>,
    pub f: Vec<f64>,
    /// DÉCROCHAGE DYNAMIQUE DE LEISHMAN–BEDDOES (`lb=True`) — ce qu'Øye n'a
    /// pas, en deux étages : (1) un RETARD DE PRESSION : le point de séparation
    /// ne répond pas à α mais à une incidence effective α_f tirée d'un C_N
    /// retardé de T_p, puis f'' relaxe vers f_st(α_f) en T_f ; (2) un
    /// TOURBILLON DE BORD D'ATTAQUE : quand le C_N retardé dépasse le seuil
    /// critique C_N1 (lu sur la polaire : f_st = 0,7), la portance que la
    /// séparation retire, C_v = C_N·(1 − K_N), K_N = ((1+√f'')/2)², est
    /// réinjectée comme portance tourbillonnaire C_N^v qui décroît en T_v
    /// pendant que le tourbillon convecte (τ_v < T_vl), puis s'éteint.
    /// La portance totale garde la FORME d'Øye — f''·CL_att + (1−f'')·CL_fs —
    /// qui coïncide avec le Kirchhoff de LB à f = f_st et rend le c81 statique
    /// EXACTEMENT à t → ∞ ; C_N^v s'y ajoute. CD et CM restent statiques
    /// (LB les retarde aussi : non modélisé, déclaré). Les états sont figés
    /// sur le pas (état de fluide), avancés en exponentielle exacte.
    pub lb: bool,
    pub lb_cn1: f64,
    pub etats_lb: Vec<EtatLb>,
    /// constantes de LB de CETTE pale (T_p, T_f, T_v, T_vl) en semi-cordes —
    /// par défaut celles de la littérature ; réglables pour MESURER la
    /// sensibilité d'un écart à un essai (S809 à 20°, 5 sept.), pas pour caler
    pub lb_const: (f64, f64, f64, f64),
    /// MASSE AJOUTÉE (`masse_ajoutee=True`) — la portance NON circulatoire de
    /// Theodorsen, L_nc = πρb²·ẇ, ẇ la dérivée de la vitesse normale de l'air
    /// vue au MILIEU de la corde (= ḧ + Vα̇ − b·a·α̈ pour un pivot à a).
    /// Absente, une section en tangage rapide perd la moitié de sa portance à
    /// k 0,76 (Kim 2017, 6 sept.). ẇ par différence arrière sur le pas —
    /// ordre 1, implicite en vitesse : l'AD porte πρb²/h dans la tangente.
    /// Pas de moment non circulatoire (déclaré). Faux par défaut.
    pub nc: bool,
    pub nc_w: Vec<f64>,
    pub nc_h: f64,
    /// INCIDENCE AU 3/4 DE CORDE (`alpha_34=True`) : la vitesse normale qui
    /// fait la portance circulatoire est celle du point à 3/4 de corde, pas
    /// celle de la ligne de référence (c/4) — c'est le b(½−a)α̇ de Theodorsen.
    /// Sans lui une section qui TANGUE ne voit pas son taux de tangage. Faux
    /// par défaut : le vol S2 et les bancs restent identiques au bit près.
    pub a34: bool,
    /// RETARD LINÉAIRE (`retard_lineaire=True`, avec `instationnaire`) : Wagner
    /// ne retarde que la partie LINÉAIRE a₀(α − α₀) ; la non-linéarité statique
    /// Δ(α) = c81(α) − a₀(α − α₀) est lue à l'incidence INSTANTANÉE. Lire toute
    /// la table à l'incidence retardée AMPLIFIE une bulle laminaire de bas Re
    /// au lieu de la lisser (Kim 2017 : C_L max −78 %). a₀ vient de la table
    /// (moindres carrés ±5°) ou de `a0_rad` quand la table n'a pas de partie
    /// linéaire — c'est le cas d'une polaire à bulle.
    pub retard_lineaire: bool,
}

/// CONTACT par PÉNALITÉ régularisée — une sphère portée par un corps contre un
/// demi-espace, loi de Hunt–Crossley plus frottement de Coulomb lissé.
///
/// ```text
///     δ = R − (n·(p − p₀))            enfoncement, nul si négatif
///     F_n = k δ^e (1 + c δ̇)          Hunt–Crossley : l'amortissement est
///                                      PROPORTIONNEL à δ, donc la force part
///                                      et revient à zéro continûment — un
///                                      amortisseur linéaire c·δ̇ colle à
///                                      l'impact et TIRE au décollement
///     F_t = −μ F_n · tanh(‖v_t‖/v_ε) · v_t/‖v_t‖
/// ```
///
/// Le choix de la pénalité plutôt que d'un multiplicateur est celui de tous
/// les codes de contact : une contrainte stricte rend le résidu discontinu à
/// l'instant du choc et fait diverger Newton. Ici tout reste C¹, donc
/// différentiable — et l'AD du noyau passe au travers.
#[derive(Clone, Debug)]
pub struct Contact {
    pub nom: String,
    pub corps: usize,
    /// premier point de l'AXE, repère corps. Une CAPSULE est un segment
    /// [p0, p1] dilaté de `rayon` ; la SPHÈRE en est le cas dégénéré p1 = p0.
    /// C'est la primitive qui couvre le plus de mécanismes réels sans
    /// maillage : tiges, barres, doigts, axes, biellettes.
    pub p0: V3,
    pub p1: V3,
    pub rayon: f64,
    /// SECOND corps porteur d'une sphère ou capsule (None = demi-espace fixe)
    pub b: Option<usize>,
    pub pb: V3,
    pub p1b: V3,
    pub rayon_b: f64,
    /// demi-espace admissible : n·(p − orig) ≥ 0 (si `b` est None)
    pub normale: V3,
    pub origine: V3,
    pub k: f64,
    /// exposant de la loi de Hertz ; `expo < 0` bascule sur la BARRIÈRE IPC
    /// (Li & al. 2020) : B(d) = −(d − d̂)² ln(d/d̂) sur 0 < d < d̂, nulle
    /// au-delà. `d_hat` est alors le seuil d'activation, `k` le facteur κ.
    pub expo: f64,
    pub d_hat: f64,
    pub c: f64,
    pub mu: f64,
    pub v_eps: f64,
    /// BOÎTE (OBB) portée par le second corps : demi-dimensions dans son
    /// repère, centrée en `pb`. Avec la sphère et la capsule, c'est la
    /// troisième primitive qui couvre l'essentiel d'un mécanisme sans
    /// maillage — un bâti, une glissière, une came plane, un carter.
    ///
    /// Le point le plus proche est un CLAMP dans le repère de la boîte, donc
    /// exact et différentiable ; les trois régimes (face, arête, coin) sont
    /// les trois façons dont le clamp mord, et ils se vérifient séparément.
    ///
    /// ⚠ NON COUVERT, déclaré : le centre de la primitive À L'INTÉRIEUR de la
    /// boîte. La distance y change de nature (il faut sortir par la face la
    /// plus proche) et une pénalité ne doit de toute façon jamais y arriver —
    /// si elle y arrive, la raideur est trop faible, et c'est ça qu'il faut
    /// corriger. Le noyau rend alors une distance nulle plutôt qu'un signe
    /// faux : mieux vaut ne pas repousser que repousser dans le mauvais sens.
    pub demi: Option<V3>,
    /// CYLINDRE FINI porté par le second corps : (axe unitaire dans son repère,
    /// rayon, demi-longueur), centré en `pb`. La primitive qui manquait à la
    /// liste — un arbre, un galet, un rouleau, un tourillon — et que ni la
    /// capsule (bouts ronds) ni la boîte ne remplacent : le point le plus
    /// proche a QUATRE régimes (flanc, fond, arête circulaire, loin) qui se
    /// vérifient séparément.
    ///
    /// ⚠ NON COUVERT, déclaré, comme pour la boîte : le centre de la primitive
    /// À L'INTÉRIEUR du cylindre — la distance y change de nature et une
    /// pénalité ne doit jamais y arriver. Le noyau rend alors une distance
    /// nulle plutôt qu'un signe faux.
    pub cylindre: Option<(V3, f64, f64)>,
    /// MAILLAGE du second corps : géométrie quelconque, index dans
    /// `Modele.maillages`. Prioritaire sur boîte, capsule et sphère.
    pub maille: Option<usize>,
    /// point le plus proche sur le maillage, repère du corps porteur. Figé sur
    /// le PAS : la recherche BVH n'est pas différentiable, mais la distance à
    /// un point fixe l'est. C'est le procédé standard du contact sur maillage
    /// — et il est licite tant que le pas ne fait pas changer de triangle,
    /// ce que le pas piloté par le contact garantit déjà près du choc.
    pub q_maille: V3,
    /// CAPSULE contre boîte, cylindre ou maillage : abscisse s ∈ [0, 1] du point
    /// de l'axe [p0, p1] le plus proche de la primitive, FIGÉE sur le pas comme
    /// `q_maille` — le résidu voit la sphère de l'axe qui s'y trouve. Exact aux
    /// bouts, d'ordre 2 à l'intérieur (la distance est stationnaire en s au
    /// minimum). 0 pour une sphère, et pour toute autre paire.
    pub s_axe: f64,
    /// CÂBLE : le contact INVERSÉ. Un contact repousse quand deux corps se
    /// rapprochent ; un câble tire quand ils s'éloignent — c'est le même
    /// élément au signe près, et il donne aussi le ressort unilatéral et la
    /// butée de translation. Le « rayon » devient la LONGUEUR au repos.
    pub cable: bool,
    /// dernier (enfoncement, force normale) — diagnostic
    pub sortie: (f64, f64),
    /// CONTACT NON LISSE : la force normale n'est plus une pénalité mais un
    /// MULTIPLICATEUR λ ≥ 0 sous complémentarité gap ≥ 0, λ·gap = 0, écrite
    /// avec la fonction de Fischer–Burmeister φ(a, b) = a + b − √(a² + b²)
    /// (semi-lisse : le même Newton la résout). L'élément `Elem::K` porte la
    /// ligne ; `row` est son indice dans λ. Le frottement reste le Coulomb
    /// lissé, appliqué sur λ (`fn_impose`, posé sur une COPIE le temps d'une
    /// évaluation).
    ///
    /// LA LIGNE EST AU NIVEAU VITESSE (Moreau–Jean). Une contrainte de
    /// POSITION sous Newmark rebondit d'elle-même : q₁ bloqué ⇒ a₁ impulsif ⇒
    /// u₁ ≈ u₀(1 − γ/β) ≈ −u₀, restitution parasite ≈ 1 — mesuré, la bille
    /// remontait plus haut que son lâcher. On impose donc, pour un contact
    /// ACTIF (gap prédit gap₀ + h·ġap₀ ≤ 0, figé sur le pas) :
    ///     ġap₁ + e·ġap₀ + gap₀/h = 0   ⊥   λ ≥ 0
    /// loi de Newton en vitesse, plus un rappel du gap qui converge en β/γ par
    /// pas (pas de dérive de pénétration). Un contact inactif porte λ = 0.
    pub nonlisse: bool,
    pub row: usize,
    pub fn_impose: Option<f64>,
    /// coefficient de restitution de Newton (0 = choc mou)
    pub restitution: f64,
    /// contact créé par l'APPARIEMENT automatique (`apparie`) : remplacé à
    /// chaque passe, toujours rangé APRÈS les contacts déclarés à la main —
    /// dont les indices (`Elem::K::ct`) doivent rester stables
    pub auto: bool,
}

/// Poutre à deux nœuds, cinématique de type Simo–Reissner.
/// Le mode historique utilise les déformations au milieu de l'élément.
/// Dans la formulation intégrée optionnelle, la rotation
/// matérielle est interpolée sur SO(3) et intégrée exactement pour relier
/// la corde à la déformation généralisée γ. κ = log(R_Aᵀ R_B)/L.
/// L'énergie est U = ½ L (γᵀ C_eff γ + κᵀ C_M κ), avec la flexibilité
/// transverse 1/C_eff = 1/GA + L²/(12 EI) issue de la condensation linéaire.
/// Cette limite retrouve la raideur nodale de Timoshenko ; la cinématique
/// reproduit l'arc de flexion pure. Ce n'est pas une solution exacte de
/// toute poutre non linéaire. Les nœuds portent les masses concentrées.
/// Forces = −∇U par duaux ; tangente = dérivée du champ (cf. `tangent`).
#[derive(Clone, Debug)]
pub struct Poutre {
    pub nom: String,
    /// Interpolation intégrée avec flexibilité condensée (option explicite).
    pub integree: bool,
    pub a: usize,
    pub b: usize,
    pub l: f64,
    /// (EA, GA_y, GA_z)
    pub cn: [f64; 3],
    /// (GJ, EI_y, EI_z)
    pub cm: [f64; 3],
    /// repère matériel au repos, vu de chaque nœud : R_mat = R_A·ra = R_B·rb
    pub ra: M3,
    pub rb: M3,
}

/// MAILLAGE TRIANGULAIRE porté par un corps — la GÉOMÉTRIE QUELCONQUE, qui
/// manquait au contact. Sphères, capsules et boîtes couvrent beaucoup ; elles
/// ne couvrent pas une pièce usinée, un carter, une came.
///
/// Le maillage vit dans le repère du corps et ne bouge jamais : c'est le
/// corps qui tourne. Le BVH se construit donc UNE FOIS, à la déclaration, et
/// reste valide pour toute la simulation — c'est ce qui rend la requête
/// logarithmique sans jamais rien reconstruire.
#[derive(Clone, Debug)]
pub struct Maillage {
    pub nom: String,
    pub corps: usize,
    pub sommets: Vec<V3>,
    pub tris: Vec<[usize; 3]>,
    /// arbre d'AABB : (min, max, premier, nb, gauche, droite) — feuille si
    /// `gauche` est usize::MAX
    pub bvh: Vec<(V3, V3, usize, usize, usize, usize)>,
    /// index des triangles, réordonné par la construction du BVH
    pub ordre: Vec<usize>,
}

impl Maillage {
    /// Construit l'arbre d'AABB par médiane sur l'axe le plus étalé.
    fn bati(
        sommets: &[V3],
        tris: &[[usize; 3]],
    ) -> (Vec<(V3, V3, usize, usize, usize, usize)>, Vec<usize>) {
        let n = tris.len();
        let mut ordre: Vec<usize> = (0..n).collect();
        let centre = |t: &[usize; 3]| (sommets[t[0]] + sommets[t[1]] + sommets[t[2]]) / 3.0;
        let mut noeuds: Vec<(V3, V3, usize, usize, usize, usize)> = Vec::new();
        // pile explicite : un arbre profond en récursion coûterait la pile
        let mut pile = vec![(0usize, n, usize::MAX, 0usize)];
        while let Some((deb, fin, parent, cote)) = pile.pop() {
            let mut lo = V3::repeat(f64::INFINITY);
            let mut hi = V3::repeat(f64::NEG_INFINITY);
            for &ti in &ordre[deb..fin] {
                for &si in &tris[ti] {
                    lo = lo.inf(&sommets[si]);
                    hi = hi.sup(&sommets[si]);
                }
            }
            let moi = noeuds.len();
            noeuds.push((lo, hi, deb, fin - deb, usize::MAX, usize::MAX));
            if parent != usize::MAX {
                if cote == 0 {
                    noeuds[parent].4 = moi;
                } else {
                    noeuds[parent].5 = moi;
                }
            }
            if fin - deb <= 4 {
                continue;
            } // feuille
            let d = hi - lo;
            let ax = if d.x >= d.y && d.x >= d.z {
                0
            } else if d.y >= d.z {
                1
            } else {
                2
            };
            let mid = (deb + fin) / 2;
            ordre[deb..fin].sort_by(|&a, &b| {
                centre(&tris[a])[ax]
                    .partial_cmp(&centre(&tris[b])[ax])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            pile.push((mid, fin, moi, 1));
            pile.push((deb, mid, moi, 0));
        }
        (noeuds, ordre)
    }

    /// Distance d'un point (repère CORPS) au maillage, et le point le plus
    /// proche. Parcours du BVH avec élagage par la meilleure distance courante
    /// — c'est l'élagage qui fait le log, pas l'arbre seul.
    pub fn plus_proche(&self, p: &V3) -> (f64, V3) {
        let mut best = (f64::INFINITY, *p);
        if self.bvh.is_empty() {
            return best;
        }
        let mut pile = vec![0usize];
        while let Some(i) = pile.pop() {
            let (lo, hi, deb, nb, g, d) = self.bvh[i];
            // distance du point à l'AABB : si elle dépasse la meilleure, tout
            // le sous-arbre est hors course
            let q = V3::new(
                p.x.clamp(lo.x, hi.x),
                p.y.clamp(lo.y, hi.y),
                p.z.clamp(lo.z, hi.z),
            );
            if (q - p).norm() >= best.0 {
                continue;
            }
            if g == usize::MAX {
                for &ti in &self.ordre[deb..deb + nb] {
                    let t = self.tris[ti];
                    let c = point_triangle(
                        p,
                        &self.sommets[t[0]],
                        &self.sommets[t[1]],
                        &self.sommets[t[2]],
                    );
                    let dd = (c - p).norm();
                    if dd < best.0 {
                        best = (dd, c);
                    }
                }
            } else {
                // ORDRE DE DESCENTE : on visite l'enfant le PLUS PROCHE en
                // premier, donc on l'empile en DERNIER. Sans ça `best` reste
                // grand longtemps et l'élagage ne coupe presque rien —
                // mesuré : exposant 0,60 en N au lieu de 0,26.
                let dd = |j: usize| -> f64 {
                    if j == usize::MAX {
                        return f64::INFINITY;
                    }
                    let (l2, h2, ..) = self.bvh[j];
                    let q2 = V3::new(
                        p.x.clamp(l2.x, h2.x),
                        p.y.clamp(l2.y, h2.y),
                        p.z.clamp(l2.z, h2.z),
                    );
                    (q2 - p).norm()
                };
                let (dg, dr) = (dd(g), dd(d));
                if dg <= dr {
                    if d != usize::MAX {
                        pile.push(d);
                    }
                    pile.push(g);
                } else {
                    pile.push(g);
                    pile.push(d);
                }
            }
        }
        best
    }
}

impl Maillage {
    /// Distance d'un SEGMENT [a, b] (repère corps) au maillage, en forme
    /// fermée triangle par triangle (`seg_tri`) : rend (distance, s sur [a, b],
    /// point du maillage). Le BVH est élagué par le test segment/AABB dilatée
    /// de la meilleure distance courante (conservatif : le cube contient la
    /// boule). C'est la requête du contact CAPSULE contre maillage — et du
    /// CCD, où la sphère balayée sur le pas est une capsule.
    pub fn plus_proche_segment(&self, a: &V3, b: &V3) -> (f64, f64, V3) {
        let mut best = (f64::INFINITY, 0.0, *a);
        if self.bvh.is_empty() {
            return best;
        }
        let d = b - a;
        let mil = 0.5 * (a + b);
        let mut pile = vec![0usize];
        while let Some(i) = pile.pop() {
            let (lo, hi, deb, nb, g, dr) = self.bvh[i];
            // slabs : le segment approche-t-il l'AABB à moins de `best` ?
            let r = best.0;
            let (mut t0, mut t1) = (0.0f64, 1.0f64);
            let mut dehors = false;
            for k in 0..3 {
                let (l, h) = (lo[k] - r, hi[k] + r);
                if d[k].abs() < 1e-300 {
                    if a[k] < l || a[k] > h {
                        dehors = true;
                    }
                    continue;
                }
                let (mut u0, mut u1) = ((l - a[k]) / d[k], (h - a[k]) / d[k]);
                if u0 > u1 {
                    std::mem::swap(&mut u0, &mut u1);
                }
                t0 = t0.max(u0);
                t1 = t1.min(u1);
            }
            if dehors || t0 > t1 {
                continue;
            }
            if g == usize::MAX {
                for &ti in &self.ordre[deb..deb + nb] {
                    let t = self.tris[ti];
                    let (dd, sx, q) = seg_tri(
                        a,
                        b,
                        &self.sommets[t[0]],
                        &self.sommets[t[1]],
                        &self.sommets[t[2]],
                    );
                    if dd < best.0 {
                        best = (dd, sx, q);
                    }
                }
            } else {
                // l'enfant le plus proche du milieu du segment d'abord (empilé
                // en dernier), pour que `best` chute vite et que l'élagage coupe
                let dd = |j: usize| -> f64 {
                    if j == usize::MAX {
                        return f64::INFINITY;
                    }
                    let (l2, h2, ..) = self.bvh[j];
                    let q2 = V3::new(
                        mil.x.clamp(l2.x, h2.x),
                        mil.y.clamp(l2.y, h2.y),
                        mil.z.clamp(l2.z, h2.z),
                    );
                    (q2 - mil).norm()
                };
                if dd(g) <= dd(dr) {
                    if dr != usize::MAX {
                        pile.push(dr);
                    }
                    pile.push(g);
                } else {
                    pile.push(g);
                    pile.push(dr);
                }
            }
        }
        best
    }

    /// CCD : la CAPSULE [a, b] de rayon `r` touche-t-elle le maillage ?
    pub fn traverse(&self, a: &V3, b: &V3, r: f64) -> bool {
        self.plus_proche_segment(a, b).0 <= r
    }
}

/// Points les plus proches de deux SEGMENTS [p1, q1] et [p2, q2] — Ericson
/// 5.1.9, forme fermée : (s, t) clampés dans [0, 1], segments dégénérés et
/// parallèles compris (repli s = 0).
pub fn seg_seg(p1: &V3, q1: &V3, p2: &V3, q2: &V3) -> (f64, f64) {
    let d1 = q1 - p1;
    let d2 = q2 - p2;
    let r = p1 - p2;
    let a = d1.dot(&d1);
    let e = d2.dot(&d2);
    let f = d2.dot(&r);
    let eps = 1e-30;
    if a <= eps && e <= eps {
        return (0.0, 0.0);
    }
    if a <= eps {
        return (0.0, (f / e).clamp(0.0, 1.0));
    }
    let c = d1.dot(&r);
    if e <= eps {
        return ((-c / a).clamp(0.0, 1.0), 0.0);
    }
    let b = d1.dot(&d2);
    let denom = a * e - b * b;
    let mut s = if denom > 1e-14 * a * e {
        ((b * f - c * e) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut t = (b * s + f) / e;
    if t < 0.0 {
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - c) / a).clamp(0.0, 1.0);
    }
    (s, t)
}

/// Distance d'un SEGMENT [a, b] à un TRIANGLE, en forme fermée (Ericson
/// 5.1.10) : le minimum est une intersection (0), ou entre le segment et une
/// des trois arêtes, ou entre un bout du segment et le triangle. Rend
/// (distance, s sur [a, b], point du triangle). Remplace la section dorée
/// (84 évaluations par triangle) du CCD du 7 sept.
pub fn seg_tri(a: &V3, b: &V3, t0: &V3, t1: &V3, t2: &V3) -> (f64, f64, V3) {
    let n = (t1 - t0).cross(&(t2 - t0));
    let da = (a - t0).dot(&n);
    let db = (b - t0).dot(&n);
    if n.norm_squared() > 0.0 && da * db <= 0.0 && (da != 0.0 || db != 0.0) {
        let s = da / (da - db);
        let p = a + (b - a) * s;
        let c0 = (t1 - t0).cross(&(p - t0)).dot(&n);
        let c1 = (t2 - t1).cross(&(p - t1)).dot(&n);
        let c2 = (t0 - t2).cross(&(p - t2)).dot(&n);
        if c0 >= 0.0 && c1 >= 0.0 && c2 >= 0.0 {
            return (0.0, s, p);
        }
    }
    let mut best = (f64::INFINITY, 0.0, *a);
    for (s, p) in [(0.0, *a), (1.0, *b)] {
        let q = point_triangle(&p, t0, t1, t2);
        let dd = (q - p).norm();
        if dd < best.0 {
            best = (dd, s, q);
        }
    }
    for (e0, e1) in [(t0, t1), (t1, t2), (t2, t0)] {
        let (s, t) = seg_seg(a, b, e0, e1);
        let p = a + (b - a) * s;
        let q = e0 + (e1 - e0) * t;
        let dd = (q - p).norm();
        if dd < best.0 {
            best = (dd, s, q);
        }
    }
    best
}

/// Minimum d'une fonction CONVEXE le long du segment [a, b], et son abscisse
/// — section dorée. La distance d'un point à un convexe (boîte, cylindre) est
/// convexe, et composée avec un mouvement affine elle l'est en t : le minimum
/// est GLOBAL, ce n'est pas un échantillonnage. 80 itérations : 0,618⁸⁰ ≈
/// 2e-17, l'intervalle est sous l'arrondi. Sert au CCD (sphère balayée) et à
/// l'abscisse `s_axe` d'une capsule contre boîte ou cylindre.
pub fn min_convexe_segment(f: impl Fn(&V3) -> f64, a: &V3, b: &V3) -> (f64, f64) {
    let d = b - a;
    let p = |t: f64| a + d * t;
    let g = 0.5 * (5f64.sqrt() - 1.0);
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    let (mut t1, mut t2) = (hi - g * (hi - lo), lo + g * (hi - lo));
    let (mut f1, mut f2) = (f(&p(t1)), f(&p(t2)));
    for _ in 0..80 {
        if f1 < f2 {
            hi = t2;
            t2 = t1;
            f2 = f1;
            t1 = hi - g * (hi - lo);
            f1 = f(&p(t1));
        } else {
            lo = t1;
            t1 = t2;
            f1 = f2;
            t2 = lo + g * (hi - lo);
            f2 = f(&p(t2));
        }
    }
    [(f1, t1), (f2, t2), (f(a), 0.0), (f(b), 1.0)]
        .into_iter()
        .fold((f64::INFINITY, 0.0), |m, c| if c.0 < m.0 { c } else { m })
}

/// Distance d'un point à une BOÎTE centrée à l'origine (repère boîte) — 0 dedans.
pub fn dist_boite(loc: &V3, demi: &V3) -> f64 {
    V3::new(
        (loc.x.abs() - demi.x).max(0.0),
        (loc.y.abs() - demi.y).max(0.0),
        (loc.z.abs() - demi.z).max(0.0),
    )
    .norm()
}

/// Distance d'un point à un CYLINDRE FINI centré à l'origine — 0 dedans.
pub fn dist_cylindre(loc: &V3, axe: &V3, rc: f64, lc: f64) -> f64 {
    let e = axe.normalize();
    let z = loc.dot(&e);
    let rho = (loc - e * z).norm();
    let dz = (z.abs() - lc).max(0.0);
    let dr = (rho - rc).max(0.0);
    (dz * dz + dr * dr).sqrt()
}

/// Point le plus proche d'un point sur un TRIANGLE — les sept régions de
/// Voronoï (Ericson, *Real-Time Collision Detection*), en forme fermée.
///
/// ⚠ C⁰ mais pas C¹ aux frontières de région, comme la capsule : la dérivée
/// est celle de la région active, correcte hors ces surfaces de mesure nulle.
pub fn point_triangle(p: &V3, a: &V3, b: &V3, c: &V3) -> V3 {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(&ap), ac.dot(&ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return *a;
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(&bp), ac.dot(&bp));
    if d3 >= 0.0 && d4 <= d3 {
        return *b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(&cp), ac.dot(&cp));
    if d6 >= 0.0 && d5 <= d6 {
        return *c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let den = 1.0 / (va + vb + vc);
    a + ab * (vb * den) + ac * (vc * den)
}

/// SUPERÉLÉMENT — la raideur d'un maillage ÉLÉMENTS FINIS quelconque,
/// importée dans le multicorps.
///
/// C'est ainsi qu'un noyau multicorps devient GÉNÉRALISTE côté flexible : on
/// ne maille pas une coque ou un solide dans le solveur multicorps (ADAMS et
/// RecurDyn ne le font pas non plus), on prend la matrice de raideur d'un
/// code EF, on la RÉDUIT — Craig–Bampton, que ce dépôt a déjà — et on
/// l'attache à des nœuds qui sont des `Corps` ordinaires. Toute la machinerie
/// existante suit : liaisons, contact, modes, adjoints.
///
/// ```text
///     f = −BᵀK·u, B = ∂u/∂q ; u = déformation relative à la pose de repos
/// ```
///
/// COROTATIONNEL : `u` est mesuré dans un repère qui SUIT le corps. Une
/// matrice de raideur constante en repère global est fausse dès qu'il y a une
/// rotation d'ensemble — elle fabrique des forces de rappel là où il n'y a
/// aucune déformation. Le repère suit la rotation du premier nœud ; c'est la
/// version la plus simple du corotationnel, valable tant que la déformation
/// RELATIVE reste petite, ce qui est l'hypothèse d'une base modale de toute
/// façon.
///
/// ⚠ NON COUVERT, déclaré : la raideur géométrique (le raidissement d'un
/// corps précontraint ou en rotation rapide). Pour ça, les poutres
/// géométriquement exactes du noyau restent le bon outil — c'est d'ailleurs
/// pourquoi elles ont été faites AVANT, et pourquoi Craig–Bampton seul manque
/// la raideur centrifuge d'une pale.
#[derive(Clone, Debug)]
pub struct Superelement {
    pub nom: String,
    /// nœuds (des corps), 6 ddl chacun dans l'ordre (translation, rotation)
    pub noeuds: Vec<usize>,
    /// K symétrique (6n × 6n), en repère de RÉFÉRENCE du superélément
    pub k: DMatrix<f64>,
    /// α non nul non pris en charge ; amortissement objectif β Bᵀ K B
    pub alpha: f64,
    pub beta: f64,
    /// pose de repos : position de chaque nœud dans le repère du nœud 0
    pub ref_p: Vec<V3>,
    pub ref_r: Vec<M3>,
}

#[derive(Clone, Debug)]
pub enum Elem {
    L(Liaison),
    D(Distance),
    E(Engrenage),
    V(Vis),
    C(Cardan),
    /// contact NON LISSE : une ligne de complémentarité sur `Modele.contacts[ct]`
    K(Unilateral),
}

#[derive(Clone, Debug)]
pub struct Unilateral {
    pub nom: String,
    pub ct: usize,
    pub a: usize,
    pub b: Option<usize>,
}

impl Elem {
    pub fn nom(&self) -> &str {
        match self {
            Elem::L(l) => &l.nom,
            Elem::D(d) => &d.nom,
            Elem::E(e) => &e.nom,
            Elem::V(v) => &v.nom,
            Elem::C(c) => &c.nom,
            Elem::K(k) => &k.nom,
        }
    }
    pub fn n(&self) -> usize {
        match self {
            Elem::L(l) => l.bt.len() + l.br.len(),
            Elem::D(_) | Elem::E(_) | Elem::V(_) | Elem::C(_) | Elem::K(_) => 1,
        }
    }
    /// Les corps sur lesquels portent les colonnes du jacobien local, par blocs de 6.
    pub fn corps(&self) -> Vec<Ref> {
        match self {
            Elem::L(l) => vec![l.a, l.b],
            Elem::D(d) => vec![d.a, d.b],
            Elem::E(e) => vec![e.a, e.b, e.c],
            Elem::V(v) => vec![v.a, v.b],
            Elem::C(c) => vec![c.a, c.b],
            Elem::K(k) => vec![Some(k.a), k.b],
        }
    }
}

const BATI_R: V3 = V3::new(0.0, 0.0, 0.0);

fn pose_de(corps: &[Corps], c: Ref) -> (V3, M3) {
    match c {
        Some(i) => (corps[i].r, corps[i].rot),
        None => (BATI_R, M3::identity()),
    }
}

fn partie_basse(corps: &[Corps], c: Ref) -> V3 {
    c.map_or(V3::zeros(), |i| corps[i].r_bas)
}

impl Liaison {
    /// Le résidu de rotation vee(skew E) = sin θ·n s'annule aussi à θ = π :
    /// solution PARASITE, hors bassin de Newton (JOURNAL, 3 sept.). Indicateur
    /// signé : > 0 dans le bassin, < 0 sur la branche retournée. Pivot (deux
    /// lignes bloquées) : E_kk sur l'axe LIBRE, insensible à la rotation
    /// propre ; encastrement : tr E. Une seule ligne ou aucune : `None` — deux
    /// axes libres rendent le critère indécidable et la parasite ne s'y voit pas.
    pub fn bassin(&self, corps: &[Corps], t: f64) -> Option<f64> {
        if self.br.len() < 2 {
            return None;
        }
        let (_, rota) = pose_de(corps, self.a);
        let (_, rotb) = pose_de(corps, self.b);
        let (ua, ub) = (rota * self.ra, rotb * self.rb);
        let cr = match &self.cible_r {
            Some((axe, loi)) => expm(&(axe * loi.valeur(t))),
            None => M3::identity(),
        };
        let e = ua.transpose() * ub * cr.transpose();
        Some(if self.br.len() == 3 {
            e.trace()
        } else {
            let k = (0..3).find(|k| !self.br.contains(k)).unwrap();
            e[(k, k)]
        })
    }

    /// ∂Φ/∂t à pose FIGÉE — ce que la ligne GGL G·u + Φ_t = 0 demande, à la pose
    /// COURANTE du Newton. Pris par différence finie à la pose du début du pas,
    /// il était en retard de O(h) sur G : sur une rotation imposée à 8 rad/s,
    /// GGL rendait Φ̇ 60× PIRE que l'index 3 (2,6e-4 contre 4,3e-6) et la tête
    /// S2 ne passait pas son second pas. Translation : −dir·θ̇ ; rotation :
    /// E = uaᵀ·ub·crᵀ, ∂E/∂t = −E·[axe]×·θ̇, φ_t = vee(½(∂E − ∂Eᵀ)).
    pub fn phi_dt(&self, corps: &[Corps], t: f64) -> DVector<f64> {
        let n = self.bt.len() + self.br.len();
        let mut out = DVector::zeros(n);
        if let Some((dir, loi)) = &self.cible_t {
            for (row, &i) in self.bt.iter().enumerate() {
                out[row] = -dir[i] * loi.derivee(t);
            }
        }
        if let Some((axe, loi)) = &self.cible_r {
            let (_, rota) = pose_de(corps, self.a);
            let (_, rotb) = pose_de(corps, self.b);
            let e = (rota * self.ra).transpose()
                * (rotb * self.rb)
                * expm(&(axe * loi.valeur(t))).transpose();
            let de = -e * skew(axe) * loi.derivee(t);
            let v = vee(&(0.5 * (de - de.transpose())));
            for (k, &i) in self.br.iter().enumerate() {
                out[self.bt.len() + k] = v[i];
            }
        }
        out
    }

    /// Φ (n) et G (n × 12) sur (δr_a, δθ_a, δr_b, δθ_b).
    pub fn phi_g(&self, corps: &[Corps], t: f64) -> (DVector<f64>, DMatrix<f64>) {
        let (ra_, rota) = pose_de(corps, self.a);
        let (rb_, rotb) = pose_de(corps, self.b);
        let ua = rota * self.ra;
        let ub = rotb * self.rb;
        let qa = rota * self.pa;
        // NON HOLONOME : le bras vers B suit le point de contact GÉOMÉTRIQUE ;
        // ce n'est PAS un point matériel figé dans B. Une roue qui roule voit
        // son point de contact changer de matière à chaque instant ; garder
        // `pb` fixe suit un point de la jante qui s'éloigne du sol — mesuré,
        // la roue décélérait et repartait en arrière.
        let db = partie_basse(corps, self.b) - partie_basse(corps, self.a);
        let qb = if self.nh {
            (ra_ + qa) - rb_ - db
        } else {
            rotb * self.pb
        };
        let dw = (rb_ + qb) - (ra_ + qa);
        let d = ua.transpose() * dw;
        let dw = dw + db;
        let cr = match &self.cible_r {
            Some((axe, loi)) => expm(&(axe * loi.valeur(t))),
            None => M3::identity(),
        };
        let ct = match &self.cible_t {
            Some((dir, loi)) => dir * loi.valeur(t),
            None => V3::zeros(),
        };
        let e = ua.transpose() * ub * cr.transpose();
        let phi_t = (d - ct) + ua.transpose() * db;
        let phi_r = vee(&(0.5 * (e - e.transpose())));
        // δd = uaᵀ([dw]× + [qa]×) δθ_a − uaᵀ δr_a + uaᵀ δr_b − uaᵀ[qb]× δθ_b
        let uat = ua.transpose();
        let gt_ra = -uat;
        let gt_ta = uat * (skew(&dw) + skew(&qa));
        let gt_rb = uat;
        let gt_tb = -uat * skew(&qb);
        // δE = [w]× E, w = uaᵀ(δθ_b − δθ_a) ; δφ = vee(skew([w]× E)) linéaire en w
        let mut dwm = M3::zeros();
        for k in 0..3 {
            let mut ek = V3::zeros();
            ek[k] = 1.0;
            let s = skew(&ek) * e;
            let col = vee(&(0.5 * (s - s.transpose())));
            dwm.set_column(k, &col);
        }
        let gr_ta = -dwm * uat;
        let gr_tb = dwm * uat;
        let n = self.bt.len() + self.br.len();
        let mut phi = DVector::zeros(n);
        let mut g = DMatrix::zeros(n, 12);
        let mut row = 0;
        for &i in &self.bt {
            phi[row] = phi_t[i];
            for c in 0..3 {
                g[(row, c)] = gt_ra[(i, c)];
                g[(row, 3 + c)] = gt_ta[(i, c)];
                g[(row, 6 + c)] = gt_rb[(i, c)];
                g[(row, 9 + c)] = gt_tb[(i, c)];
            }
            row += 1;
        }
        for &i in &self.br {
            phi[row] = phi_r[i];
            for c in 0..3 {
                g[(row, 3 + c)] = gr_ta[(i, c)];
                g[(row, 9 + c)] = gr_tb[(i, c)];
            }
            row += 1;
        }
        (phi, g)
    }
}

impl Distance {
    pub fn phi_g(&self, corps: &[Corps]) -> (DVector<f64>, DMatrix<f64>) {
        let (ra_, rota) = pose_de(corps, self.a);
        let (rb_, rotb) = pose_de(corps, self.b);
        let qa = rota * self.pa;
        let qb = rotb * self.pb;
        let dw = (rb_ + qb) - (ra_ + qa);
        let db = partie_basse(corps, self.b) - partie_basse(corps, self.a);
        let residu = tangent::distance_ecart::<f64>(dw.into(), db.into(), self.l);
        let dw = dw + db;
        let l = dw.norm();
        let nrm = dw / l;
        let mut g = DMatrix::zeros(1, 12);
        let ga = nrm.transpose() * skew(&qa);
        let gb = -nrm.transpose() * skew(&qb);
        for c in 0..3 {
            g[(0, c)] = -nrm[c];
            g[(0, 3 + c)] = ga[c];
            g[(0, 6 + c)] = nrm[c];
            g[(0, 9 + c)] = gb[c];
        }
        (DVector::from_element(1, residu), g)
    }
}

// ── modèle ───────────────────────────────────────────────────────────────────
#[derive(Clone, Debug)]
pub struct Modele {
    pub corps: Vec<Corps>,
    pub elems: Vec<Elem>,
    pub g: V3,
    pub t: f64,
    pub lam: DVector<f64>,
    /// efforts extérieurs constants au CdM, en repère monde : (corps, F, M)
    pub efforts: Vec<(usize, V3, V3)>,
    /// couples à loi (servos, gouverneur, ressorts)
    pub couples: Vec<Couple>,
    pub pales: Vec<Pale>,
    pub contacts: Vec<Contact>,
    /// sphères candidates à l'APPARIEMENT automatique : (corps, centre, rayon)
    pub spheres: Vec<(usize, V3, f64)>,
    /// loi commune des contacts découverts : (k, expo, c, mu, v_eps, marge, d_hat)
    pub loi_contact: Option<(f64, f64, f64, f64, f64, f64, f64)>,
    /// paires actives au dernier appariement, et coût de la dernière passe
    pub n_paires: usize,
    pub t_appar: f64,
    pub inflows: Vec<Inflow>,
    pub poutres: Vec<Poutre>,
    pub supers: Vec<Superelement>,
    pub maillages: Vec<Maillage>,
    /// vitesse de l'air ambiant en repère MONDE : une section voit −vent en plus
    /// de sa vitesse propre. C'est ce qui donne μ (le rotor de banc est encastré,
    /// « avancer à V » s'écrit « vent de face à V »).
    pub vent: V3,
    /// itérations de Newton cumulées, et repli SVD (singularité) cumulés
    pub newton_total: usize,
    pub svd_total: usize,
    pub jacobien_total: usize,
    /// lignes de contrainte ACTIVES : les redondantes STRUCTURELLES (rang de G
    /// au départ) sont neutralisées — leur λ est fixé à 0, le mouvement est le
    /// même (mesuré : quatre-barres redondant = non redondant à 1e-14 m)
    pub actif: Vec<bool>,
    /// chronos cumulés (s) : jacobien, résolution, résidus — pour attribuer le temps
    pub t_jac: f64,
    pub t_sol: f64,
    pub t_res: f64,
    /// total passé dans la boucle de pas, et fin de pas (actualise, contrôles,
    /// projections) — sans eux, 82 % du temps n'était chronométré nulle part
    pub t_pas: f64,
    pub t_fin: f64,
    /// part du jacobien passée à construire la tangente de exp (les J_l) — le
    /// seul chiffre qui tranche le refus de Cayley ; borne basse, cf.
    /// `tangent::jacobien_ad`
    pub t_exp: f64,
    /// structure symbolique du LU creux, gardée d'un jacobien à l'autre
    pub symb: Option<std::sync::Arc<faer::sparse::linalg::lu::SymbolicLu<usize>>>,
    /// empreinte du MOTIF qui a servi à cette symbolique.
    ///
    /// ⚠ La symbolique gardée n'est valide que pour SON motif. Si le motif
    /// change — un terme qui n'est plus poussé, une ligne neutralisée — faer
    /// déborde ses tableaux et PANIQUE au travers de PyO3. Mesuré le
    /// 3 sept. sur un treillis 20×20 à boucles fermées (7 480 inconnues) :
    /// « range end index 95682 out of range for slice of length 95634 ».
    /// Le défaut est ANTÉRIEUR et ne se voyait qu'à cette taille.
    pub symb_motif: u64,
    pub h_dernier: f64,
    /// DÉTECTION DE COLLISION CONTINUE : refuse un pas dont la TRAJECTOIRE
    /// traverse un contact, même si ses deux extrémités sont saines. C'est la
    /// garantie topologique d'IPC — sans elle, un pas assez grand fait passer
    /// une sphère de l'autre côté sans que rien ne le voie.
    pub ccd: bool,
    /// PAS PILOTÉ PAR LE CONTACT : (h_libre, h_contact, marge). Hors contact
    /// et loin de lui, on avance à `h_libre` ; dès qu'une paire est à moins de
    /// `marge` de se toucher, on passe à `h_contact`. C'est la décomposition
    /// en temps sous sa forme la plus simple — le pas suit l'ÉCHELLE ACTIVE,
    /// et non la plus raide du modèle. Contrairement au résidu de demi-pas,
    /// l'estimateur est ici EXACT et gratuit : on connaît la distance.
    pub pas_contact: Option<(f64, f64, f64)>,
    /// STABILISATION GGL (Gear–Gupta–Leimkuhler, index 2) : m multiplicateurs
    /// ζ de plus par pas, la pose reçoit βh·G₀ᵀζ et la ligne G·u + Φ_t = 0
    /// entre au système. Φ ET Φ̇ sont tenus à la précision machine, l'ordre 2
    /// du schéma est conservé (Arnold–Brüls 2007) — là où la projection de
    /// vitesse après coup le cassait. Coût : m inconnues et m colonnes de
    /// jacobien par différences finies. Opt-in : l'index 3 direct reste le
    /// défaut, c'est lui qui tient l'échelle.
    pub ggl: bool,
    /// Modification géométrique expérimentale : 0 = schéma historique.
    /// La relation cinématique implicite est résolue dans ||theta|| < pi.
    pub sigma_lie: f64,
    /// échelle de λ dans Fischer–Burmeister : ech_l / (1 + ‖f‖), posée par pas
    pub s_fb: f64,
    /// contacts non lisses, par pas : (actif, e·ġap₀ + gap₀/h), indexé comme `contacts`
    pub k_pas: Vec<(bool, f64)>,
    /// PROJECTION SUR LA VARIÉTÉ DU MOMENT (Noether discret par projection)
    /// Tolérance sur ‖Φ‖ après chaque pas accepté (0 = contrôle désactivé).
    /// Ce n'est pas un réglage de précision : c'est le garde-fou qui distingue
    /// une redondance LICITE (Φ reste nul) d'une contradiction entre liaisons.
    pub tol_phi: f64,
    pub proj_moment: bool,
    /// moment cinétique de référence, capturé au premier pas
    pub l_ref: Option<V3>,
    /// PROJECTION SUR L'ÉNERGIE (schéma énergie-moment par projection)
    pub proj_energie: bool,
    /// énergie de référence, capturée au premier pas
    pub e_ref: Option<f64>,
    /// tolérance RELATIVE du résidu de demi-pas ; None = pas FIXE (défaut).
    /// Voir `simule` : c'est le contrôle de pas de Hibbitt & Karlsson (1979).
    pub adapt: Option<f64>,
    /// pas rejoués par le contrôle adaptatif au dernier `simule`
    pub n_rejeu: usize,
    /// diagnostic du contrôle de pas : (max de R_half, max du seuil) au dernier run
    pub r_half_max: f64,
    /// max de |Φ| à mi-pas, relatif à la taille du modèle (adaptatif)
    pub phi_half_max: f64,
    pub seuil_max: f64,
    /// bornes du pas adaptatif, en multiples du pas demandé
    pub adapt_bornes: (f64, f64),
    /// accumulateur du torseur aérodynamique TOTAL réduit à l'origine, armé par
    /// `moyenne_aero(t0)` : (t0, ∫F dt, ∫M dt, ∫dt). C'est la sortie que trime
    /// un banc de rotor — une moyenne sur un tour, pas une valeur instantanée.
    pub moy_aero: Option<(f64, V3, V3, f64)>,
    /// HISTORIQUE DU SCHÉMA (opt-in) : (u̇, a, λ) à l'état initial puis après
    /// chaque pas accepté. C'est ce qu'un adjoint EN TEMPS doit connaître pour
    /// remonter — les frames de `simule` ne portent que (q, u, λ), pas les deux
    /// accélérations que l'α-généralisé traîne d'un pas à l'autre.
    pub hist_schema: Option<Vec<(Vec<f64>, Vec<f64>, Vec<f64>)>>,
    /// MULTI-RYTHME : corps GELÉS pendant un sous-pas (posés de l'extérieur,
    /// u̇ = 0 pour Newton, ni avancés ni relus), et éléments dont toutes les
    /// lignes sont gelées avec eux. Vides hors multi-rythme.
    pub gel: Vec<bool>,
    pub gel_elem: Vec<bool>,
    /// ÉTAT DU SCHÉMA à passer d'un `simule` au suivant : (u̇, a, λ). Sans lui,
    /// chaque appel REDÉMARRE sur l'accélération consistante — et l'α-généralisé
    /// redémarré à chaque pas n'est plus que d'ORDRE 1 (a₀ = u̇₀ casse le
    /// retard algorithmique qui fait l'ordre 2 : mesuré, 4 cm d'écart au bout
    /// d'une chaîne de 2 m en 0,3 s sous multi-rythme). `simule` consomme
    /// `schema_init` s'il est posé et dépose `schema_fin` en sortant.
    pub schema_init: Option<(DVector<f64>, DVector<f64>, DVector<f64>)>,
    pub schema_fin: Option<(DVector<f64>, DVector<f64>, DVector<f64>)>,
}

impl Corps {
    fn pose_precise(&self) -> Pose {
        (self.r, self.rot, self.r_bas)
    }
    fn deplace_depuis(&mut self, r: V3, bas: V3, dr: V3) {
        (self.r, self.r_bas) = numerique::deplace_compense(r, bas, dr);
    }
    fn deplace(&mut self, dr: V3) {
        self.deplace_depuis(self.r, self.r_bas, dr);
    }
    /// Différence des positions, sans perdre les petits termes au voisinage
    /// de grandes coordonnées absolues. Les deux parties restent séparées.
    fn difference_precise(&self, a: &Corps) -> (V3, V3) {
        let mut haut = V3::zeros();
        let mut bas = V3::zeros();
        for i in 0..3 {
            let (s, e) = numerique::deux_sommes(self.r[i], -a.r[i]);
            (haut[i], bas[i]) = numerique::deux_sommes(s, e + (self.r_bas[i] - a.r_bas[i]));
        }
        (haut, bas)
    }
    /// composante c de (v, ω)
    pub fn u6(&self, c: usize) -> f64 {
        if c < 3 {
            self.v[c]
        } else {
            self.w[c - 3]
        }
    }
}

/// Position haute, orientation, position basse ; snapshot interne complet de pose.
type Pose = (V3, M3, V3);

impl Modele {
    pub fn new(g: V3) -> Self {
        Modele {
            corps: vec![],
            elems: vec![],
            g,
            t: 0.0,
            lam: DVector::zeros(0),
            efforts: vec![],
            couples: vec![],
            pales: vec![],
            contacts: vec![],
            spheres: vec![],
            loi_contact: None,
            n_paires: 0,
            t_appar: 0.0,
            inflows: vec![],
            poutres: vec![],
            supers: vec![],
            maillages: vec![],
            vent: V3::zeros(),
            newton_total: 0,
            svd_total: 0,
            jacobien_total: 0,
            actif: vec![],
            t_jac: 0.0,
            t_sol: 0.0,
            t_res: 0.0,
            t_pas: 0.0,
            t_fin: 0.0,
            t_exp: 0.0,
            symb: None,
            symb_motif: 0,
            h_dernier: 0.0,
            ccd: false,
            pas_contact: None,
            tol_phi: 1e-6,
            ggl: false,
            sigma_lie: 0.0,
            s_fb: 1.0,
            k_pas: vec![],
            proj_moment: false,
            l_ref: None,
            proj_energie: false,
            e_ref: None,
            adapt: None,
            n_rejeu: 0,
            r_half_max: 0.0,
            phi_half_max: 0.0,
            seuil_max: 0.0,
            adapt_bornes: (1.0 / 64.0, 16.0),
            moy_aero: None,
            hist_schema: None,
            gel: vec![],
            gel_elem: vec![],
            schema_init: None,
            schema_fin: None,
        }
    }

    fn gele(&self, i: usize) -> bool {
        self.gel.get(i).copied().unwrap_or(false)
    }

    pub fn n(&self) -> usize {
        6 * self.corps.len()
    }
    /// ÉCHELLE DE LONGUEUR du modèle : la plus grande distance entre corps
    /// (plancher 1). Toute borne exprimée en LONGUEUR doit passer par elle —
    /// une constante en mètres écrite en dur fait qu'un mécanisme coté en
    /// micromètres n'est plus le même mécanisme. Mesuré le 3 sept. par le banc
    /// d'invariance d'unités : à l'échelle MEMS, le clamp d'accélération
    /// initiale (0,05 m) décalait la réponse de 3,7 % en silence.
    pub fn echelle_l(&self) -> f64 {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for c in &self.corps {
            for k in 0..3 {
                lo[k] = lo[k].min(c.r[k]);
                hi[k] = hi[k].max(c.r[k]);
            }
        }
        let d = (0..3)
            .map(|k| (hi[k] - lo[k]).max(0.0))
            .fold(0.0_f64, |a, b| a + b * b)
            .sqrt();
        if d.is_finite() && d > 0.0 {
            d.max(1e-300)
        } else {
            1.0
        }
    }

    pub fn m(&self) -> usize {
        self.elems.iter().map(|e| e.n()).sum()
    }

    fn u(&self) -> DVector<f64> {
        let mut u = DVector::zeros(self.n());
        for (i, c) in self.corps.iter().enumerate() {
            for k in 0..3 {
                u[6 * i + k] = c.v[k];
                u[6 * i + 3 + k] = c.w[k];
            }
        }
        u
    }
    fn set_u(&mut self, u: &DVector<f64>) {
        for (i, c) in self.corps.iter_mut().enumerate() {
            if self.gel.get(i).copied().unwrap_or(false) {
                continue;
            }
            c.v = V3::new(u[6 * i], u[6 * i + 1], u[6 * i + 2]);
            c.w = V3::new(u[6 * i + 3], u[6 * i + 4], u[6 * i + 5]);
        }
    }
    fn poses(&self) -> Vec<Pose> {
        self.corps.iter().map(Corps::pose_precise).collect()
    }
    fn set_poses(&mut self, p: &[Pose]) {
        for (c, q) in self.corps.iter_mut().zip(p) {
            c.r = q.0;
            c.r_bas = q.2;
            c.rot = q.1;
        }
    }
    /// q_{n+1} = q_n ⊕ h·dq : translation vectorielle, rotation sur le groupe.
    fn avance(&mut self, poses0: &[Pose], dq: &DVector<f64>, h: f64) {
        for (i, c) in self.corps.iter_mut().enumerate() {
            if self.gel.get(i).copied().unwrap_or(false) {
                continue;
            }
            let dr = V3::new(dq[6 * i], dq[6 * i + 1], dq[6 * i + 2]);
            let dth = V3::new(dq[6 * i + 3], dq[6 * i + 4], dq[6 * i + 5]);
            c.deplace_depuis(poses0[i].0, poses0[i].2, h * dr);
            c.rot = expm(&(h * dth)) * poses0[i].1;
        }
    }

    /// Échelle PHYSIQUE du bloc dynamique du résidu : ‖M·u̇‖ + ‖ω × J_s ω‖.
    ///
    /// Sans elle, le critère de Newton est mis à l'échelle sur les seules
    /// forces EXTÉRIEURES — et un système dont toute la dynamique est
    /// inertielle (un rotor : pas de pesanteur, aucun effort extérieur au
    /// départ) a `forces()` ≈ 0. Le critère devient alors ABSOLU à 1e-12 sur
    /// des termes qui valent 2e4, et Newton « ne converge pas » sur un résidu
    /// déjà relatif à 1e-10. Mesuré le 3 sept. sur la pale battante à
    /// Ω = 109 rad/s : convergent à h, refusé à h/2.
    fn echelle_inertie(&self, udot: &DVector<f64>) -> f64 {
        let mut s = 0.0;
        for (i, c) in self.corps.iter().enumerate() {
            let k = 6 * i;
            let js = c.rot * c.j * c.rot.transpose();
            let wd = V3::new(udot[k + 3], udot[k + 4], udot[k + 5]);
            let vd = V3::new(udot[k], udot[k + 1], udot[k + 2]);
            s += (c.m * vd).norm() + (js * wd).norm() + c.w.cross(&(js * c.w)).norm();
        }
        s
    }

    fn forces(&self) -> DVector<f64> {
        self.forces_a(self.t)
    }

    fn forces_a(&self, t: f64) -> DVector<f64> {
        let mut f = DVector::zeros(self.n());
        for (i, c) in self.corps.iter().enumerate() {
            let js = c.rot * c.j * c.rot.transpose();
            let fg = c.m * self.g;
            let gyr = -c.w.cross(&(js * c.w));
            for k in 0..3 {
                f[6 * i + k] = fg[k];
                f[6 * i + 3 + k] = gyr[k];
            }
        }
        for (i, fo, mo) in &self.efforts {
            for k in 0..3 {
                f[6 * i + k] += fo[k];
                f[6 * i + 3 + k] += mo[k];
            }
        }
        for c in &self.couples {
            let tau = c.valeur(&self.corps, t);
            if let Some(ib) = c.b {
                for k in 0..3 {
                    f[6 * ib + 3 + k] += tau[k];
                }
            }
            if let Some(ia) = c.a {
                for k in 0..3 {
                    f[6 * ia + 3 + k] -= tau[k];
                }
            }
        }
        for p in &self.pales {
            let (fo, mo, _) = p.valeur(&self.corps, &self.inflows, self.vent);
            for k in 0..3 {
                f[6 * p.corps + k] += fo[k];
                f[6 * p.corps + 3 + k] += mo[k];
            }
        }
        for b in &self.poutres {
            let fl = b.forces(&self.corps);
            for k in 0..6 {
                f[6 * b.a + k] += fl[k];
                f[6 * b.b + k] += fl[6 + k];
            }
        }
        for se in &self.supers {
            let fl = se.efforts(&self.corps, se.beta);
            for (i, &ni) in se.noeuds.iter().enumerate() {
                for k in 0..6 {
                    f[6 * ni + k] += fl[6 * i + k];
                }
            }
        }
        for ct in &self.contacts {
            // non lisse : le frottement porte sur le DERNIER λ (la normale est
            // Gᵀλ, hors des forces) — sans quoi l'accélération lisse du
            // redémarrage après impact ignorait Coulomb (mesuré : +7,9 % sur
            // l'accélération d'une bille qui glisse)
            let ct_nl;
            let ct = if ct.nonlisse {
                ct_nl = Contact {
                    fn_impose: Some(self.lam.get(ct.row).copied().unwrap_or(0.0).max(0.0)),
                    ..ct.clone()
                };
                &ct_nl
            } else {
                ct
            };
            let (fo, mo, fb, mb, _) = ct.valeur(&self.corps);
            for k in 0..3 {
                f[6 * ct.corps + k] += fo[k];
                f[6 * ct.corps + 3 + k] += mo[k];
            }
            if let Some(ib) = ct.b {
                for k in 0..3 {
                    f[6 * ib + k] += fb[k];
                    f[6 * ib + 3 + k] += mb[k];
                }
            }
        }
        f
    }

    /// G·u SANS assembler G — le produit se fait élément par élément.
    ///
    /// Même motif que `phi_seul` : `phi_g` alloue m × n en dense. `acc_init`
    /// l'appelait DEUX fois pour n'en tirer que deux produits matrice-vecteur,
    /// et `detecte_redondance` une fois avant même de tester son seuil de
    /// saut. Mesuré le 3 sept. : trois assemblages, ~5 Go et ~4 s
    /// d'initialisation à 30 000 inconnues, à chaque appel de `simule`.
    pub fn g_u(&self, t: f64, u: &DVector<f64>) -> Result<DVector<f64>, String> {
        let m = self.m();
        let mut out = DVector::zeros(m);
        let mut row = 0;
        for e in &self.elems {
            let (_, ge) = match e {
                Elem::L(l) => l.phi_g(&self.corps, t),
                Elem::D(d) => d.phi_g(&self.corps),
                Elem::E(en) => en.phi_g(&self.corps)?,
                Elem::V(vi) => vi.phi_g(&self.corps)?,
                Elem::C(c) => c.phi_g(&self.corps),
                Elem::K(_) => (DVector::zeros(1), DMatrix::zeros(1, 12)),
            };
            let ne = e.n();
            for i in 0..ne {
                let mut s = 0.0;
                for (kb, cb) in e.corps().iter().enumerate() {
                    if let Some(ib) = cb {
                        for c in 0..6 {
                            s += ge[(i, 6 * kb + c)] * u[6 * ib + c];
                        }
                    }
                }
                out[row + i] = s;
            }
            row += ne;
        }
        Ok(out)
    }

    /// max |Φ| sur les seules contraintes HOLONOMES.
    ///
    /// Une contrainte non holonome n'a pas de Φ à annuler : le `phi` de sa
    /// liaison mesure la position du point de contact, qui DÉRIVE
    /// légitimement — c'est même tout l'intérêt, la roue avance. Le contrôle
    /// de fin de pas doit donc la sauter, sinon il refuse le roulement.
    pub fn phi_holonome_max(&self, t: f64) -> Result<f64, String> {
        let mut e: f64 = 0.0;
        for (ie, el) in self.elems.iter().enumerate() {
            if matches!(el, Elem::L(l) if l.nh) || self.gel_elem.get(ie).copied().unwrap_or(false) {
                continue;
            }
            let (p, _) = match el {
                Elem::L(l) => l.phi_g(&self.corps, t),
                Elem::D(d) => d.phi_g(&self.corps),
                Elem::E(en) => en.phi_g(&self.corps)?,
                Elem::V(vi) => vi.phi_g(&self.corps)?,
                Elem::C(c) => c.phi_g(&self.corps),
                Elem::K(_) => (DVector::zeros(1), DMatrix::zeros(1, 12)),
            };
            e = e.max(p.amax());
        }
        Ok(e)
    }

    /// quelles lignes de λ sont des contacts non lisses (Fischer–Burmeister)
    fn lignes_k(&self) -> Vec<bool> {
        self.elems
            .iter()
            .flat_map(|e| std::iter::repeat_n(matches!(e, Elem::K(_)), e.n()))
            .collect()
    }
    /// … et parmi elles, celles dont le contact est INACTIF sur ce pas
    fn lignes_k_inactives(&self) -> Vec<bool> {
        self.elems
            .iter()
            .flat_map(|e| {
                let k = match e {
                    Elem::K(k) => !matches!(self.k_pas.get(k.ct), Some((true, _))),
                    _ => false,
                };
                std::iter::repeat_n(k, e.n())
            })
            .collect()
    }

    /// G en (ligne, colonne, valeur), lignes ACTIVES seulement — élément par
    /// élément, jamais densifiée.
    fn g_creuse(&self, t: f64) -> Result<Vec<(usize, usize, f64)>, String> {
        self.g_creuse_pour_temoin(t, false)
    }

    fn g_creuse_pour_temoin(
        &self,
        t: f64,
        toutes: bool,
    ) -> Result<Vec<(usize, usize, f64)>, String> {
        let mut out = Vec::with_capacity(12 * self.m());
        let mut row = 0;
        for e in &self.elems {
            let (_, ge) = match e {
                Elem::L(l) => l.phi_g(&self.corps, t),
                Elem::D(d) => d.phi_g(&self.corps),
                Elem::E(en) => en.phi_g(&self.corps)?,
                Elem::V(vi) => vi.phi_g(&self.corps)?,
                Elem::C(c) => c.phi_g(&self.corps),
                Elem::K(k) => match self.k_pas.get(k.ct) {
                    Some((true, _)) => self.contacts[k.ct].phi_g(&self.corps),
                    _ => (DVector::zeros(1), DMatrix::zeros(1, 12)),
                },
            };
            let ne = e.n();
            for i in 0..ne {
                if !toutes && !self.actif.get(row + i).copied().unwrap_or(true) {
                    continue;
                }
                for (kb, cb) in e.corps().iter().enumerate() {
                    if let Some(ib) = cb {
                        for c in 0..6 {
                            let v = ge[(i, 6 * kb + c)];
                            if v != 0.0 {
                                out.push((row + i, 6 * ib + c, v));
                            }
                        }
                    }
                }
            }
            row += ne;
        }
        Ok(out)
    }

    /// Refuse un état où une liaison est sur la branche retournée de son
    /// résidu de rotation (`Liaison::bassin`) — un pas ou un assemblage qui y
    /// « converge » a satisfait Φ = 0 sur la mauvaise solution.
    pub fn verifie_bassin(&self, t: f64) -> Result<(), String> {
        for (ie, e) in self.elems.iter().enumerate() {
            if self.gel_elem.get(ie).copied().unwrap_or(false) {
                continue;
            }
            if let Elem::L(l) = e {
                if let Some(b) = l.bassin(&self.corps, t) {
                    if b < 0.0 {
                        return Err(format!(
                            "liaison « {} » : rotation relative au-delà de 90° sur ses axes bloqués \
                             (indicateur {b:.3}) — solution parasite θ = π du résidu vee(skew E)",
                            l.nom
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Φ SEUL, sans assembler G.
    ///
    /// ⚠ `phi_g` alloue une matrice DENSE m × n — 1,7 Go et 1,3 s à 30 000
    /// inconnues. Tout ce qui ne veut que Φ doit passer par ici : le contrôle
    /// de fin de pas, le diagnostic, l'assemblage. Mesuré le 3 sept. après
    /// avoir introduit MOI-MÊME la régression : le contrôle `tol_phi` posé
    /// le matin appelait `phi_g` à chaque pas accepté, ce qui représentait
    /// 97 % du temps par pas et toute la mémoire — et j'ai mesuré cette
    /// régression comme si c'était une propriété du solveur.
    pub fn phi_seul(&self, t: f64) -> Result<DVector<f64>, String> {
        let m = self.m();
        let mut phi = DVector::zeros(m);
        let mut row = 0;
        for e in &self.elems {
            let (p, _) = match e {
                Elem::L(l) => l.phi_g(&self.corps, t),
                Elem::D(d) => d.phi_g(&self.corps),
                Elem::E(en) => en.phi_g(&self.corps)?,
                Elem::V(vi) => vi.phi_g(&self.corps)?,
                Elem::C(c) => c.phi_g(&self.corps),
                Elem::K(_) => (DVector::zeros(1), DMatrix::zeros(1, 12)),
            };
            let ne = e.n();
            phi.rows_mut(row, ne).copy_from(&p.rows(0, ne));
            row += ne;
        }
        Ok(phi)
    }

    pub fn phi_g(&self, t: f64) -> Result<(DVector<f64>, DMatrix<f64>), String> {
        let (n, m) = (self.n(), self.m());
        let mut phi = DVector::zeros(m);
        let mut g = DMatrix::zeros(m, n);
        let mut row = 0;
        for e in &self.elems {
            let (p, ge) = match e {
                Elem::L(l) => l.phi_g(&self.corps, t),
                Elem::D(d) => d.phi_g(&self.corps),
                Elem::E(en) => en.phi_g(&self.corps)?,
                Elem::V(vi) => vi.phi_g(&self.corps)?,
                Elem::C(c) => c.phi_g(&self.corps),
                Elem::K(_) => (DVector::zeros(1), DMatrix::zeros(1, 12)),
            };
            let ne = e.n();
            for i in 0..ne {
                phi[row + i] = p[i];
                for (kb, cb) in e.corps().iter().enumerate() {
                    if let Some(ib) = cb {
                        for c in 0..6 {
                            g[(row + i, 6 * ib + c)] += ge[(i, 6 * kb + c)];
                        }
                    }
                }
            }
            row += ne;
        }
        Ok((phi, g))
    }

    /// Après un pas accepté : les engrenages mémorisent leurs angles continus.
    fn actualise(&mut self) -> Result<(), String> {
        let corps = &self.corps;
        for e in &mut self.elems {
            if let Elem::V(vi) = e {
                if let Ok(th) = vi.angle(corps) {
                    vi.prev = th;
                }
            }
            if let Elem::E(en) = e {
                en.prev = en.angles(corps)?;
                // REBASAGE : Φ = θ_a − r·θ_b est invariant sous (θ_a − 2πk·r, θ_b − 2πk),
                // et un angle DÉROULÉ n'a que la précision de son ulp — à 65 536 rad
                // (le pignon S2 à 50,8 s de vol) l'ulp vaut 1,5e-11, la ligne mise
                // à l'échelle 1/(βh²) ne descend plus sous la tolérance de Newton, et
                // le vol casse (mesuré : résidu 5,25e-3 à h = 1e-4, 1,31e-1 à 2e-5).
                // On ramène θ_b dans un tour, θ_a suit ; rien de physique ne bouge.
                // ponytail: la vis-écrou garde son angle déroulé (sa translation
                // porte la même origine) — même mur à ~10 000 tours de vis.
                let k = (en.prev.1 / std::f64::consts::TAU).round();
                // Le déroulement n'autorise que des décalages de TOURS ENTIERS
                // sur chacun des arbres. Pour un rapport fractionnaire, 2πkr
                // n'est pas un tour : le pas suivant voyait un faux saut à π.
                // On conserve alors les angles déroulés, sans changer le repère.
                if k != 0.0 && en.rapport.fract() == 0.0 {
                    en.prev.1 -= std::f64::consts::TAU * k;
                    en.prev.0 -= std::f64::consts::TAU * k * en.rapport;
                }
            }
        }
        for c in &mut self.couples {
            c.prev = c.angle(corps);
        }
        // inflow : la poussée du pas accepté nourrit v_i au pas suivant (retard τ)
        let mut pousse = vec![0.0; self.inflows.len()];
        let inflows = &self.inflows;
        let vent = self.vent;
        let (mut f_aero, mut m_aero) = (V3::zeros(), V3::zeros());
        let mut m_rt = vec![(0.0f64, 0.0f64); self.inflows.len()];
        let h_p = self.h_dernier;
        for p in &mut self.pales {
            let (fo, mo_, tq) = p.valeur(corps, inflows, vent);
            p.sortie = tq;
            if p.instat && h_p > 0.0 {
                // ż = −b z + α, en temps SEMI-CORDE : ds = 2|V| dt / c. C'est
                // ce changement d'échelle qui fait que le retard est une
                // distance parcourue en cordes, non un temps — une section
                // rapide « oublie » son sillage plus vite.
                let st = p.stations(corps, inflows, vent);
                for (k, (al, v, _)) in st.iter().enumerate() {
                    let ds = 2.0 * v.abs() * h_p / p.corde;
                    let (z1, z2) = p.z[k];
                    // exponentielle exacte sur le pas : stable quel que soit ds
                    let (e1, e2) = ((-W_B1 * ds).exp(), (-W_B2 * ds).exp());
                    p.z[k] = (
                        z1 * e1 + al / W_B1 * (1.0 - e1),
                        z2 * e2 + al / W_B2 * (1.0 - e2),
                    );
                }
            }
            if p.lb && h_p > 0.0 {
                // LEISHMAN–BEDDOES, en temps semi-corde ds = 2|V|dt/c, tout en
                // exponentielle exacte : retard de pression → α_f → f' = f_st(α_f)
                // → f'' (retard T_f) → tourbillon (seuil C_N1, convection T_vl,
                // décroissance T_v). Les états sont figés sur le pas suivant.
                let st = p.stations(corps, inflows, vent);
                let (t_p, t_f, t_v, t_vl) = p.lb_const;
                for (k, (al, v, _)) in st.iter().enumerate() {
                    let ds = 2.0 * v.abs() * h_p / p.corde;
                    let e = &mut p.etats_lb[k];
                    let cn_pot = p.a0 * (al.to_degrees() - p.alpha0);
                    e.dp =
                        e.dp * (-ds / t_p).exp() + (cn_pot - e.cn_prev) * (-0.5 * ds / t_p).exp();
                    e.cn_prev = cn_pot;
                    let cn_lag = cn_pot - e.dp;
                    let alpha_f = cn_lag / p.a0 + p.alpha0;
                    let f_p = tangent::interp::<f64>(&p.polaire.alpha, &p.fst, alpha_f);
                    e.ff = f_p + (e.ff - f_p) * (-ds / t_f).exp();
                    let kn = 0.25 * (1.0 + e.ff.max(0.0).sqrt()).powi(2);
                    let cv = cn_pot * (1.0 - kn);
                    // le tourbillon n'existe qu'au-delà du seuil critique ; il
                    // convecte pendant T_vl semi-cordes puis ne nourrit plus.
                    // Son horloge repart sous le seuil, ET au début d'une
                    // nouvelle MONTÉE (dα/dt qui repasse positif) : une boucle
                    // entièrement au-dessus du seuil lâche un tourbillon par
                    // cycle, comme la mesure le montre (S809 à 20° ± 5,5°).
                    let monte = *al > e.al_prev;
                    if monte && e.descend {
                        e.tau_v = 0.0;
                    }
                    e.descend = !monte;
                    e.al_prev = *al;
                    if cn_lag.abs() > p.lb_cn1 {
                        e.tau_v += ds;
                    } else {
                        e.tau_v = 0.0;
                    }
                    let nourrit = e.tau_v > 0.0 && e.tau_v < t_vl;
                    e.cnv = e.cnv * (-ds / t_v).exp()
                        + if nourrit {
                            (cv - e.cv_prev) * (-0.5 * ds / t_v).exp()
                        } else {
                            0.0
                        };
                    e.cv_prev = cv;
                }
            }
            if p.nc && h_p > 0.0 {
                // masse ajoutée : w au milieu de corde mémorisé, avec le pas qui
                // servira à la différence arrière au pas SUIVANT (le premier pas
                // n'a pas de terme : nc_h = 0, et c'est voulu)
                let st = p.stations(corps, inflows, vent);
                for (k, (_, _, wm)) in st.iter().enumerate() {
                    p.nc_w[k] = *wm;
                }
                p.nc_h = h_p;
            }
            if p.oye && h_p > 0.0 {
                // ḟ = (f_st(α) − f)/τ, τ = 4c/|V| (Øye) — exponentielle exacte sur le pas
                let st = p.stations(corps, inflows, vent);
                for (k, (al, v, _)) in st.iter().enumerate() {
                    let fs = tangent::interp::<f64>(&p.polaire.alpha, &p.fst, al.to_degrees());
                    let tau = 4.0 * p.corde / v.abs().max(1e-3);
                    p.f[k] = fs + (p.f[k] - fs) * (-h_p / tau).exp();
                }
            }
            if let Some(i) = p.inflow {
                pousse[i] += tq.0;
                // MOMENTS 1/rev autour du CENTRE DU DISQUE, projetés sur (e1, e2).
                // Ce sont eux qui alimentent λ₁s et λ₁c : sans eux Pitt–Peters
                // n'aurait que son état uniforme, c'est-à-dire rien de plus que
                // le modèle qu'il remplace.
                let inf = &inflows[i];
                if inf.pitt {
                    let bras = corps[p.corps].r - inf.centre;
                    let m_c = mo_ + bras.cross(&fo);
                    let e2 = inf.axe.cross(&inf.e1);
                    m_rt[i].0 += m_c.dot(&inf.e1);
                    m_rt[i].1 += m_c.dot(&e2);
                }
            }
            f_aero += fo;
            m_aero += mo_ + corps[p.corps].r.cross(&fo); // réduit à l'origine du monde
        }
        if let Some((t0, sf, sm, sh)) = self.moy_aero.as_mut() {
            if self.t >= *t0 {
                let h = self.h_dernier;
                *sf += f_aero * h;
                *sm += m_aero * h;
                *sh += h;
            }
        }
        self.apparie();
        let corps = &self.corps;
        // POINT LE PLUS PROCHE sur les maillages (une requête BVH par contact
        // et par pas) et, pour une CAPSULE contre boîte, cylindre ou maillage,
        // l'abscisse `s_axe` du point de l'axe le plus proche. C'est la seule
        // partie non différentiable du contact, et elle est isolée ici — le
        // résidu, lui, ne voit qu'un point fixe et une sphère de l'axe.
        for ci in 0..self.contacts.len() {
            let (p0, p1, b, pb, demi, cyl, maille, ic) = {
                let c = &self.contacts[ci];
                (c.p0, c.p1, c.b, c.pb, c.demi, c.cylindre, c.maille, c.corps)
            };
            let capsule = p1 != p0 && (demi.is_some() || cyl.is_some() || maille.is_some());
            let ca = &corps[ic];
            if let Some(mi) = maille {
                let m = &self.maillages[mi];
                let cb = &corps[m.corps];
                let l0 = cb.rot.transpose() * (ca.r + ca.rot * p0 - cb.r);
                if capsule {
                    let l1 = cb.rot.transpose() * (ca.r + ca.rot * p1 - cb.r);
                    let (_, sx, q) = m.plus_proche_segment(&l0, &l1);
                    self.contacts[ci].s_axe = sx;
                    self.contacts[ci].q_maille = q;
                } else {
                    let (_, q) = m.plus_proche(&l0);
                    self.contacts[ci].q_maille = q;
                }
            } else if capsule {
                let cb =
                    &corps[b.expect("boîte ou cylindre sans porteur : refusé à la déclaration")];
                let l0 = cb.rot.transpose() * (ca.r + ca.rot * p0 - cb.r) - pb;
                let l1 = cb.rot.transpose() * (ca.r + ca.rot * p1 - cb.r) - pb;
                let sx = if let Some(dm) = demi {
                    min_convexe_segment(|q| dist_boite(q, &dm), &l0, &l1).1
                } else {
                    let (axe, rc, lc) = cyl.unwrap();
                    min_convexe_segment(|q| dist_cylindre(q, &axe, rc, lc), &l0, &l1).1
                };
                self.contacts[ci].s_axe = sx;
            }
        }
        for ct in &mut self.contacts {
            let (_, _, _, _, d) = ct.valeur(corps);
            ct.sortie = d;
        }
        for (i, inf) in self.inflows.iter_mut().enumerate() {
            inf.poussee = pousse[i];
            if inf.impose {
                continue; // imposé de l'extérieur : rien n'évolue
            }
            if inf.pitt {
                inf.m_roul = m_rt[i].0;
                inf.m_tang = m_rt[i].1;
                pitt_peters(inf, self.vent, self.h_dernier);
                continue;
            }
            // vitesse de l'air VUE DU DISQUE, décomposée sur l'axe et dans le plan
            let vz = vent.dot(&inf.axe);
            let vp = (vent - inf.axe * vz).norm();
            let t = pousse[i].max(0.0);
            let k = 2.0 * inf.rho * inf.aire;
            let mut v = inf.v_i.max(1e-6);
            for _ in 0..40 {
                let u = ((vp * vp) + (vz - v) * (vz - v)).sqrt().max(1e-9);
                v = 0.5 * v + 0.5 * t / (k * u);
            }
            let h = self.h_dernier;
            inf.v_i += (h / inf.tau).min(1.0) * (v - inf.v_i);
        }
        Ok(())
    }

    /// La trajectoire du pas traverse-t-elle un contact ?
    ///
    /// Pour nos primitives, le test est ANALYTIQUE et il n'a rien d'une
    /// heuristique : les centres se déplacent linéairement sur le pas, donc
    /// · sphère/plan à normale fixe : la distance est AFFINE en t, ses deux
    ///   extrémités positives suffisent ;
    /// · sphère/sphère : la distance entre centres est |Δc₀ + t·Δv|, dont le
    ///   minimum sur [0,1] se calcule en fermé — c'est lui qu'il faut tester,
    ///   et non les extrémités, sinon une sphère traverse l'autre de part en
    ///   part sans qu'aucun instant échantillonné ne la voie ;
    /// · sphère contre boîte, cylindre ou maillage (7 sept.) : la sphère
    ///   balayée est une CAPSULE [pa0, pa1] de rayon r, et la distance d'un
    ///   point à un convexe est convexe le long d'un segment — son minimum est
    ///   global (`min_convexe_segment`), pas un échantillon ; le maillage se
    ///   parcourt triangle par triangle sous le BVH dilaté de r.
    fn traverse(&self, poses0: &[Pose]) -> bool {
        for ct in &self.contacts {
            let pa0 = poses0[ct.corps].0 + poses0[ct.corps].1 * ct.p0;
            let ca = &self.corps[ct.corps];
            let pa1 = ca.r + ca.rot * ct.p0;
            match ct.b {
                Some(ib) if ct.demi.is_some() || ct.cylindre.is_some() || ct.maille.is_some() => {
                    // boîte, cylindre, maillage : la sphère balayée sur le pas
                    // est une CAPSULE [pa0, pa1] de rayon r, et elle traverse
                    // si la distance de son segment à la primitive est ≤ r.
                    // Le segment est pris dans le repère du PORTEUR, chaque
                    // bout à sa pose du moment — exact si le porteur ne
                    // tourne pas sur le pas, approché sinon (même limite que
                    // le plan porté, déclarée). Jusqu'au 5 sept. ces contacts
                    // tombaient dans la branche « plan » avec la normale par
                    // défaut ; du 5 au 7, ils n'avaient pas de CCD du tout.
                    let (rb0, qb0, _) = poses0[ib];
                    let cb = &self.corps[ib];
                    // capsule : trois sphères balayées — les deux bouts de l'axe
                    // et le point figé sur le pas (`s_axe`). Déclaré approché :
                    // une tige fine peut enjamber un petit obstacle entre eux.
                    let pts: Vec<V3> = if ct.p1 != ct.p0 {
                        vec![ct.p0, ct.p1, ct.p0 + (ct.p1 - ct.p0) * ct.s_axe]
                    } else {
                        vec![ct.p0]
                    };
                    for p in &pts {
                        let pa0 = poses0[ct.corps].0 + poses0[ct.corps].1 * p;
                        let pa1 = ca.r + ca.rot * p;
                        let l0 = qb0.transpose() * (pa0 - rb0);
                        let l1 = cb.rot.transpose() * (pa1 - cb.r);
                        let touche = if let Some(mi) = ct.maille {
                            self.maillages[mi].traverse(&l0, &l1, ct.rayon)
                        } else if let Some(demi) = ct.demi {
                            min_convexe_segment(|q| dist_boite(&(q - ct.pb), &demi), &l0, &l1).0
                                <= ct.rayon
                        } else {
                            let (axe, rc, lc) = ct.cylindre.unwrap();
                            min_convexe_segment(
                                |q| dist_cylindre(&(q - ct.pb), &axe, rc, lc),
                                &l0,
                                &l1,
                            )
                            .0 <= ct.rayon
                        };
                        if touche {
                            return true;
                        }
                    }
                }
                Some(ib) if ct.rayon_b > 0.0 => {
                    let pb0 = poses0[ib].0 + poses0[ib].1 * ct.pb;
                    let cb = &self.corps[ib];
                    let pb1 = cb.r + cb.rot * ct.pb;
                    let d0 = pa0 - pb0;
                    let dv = (pa1 - pb1) - d0;
                    let n2 = dv.norm_squared();
                    let t = if n2 > 1e-30 {
                        (-d0.dot(&dv) / n2).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    if (d0 + dv * t).norm() <= ct.rayon + ct.rayon_b {
                        return true;
                    }
                }
                _ => {
                    // plan fixe ou porté : la distance est affine en t si le
                    // porteur ne tourne pas ; on teste les deux bouts, et on
                    // DÉCLARE la limite (une rotation rapide du porteur peut
                    // rendre la distance non affine — non couvert).
                    let d1 = ct.rayon - (pa1 - ct.origine).dot(&ct.normale);
                    let d0 = ct.rayon - (pa0 - ct.origine).dot(&ct.normale);
                    if d0 >= 0.0 || d1 >= 0.0 {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// APPARIEMENT automatique des sphères, par grille de hachage spatial.
    ///
    /// Le naïf teste les N(N−1)/2 paires ; la grille ne teste que les voisines
    /// des 27 cellules adjacentes, ce qui est O(N) tant que la densité reste
    /// bornée. C'est la version « broad phase » du problème — la même fin
    /// qu'une hiérarchie de volumes englobants, avec une structure plate qui
    /// se reconstruit à chaque pas au lieu de se mettre à jour.
    ///
    /// ⚠ La liste des contacts change la STRUCTURE du jacobien : l'analyse
    /// symbolique du LU creux est invalidée dès qu'une paire apparaît ou
    /// disparaît, sans quoi la factorisation suivante écrit hors motif.
    fn apparie(&mut self) {
        let Some((k, expo, c, mu, v_eps, marge, d_hat)) = self.loi_contact else {
            return;
        };
        if self.spheres.is_empty() {
            return;
        }
        let t0 = std::time::Instant::now();
        let mut pos: Vec<(usize, V3, f64)> = Vec::with_capacity(self.spheres.len());
        let mut r_max = 0.0f64;
        for (ic, p0, r) in &self.spheres {
            let cc = &self.corps[*ic];
            pos.push((*ic, cc.r + cc.rot * p0, *r));
            r_max = r_max.max(*r);
        }
        let pas = (2.0 * r_max + marge).max(1e-9);
        let mut grille: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
            std::collections::HashMap::new();
        let cell = |p: &V3, pas: f64| {
            (
                (p.x / pas).floor() as i64,
                (p.y / pas).floor() as i64,
                (p.z / pas).floor() as i64,
            )
        };
        for (i, (_, p, _)) in pos.iter().enumerate() {
            grille.entry(cell(p, pas)).or_default().push(i);
        }
        let mut paires: Vec<(usize, usize)> = vec![];
        for (i, (ic, p, r)) in pos.iter().enumerate() {
            let (cx, cy, cz) = cell(p, pas);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(v) = grille.get(&(cx + dx, cy + dy, cz + dz)) else {
                            continue;
                        };
                        for &j in v {
                            if j <= i {
                                continue;
                            }
                            let (jc, pj, rj) = &pos[j];
                            if jc == ic {
                                continue;
                            } // même corps : pas de contact
                            if (p - pj).norm() < r + rj + marge {
                                paires.push((i, j));
                            }
                        }
                    }
                }
            }
        }
        // la structure du jacobien ne change que si la liste des PAIRES change.
        // Seuls les contacts AUTOMATIQUES (en queue de liste) sont comparés et
        // remplacés : ceux déclarés à la main restent. Jusqu'au 6 sept. la liste
        // ENTIÈRE était écrasée — un contact manuel disparaissait au premier pas
        // dès que l'appariement était actif, sans un mot.
        let n_man = self
            .contacts
            .iter()
            .position(|c| c.auto)
            .unwrap_or(self.contacts.len());
        let autos = &self.contacts[n_man..];
        let change = paires.len() != autos.len()
            || paires.iter().zip(autos.iter()).any(|((i, j), ct)| {
                let (ia, pa, ra) = self.spheres[*i];
                let (ib, pb, rb) = self.spheres[*j];
                ct.corps != ia
                    || ct.b != Some(ib)
                    || ct.p0 != pa
                    || ct.pb != pb
                    || ct.rayon != ra
                    || ct.rayon_b != rb
                    || ct.k != k
                    || ct.expo != expo
                    || ct.c != c
                    || ct.mu != mu
                    || ct.v_eps != v_eps
                    || ct.d_hat != d_hat
            });
        if change {
            self.contacts.truncate(n_man);
            for (i, j) in &paires {
                let (ia, pa, ra) = self.spheres[*i];
                let (ib, pb, rb) = self.spheres[*j];
                self.contacts.push(Contact {
                    nom: format!("auto{i}-{j}"),
                    corps: ia,
                    p0: pa,
                    p1: pa,
                    rayon: ra,
                    b: Some(ib),
                    pb,
                    p1b: pb,
                    rayon_b: rb,
                    normale: V3::z(),
                    origine: V3::zeros(),
                    k,
                    expo,
                    d_hat,
                    c,
                    mu,
                    v_eps,
                    demi: None,
                    cylindre: None,
                    maille: None,
                    q_maille: V3::zeros(),
                    s_axe: 0.0,
                    cable: false,
                    sortie: (0.0, 0.0),
                    nonlisse: false,
                    row: 0,
                    fn_impose: None,
                    restitution: 0.0,
                    auto: true,
                });
            }
            // le MOTIF du jacobien change (blocs de contact) : la symbolique
            // creuse est à refaire. Le masque `actif`, lui, ne concerne que
            // les lignes de CONTRAINTE, et une paire de pénalité n'en ajoute
            // aucune. Le vider ici — jusqu'au 6 sept. — laissait `jacobien_ad`
            // l'indexer à vide au pas suivant : panique dès qu'un modèle
            // apparié portait aussi une liaison.
            self.symb = None;
        }
        self.n_paires = paires.len();
        self.t_appar += t0.elapsed().as_secs_f64();
    }

    /// Moment cinétique total autour de l'origine du monde.
    pub fn moment(&self) -> V3 {
        self.corps.iter().fold(V3::zeros(), |s, c| {
            s + c.m * c.r.cross(&c.v) + (c.rot * c.j * c.rot.transpose()) * c.w
        })
    }

    /// Le moment cinétique EST-IL un invariant de ce modèle ? Le théorème de
    /// Noether ne le donne que si le lagrangien est invariant par l'action
    /// GLOBALE de SO(3) : aucun effort extérieur, aucune liaison au bâti.
    /// Projeter en dehors de ce domaine fabriquerait une conservation FAUSSE —
    /// c'est le seul endroit où cette brique peut mentir, donc elle refuse.
    pub fn moment_invariant(&self) -> Option<String> {
        if self.g.norm() > 0.0 {
            return Some("gravité non nulle".into());
        }
        if !self.efforts.is_empty() {
            return Some("efforts extérieurs imposés".into());
        }
        if !self.pales.is_empty() {
            return Some("forces aérodynamiques".into());
        }
        for e in &self.elems {
            let (a, b) = match e {
                Elem::L(l) => (l.a, l.b),
                Elem::D(d) => (d.a, d.b),
                Elem::E(g) => (g.a, g.b),
                Elem::V(v) => (v.a, v.b),
                Elem::C(c) => (c.a, c.b),
                Elem::K(k) => (Some(k.a), k.b),
            };
            if a.is_none() || b.is_none() {
                return Some(format!("liaison « {} » au bâti", e.nom()));
            }
        }
        for c in &self.couples {
            if c.a.is_none() || c.b.is_none() {
                return Some(format!("couple « {} » réagissant au bâti", c.nom));
            }
        }
        for c in &self.contacts {
            if c.b.is_none() {
                return Some("contact contre un plan fixe".into());
            }
        }
        None
    }

    /// L'ÉNERGIE est-elle un invariant de ce modèle ? Elle l'est pour un
    /// système à liaisons parfaites et SCLÉRONOMES dont toutes les forces
    /// dérivent d'un potentiel compté par `energie()` — c'est-à-dire la
    /// gravité, et elle seule. Tout le reste ou bien dissipe (contact), ou
    /// bien injecte (servo, gouverneur, aéro, cible mobile), ou bien stocke
    /// une énergie que `energie()` ne compte pas (poutres).
    pub fn energie_invariante(&self) -> Option<String> {
        // MESURÉ, ET C'EST UN REFUS DE DOMAINE, pas de principe. Sans
        // contrainte, ∂E/∂u = (Mu)ᵀ donne δu = μ·u — COLINÉAIRE à la vitesse,
        // donc une reparamétrisation infinitésimale du temps : l'ordre est
        // préservé exactement (2,00 · 2,00 · 2,02, erreurs à 1 % de celles
        // sans projection). Dès qu'une contrainte existe, `G δu = 0` interdit
        // cette colinéarité, la correction devient TRANSVERSE à la
        // trajectoire, et l'erreur cesse de converger : plancher 1,9e-5,
        // ordre ~0 sur un pendule sphérique. Mesuré dans les deux sens.
        if !self.elems.is_empty() {
            return Some(
                "liaisons présentes : la correction d'énergie devient transverse \
                         (mesuré : ordre 2,00 → 0 sur un pendule sphérique)"
                    .into(),
            );
        }
        if !self.supers.is_empty() || self.loi_contact.is_some() {
            return Some("superéléments ou appariement : potentiel non compté".into());
        }
        if !self.poutres.is_empty() {
            return Some("énergie élastique des poutres non comptée".into());
        }
        if !self.couples.is_empty() {
            return Some("couples à loi (servo, gouverneur, ressort)".into());
        }
        if !self.pales.is_empty() {
            return Some("forces aérodynamiques".into());
        }
        if !self.contacts.is_empty() {
            return Some("contacts (dissipatifs)".into());
        }
        if !self.efforts.is_empty() {
            return Some("efforts extérieurs imposés".into());
        }
        for e in &self.elems {
            if let Elem::L(l) = e {
                if l.cible_t.is_some() || l.cible_r.is_some() {
                    return Some(format!("liaison « {} » à cible imposée (rhéonome)", l.nom));
                }
            }
        }
        None
    }

    /// PROJECTION SUR LA VARIÉTÉ DES INVARIANTS — Noether discret, par projection.
    ///
    /// L'α-généralisé n'est ni symplectique ni conservatif : il dissipe, et
    /// la dissipation agit AUSSI le long des directions engendrées par
    /// l'action du groupe de symétrie. Le moment `L = Σ mᵢrᵢ×vᵢ + Jˢᵢωᵢ` y
    /// perd une part SÉCULAIRE (mesurée : ×253 entre ρ∞ 0,9 et 0,5), et le
    /// bouton ρ∞ ne la répare pas — au-delà de 0,95 il dégrade les systèmes
    /// contraints, et Newton lâche à 0,999 (`bancs.dissipation`).
    ///
    /// On corrige donc la seule chose qui doit l'être, avec la plus petite
    /// correction possible au sens de l'énergie cinétique :
    ///
    /// ```text
    ///     min ½ δuᵀ M δu   s.c.   A δu = e   et   G δu = 0
    /// ```
    ///
    /// où `A = ∂L/∂u` (3 × n) et `e` est l'erreur de moment. La contrainte
    /// `G δu = 0` garde la correction ADMISSIBLE — sans elle on réparerait un
    /// invariant en cassant les liaisons. L'élimination du multiplicateur des
    /// liaisons laisse un système **3 × 3** :
    ///
    /// ```text
    ///     [A M⁻¹Aᵀ − A M⁻¹Gᵀ (G M⁻¹Gᵀ)⁻¹ G M⁻¹Aᵀ] μ = e ,  δu = M⁻¹(Aᵀμ + Gᵀν)
    /// ```
    ///
    /// POURQUOI ÇA NE CASSE PAS L'ORDRE, là où `projette_vitesse` le casse :
    /// la correction est de la taille de l'ERREUR LOCALE, O(h^{p+1}), et non
    /// de la taille de la solution. C'est l'hypothèse du théorème de
    /// préservation d'ordre des méthodes de projection (Hairer, Lubich &
    /// Wanner, *Geometric Numerical Integration*, IV.4). La projection de Φ̇,
    /// elle, corrige un O(h²) qui n'est pas une erreur mais une partie de la
    /// solution index-3 : d'où sa chute d'ordre 2,00 → 0,8. La même raison
    /// autorise à ne pas retoucher l'accélération algorithmique `a` : son
    /// effet sur q au pas suivant serait O(h^{p+2}).
    fn projette_invariants(&mut self, t: f64) -> Result<(), String> {
        let n = self.n();
        // On empile les invariants demandés : 3 lignes pour le moment, 1 pour
        // l'énergie. Le système final est (nl × nl) — au plus 4 × 4, quel que
        // soit le nombre de corps.
        let mut cibles: Vec<f64> = Vec::with_capacity(4);
        let mut ech = 0.0f64;
        if let Some(l0) = self.l_ref {
            let d = l0 - self.moment();
            cibles.extend([d.x, d.y, d.z]);
            ech = ech.max(1.0 + l0.norm());
        }
        if let Some(e0) = self.e_ref {
            cibles.push(e0 - self.energie());
            ech = ech.max(1.0 + e0.abs());
        }
        if cibles.is_empty() {
            return Ok(());
        }
        let nl = cibles.len();
        let e = DVector::from_vec(cibles);
        if e.norm() <= 1e-15 * ech {
            return Ok(());
        }
        // A = ∂L/∂u : ∂(m r×v)/∂v = m [r]ₓ ; ∂(Jˢω)/∂ω = Jˢ
        let lm = if self.l_ref.is_some() { 3 } else { 0 };
        let mut a = DMatrix::<f64>::zeros(nl, n);
        let mut minv = DVector::<f64>::zeros(n); // M⁻¹ de la partie translation
        let mut js_inv = Vec::with_capacity(self.corps.len());
        for (i, c) in self.corps.iter().enumerate() {
            let rx = skew(&c.r) * c.m;
            let js = c.rot * c.j * c.rot.transpose();
            for l in 0..3 {
                if lm == 3 {
                    for k in 0..3 {
                        a[(l, 6 * i + k)] = rx[(l, k)];
                        a[(l, 6 * i + 3 + k)] = js[(l, k)];
                    }
                }
                minv[6 * i + l] = 1.0 / c.m;
            }
            if self.e_ref.is_some() {
                // ∂E/∂u = (M u)ᵀ — l'énergie potentielle ne dépend que de q
                let jw = js * c.w;
                for k in 0..3 {
                    a[(lm, 6 * i + k)] = c.m * c.v[k];
                    a[(lm, 6 * i + 3 + k)] = jw[k];
                }
            }
            js_inv.push(match js.try_inverse() {
                Some(v) => v,
                None => return Ok(()),
            });
        }
        // W = M⁻¹Aᵀ  (n × 3)
        let mut w = DMatrix::<f64>::zeros(n, nl);
        for i in 0..self.corps.len() {
            for col in 0..nl {
                for k in 0..3 {
                    w[(6 * i + k, col)] = minv[6 * i + k] * a[(col, 6 * i + k)];
                }
                let mut b = V3::zeros();
                for k in 0..3 {
                    b[k] = a[(col, 6 * i + 3 + k)];
                }
                let bi = js_inv[i] * b;
                for k in 0..3 {
                    w[(6 * i + 3 + k, col)] = bi[k];
                }
            }
        }
        let (_, g) = self.phi_g(t)?;
        let mut s3 = &a * &w; // A M⁻¹Aᵀ  (nl × nl)
        let mut nu: Option<DMatrix<f64>> = None;
        if g.nrows() > 0 {
            let gw = &g * &w; // G M⁻¹Aᵀ  (m × nl)
            let mut gmg = DMatrix::<f64>::zeros(g.nrows(), g.nrows());
            for i in 0..g.nrows() {
                // G M⁻¹Gᵀ
                for j in 0..g.nrows() {
                    let mut acc = 0.0;
                    for kk in 0..self.corps.len() {
                        for k in 0..3 {
                            acc += g[(i, 6 * kk + k)] * minv[6 * kk + k] * g[(j, 6 * kk + k)];
                        }
                        let (mut bi, mut bj) = (V3::zeros(), V3::zeros());
                        for k in 0..3 {
                            bi[k] = g[(i, 6 * kk + 3 + k)];
                            bj[k] = g[(j, 6 * kk + 3 + k)];
                        }
                        acc += bi.dot(&(js_inv[kk] * bj));
                    }
                    gmg[(i, j)] = acc;
                }
            }
            // contraintes possiblement redondantes (quatre-barres) : SVD en repli
            let y = match gmg.clone().cholesky() {
                Some(ch) => ch.solve(&gw),
                None => match svd_sure(gmg)
                    .and_then(|s| s.solve(&gw, 1e-10).map_err(|e| e.to_string()))
                {
                    Ok(v) => v,
                    Err(_) => return Ok(()),
                },
            };
            s3 -= gw.transpose() * &y;
            nu = Some(-y);
        }
        let mu = match s3.clone().lu().solve(&e) {
            Some(v) => v,
            None => return Ok(()),
        };
        let mut du = &w * &mu;
        if let Some(y) = nu {
            // + M⁻¹Gᵀν
            let v = &y * &mu;
            for i in 0..self.corps.len() {
                for k in 0..3 {
                    let mut acc = 0.0;
                    for j in 0..g.nrows() {
                        acc += g[(j, 6 * i + k)] * v[j];
                    }
                    du[6 * i + k] += minv[6 * i + k] * acc;
                }
                let mut b = V3::zeros();
                for k in 0..3 {
                    let mut acc = 0.0;
                    for j in 0..g.nrows() {
                        acc += g[(j, 6 * i + 3 + k)] * v[j];
                    }
                    b[k] = acc;
                }
                let bi = js_inv[i] * b;
                for k in 0..3 {
                    du[6 * i + 3 + k] += bi[k];
                }
            }
        }
        let un = self.u() + du;
        self.set_u(&un);
        Ok(())
    }

    /// K = -∂f/∂q + ∂(Gᵀλ)/∂q, assemblée élément par élément.
    /// Les dérivées locales suivent le champ à états internes figés ; les
    /// poutres différencient leurs efforts analytiques par duaux d'ordre un.
    pub(crate) fn raideur(&self, t: f64) -> Result<DMatrix<f64>, String> {
        let k = analyse::dense(&self.raideur_locale(t)?, self.n());
        fini(k.as_slice(), "raideur tangente")?;
        Ok(k)
    }

    /// C = -∂f/∂u, par dérivées locales du même champ de forces.
    pub(crate) fn amortissement(&self, t: f64) -> Result<DMatrix<f64>, String> {
        let (_, c) = self.tangentes_forces_locales(t, false, true)?;
        let c = analyse::dense(&c, self.n());
        fini(c.as_slice(), "amortissement tangent")?;
        Ok(c)
    }

    /// Base du noyau de G : les mouvements admissibles (n × p).
    pub(crate) fn base_admissible(&self, t: f64) -> Result<DMatrix<f64>, String> {
        fini(&[t], "temps d'analyse")?;
        let (n, m) = (self.n(), self.m());
        if m == 0 {
            return Ok(DMatrix::identity(n, n));
        }
        let (_, g) = self.phi_g(t)?;
        numerique::noyau(g)
    }

    pub(crate) fn masse_spatiale(&self) -> DMatrix<f64> {
        let n = self.n();
        let mut mm = DMatrix::zeros(n, n);
        for (i, c) in self.corps.iter().enumerate() {
            let js = c.rot * c.j * c.rot.transpose();
            for a in 0..3 {
                mm[(6 * i + a, 6 * i + a)] = c.m;
                for b in 0..3 {
                    mm[(6 * i + 3 + a, 6 * i + 3 + b)] = js[(a, b)];
                }
            }
        }
        mm
    }

    /// Racines de la linéarisation existante, sans filtrer les racines réelles.
    /// Le domaine historique de modes_complexes est conservé par ce chemin
    /// interne ; spectre applique en plus son domaine mécanique explicite.
    fn racines_lineaires(&self, t: f64) -> Result<Vec<(f64, f64)>, String> {
        fini(&[t], "temps d'analyse")?;
        self.verifie_etat_fini()?;
        if self.n() == 0 {
            return Ok(vec![]);
        }
        let z = self.base_admissible(t)?;
        let p = z.ncols();
        if p == 0 {
            return Ok(vec![]);
        }
        let n = self.n();
        // CINÉMATIQUE EXACTE DE LA PERTURBATION AUTOUR D'UN ÉTAT TOURNANT (5 sept.).
        // Le refus « un corps contraint tourne » tenait depuis le 4 sept., avec
        // pour cause déclarée une « fuite de transport » ; la cause réelle est
        // que le système d'état supposait δq̇ = δu. Pour R = exp(δθ)·R₀(t) qui
        // tourne, ω = ω₀ + δθ̇ + δθ×ω₀, donc δω = δθ̇ − ω₀×δθ ; et une base Z(q)
        // qui suit le mouvement a une dérivée Ż. D'où :
        //     δq = Z ξ ,   δu = Z ξ̇ + S ξ ,   S = Ż − T_L Z ,   T_L = diag(0, [ω₀]×)
        // et le résidu projeté Zᵀ(M δu̇ + C δu + K δq) donne
        //     M̃ = ZᵀMZ ,  C̃ = Zᵀ(MS + CZ) ,  K̃ = Zᵀ(CS + KZ).
        // Ż se prend par différence finie le long de u₀, la base à q ⊕ h·u₀
        // étant ALIGNÉE sur celle de q (Z_p orthonormale ⇒ Q = Z_pᵀZ).
        // Mesuré sur la toupie rapide (Ja 2e-3, Jt 1e-3, 300 rad/s) : 0,268524 et
        // 8,412654 Hz, les deux racines EXACTES du trinôme Jt'Ω² − Jaω₃Ω + mgl,
        // σ 1e-10 — là où « Zᵀ[ω]×Z » (transport projeté) rendait 0,004 et 51,6.
        // À ω₀ = 0 tout ceci vaut zéro et la formule d'origine est retrouvée.
        let u0 = self.u();
        let ech = self.echelle_l().max(1e-9);
        let umax = u0.amax();
        let zd = if umax > 1e-12 {
            let h = 1e-6 * ech / umax;
            let mut me = self.clone();
            me.avance(&self.poses(), &u0, h);
            let zp = me.base_admissible(t)?;
            if zp.ncols() != p {
                return Err("modes_complexes : la dimension de l'espace admissible change le long du mouvement".into());
            }
            let q = zp.transpose() * &z;
            (&zp * q - &z) / h
        } else {
            DMatrix::zeros(n, p)
        };
        let mut tl = DMatrix::<f64>::zeros(n, n);
        for (i, c) in self.corps.iter().enumerate() {
            let sk = skew(&c.w);
            for a_ in 0..3 {
                for b_ in 0..3 {
                    tl[(6 * i + 3 + a_, 6 * i + 3 + b_)] = sk[(a_, b_)];
                }
            }
        }
        let s_ = &zd - &tl * &z;
        let (k, c, m) = (
            self.raideur(t)?,
            self.amortissement(t)?,
            self.masse_spatiale(),
        );
        let mm = z.transpose() * &m * &z;
        let cm = z.transpose() * (&m * &s_ + &c * &z);
        let km = z.transpose() * (&c * &s_ + &k * &z);
        fini(mm.as_slice(), "masse réduite")?;
        let ch = mm.cholesky().ok_or("masse réduite non définie positive")?;
        let mut a = DMatrix::zeros(2 * p, 2 * p);
        for i in 0..p {
            a[(i, p + i)] = 1.0;
        }
        let mk = -ch.solve(&km);
        let mc = -ch.solve(&cm);
        for i in 0..p {
            for j in 0..p {
                a[(p + i, j)] = mk[(i, j)];
                a[(p + i, p + j)] = mc[(i, j)];
            }
        }
        fini(a.as_slice(), "matrice d'état modale")?;
        let schur = nalgebra::linalg::Schur::try_new(a, f64::EPSILON, 2000 * p)
            .ok_or("modes complexes : plafond d'itérations atteint")?;
        let valeurs = schur.complex_eigenvalues();
        if valeurs
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err("modes complexes : valeurs propres non finies".into());
        }
        Ok(valeurs.iter().map(|l| (l.re, l.im)).collect())
    }

    fn limites_domaine_spectre(&self) -> Vec<&'static str> {
        let mut raisons = Vec::new();
        if !self.contacts.is_empty() {
            raisons.push("contacts présents : linéarisation complète non prise en charge");
        }
        if !self.spheres.is_empty() || self.loi_contact.is_some() {
            raisons.push("appariement automatique de contacts présent");
        }
        if !self.pales.is_empty() || !self.inflows.is_empty() {
            raisons.push(
                "éléments aérodynamiques présents : états de fluide hors du spectre mécanique",
            );
        }
        raisons
    }

    /// Spectre mécanique local complet : (partie réelle, partie imaginaire)
    /// en s⁻¹, multiplicités conservées, tri réel décroissant puis imaginaire
    /// croissant. Contacts et aérodynamique refusés ; pas de conclusion sur
    /// la stabilité globale d'une trajectoire ou d'un système périodique.
    pub fn spectre(&self, t: f64) -> Result<Vec<(f64, f64)>, String> {
        fini(&[t], "temps d'analyse")?;
        self.verifie_etat_fini()?;
        let raisons = self.limites_domaine_spectre();
        if !raisons.is_empty() {
            return Err(format!("spectre hors domaine : {}", raisons.join(" ; ")));
        }
        let mut valeurs = self.racines_lineaires(t)?;
        valeurs.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.total_cmp(&b.1)));
        Ok(valeurs)
    }

    /// Modes oscillants : (fréquence Hz, amortissement réduit ζ, σ), triés
    /// par fréquence. Les racines purement réelles restent exclues pour
    /// compatibilité ; utiliser spectre pour la linéarisation mécanique complète.
    pub fn modes_complexes(&self, t: f64, combien: usize) -> Result<Vec<(f64, f64, f64)>, String> {
        fini(&[t], "temps d'analyse")?;
        self.verifie_etat_fini()?;
        if combien == 0 {
            return Ok(vec![]);
        }
        let mut out: Vec<(f64, f64, f64)> = self
            .racines_lineaires(t)?
            .into_iter()
            .filter(|&(_, im)| im > 1e-9)
            .map(|(re, im)| {
                let module = re.hypot(im);
                (
                    im / (2.0 * std::f64::consts::PI),
                    -re / module.max(1e-300),
                    re,
                )
            })
            .collect();
        out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        out.truncate(combien);
        Ok(out)
    }

    pub fn modes(&self, t: f64, combien: usize) -> Result<Vec<(f64, Vec<f64>)>, String> {
        fini(&[t], "temps d'analyse")?;
        self.verifie_etat_fini()?;
        if combien == 0 || self.n() == 0 {
            return Ok(vec![]);
        }
        let k = self.raideur(t)?;
        let ks = (&k + k.transpose()) * 0.5;
        let z = self.base_admissible(t)?;
        if z.ncols() == 0 {
            return Ok(vec![]);
        }
        let kr = z.transpose() * &ks * &z;
        let mr = z.transpose() * self.masse_spatiale() * &z;
        // Mr SPD : Cholesky, puis problème symétrique standard L⁻¹KrL⁻ᵀ
        let ch = mr
            .clone()
            .cholesky()
            .ok_or("masse réduite non définie positive")?;
        let li = ch.l().try_inverse().ok_or("Cholesky non inversible")?;
        let a = &li * kr * li.transpose();
        let a = (&a + a.transpose()) * 0.5;
        let eig = numerique::propres_sym(a)?;
        let mut idx: Vec<usize> = (0..eig.eigenvalues.len()).collect();
        idx.sort_by(|&i, &j| eig.eigenvalues[i].partial_cmp(&eig.eigenvalues[j]).unwrap());
        let mut out = vec![];
        for &i in idx.iter() {
            let lam = eig.eigenvalues[i];
            if lam <= 1e-6 {
                continue;
            } // modes rigides et numériques
            let v = &z * (li.transpose() * eig.eigenvectors.column(i));
            out.push((
                lam.sqrt() / (2.0 * std::f64::consts::PI),
                v.iter().copied().collect(),
            ));
            if out.len() >= combien {
                break;
            }
        }
        Ok(out)
    }

    /// Multiplicateurs de chaque élément, par nom.
    ///
    /// `Err` tant qu'aucun λ n'a été calculé : `self.lam` n'est dimensionné
    /// que par un solveur, et l'indexer avant faisait PANIQUER le noyau à
    /// travers la frontière PyO3 (« Matrix index out of bounds »). Un panic
    /// est le pire message qu'un utilisateur puisse recevoir : il ne dit pas
    /// quoi faire. (Mesuré le 3 sept. en appelant `reactions()` juste après
    /// `assemble()`, ce que fait quiconque veut les efforts d'un assemblage.)
    pub fn reactions(&self) -> Result<Vec<(String, Vec<f64>)>, String> {
        let m: usize = self.elems.iter().map(|e| e.n()).sum();
        if self.lam.len() < m {
            return Err(format!(
                "reactions() : aucun multiplicateur calculé ({} attendus). \
                 Lancer `simule` ou `statique` d'abord — un assemblage résout \
                 la POSE, pas les efforts.",
                m
            ));
        }
        let mut out = vec![];
        let mut row = 0;
        for e in &self.elems {
            out.push((
                e.nom().to_string(),
                (0..e.n()).map(|i| self.lam[row + i]).collect(),
            ));
            row += e.n();
        }
        Ok(out)
    }

    pub fn energie(&self) -> f64 {
        let mut t = 0.0;
        let mut v = 0.0;
        for c in &self.corps {
            let js = c.rot * c.j * c.rot.transpose();
            t += 0.5 * c.m * c.v.dot(&c.v) + 0.5 * c.w.dot(&(js * c.w));
            v -= c.m * self.g.dot(&c.r);
        }
        t + v
    }

    /// Moindres carrés par SVD, sans panique ni blocage : une entrée non finie
    /// rend un vecteur NaN tout de suite — les appelants jugent la finitude.
    fn lstsq(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
        svd_sure(a.clone())
            .and_then(|s| s.solve(b, 1e-13).map_err(|e| e.to_string()))
            .unwrap_or_else(|_| DVector::from_element(a.ncols(), f64::NAN))
    }

    /// Résolution avec une factorisation GARDÉE : LU à pivotage partiel de
    /// faer, calculé UNE fois par jacobien (le jacobien est gardé sur le pas),
    /// puis des descentes-remontées. Mesuré sur la tête S2 : la résolution
    /// prenait 81 % du pas quand nalgebra refactorisait à chaque itération.
    /// Si la solution ne vérifie pas le système (singularité de configuration),
    /// repli SVD de norme minimale. Rend (x, repli_svd).
    /// `controle` : vérifier Jx − b contre les triplets. Ce contrôle juge la
    /// FACTORISATION (singulière ou presque), pas le second membre : une fois
    /// par factorisation suffit — la tête S2 fait dix descentes par jacobien,
    /// et le contrôle à chaque descente coûtait 8 % du pas (mesuré).
    fn resout(
        f: &Facto,
        trip: &Triplets,
        dim: usize,
        b: &DVector<f64>,
        controle: bool,
    ) -> (DVector<f64>, bool) {
        let t0 = std::time::Instant::now();
        let x = f.solve(b);
        T_SV.fetch_add(
            t0.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
        let fini = x.iter().all(|v| v.is_finite());
        if fini && !controle {
            return (x, false);
        }
        let t0 = std::time::Instant::now();
        let res = (trip_mv(trip, dim, &x) - b).norm();
        T_MV.fetch_add(
            t0.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
        if fini && res <= 1e-9 * (1.0 + b.norm()) {
            return (x, false);
        }
        let mut a = DMatrix::zeros(dim, dim);
        for t in trip {
            a[(t.row, t.col)] += t.val;
        }
        (Self::lstsq(&a, b), true)
    }

    /// Redondance STRUCTURELLE : QR à pivotage de colonnes de Gᵀ au départ ; les
    /// lignes de G qui ne sont pas dans une base des lignes sont déclarées
    /// inactives (λ = 0). Un mécanisme plan modélisé en 3D en a ; sans ça le
    /// système de Newton est singulier et chaque itération retombe sur une
    /// SVD (mesuré : 0,4 ms/pas au lieu de 0,05 sur le quatre-barres).
    fn detecte_redondance(&mut self, t: f64) -> Result<usize, String> {
        let m = self.m();
        self.actif = vec![true; m];
        if m == 0 {
            return Ok(0);
        }
        // ponytail: QR dense à pivotage, O(n·m²) — au-delà de 2 000 contraintes on
        // suppose le rang plein (les grands systèmes de ce noyau sont des chaînes et
        // des maillages, pas des mécanismes plans) ; une QR creuse à pivotage le
        // remplacera quand un banc redondant de cette taille existera.
        // SEUIL MESURÉ, pas posé. Le QR pivoté est DENSE et coûte O(m²n) : à
        // 1 800 contraintes (600 maillons) il pèse 5,9 s, soit 58,9 ms par pas
        // amortis sur cent pas — SEPT FOIS le coût du pas lui-même (8,2 ms).
        // À 3 600 contraintes le seuil le sautait déjà et le coût tombait à
        // 5,2. On le ramène là où il reste sous le pas : 900 contraintes coûtent
        // (900/1800)³ = un huitième, soit ~0,7 s. Au-delà, la redondance n'est
        // plus détectée A PRIORI — c'est le repli SVD qui la traite, plus cher
        // par pas mais correct, et `stats()` le compte.
        // ⚠ LE SEUIL SE TESTE AVANT D'ASSEMBLER. Jusqu'au 3 sept. `phi_g`
        // était appelé en tête de fonction : on payait 1,7 Go et 1,3 s
        // d'assemblage DENSE à 30 000 inconnues... pour retourner 0 juste
        // après. Un garde-fou qui coûte plus cher que ce qu'il garde.
        if m > REDONDANCE_MAX {
            return Ok(0);
        }
        let (_, g) = self.phi_g(t)?;
        let qr = g.transpose().col_piv_qr();
        let r = qr.r();
        let diag_max = (0..r.nrows().min(r.ncols()))
            .map(|i| r[(i, i)].abs())
            .fold(0.0, f64::max);
        let perm = qr.p();
        let mut rang = 0;
        for i in 0..r.nrows().min(r.ncols()) {
            if r[(i, i)].abs() > 1e-10 * diag_max.max(1e-300) {
                rang += 1;
            } else {
                break;
            }
        }
        // les colonnes de Gᵀ (= lignes de G) au-delà du rang, dans l'ordre de pivotage
        // nalgebra applique la permutation aux colonnes : on la lit sur les indices
        let mut ordre = DMatrix::<f64>::from_fn(1, m, |_, j| j as f64);
        perm.permute_columns(&mut ordre);
        for j in rang..m {
            self.actif[ordre[(0, j)] as usize] = false;
        }
        // une ligne de contact non lisse n'est jamais « redondante » : son G
        // est nul pour le QR (contact absent des analyses), elle vit par FB
        for (i, k) in self.lignes_k().iter().enumerate() {
            if *k {
                self.actif[i] = true;
            }
        }
        Ok(m - rang)
    }

    /// LE MASQUE DE REDONDANCE N'EST PAS FIGÉ. Détecté au départ, il peut être
    /// faux plus loin : un mécanisme parti d'un point de bifurcation
    /// (parallélogramme à plat) y a UNE ligne de moins qu'ailleurs, et la ligne
    /// neutralisée là devient structurelle dès qu'il en sort — Φ dérive, ou
    /// Newton cale sur un G quasi singulier. On re-détecte à la pose courante ;
    /// si le masque change, l'état d'avant le pas est RAMENÉ sur Φ = 0 et
    /// Φ̇ = 0 sous le nouveau masque (la dérive accumulée, < tol_phi, est
    /// effacée d'un coup au lieu d'être laissée à Newton) et rend `true` : le
    /// pas se rejoue.
    fn masque_rejoue(
        &mut self,
        t_detect: f64,
        t0: f64,
        poses0: &[Pose],
        u0: &DVector<f64>,
    ) -> Result<bool, String> {
        if self.m() > REDONDANCE_MAX {
            return Ok(false);
        }
        let avant = self.actif.clone();
        self.detecte_redondance(t_detect)?;
        let mut row = 0;
        for (ie, e) in self.elems.iter().enumerate() {
            if self.gel_elem.get(ie).copied().unwrap_or(false) {
                self.actif[row..row + e.n()].fill(false);
            }
            row += e.n();
        }
        if self.actif == avant {
            return Ok(false);
        }
        self.set_poses(poses0);
        self.set_u(u0);
        self.assemble(t0, 1e-12, 30)?;
        self.projette_vitesse(t0, 1e-12)?;
        Ok(true)
    }

    /// u̇₀ consistant : [M Gᵀ; G 0][u̇; λ] = [f; −γ], avec γ obtenu par AD
    /// directionnelle (géométrie et commandes). LU contrôlé, puis projection
    /// de Gauss par composante si le système est singulier. Les lignes
    /// explicitement désactivées portent λ = 0 ; les redondances restantes
    /// sont conservées dans la projection de norme minimale des réactions.
    fn acc_init(&mut self, t: f64) -> Result<(DVector<f64>, DVector<f64>), String> {
        self.acc_init_temoin(t, None)
    }

    fn acc_init_temoin(
        &mut self,
        t: f64,
        temoin: Option<&mut Option<initialisation::TemoinLineaire>>,
    ) -> Result<(DVector<f64>, DVector<f64>), String> {
        let (n, m) = (self.n(), self.m());
        let gdot_u = self.biais_acceleration(t)?;
        let mut trip: Triplets = Vec::with_capacity(36 * self.corps.len() + 12 * m);
        let tp = |r: usize, c: usize, v: f64| faer::sparse::Triplet {
            row: r,
            col: c,
            val: v,
        };
        for (i, c) in self.corps.iter().enumerate() {
            let k = 6 * i;
            let js = c.rot * c.j * c.rot.transpose();
            for a in 0..3 {
                trip.push(tp(k + a, k + a, c.m));
            }
            for a in 0..3 {
                for b in 0..3 {
                    trip.push(tp(k + 3 + a, k + 3 + b, js[(a, b)]));
                }
            }
        }
        // G en triplets, élément par élément (la boucle dense coûtait m·n visites)
        let g_temoin = {
            let capture = temoin.is_some();
            let g = self.g_creuse_pour_temoin(t, capture)?;
            for &(i, j, v) in &g {
                if !capture || self.actif.get(i).copied().unwrap_or(true) {
                    trip.push(tp(n + i, j, v));
                    trip.push(tp(j, n + i, v));
                }
            }
            // Sans témoin, libérer G avant la factorisation, comme le faisait
            // la boucle d'origine ; ne pas prolonger son allocation temporaire.
            if capture {
                g
            } else {
                Vec::new()
            }
        };
        let lki = self.lignes_k_inactives();
        for i in 0..m {
            if !self.actif.get(i).copied().unwrap_or(true) || lki[i] {
                trip.push(tp(n + i, n + i, 1.0));
            }
        }
        let mut rhs = DVector::zeros(n + m);
        rhs.rows_mut(0, n).copy_from(&self.forces_a(t));
        for i in 0..m {
            rhs[n + i] = if self.actif.get(i).copied().unwrap_or(true) {
                -gdot_u[i]
            } else {
                0.0
            };
        }
        let indices = Indices::du_modele(self);
        let mut vitesse_contact = vec![0.; m];
        let mut echelle_contact = vec![0.; m];
        if self.k_pas.iter().any(|&(actif, _)| actif) {
            for t in trip.iter().filter(|t| t.row >= n && t.col < n) {
                let v = t.val * self.corps[t.col / 6].u6(t.col % 6);
                vitesse_contact[t.row - n] += v;
                echelle_contact[t.row - n] += v.abs();
            }
        }
        let contacts: Vec<_> = self
            .contacts
            .iter()
            .enumerate()
            .filter(|(i, ct)| {
                ct.nonlisse
                    && self.actif.get(ct.row).copied().unwrap_or(true)
                    && matches!(self.k_pas.get(*i), Some((true, _)))
            })
            .filter_map(|(_, ct)| {
                let row = n + ct.row;
                let local = match &indices {
                    Some(i) => i.global.binary_search(&row).ok()?,
                    None => row,
                };
                // δ = -gap. Une vitesse δ̇ < 0 indique un décollement ;
                // le redémarrage lisse ne doit pas imposer δ̈ = 0.
                Some((
                    local,
                    vitesse_contact[ct.row] >= -64. * f64::EPSILON * echelle_contact[ct.row],
                ))
            })
            .collect();
        let mut impose = rhs.clone();
        if indices.is_some() {
            // Les corps gelés n'ont aucune contrainte active. Leurs blocs de
            // masse sont indépendants : les résoudre préserve exactement
            // l'initialisation globale, y compris ses composantes non nulles.
            for (i, c) in self.corps.iter().enumerate() {
                if !self.gele(i) {
                    continue;
                }
                let k = 6 * i;
                impose.rows_mut(k, 3).scale_mut(1.0 / c.m);
                let js = c.rot * c.j * c.rot.transpose();
                let b = V3::new(rhs[k + 3], rhs[k + 4], rhs[k + 5]);
                let w = js.lu().solve(&b).ok_or("masse gelée singulière")?;
                impose.rows_mut(k + 3, 3).copy_from(&w);
            }
        }
        let (reduction, trip) = if let Some(indices) = indices {
            let (r, trip) = Reduction::nouvelle(indices, trip);
            (Some(r), trip)
        } else {
            (None, trip)
        };
        let rhs = reduction
            .as_ref()
            .map_or_else(|| rhs.clone(), |r| r.rhs(&rhs, &impose));
        let physiques = reduction
            .as_ref()
            .map_or(n, |r| r.indices.global.partition_point(|&i| i < n));
        if temoin.is_some() && (reduction.is_some() || !contacts.is_empty()) {
            return Err(
                "certificat initialisation : réduction ou contact unilatéral hors périmètre".into(),
            );
        }
        let x = initialisation::complementarite(&trip, physiques, &rhs, &contacts)?;
        if let Some(out) = temoin {
            *out = Some(initialisation::TemoinLineaire {
                a: trip,
                b: rhs,
                x: x.clone(),
                g_complet: g_temoin,
                c_complet: -gdot_u,
                actives: (0..m)
                    .map(|i| self.actif.get(i).copied().unwrap_or(true))
                    .collect(),
            });
        }
        let x = reduction
            .as_ref()
            .map_or_else(|| x.clone(), |r| r.etend(&x, impose));
        if x.iter().any(|v| !v.is_finite()) {
            return Err(
                "accélération initiale non finie : masse, inertie, gravité ou force mal posée"
                    .into(),
            );
        }
        Ok((x.rows(0, n).into_owned(), x.rows(n, m).into_owned()))
    }

    /// SCHÉMA ÉNERGIE-MOMENT DANS LA FORMULATION (Simo–Wong 1991, forme de
    /// Cayley au point milieu) — ce que la PROJECTION ne peut pas faire.
    ///
    /// `projette_invariants` conserve L *ou* E, jamais les deux : pour un
    /// corps rigide ∂E/∂ω = ωᵀ·∂L/∂ω, la pile est de rang 3 sur 4 et corriger
    /// les seules vitesses ne peut pas tenir les deux (mesuré, |L| ×10⁴ pire).
    /// Ici les deux sont conservés EXACTEMENT, par construction du pas :
    ///
    /// ```text
    ///   corps, moment cinétique matériel Π = J ω_b :
    ///     (Π₁ − Π₀)/h = Π̄ × J⁻¹Π̄ ,   Π̄ = (Π₀ + Π₁)/2 ,  ω̄ = J⁻¹Π̄      (Newton 3×3)
    ///     R₁ = R₀ · cay(h ω̄)  ,   cay(x̂) = (I − x̂/2)⁻¹ (I + x̂/2)
    ///   translation, potentiel linéaire (gravité) :
    ///     v₁ = v₀ + h g ,   r₁ = r₀ + h (v₀ + v₁)/2
    /// ```
    ///
    /// Pourquoi c'est exact, et pourquoi c'est CAYLEY et non exp : l'énergie de
    /// rotation ½ΠᵀJ⁻¹Π est quadratique, donc son gradient discret (Gonzalez)
    /// EST la valeur au point milieu J⁻¹Π̄ — multiplier l'équation par ω̄ donne
    /// E₁ − E₀ = 0 identiquement ; et l'équation du moment s'écrit
    /// (I + hω̄̂/2)Π₁ = (I − hω̄̂/2)Π₀, soit Π₁ = cay(−hω̄̂)Π₀ : avec R₁ = R₀cay(hω̄̂)
    /// le moment SPATIAL R₁Π₁ = R₀Π₀ est conservé en arithmétique exacte. Avec exp à la
    /// place de cay, cette identité n'est vraie qu'à O(h³). En translation le
    /// point milieu est exact pour un potentiel linéaire. Ordre 2.
    ///
    /// DOMAINE, refusé sinon (`energie_invariante`) : corps libres sous la seule
    /// gravité. Les liaisons demandent le gradient discret de Φ (Gonzalez 1999,
    /// Betsch–Steinmann 2001), les forces non potentielles n'ont pas
    /// d'invariant à tenir — non fait, nommé.
    pub fn simule_em<F: FnMut(f64, &Modele)>(
        &mut self,
        t_end: f64,
        h: f64,
        mut sortie: F,
    ) -> Result<usize, String> {
        reglages(self.t, t_end, h, 1e-12, 30)?;
        self.verifie_etat_fini()?;
        self.verifie_options()?;
        let mut sauve = Sauvegarde::new(self);
        let res = self.simule_em_interne(t_end, h, |t, m, _, initial| {
            if !initial {
                sauve.actualise(m);
                sortie(t, m);
            }
        });
        if res.is_err() {
            sauve.restaure(self);
        }
        res
    }

    fn simule_em_interne<F: FnMut(f64, &Modele, &[V3], bool)>(
        &mut self,
        t_end: f64,
        h: f64,
        mut sortie: F,
    ) -> Result<usize, String> {
        if let Some(c) = self.energie_invariante() {
            return Err(format!("énergie-moment : hors domaine ({c})"));
        }
        if !self.elems.is_empty() {
            return Err(
                "énergie-moment : liaisons présentes (gradient discret de Φ non fait)".into(),
            );
        }
        let mut npas = 0;
        // LE MOMENT MATÉRIEL EST PORTÉ D'UN PAS À L'AUTRE, pas relu de (R, ω).
        // Relu, il fait l'aller-retour J·Rᵀ·R·J⁻¹·Π à chaque pas, et R n'est
        // orthonormale qu'à N·ε près : la dérive de E et L montait en N²
        // (mesuré : 4e-13 à 2 000 pas, 7,5e-11 à 20 000, 5,6e-9 à 200 000).
        // Porté, Π₁ = cay(−hω̄̂)Π₀ évite cet aller-retour ; R est
        // re-orthonormalisée par projection polaire. Le Newton et les produits
        // flottants gardent une erreur d'arrondi, mesurée sur les invariants.
        let mut pis: Vec<V3> = self
            .corps
            .iter()
            .map(|c| c.j * (c.rot.transpose() * c.w))
            .collect();
        sortie(self.t, self, &pis, true);
        while self.t < t_end {
            let t0 = self.t;
            let t1 = numerique::fin_pas(t0, t_end, h, true)?;
            let hh = t1 - t0;
            for (ic, c) in self.corps.iter_mut().enumerate() {
                // translation : point milieu, exact pour un potentiel linéaire
                let v1 = c.v + hh * self.g;
                c.deplace(0.5 * hh * (c.v + v1));
                c.v = v1;
                // rotation : Newton sur Π₁ (repère matériel)
                let jinv = match c.j.try_inverse() {
                    Some(v) => v,
                    None => return Err("énergie-moment : inertie singulière".into()),
                };
                let p0 = pis[ic];
                let mut p1 = p0;
                let mut ok = false;
                let mut f_prec = f64::INFINITY;
                let mut meilleur = p0;
                for _ in 0..30 {
                    let pm = 0.5 * (p0 + p1);
                    let f = p1 - p0 - hh * pm.cross(&(jinv * pm));
                    // Résidu RELATIF au moment, sans plancher dimensionné 1.
                    // Un seuil absolu 1e-17 acceptait p1=p0 dès le premier
                    // tour pour les petites inerties : le solide tournait
                    // alors à vitesse matérielle constante, hors axes propres.
                    // La norme infinie évite le carré d'un petit/grand moment.
                    let scale = p0.amax().max(p1.amax());
                    if !f.iter().all(|x| x.is_finite()) || !scale.is_finite() {
                        break;
                    }
                    let fn_ = if scale == 0.0 {
                        f.amax()
                    } else {
                        f.amax() / scale
                    };
                    if !fn_.is_finite() {
                        break;
                    }
                    if fn_ <= f64::EPSILON {
                        ok = true;
                        break;
                    }
                    if fn_ >= f_prec {
                        if f_prec <= 32.0 * f64::EPSILON {
                            p1 = meilleur;
                            ok = true;
                        }
                        break;
                    }
                    meilleur = p1;
                    f_prec = fn_;
                    // ∂f/∂Π₁ = I − (h/2)(−[J⁻¹Π̄]× + [Π̄]× J⁻¹)
                    let jf = M3::identity() - 0.5 * hh * (-skew(&(jinv * pm)) + skew(&pm) * jinv);
                    let dx = match jf.try_inverse() {
                        Some(ji) => ji * f,
                        None => return Err("énergie-moment : Newton singulier".into()),
                    };
                    p1 -= dx;
                }
                if !ok {
                    return Err(format!(
                        "énergie-moment : Newton ne converge pas à t={t1:.5}"
                    ));
                }
                let pm = 0.5 * (p0 + p1);
                let x = skew(&(hh * (jinv * pm)));
                let cay = (M3::identity() - 0.5 * x)
                    .try_inverse()
                    .ok_or("énergie-moment : Cayley singulier")?
                    * (M3::identity() + 0.5 * x);
                c.rot *= cay;
                // Projection polaire, sans donner priorité à une colonne.
                c.rot = rotation::projette(&c.rot)?;
                pis[ic] = p1;
                c.w = c.rot * (jinv * p1);
            }
            self.verifie_etat_fini()?;
            if !self.energie().is_finite() {
                return Err("énergie-moment : énergie non finie".into());
            }
            self.t = t1;
            self.h_dernier = hh;
            npas += 1;
            sortie(t1, self, &pis, false);
        }
        Ok(npas)
    }

    /// AUDIT DU JACOBIEN : le jacobien AD de Newton contre une différence
    /// finie centrée du MÊME résidu, au pas `h`, à l'état courant, bloc par
    /// bloc — (dynamique | contraintes | GGL) × (u̇ | λ | ζ). C'est
    /// l'instrument qui aurait attrapé le double compte gyroscopique du
    /// 5 sept. s'il avait porté sur `k_c_m_z` ; ici il porte sur ce qui fait
    /// converger Newton. Ce qu'il rend : par bloc, l'écart max relatif au max
    /// du bloc de référence, et les dix pires coefficients (i, j, AD, DF).
    /// Un terme LAISSÉ DE CÔTÉ sciemment (raideur de contact λ·∂G/∂q,
    /// (∂G/∂q·u) des lignes en vitesse) s'y lit comme un écart : c'est le but.
    #[allow(clippy::type_complexity)]
    pub fn audit_jacobien(
        &mut self,
        h: f64,
        rho: f64,
        ggl: bool,
    ) -> Result<(Vec<(String, f64, f64)>, Vec<(usize, usize, f64, f64)>), String> {
        reglages(self.t, self.t + h, h, 1e-12, 1)?;
        if !(0.0..1.0).contains(&rho) {
            return Err("rho hors [0,1[".into());
        }
        self.verifie_options()?;
        self.verifie_etat_fini()?;
        let sauve = Sauvegarde::new(self);
        let ggl0 = self.ggl;
        let r = self.audit_jacobien_interne(h, rho, ggl);
        sauve.restaure(self);
        self.ggl = ggl0;
        r
    }
    fn audit_jacobien_interne(
        &mut self,
        h: f64,
        rho: f64,
        ggl: bool,
    ) -> Result<(Vec<(String, f64, f64)>, Vec<(usize, usize, f64, f64)>), String> {
        let am = (2.0 * rho - 1.0) / (rho + 1.0);
        let af = rho / (rho + 1.0);
        let gam = 0.5 + af - am;
        let bet = 0.25 * (gam + 0.5) * (gam + 0.5);
        self.ggl = ggl;
        let (n, m) = (self.n(), self.m());
        if self.actif.len() != m {
            self.detecte_redondance(self.t)?;
        }
        let t0 = self.t;
        let t1 = t0 + h;
        let (udot, lam) = self.acc_init(t0)?;
        let nz = if ggl { m } else { 0 };
        let dim = n + m + nz;
        let ech_l = self.echelle_l();
        // activation des contacts non lisses, comme dans `simule`
        let u_now = self.u();
        self.k_pas = self
            .contacts
            .iter()
            .map(|ct| {
                if !ct.nonlisse {
                    return (false, 0.0);
                }
                let (phi, g) = ct.phi_g(&self.corps);
                let gap0 = -phi[0];
                let mut gdot0 = 0.0;
                for (kb, cb) in [Some(ct.corps), ct.b].iter().enumerate() {
                    if let Some(ib) = cb {
                        for c in 0..6 {
                            gdot0 -= g[(0, 6 * kb + c)] * u_now[6 * ib + c];
                        }
                    }
                }
                (
                    gap0 + h * gdot0 <= 0.0 || gap0 <= 1e-6 * ech_l.max(1.0),
                    ct.restitution * gdot0 + 0.2 * gap0.min(0.0) / h,
                )
            })
            .collect();
        let gglp = if ggl {
            let d = 1e-7_f64.max(1e-9 * t1.abs());
            let phi_t = (self.phi_seul(t1 + d)? - self.phi_seul(t1)?) / d;
            Some((self.g_creuse(t0)?, phi_t))
        } else {
            None
        };
        let pas = Pas {
            am,
            af,
            gam,
            bet,
            sigma: self.sigma_lie,
            hh: h,
            t1,
            u0: self.u(),
            poses0: self.poses(),
            udot: udot.clone(),
            a: udot.clone(),
            ggl: gglp,
        };
        let fnorm = self.forces().norm() + self.echelle_inertie(&udot);
        self.s_fb = ech_l / (1.0 + fnorm);
        let mut x = DVector::zeros(dim);
        x.rows_mut(0, n).copy_from(&udot);
        x.rows_mut(n, m).copy_from(&lam);
        let r0 = pas.residu(self, &x);
        let (trip, _) = self.jacobien_ad(&pas, &x, &pas.dtheta(&x, n))?;
        let mut jad = DMatrix::<f64>::zeros(dim, dim);
        for t in &trip {
            jad[(t.row, t.col)] += t.val;
        }
        let mut jdf = DMatrix::<f64>::zeros(dim, dim);
        // le pas de la différence finie se cale sur ce qu'il DÉPLACE : u̇ et ζ
        // n'agissent sur Φ que par la pose, à βh²c et βh² près — 1e-6 sur u̇
        // déplaçait la pose de 1e-13 et la colonne de contrainte était du bruit
        // à 1e-3 (mesuré). Même leçon que les colonnes ζ de GGL du 5 sept.
        let kq = bet * h * h * (1.0 - af) / (1.0 - am);
        let el = ech_l.max(1.0);
        for j in 0..dim {
            let e = if j < n {
                1e-7 * el / kq
            } else if j < n + m {
                1e-6 * x[j].abs().max(1.0)
            } else {
                1e-7 * el / (bet * h * h)
            };
            let mut xp = x.clone();
            xp[j] += e;
            let rp = pas.residu(self, &xp);
            let mut xm = x.clone();
            xm[j] -= e;
            let rm = pas.residu(self, &xm);
            for i in 0..dim {
                jdf[(i, j)] = (rp[i] - rm[i]) / (2.0 * e);
            }
        }
        let _ = r0;
        self.set_poses(&pas.poses0);
        self.set_u(&pas.u0);
        self.k_pas = vec![];
        let blocs_l = [("dyn", 0, n), ("con", n, n + m), ("ggl", n + m, dim)];
        let blocs_c = [("udot", 0, n), ("lam", n, n + m), ("zeta", n + m, dim)];
        let mut out = vec![];
        let mut pires: Vec<(usize, usize, f64, f64)> = vec![];
        for (nl, i0, i1) in blocs_l {
            for (nc, j0, j1) in blocs_c {
                if i1 <= i0 || j1 <= j0 {
                    continue;
                }
                let mut ref_max = 0.0f64;
                let mut err = 0.0f64;
                for i in i0..i1 {
                    for j in j0..j1 {
                        ref_max = ref_max.max(jdf[(i, j)].abs()).max(jad[(i, j)].abs());
                    }
                }
                for i in i0..i1 {
                    for j in j0..j1 {
                        let d = (jad[(i, j)] - jdf[(i, j)]).abs();
                        err = err.max(d);
                        if d > 1e-4 * ref_max.max(1e-300) {
                            pires.push((i, j, jad[(i, j)], jdf[(i, j)]));
                        }
                    }
                }
                out.push((
                    format!("{nl}/{nc}"),
                    if ref_max > 0.0 { err / ref_max } else { 0.0 },
                    ref_max,
                ));
            }
        }
        pires.sort_by(|a, b| ((b.2 - b.3).abs()).partial_cmp(&(a.2 - a.3).abs()).unwrap());
        pires.truncate(12);
        Ok((out, pires))
    }

    /// SOUS-CYCLAGE MULTI-RYTHME — deux partitions de corps, deux pas.
    ///
    /// Les corps `rapides` avancent à h/k, les autres à h ; le couplage entre
    /// partitions passe PAR LES FORCES SEULEMENT (poutres, couples, contacts
    /// par pénalité) — une liaison qui enjambe est REFUSÉE, parce que sa
    /// contrainte devrait être tenue aux deux rythmes à la fois. Ordre des
    /// opérations, « slowest-first » (Gear–Wells 1984) :
    ///
    /// ```text
    ///   1. pas LENT [t₀, t₀+h] : les corps rapides sont gelés sur leur
    ///      extrapolation de pose à vitesse constante depuis t₀ ;
    ///   2. k pas RAPIDES : les corps lents sont gelés sur l'INTERPOLATION de
    ///      leur trajectoire entre t₀ et t₀+h (géodésique sur SO(3), linéaire
    ///      en r et en u) ;
    ///   3. à t₀+h, tout le monde est au même instant.
    /// ```
    ///
    /// Chaque partition conserve ses accélérations, ses multiplicateurs et
    /// son cache de jacobien. Newton et l'accélération initiale factorisent
    /// uniquement les corps mobiles et leurs contraintes actives. Les lignes
    /// éliminées sont résolues explicitement et reportées dans le second membre.
    /// Le cache reste soumis au contrôle de descente de Newton ; les états et
    /// sorties restent globaux. Les réactions lentes sont restituées au
    /// macro-pas depuis leur propre schéma. Pas de correcteur de couplage.
    ///
    /// DOMAINE : pas d'aéro (les états de fluide avanceraient aux deux
    /// rythmes), pas de contact non lisse, pas de GGL ni d'adaptatif.
    #[allow(clippy::too_many_arguments)]
    pub fn simule_multi<F: FnMut(f64, &Modele)>(
        &mut self,
        t_end: f64,
        h: f64,
        k: usize,
        rapides: &[bool],
        rho: f64,
        tol: f64,
        newton_max: usize,
        mut sortie: F,
    ) -> Result<usize, String> {
        reglages(self.t, t_end, h, tol, newton_max)?;
        if !(0.0..1.0).contains(&rho) {
            return Err("rho hors [0,1[".into());
        }
        self.verifie_etat_fini()?;
        self.verifie_options()?;
        let mut sauve = Sauvegarde::new(self);
        let res = self.simule_multi_interne(t_end, h, k, rapides, rho, tol, newton_max, |t, m| {
            sauve.actualise(m);
            sortie(t, m);
        });
        if res.is_err() {
            sauve.restaure(self);
        }
        res
    }

    fn simule_multi_interne<F: FnMut(f64, &Modele)>(
        &mut self,
        t_end: f64,
        h: f64,
        k: usize,
        rapides: &[bool],
        rho: f64,
        tol: f64,
        newton_max: usize,
        mut sortie: F,
    ) -> Result<usize, String> {
        let nb = self.corps.len();
        if rapides.len() != nb || k == 0 {
            return Err("multi-rythme : partition ou k invalides".into());
        }
        if self.sigma_lie != 0.0 {
            return Err("multi-rythme : cinématique sigma non prise en charge".into());
        }
        if !self.pales.is_empty() || !self.inflows.is_empty() {
            return Err(
                "multi-rythme : pas d'aérodynamique (états de fluide à deux rythmes)".into(),
            );
        }
        if self.contacts.iter().any(|c| c.nonlisse) || self.ggl || self.adapt.is_some() {
            return Err("multi-rythme : hors domaine (contact non lisse, GGL ou adaptatif)".into());
        }
        // partition des éléments : tout d'un côté, ou refus
        let mut elem_rapide = Vec::with_capacity(self.elems.len());
        for e in &self.elems {
            let cs: Vec<usize> = e.corps().into_iter().flatten().collect();
            let nr = cs.iter().filter(|&&i| rapides[i]).count();
            if nr != 0 && nr != cs.len() {
                return Err(format!(
                    "multi-rythme : la liaison « {} » enjambe les deux partitions — \
                     le couplage ne passe que par les forces",
                    e.nom()
                ));
            }
            elem_rapide.push(nr > 0);
        }
        let m = self.m();
        if self.actif.len() != m {
            self.detecte_redondance(self.t)?;
        }
        let mut actif0 = self.actif.clone();
        let lignes_de = |elems: &[Elem], quels: &dyn Fn(usize) -> bool| -> Vec<bool> {
            let mut v = Vec::with_capacity(m);
            for (ie, e) in elems.iter().enumerate() {
                for _ in 0..e.n() {
                    v.push(quels(ie));
                }
            }
            v
        };
        let lignes_rapides = lignes_de(&self.elems, &|ie| elem_rapide[ie]);
        let lents: Vec<bool> = rapides.iter().map(|r| !r).collect();
        let elem_lent: Vec<bool> = elem_rapide.iter().map(|r| !r).collect();
        // l'état du schéma (u̇, a, λ) de CHAQUE partition survit d'un sous-pas
        // à l'autre : c'est ce qui garde l'ordre 2 (cf. `schema_init`)
        let mut sch_lent: Option<(DVector<f64>, DVector<f64>, DVector<f64>)> = None;
        let mut sch_rap: Option<(DVector<f64>, DVector<f64>, DVector<f64>)> = None;
        let mut jac_lent = None;
        let mut jac_rapide = None;
        let mut npas = 0;
        while self.t < t_end {
            let t0 = self.t;
            let t1 = numerique::fin_pas(t0, t_end, h, true)?;
            let hh = t1 - t0;
            let (q0, u0) = (self.poses(), self.u());
            // 1. pas LENT, rapides gelés sur leur EXTRAPOLATION d'ordre 1
            // (q₀ ⊕ h·u₀) : gelés à t₀ (ordre 0), le bout d'une chaîne de 2 m
            // sortait 2 cm à côté de la référence en 0,3 s — le montage
            // souple voyait en permanence un partenaire en retard d'un pas.
            for i in 0..nb {
                if rapides[i] {
                    let c = &mut self.corps[i];
                    c.deplace(hh * c.v);
                    c.rot = expm(&(hh * c.w)) * c.rot;
                }
            }
            self.gel = rapides.to_vec();
            self.gel_elem = elem_rapide.clone();
            let gl = lignes_de(&self.elems, &|ie| elem_rapide[ie]);
            for i in 0..m {
                self.actif[i] = actif0[i] && !gl[i];
            }
            self.schema_init = sch_lent.take();
            let r = self.simule_interne(t1, hh, rho, tol, newton_max, &mut jac_lent, |_, _| {});
            r?;
            for i in 0..m {
                if !lignes_rapides[i] {
                    actif0[i] = self.actif[i];
                }
            }
            self.actif.clone_from(&actif0);
            sch_lent = self.schema_fin.take();
            let (q1, u1) = (self.poses(), self.u());
            // 2. pas RAPIDES, lents interpolés entre t₀ et t₁ ; tout le monde
            // repart de t₀ (les rapides de leur vrai état, pas de l'extrapolé)
            self.gel = vec![];
            self.set_poses(&q0);
            self.set_u(&u0);
            self.t = t0;
            self.gel = lents.clone();
            self.gel_elem = elem_lent.clone();
            let gl = lignes_de(&self.elems, &|ie| elem_lent[ie]);
            for i in 0..m {
                self.actif[i] = actif0[i] && !gl[i];
            }
            let rotations_lentes: Vec<_> = (0..nb)
                .map(|i| {
                    if rapides[i] {
                        V3::zeros()
                    } else {
                        log_so3(&(q1[i].1 * q0[i].1.transpose()))
                    }
                })
                .collect();
            let hk = hh / k as f64;
            let mut res = Ok(0);
            for j in 1..=k {
                let s = j as f64 / k as f64;
                let tj = if j == k { t1 } else { t0 + s * hh };
                for i in 0..nb {
                    if rapides[i] {
                        continue;
                    }
                    let th = rotations_lentes[i];
                    let (d, db) = numerique::deplace_compense(q1[i].0, q1[i].2, -q0[i].0);
                    self.corps[i].deplace_depuis(q0[i].0, q0[i].2, s * (d + (db - q0[i].2)));
                    self.corps[i].rot = expm(&(s * th)) * q0[i].1;
                    for c in 0..3 {
                        self.corps[i].v[c] = u0[6 * i + c] + s * (u1[6 * i + c] - u0[6 * i + c]);
                        self.corps[i].w[c] =
                            u0[6 * i + 3 + c] + s * (u1[6 * i + 3 + c] - u0[6 * i + 3 + c]);
                    }
                }
                self.schema_init = sch_rap.take();
                res = self.simule_interne(tj, hk, rho, tol, newton_max, &mut jac_rapide, |_, _| {});
                if res.is_err() {
                    break;
                }
                sch_rap = self.schema_fin.take();
            }
            res?;
            for i in 0..m {
                if lignes_rapides[i] {
                    actif0[i] = self.actif[i];
                }
                // Le dernier micro-pas porte zéro sur les lignes lentes :
                // restituer les réactions de LEUR pas, pas ces zéros de gel.
                let schema = if lignes_rapides[i] {
                    &sch_rap
                } else {
                    &sch_lent
                };
                self.lam[i] = schema.as_ref().expect("partition intégrée").2[i];
            }
            self.actif.clone_from(&actif0);
            self.gel.clear();
            self.gel_elem.clear();
            self.t = t1;
            npas += 1;
            sortie(t1, self);
        }
        Ok(npas)
    }

    /// Intègre de `self.t` à `t_end` au pas `h` (α-généralisé, index 3, Lie).
    /// `sortie` reçoit (t, modèle) après chaque pas. Rend le nombre de pas.
    ///
    /// MULTITHREAD PAR DÉFAUT : les blocs d'éléments du jacobien AD se
    /// calculent en parallèle (rayon, `RAYON_NUM_THREADS` pour borner).
    pub fn simule<F: FnMut(f64, &Modele)>(
        &mut self,
        t_end: f64,
        h: f64,
        rho: f64,
        tol: f64,
        newton_max: usize,
        mut sortie: F,
    ) -> Result<usize, String> {
        reglages(self.t, t_end, h, tol, newton_max)?;
        self.verifie_etat_fini()?;
        self.verifie_options()?;
        let mut sauve = Sauvegarde::new(self);
        let res = self.simule_interne(t_end, h, rho, tol, newton_max, &mut None, |t, m| {
            sauve.actualise(m);
            sortie(t, m);
        });
        if res.is_err() {
            sauve.restaure(self);
        }
        res
    }

    fn simule_interne<F: FnMut(f64, &Modele)>(
        &mut self,
        t_end: f64,
        h: f64,
        rho: f64,
        tol: f64,
        newton_max: usize,
        jac_garde: &mut Option<JacGarde>,
        mut sortie: F,
    ) -> Result<usize, String> {
        if !(h > 0.0) || !h.is_finite() {
            return Err("le pas doit être fini et strictement positif".into());
        }
        if !self.t.is_finite() || !t_end.is_finite() || t_end < self.t {
            return Err("temps courant et final finis, fin ≥ temps courant requis".into());
        }
        // ρ∞ = 1 (trapèze, aucune dissipation) est INSTABLE sur une DAE d'index 3 :
        // mesuré le 3 sept. — pendule à 30°, h = T/400, Newton diverge à t = 16 s
        // (résidu 1,9e5) après des oscillations croissantes du multiplicateur.
        // C'est le résultat classique (Cardona–Géradin 1989) qui a motivé HHT-α
        // puis l'α-généralisé : il FAUT de la dissipation haute fréquence.
        if !(0.0..1.0).contains(&rho) {
            return Err(format!(
                "rho = {rho} : ρ∞ doit être dans [0, 1[ — ρ∞ = 1 est instable en index 3"
            ));
        }
        let am = (2.0 * rho - 1.0) / (rho + 1.0);
        let af = rho / (rho + 1.0);
        let gam = 0.5 + af - am;
        let bet = 0.25 * (gam + 0.5) * (gam + 0.5);
        let (n, m) = (self.n(), self.m());
        if self.actif.len() != m {
            self.detecte_redondance(self.t)?;
        }
        let (mut udot, mut lam, a_init) = match self.schema_init.take() {
            Some((u, a, l)) if u.len() == n && a.len() == n && l.len() == m => (u, l, Some(a)),
            _ => {
                let (u, l) = self.acc_init(self.t)?;
                (u, l, None)
            }
        };
        let nz = if self.ggl { m } else { 0 };
        let dim = n + m + nz;
        let mut ze = DVector::zeros(nz);
        let ech_l = self.echelle_l();
        // L'accélération physique et l'accélération algorithmique initiale
        // doivent être consistantes. Borner udot en fonction de h modifiait
        // la chute libre et l'énergie d'un pendule en silence. Les difficultés
        // de Newton se traitent par les reprises et subdivisions du pas,
        // sans altérer l'état initial qui définit le problème mécanique.
        let mut a = a_init.unwrap_or_else(|| udot.clone());
        if let Some(h) = self.hist_schema.as_mut() {
            h.push((
                udot.as_slice().to_vec(),
                a.as_slice().to_vec(),
                lam.as_slice().to_vec(),
            ));
        }
        let mut npas = 0;
        // PAS ADAPTATIF (opt-in) : résidu de demi-pas, Hibbitt & Karlsson 1979.
        // L'équilibre est résolu à t+h, ce qui ne dit RIEN de l'équilibre à
        // l'intérieur du pas. On interpole l'état à t+h/2 en supposant
        // l'accélération LINÉAIRE sur le pas — la base même de Newmark — et on
        // y évalue le résidu d'équilibre. Grand devant les forces du problème,
        // le pas est trop long : on le REJOUE à h/2. Petit, on allonge le
        // suivant. Le pas fixe reste le DÉFAUT (`adapt = None`), comme chez
        // ABAQUS qui offre les deux.
        let mut h_cur = h;
        let (h_min, h_max) = (h * self.adapt_bornes.0, h * self.adapt_bornes.1);
        let mut n_rejeu = 0usize;
        self.r_half_max = 0.0;
        self.phi_half_max = 0.0;
        self.seuil_max = 0.0;
        // NEWTON MODIFIÉ D'UN PAS À L'AUTRE : le jacobien factorisé est GARDÉ
        // au-delà du pas, avec le même juge qu'à l'intérieur du pas (le résidu
        // doit chuter d'au moins 10× par itération, sinon rafraîchi). Mesuré le
        // 7 sept. : jacobien + factorisation pesaient 81 % du pas sur Princeton
        // (poutres, tangente par 24 passages de duaux) et 84 % sur la tête S2
        // (LU dense 143). Une itération de plus coûte un résidu et une descente,
        // un jacobien en coûte dix. Le jacobien porte h et dim : il tombe si
        // l'un des deux change.
        // LE JACOBIEN TOURNE AVEC LES CORPS. Les lignes de Φ sont matérielles
        // (d = uaᵀ·dw, phi_r dans le repère de a) : elles ne tournent pas ; les
        // colonnes (δr, δθ, ẇ spatiaux) et les lignes dynamiques d'un corps tournent
        // avec lui, Qᵢ = Rᵢ·Rᵢ⁰ᵀ depuis la pose d'assemblage. Donc J(q) ≈ T_r·J₀·T_cᵀ,
        // EXACT pour tout élément entre corps co-rotatifs ou vers le bâti, approché
        // sur une liaison entre un corps qui tourne et un qui ne tourne pas — et
        // J⁻¹ = T_c·J₀⁻¹·T_rᵀ : le LU gardé sert tel quel, on tourne b et x. Sur une
        // tête de rotor le jacobien mourait tous les deux pas (1,9° de rotation par
        // pas) ; la règle des 10× reste seule juge de sa péremption.
        let mut jm_garde = jac_garde.take();
        // PRÉDICTEUR D'ORDRE 1 : accélération et multiplicateurs du pas d'AVANT
        // (a_{n-1}, λ_{n-1}), pour partir de 2a_n − a_{n-1} au lieu de a_n. Sur
        // la tête S2 (rotor établi, h·ω ≈ 0,03) le résidu de départ tombe d'autant
        // et Newton en fait moins ; réinitialisé sur impact et si h change.
        let mut pred_prec: Option<(DVector<f64>, DVector<f64>, f64)> = None;
        // RECUL EN PAS FIXE : un Newton qui lâche sur un événement (butée du KUKA
        // à h 2e-3 — mesuré, l'ordre 0 lâchait à 1,8, 2,1 et 2,2e-3 et ne passait
        // qu'à 2,0 par chance) se rejoue au demi-pas jusqu'au point de grille
        // manqué, puis le pas nominal reprend. Ce qui échouait passe ; le reste
        // est inchangé au bit. Plancher h/64, sinon l'erreur d'origine.
        let debug = std::env::var("VINKULUM_DEBUG").is_ok();
        let pred0 = std::env::var("VINKULUM_PRED0").is_ok();
        let diag_j = std::env::var("VINKULUM_DIAG_J").is_ok();
        let mut t_grille: Option<f64> = None;
        let mut t_contact = f64::NEG_INFINITY;
        let mut rejoues_au_pas = 0usize;
        let mut t_tentative = f64::NEG_INFINITY;
        while self.t < t_end {
            let tp0 = std::time::Instant::now();
            let t0 = self.t;
            if self.t == t_tentative {
                rejoues_au_pas += 1;
            } else {
                t_tentative = self.t;
                rejoues_au_pas = 0;
            }
            if rejoues_au_pas > 128 {
                return Err("nombre maximal de rejeux du pas atteint".into());
            }
            // PAS PILOTÉ PAR LE CONTACT : la distance à la prochaine collision
            // est CONNUE, donc l'estimateur n'a rien d'empirique. On resserre
            // avant le choc, on relâche après — le pas suit l'échelle active.
            if let Some((h_l, h_c, marge)) = self.pas_contact {
                let proche = self.contacts.iter().any(|ct| {
                    let (_, _, _, _, (d, _)) = ct.valeur(&self.corps);
                    d > -marge
                });
                if t_contact != t0 {
                    h_cur = if proche { h_c } else { h_l };
                    t_contact = t0;
                }
            }
            // Le résidu de fin ne fait PAS un pas : `t_end / h` tombe rarement
            // juste, et un dernier pas de quelques 1e-12 s rend le système
            // (termes en 1/h²) si mal conditionné que Newton diverge — mesuré
            // le 3 sept. sur un pendule à h = 1,25e-5 : « ne converge pas à
            // t=2.00000, résidu 4,5e7 », sans rien à voir avec la physique.
            // On l'ABSORBE : le dernier pas vaut entre h et 1,5 h. En
            // adaptatif le pas est rejoué, donc on se contente de borner.
            let t1 = numerique::fin_pas(t0, t_end, h_cur, self.adapt.is_none())?;
            // activation des contacts non lisses, figée sur le pas
            let u_now = self.u();
            self.k_pas = self
                .contacts
                .iter()
                .map(|ct| {
                    if !ct.nonlisse {
                        return (false, 0.0);
                    }
                    let (phi, g) = ct.phi_g(&self.corps);
                    let gap0 = -phi[0];
                    let mut gdot0 = 0.0;
                    for (kb, cb) in [Some(ct.corps), ct.b].iter().enumerate() {
                        if let Some(ib) = cb {
                            for c in 0..6 {
                                gdot0 -= g[(0, 6 * kb + c)] * u_now[6 * ib + c];
                            }
                        }
                    }
                    // ponytail: rappel de pénétration κ = 0,2 (l'ERP d'ODE) — converge
                    // en κ·β/γ par pas sans dépasser ; la loi d'impact reste pure
                    // ACTIF si le gap prédit ferme OU si le contact est déjà fermé à
                    // une tolérance près : à gap = 0⁺ exactement, l'activation par
                    // le seul gap prédit BATTAIT (inactif → chute de gh²/2 → impulsion
                    // ×3 au pas suivant — mesuré sur la pente, +7,8 % sur a). Rester
                    // actif ne coûte rien : λ ≥ 0 laisse partir librement.
                    (
                        gap0 + (t1 - t0) * gdot0 <= 0.0 || gap0 <= 1e-6 * ech_l.max(1.0),
                        ct.restitution * gdot0 + 0.2 * gap0.min(0.0) / (t1 - t0),
                    )
                })
                .collect();
            let ggl = if self.ggl {
                let d = 1e-7_f64.max(1e-9 * t1.abs());
                let phi_t = (self.phi_seul(t1 + d)? - self.phi_seul(t1)?) / d;
                Some((self.g_creuse(t0)?, phi_t))
            } else {
                None
            };
            let pas = Pas {
                am,
                af,
                gam,
                bet,
                sigma: self.sigma_lie,
                hh: t1 - t0,
                t1,
                u0: self.u(),
                poses0: self.poses(),
                udot: udot.clone(),
                a: a.clone(),
                ggl,
            };
            let mut x = DVector::zeros(dim);
            x.rows_mut(0, n).copy_from(&udot);
            x.rows_mut(n, m).copy_from(&lam);
            x.rows_mut(n + m, nz).copy_from(&ze);
            let mut extrapole = false;
            if let Some((up, lp, hp)) = pred_prec.as_ref() {
                // mesuré : S2 (143) Newton 4,94 → 3,94 par pas ; Kapitza (18,
                // consigne en table à cassures) +6 % — sous le seuil on garde
                // l'ordre 0, comme pour le jacobien hérité
                // `VINKULUM_PRED0=1` rend le prédicteur d'ordre 0 partout (A/B, robustesse)
                if dim >= JAC_HERITE_MIN
                    && *hp == pas.hh
                    && up.len() == n
                    && lp.len() == m
                    && !pred0
                {
                    x.rows_mut(0, n).copy_from(&(2.0 * &udot - up));
                    x.rows_mut(n, m).copy_from(&(2.0 * &lam - lp));
                    extrapole = true;
                }
            }
            let (ud0, lam0, a0) = (udot.clone(), lam.clone(), a.clone());
            // PRÉDICTEUR BORNÉ : l'itéré de départ extrapole l'accélération du pas
            // précédent ; sur un nœud d'inertie minuscule sous un moment (poutre),
            // ça prédit 12 rad de rotation et Newton part hors de son bassin
            // (mesuré : divergence au premier pas, OK à h/4). On ramène la ROTATION
            // PRÉDITE sous 0,3 rad en réduisant u̇ de départ — l'itéré, pas la solution.
            {
                let dth = pas.dtheta(&x, n);
                let pire = dth.iter().map(|v| v.norm()).fold(0.0, f64::max);
                if pire > seuil(0.3) {
                    let s = seuil(0.3) / pire;
                    for i in 0..n {
                        x[i] *= s;
                    }
                }
            }
            // deux échelles, et les confondre est un piège mesuré : le critère
            // de NEWTON se cale sur tout ce qui pèse dans le résidu, inertie
            // comprise (sans quoi un rotor, dont c'est toute la dynamique, a une
            // échelle nulle) ; le contrôle de PAS, lui, se cale sur les forces
            // RÉELLES du problème, comme le dit ABAQUS (« P is a typical
            // magnitude of real forces »). Y mettre l'inertie en fait un
            // cliquet : plus le pas est long, plus l'accélération numérique est
            // grande, plus le seuil monte — mesuré, le seuil valait 1223 pour
            // des forces de 70, et le pas montait à 0,6 point par période.
            let f_reel = self.forces().norm();
            let fnorm = f_reel + self.echelle_inertie(&udot);
            if !f_reel.is_finite() || !fnorm.is_finite() {
                return Err("norme des forces ou de l'inertie non finie".into());
            }
            self.s_fb = ech_l / (1.0 + fnorm);
            let tr0 = std::time::Instant::now();
            let mut r = pas.residu(self, &x);
            self.t_res += tr0.elapsed().as_secs_f64();
            let mut converge = false;
            // TOLÉRANCE 1e-12 PAR DÉFAUT, mesurée : à tol 1e-9 l'ordre de convergence
            // du multiplicateur tombe à 1,77 à h = T/800 (l'erreur de Newton passe
            // devant l'erreur du schéma) ; à 1e-13 il vaut 2,00. Le jacobien par
            // différences finies l'atteint (pas 1e-6, résidu quadratique en dessous).
            // CRITÈRE PAR BLOCS, mesuré le 3 sept. : les lignes de contrainte sont
            // mises à l'échelle 1/(β·h²) — à h = 1e-4 c'est 4e8, et un Φ au plancher
            // d'arrondi (1e-17 m) pèse 4e-9 dans la norme globale : Newton « ne
            // convergeait pas » à 2,5e-9 sur un résidu déjà nul. On juge la
            // dynamique en relatif aux forces, et Φ en ABSOLU, dans ses unités.
            let phi_scale = pas.bet * pas.hh * pas.hh;
            // PLANCHER D'ARRONDI des forces internes : une poutre de raideur EA
            // calcule γ comme une différence d'O(1), donc f = EA·γ porte EA·ε de
            // bruit (7e6 N/m → 1e-9 N). Le critère en résidu doit le savoir, sinon
            // Newton « ne converge pas » sur un résidu déjà au plancher (mesuré).
            let plancher = 1e-12
                * self
                    .poutres
                    .iter()
                    .map(|b| {
                        b.cn.iter().cloned().fold(0.0, f64::max)
                            + b.cm.iter().cloned().fold(0.0, f64::max) / (b.l * b.l)
                    })
                    .sum::<f64>();
            let tol_dyn = tol * (1.0 + fnorm) + plancher;
            // NEWTON SIMPLIFIÉ : le jacobien se calcule à la première itération
            // et se GARDE tant que le résidu chute d'au moins 10× par itération
            // (au point prédit il est déjà à O(h) du jacobien convergé) ; sinon
            // il se rafraîchit. Deux itérations, un seul jacobien — mesuré.
            // le jacobien gardé du pas précédent reste licite tant que le pas et la
            // dimension du système n'ont pas bougé (contacts activés, masque rejoué)
            if dim < JAC_HERITE_MIN
                || jm_garde.as_ref().is_some_and(|g| {
                    // L'interpolation des dates des micro-pas produit des
                    // différences d'arrondi dans h. Le pas intégré reste exact ;
                    // seul le cache tolère 1e-10 relatif, sous le juge de Newton.
                    (g.h != pas.hh
                        && (self.gel.is_empty()
                            || (g.h - pas.hh).abs() > 1e-10 * g.h.abs().max(pas.hh.abs())))
                        || g.dim != dim
                        || g.reduction.as_ref().map(|r| &r.indices)
                            != Indices::du_modele(self).as_ref()
                })
            {
                jm_garde = None;
            }
            let mut r_prec = f64::INFINITY;
            let mut erreur_jac: Option<String> = None;
            // jacobien calculé DANS ce pas (frais) ou hérité du pas précédent
            let mut frais = false;
            // un ACCROC (rebroussement de Newton) dit que le pas n'est pas lisse :
            // le suivant ne s'extrapole pas (KUKA sur butée, mesuré : l'extrapolation
            // à travers l'événement change de branche). Un pas hérité refusé n'en
            // est pas un — c'est un jacobien périmé, pas un événement (S2 : compter
            // les refus coûtait la moitié du gain, 0,172 → 0,185 ms/pas).
            let mut accroc = false;
            for it in 0..newton_max {
                if !r.iter().all(|v| v.is_finite()) {
                    break;
                }
                let r_dyn = r.rows(0, n).norm();
                let r_phi = r.rows(n, m).norm() * phi_scale;
                // GGL : G·u en ABSOLU dans ses unités (vitesse), relatif à ‖u₀‖
                let r_gu = r.rows(n + m, nz).norm() * pas.gam * pas.hh / (1.0 + pas.u0.norm());
                if debug {
                    let (imax, vmax) = r.iter().enumerate().fold((0, 0.0), |a, (i, v)| {
                        if v.abs() > a.1 {
                            (i, v.abs())
                        } else {
                            a
                        }
                    });
                    eprintln!(
                        "  t={:.6} it={} r_dyn={:.3e} r_phi={:.3e} (tol_dyn {:.1e}) frais={} max|r| {:.2e} @{}{}",
                        t1,
                        it,
                        r_dyn,
                        r_phi,
                        tol_dyn,
                        frais,
                        vmax,
                        imax,
                        if imax < n { format!(" corps {} ddl {}", imax / 6, imax % 6) } else { format!(" contrainte {}", imax - n) }
                    );
                }
                let r_tot = r.norm();
                let sous_tol = r_dyn < tol_dyn && r_phi < 1e-11 && r_gu < 1e-11;
                // à moins de 100× du seuil, un jacobien hérité qui descend encore
                // finit en une ou deux descentes : on ne lui demande plus les 10×
                // (mesuré sur la tête S2 : refusé à 3,3× du seuil, il coûtait un
                // jacobien tous les deux pas)
                let pres = r_dyn < seuil(100.0) * tol_dyn && r_phi < 1e-9 && r_gu < 1e-9;
                // Sous un jacobien HÉRITÉ la convergence est linéaire : le dernier
                // itéré s'arrête AU seuil, là où le quadratique le dépassait de
                // plusieurs décades — mesuré : invariance de repère 3e-10 m au lieu
                // de 5e-13, DF de sensibilité modale à 3,5 % au lieu de 0,6. Sous le
                // seuil on continue donc de DESCENDRE avec le même jacobien (une
                // descente coûte un résidu, pas un jacobien) jusqu'à ce qu'il ne
                // rende plus 2× — le plancher d'arrondi, celui que le quadratique
                // atteint en dépassant — ou six décades sous le seuil.
                if sous_tol && (frais || r_tot > 0.5 * r_prec || r_dyn < 1e-6 * tol * (1.0 + fnorm))
                {
                    converge = true;
                    break;
                }
                // STAGNATION : le résidu ne descend plus (rapport > 0,9 deux fois
                // de suite) alors qu'un jacobien frais a déjà été essayé — on est
                // au plancher d'arrondi du modèle, pas en divergence. Mesuré sur
                // une console à 4 éléments : r stagne à 2e-9 pendant 25 itérations.
                // RETIRÉ PUIS REMIS le 4 sept. : sans lui, `verification` casse
                // (« Newton ne converge pas à t=2.36332, résidu 1.14e-7 ») —
                // le critère sur l'incrément ne le couvre PAS.
                if it > 2 && r_tot > 0.9 * r_prec && r_tot < seuil(1e-6) * (1.0 + fnorm) {
                    converge = true;
                    break;
                }
                // sous le seuil, un jacobien frais coûterait dix descentes pour
                // gagner ce qu'une descente héritée rend presque gratuitement
                // près du seuil, r_tot est porté par les lignes de Φ à l'échelle
                // 1/(βh²) — de l'arrondi : il ne juge plus la péremption
                let rafraichir =
                    jm_garde.is_none() || (!sous_tol && !pres && r_tot > seuil(0.1) * r_prec);
                r_prec = r_tot;
                if rafraichir {
                    let tj0 = std::time::Instant::now();
                    // le modèle est posé à l'état de x par le dernier `residu`
                    let (trip, t_exp) = match self.jacobien_ad(&pas, &x, &pas.dtheta(&x, n)) {
                        Ok(v) => v,
                        Err(e) => {
                            // l'itéré est hors du domaine d'un élément (angle
                            // d'engrenage ou de vis sauté) : le résidu y vaut
                            // déjà 1e30 et le jacobien n'existe pas. On sort
                            // du Newton comme d'une non-convergence.
                            erreur_jac = Some(e);
                            break;
                        }
                    };
                    self.t_exp += t_exp;
                    if diag_j && t1 > 3.0 {
                        if let Some(g) = jm_garde.as_ref() {
                            static UNE: std::sync::Once = std::sync::Once::new();
                            UNE.call_once(|| self.diag_derive_j(g, &trip, dim, n, t1));
                        }
                    }
                    // les colonnes ζ de GGL sont dans `trip`, EXACTES (`jacobien_ad`) :
                    // les différences finies par colonne — une copie du modèle
                    // chacune — sont parties le 5 sept. au soir avec le point mort
                    // de `srscm` qu'elles ne traversaient pas.
                    self.jacobien_total += 1;
                    frais = true;
                    self.t_jac += tj0.elapsed().as_secs_f64();
                    let tf0 = std::time::Instant::now();
                    let _fac = ChronoFac(tf0);
                    let (reduction, trip) = if let Some(indices) = Indices::du_modele(self) {
                        let (r, trip) = Reduction::nouvelle(indices, trip);
                        (Some(r), trip)
                    } else {
                        (None, trip)
                    };
                    let dim_fact = reduction.as_ref().map_or(dim, |r| r.indices.dim());
                    let f = if dim_fact == 0 {
                        None
                    } else {
                        Some(if dim_fact <= DENSE_MAX {
                            Facto::dense_seule(&trip, dim_fact)
                        } else {
                            {
                                // le motif a-t-il bougé depuis la symbolique gardée ?
                                use std::hash::{Hash, Hasher};
                                let th0 = std::time::Instant::now();
                                let mut hs = std::hash::DefaultHasher::new();
                                dim_fact.hash(&mut hs);
                                reduction.as_ref().map(|r| &r.indices).hash(&mut hs);
                                for t in &trip {
                                    (t.row, t.col).hash(&mut hs);
                                }
                                let h = hs.finish();
                                T_HASH.fetch_add(
                                    th0.elapsed().as_nanos() as u64,
                                    std::sync::atomic::Ordering::Relaxed,
                                );
                                NNZ.store(trip.len() as u64, std::sync::atomic::Ordering::Relaxed);
                                if h != self.symb_motif {
                                    self.symb = None;
                                    self.symb_motif = h;
                                }
                                Facto::creux(&trip, dim_fact, &mut self.symb)?
                            }
                        })
                    };
                    self.t_sol += tf0.elapsed().as_secs_f64();
                    jm_garde = Some(JacGarde::neuf(trip, f, pas.hh, dim, &self.corps, reduction));
                }
                let ts0 = std::time::Instant::now();
                let (dx, svd) =
                    jm_garde
                        .as_mut()
                        .unwrap()
                        .resout(&self.corps, &self.elems, dim, &r);
                self.t_sol += ts0.elapsed().as_secs_f64();
                if svd {
                    self.svd_total += 1;
                }
                let dx_rel = dx.norm() / (1.0 + x.norm());
                self.newton_total += 1;
                // NEWTON AMORTI (rebroussement) : sur une poutre à nœuds d'inertie
                // minuscule le prédicteur est loin (résidu 1e6) et le pas plein sort
                // du bassin — mesuré, divergence en 25 itérations. Si le résidu
                // MONTE, on divise le pas par 2, jusqu'à 6 fois ; le jacobien est
                // rafraîchi à l'itération suivante (le résidu n'a pas chuté de 10×).
                let r0n = r.norm();
                let tr0 = std::time::Instant::now();
                let mut alpha = 1.0;
                let mut k_bt = 0;
                let mut refuse = false;
                loop {
                    let x_essai = &x - alpha * &dx;
                    let r_essai = pas.residu(self, &x_essai);
                    let fini = r_essai.iter().all(|v| v.is_finite());
                    if !frais {
                        // JACOBIEN HÉRITÉ : son pas n'est accepté que s'il fait les
                        // 10× d'un Newton sain ; sinon on RESTE au point courant et on
                        // le rafraîchit — ni rebroussement (halver un pas faux ne le
                        // rend pas juste) ni pas forcé. Mesuré le 7 sept. sur le
                        // quatre-barres flexible à h = 4e-3 : un premier pas hérité
                        // accepté à 2× de chute sortait du bassin et Newton lâchait à
                        // t = 0,016 ; refusé, la trajectoire de Newton est celle du
                        // jacobien frais partout où l'hérité n'est pas aussi bon.
                        // Sous le seuil on descend vers le plancher : toute descente
                        // est prise, et un pas qui ne descend plus EST le plancher.
                        let rn = r_essai.norm();
                        if fini
                            && if sous_tol || pres {
                                rn < r0n
                            } else {
                                rn <= seuil(0.1) * r0n
                            }
                        {
                            x = x_essai;
                            r = r_essai;
                        } else if sous_tol || pres {
                            // PLANCHER, pas péremption : mesuré sur la tête S2, le taux
                            // de descente du J tourné ne bouge pas avec son âge (1e-7 à
                            // l'âge 1, 2e-8 à l'âge 11), et le refus tombait toujours à
                            // moins de 100× du seuil, sur des lignes de Φ à l'arrondi.
                            converge = true;
                            refuse = true;
                        } else {
                            jm_garde = None;
                            refuse = true;
                            // un pas hérité qui fait MONTER le résidu est le signe
                            // d'un événement (butée du KUKA, mesuré) ; qui descend
                            // de moins de 10× n'est qu'un jacobien périmé
                            if !(fini && rn <= r0n) {
                                accroc = true;
                            }
                        }
                        break;
                    }
                    if fini && (k_bt >= 6 || r_essai.norm() <= r0n) {
                        x = x_essai;
                        r = r_essai;
                        break;
                    }
                    if k_bt >= 6 {
                        pas.residu(self, &x);
                        accroc = true;
                        break;
                    }
                    alpha *= 0.5;
                    k_bt += 1;
                    accroc = true;
                }
                self.t_res += tr0.elapsed().as_secs_f64();
                if refuse {
                    // le modèle est reposé à x par ce dernier `residu`
                    pas.residu(self, &x);
                    if converge {
                        break;
                    }
                    continue;
                }
                // CRITÈRE SUR L'INCRÉMENT, mesuré le 4 sept. sur une poutre : les
                // forces internes EA·γ portent un plancher d'arrondi EA·ε ≈ 1e-9 N
                // que le critère en résidu (1e-12·(1 + |f|)) ne peut pas franchir —
                // Newton « ne convergeait pas » à 3e-9. Quand le pas de Newton
                // devient négligeable devant x, la solution est au plancher : convergé.
                if dx_rel < seuil(1e-11)
                    && r.iter().all(|v| v.is_finite())
                    && r.rows(0, n).norm() <= tol_dyn.max(1e-6 * (1.0 + fnorm))
                    && r.rows(n, m).norm() * phi_scale < 1e-9
                    && r.rows(n + m, nz).norm() * pas.gam * pas.hh / (1.0 + pas.u0.norm()) < 1e-9
                {
                    converge = true;
                    break;
                }
            }
            if !converge && r.iter().all(|v| v.is_finite()) {
                converge = r.rows(0, n).norm() < tol_dyn
                    && r.rows(n, m).norm() * phi_scale < 1e-11
                    && r.rows(n + m, nz).norm() * pas.gam * pas.hh / (1.0 + pas.u0.norm()) < 1e-11;
            }
            if !converge && extrapole {
                // le départ extrapolé a sorti Newton de son bassin (quatre-barres
                // flexible à h 4e-3, mesuré) : le pas se rejoue à l'ordre 0, et
                // c'est seulement là que l'échec est un échec
                pred_prec = None;
                self.set_poses(&pas.poses0);
                self.set_u(&pas.u0);
                n_rejeu += 1;
                continue;
            }
            if !converge {
                if let Some(e) = erreur_jac {
                    // PAS de `masque_rejoue` ici : il re-détecte la redondance à
                    // la pose courante, celle-là même où Φ n'est pas défini.
                    // En adaptatif, c'est un pas trop long : on le rejoue plus
                    // court ; sinon on rend l'erreur de l'élément, qui dit quoi
                    // faire (réduire le pas).
                    if self.adapt.is_some() && h_cur > h_min * 1.5 {
                        self.set_poses(&pas.poses0);
                        self.set_u(&pas.u0);
                        h_cur *= 0.5;
                        n_rejeu += 1;
                        continue;
                    }
                    return Err(format!("Newton à t={t1:.5} : jacobien impossible — {e}"));
                }
                // un masque de redondance périmé (G quasi singulier au point de
                // départ) fait caler Newton : on le re-détecte à l'itéré courant
                if self.masque_rejoue(t1, t0, &pas.poses0, &pas.u0)? {
                    jm_garde = None;
                    n_rejeu += 1;
                    continue;
                }
                // en adaptatif, un Newton qui ne converge pas n'est pas un
                // échec : c'est un pas trop long. On le rejoue plus court.
                if self.adapt.is_some() && h_cur > h_min * 1.5 {
                    self.set_poses(&pas.poses0);
                    self.set_u(&pas.u0);
                    h_cur *= 0.5;
                    n_rejeu += 1;
                    continue;
                }
                if self.adapt.is_none() && h_cur > h / 64.0 {
                    self.set_poses(&pas.poses0);
                    self.set_u(&pas.u0);
                    t_grille.get_or_insert(t1);
                    h_cur *= 0.5;
                    n_rejeu += 1;
                    continue;
                }
                if let Err(cause) = pas.a1_dq(&x) {
                    return Err(format!("Newton à t={t1:.5} : {cause}"));
                }
                return Err(format!(
                    "Newton ne converge pas à t={:.5} (résidu {:.2e})",
                    t1,
                    r.norm()
                ));
            }
            let (ud, la, a1) = pas.etat(self, &x)?;
            if let Some(tol_rel) = self.adapt {
                // état à t+h/2 par intégration EXACTE d'une accélération
                // linéaire entre les deux accélérations VRAIES du pas :
                //   ü(h/2)  = (a₀ + a₁)/2
                //   u̇(h/2)  = u̇₀ + (h/8)(3a₀ + a₁)
                //   Δq(h/2) = (h/2)u̇₀ + h²(5a₀ + a₁)/48
                let hh = pas.hh;
                let poses1 = self.poses();
                let u1 = self.u();
                let ah = 0.5 * (&pas.udot + &ud);
                let uh = &pas.u0 + (hh / 8.0) * (3.0 * &pas.udot + &ud);
                let dq = 0.5 * &pas.u0 + (hh / 48.0) * (5.0 * &pas.udot + &ud);
                self.avance(&pas.poses0, &dq, hh);
                self.set_u(&uh);
                let lh = 0.5 * (&lam + &la);
                let rh = self.residu_pose(&ah, &lh, t0 + 0.5 * hh, 1.0, 1.0, None);
                let r_half = (0..n).fold(0.0f64, |acc, i| acc.max(rh[i].abs()));
                // et la CONTRAINTE à mi-pas : Φ(q_{h/2}) en mètres, relatif à la
                // taille du modèle — l'index 3 ne la tient qu'aux extrémités du pas,
                // sa violation au milieu est l'erreur locale que λ paie
                let lk = self.lignes_k();
                let phi_half = (n..n + m)
                    .filter(|&i| self.actif[i - n] && !lk[i - n])
                    .fold(0.0f64, |acc, i| acc.max(rh[i].abs()))
                    / ech_l.max(1.0);
                self.set_poses(&poses1);
                self.set_u(&u1);
                let seuil = tol_rel * (1.0 + f_reel);
                self.r_half_max = self.r_half_max.max(r_half);
                self.seuil_max = self.seuil_max.max(seuil);
                self.phi_half_max = self.phi_half_max.max(phi_half);
                if (r_half > seuil || phi_half > tol_rel) && h_cur > h_min * 1.5 {
                    self.set_poses(&pas.poses0);
                    self.set_u(&pas.u0);
                    h_cur = (h_cur * 0.5).max(h_min);
                    n_rejeu += 1;
                    continue;
                }
                // bien en dessous : on allonge, mais doucement — un pas qui
                // change injecte du bruit haute fréquence, que la dissipation
                // du schéma doit avoir le temps d'absorber (d'où ρ∞ < 1)
                if r_half < 0.2 * seuil && phi_half < 0.2 * tol_rel {
                    h_cur = (h_cur * 1.25).min(h_max);
                }
            }
            // CCD : la trajectoire du pas traverse-t-elle ? Si oui, on la
            // REJOUE plus court. C'est la garantie topologique d'IPC portée au
            // niveau du pas — la barrière rend l'approche coûteuse, le CCD
            // rend la traversée impossible.
            if self.ccd && self.traverse(&pas.poses0) {
                if h_cur > h * 1.0 / 256.0 {
                    self.set_poses(&pas.poses0);
                    self.set_u(&pas.u0);
                    h_cur *= 0.5;
                    n_rejeu += 1;
                    continue;
                }
                return Err(format!("CCD : traversée à t={t1:.6} même au pas minimal"));
            }
            pred_prec = if accroc {
                None
            } else {
                Some((ud0.clone(), lam0.clone(), pas.hh))
            };
            if let Some(tg) = t_grille {
                // retour sur la grille : on reste au pas réduit jusqu'au point
                // manqué, puis le pas nominal reprend
                if t1 >= tg - 1e-12 {
                    t_grille = None;
                    h_cur = h;
                } else {
                    h_cur = h_cur.min(tg - t1);
                }
            }
            udot = ud;
            lam = la;
            a = a1;
            ze = x.rows(n + m, nz).into_owned();
            if self.k_pas.iter().any(|(k, _)| *k) {
                pred_prec = None;
                // IMPACT : l'accélération du pas est impulsive et l'α-généralisé
                // en garde la mémoire ((1−γ)h·a au pas suivant) — mesuré : une
                // bille à e = 0 repartait à 0,8 fois sa vitesse d'arrivée. On
                // REDÉMARRE le schéma sur l'accélération LISSE, consistante avec
                // les contacts actifs (Chen–Acary–Virlez–Brüls 2013 séparent de
                // même la partie impulsive). Le pas de contact est d'ordre 1.
                let (ud_s, _) = self.acc_init(t1)?;
                udot = ud_s.clone();
                a = ud_s;
            }
            self.verifie_bassin(t1)?;
            self.lam = lam.clone();
            self.t = t1;
            // Φ̇ reste violé à O(h²) : c'est l'index 3 direct. La projection de
            // vitesse après le pas a été essayée et RETIRÉE (JOURNAL, 4 sept.) :
            // elle annule Φ̇ mais fait tomber l'ordre de 2,00 à ~0,8. Le remède
            // correct est GGL — des multiplicateurs dans la FORMULATION.
            if self.proj_moment || self.proj_energie {
                self.projette_invariants(t1)?;
            }
            let tf1 = std::time::Instant::now();
            // LES CONTRAINTES NE SE VIOLENT PAS EN SILENCE. Le noyau NEUTRALISE
            // les lignes de G redondantes — c'est ce qui fait passer un
            // quatre-barres surcontraint ou un Bennett, et c'est un acquis. Mais
            // deux contraintes qui se CONTREDISENT (un encastrement et une
            // rotation imposée sur le même axe) sont, elles aussi, redondantes :
            // le solveur en satisfait une, laisse Φ dériver, et rendait
            // jusqu'ici un résultat entièrement faux sans un mot.
            //
            // Trouvé le 3 sept. en montant l'arbre tournant du jeu de tests
            // MBDyn : l'arbre ne tournait pas du tout, Φ valait 0,96, et
            // `simule` rendait normalement. La différence entre une redondance
            // LICITE et une redondance FAUTIVE se lit sur Φ, pas sur le rang.
            if self.tol_phi > 0.0 {
                let e = self.phi_holonome_max(t1)?;
                // relatif à la taille du modèle : Φ mêle longueurs et angles,
                // et une tolérance en mètres n'a pas de sens sur un mécanisme
                // coté en micromètres (ni sur un pont de 200 m).
                let seuil = self.tol_phi * ech_l.max(1.0);
                if e > seuil && self.masque_rejoue(t1, t0, &pas.poses0, &pas.u0)? {
                    udot = ud0;
                    self.lam = lam0.clone();
                    lam = lam0;
                    a = a0;
                    self.t = t0;
                    n_rejeu += 1;
                    continue;
                }
                if e > seuil {
                    let n_red = self.actif.iter().filter(|a| !**a).count();
                    return Err(format!(
                        "t = {t1:.6} : ‖Φ‖ = {e:.3e} au-dessus de la tolérance \
                         {:.1e} — les contraintes ne sont plus satisfaites{}",
                        seuil,
                        if n_red > 0 {
                            format!(
                                " ({n_red} ligne(s) neutralisée(s) comme redondante(s) : \
                                 deux liaisons se contredisent-elles ?)"
                            )
                        } else {
                            String::new()
                        }
                    ));
                }
            }
            self.verifie_etat_fini()?;
            self.h_dernier = pas.hh;
            self.actualise()?;
            self.verifie_etat_fini()?;
            if let Some(hs) = self.hist_schema.as_mut() {
                hs.push((
                    udot.as_slice().to_vec(),
                    a.as_slice().to_vec(),
                    lam.as_slice().to_vec(),
                ));
            }
            self.schema_fin = Some((udot.clone(), a.clone(), lam.clone()));
            self.t_fin += tf1.elapsed().as_secs_f64();
            self.t_pas += tp0.elapsed().as_secs_f64();
            npas += 1;
            sortie(t1, self);
        }
        self.n_rejeu = n_rejeu;
        *jac_garde = jm_garde;
        self.schema_fin = Some((udot, a, lam));
        Ok(npas)
    }
}

/// SVD QUI NE BOUCLE PAS. `Matrix::svd` de nalgebra itère SANS plafond, et sur
/// un NaN il ne converge jamais — trouvé par sondage le 7 sept. (masse NaN,
/// gravité NaN, normale nulle dont `normalize()` rend NaN) : `acc_init` bloqué
/// dans la SVD, aucune trace, un processus à tuer. Une entrée non finie est
/// refusée AVANT ; le calcul borné de faer signale la non-convergence.
pub fn svd_sure(
    a: DMatrix<f64>,
) -> Result<nalgebra::linalg::SVD<f64, nalgebra::Dyn, nalgebra::Dyn>, String> {
    if a.iter().any(|v| !v.is_finite()) {
        return Err(
            "matrice non finie (NaN ou ∞) : masse, inertie, gravité, normale ou force mal posée"
                .into(),
        );
    }
    svd::decompose(&a)
}

/// Le jacobien de Newton en TRIPLETS (i, j, v) — assemblé par éléments, jamais
/// densifié au-delà de `DENSE_MAX` inconnues.
pub type Triplets = Vec<faer::sparse::Triplet<usize, usize, f64>>;
/// en dessous : LU dense (plus rapide qu'un creux sur une petite matrice)
/// En dessous de cette dimension Newton reste celui d'origine — jacobien frais à
/// chaque pas, prédicteur d'ordre 0 : le jacobien y est bon marché, et le Newton
/// quadratique y achète, en dépassant le seuil de plusieurs décades, l'invariance
/// au plancher (banc `invariances`, 30 inconnues : 5e-13 m ; sous jacobien hérité
/// 1e-12, deux fois moins). Au-dessus : jacobien gardé d'un pas à l'autre et
/// prédicteur 2a_n − a_{n−1}. Même seuil que le parallélisme des éléments : sur
/// un petit modèle, rien ne vaut d'être économisé.
pub const JAC_HERITE_MIN: usize = 48;
pub const DENSE_MAX: usize = 160;
/// au-delà : le LU creux est autorisé à paralléliser (rayon)
pub const CREUX_PAR: usize = 2000;
/// au-delà, la détection A PRIORI des contraintes redondantes est sautée : son
/// QR pivoté est dense et coûte O(m²n) — mesuré 5,9 s à m = 1 800, contre
/// 8,2 ms pour un pas. Le repli SVD prend alors le relais.
pub const REDONDANCE_MAX: usize = 900;

/// Factorisation LU gardée avec le jacobien : dense (pivotage partiel) ou
/// creuse (structure symbolique réutilisée d'un jacobien à l'autre — le motif
/// ne change pas, seuls les coefficients).
/// chronos FINS de la résolution (ns cumulés, lus par `chronos()` sous
/// `VINKULUM_CHRONO`) : factorisation, descente, contrôle Jx − b
pub static T_FAC: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static T_SV: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static T_MV: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static T_HASH: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static NNZ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Le jacobien gardé d'un pas à l'autre : J₀ factorisé, ses triplets (contrôle
/// Jx − b une fois), le pas et la dimension qui l'invalident, et les rotations
/// des corps à l'assemblage — pour le TOURNER au lieu de le refaire.
pub struct JacGarde {
    trip: Triplets,
    f: Option<Facto>,
    reduction: Option<Reduction>,
    h: f64,
    dim: usize,
    verifie: bool,
    rot0: Vec<M3>,
    tourne: bool,
    chrono: bool,
}

impl JacGarde {
    fn neuf(
        trip: Triplets,
        f: Option<Facto>,
        h: f64,
        dim: usize,
        corps: &[Corps],
        reduction: Option<Reduction>,
    ) -> Self {
        JacGarde {
            trip,
            f,
            reduction,
            h,
            dim,
            verifie: false,
            rot0: corps.iter().map(|c| c.rot).collect(),
            tourne: std::env::var("VINKULUM_TOURNE0").is_err(),
            chrono: std::env::var("VINKULUM_CHRONO").is_ok(),
        }
    }
    /// x = T_c·J₀⁻¹·T_rᵀ·b — les 6 lignes et 6 colonnes de chaque corps tournées
    /// de Qᵢ = Rᵢ·Rᵢ⁰ᵀ ; lignes de contrainte (matérielles) et de GGL inchangées.
    /// `VINKULUM_TOURNE0=1` désactive la rotation (A/B).
    fn resout(
        &mut self,
        corps: &[Corps],
        elems: &[Elem],
        dim: usize,
        b: &DVector<f64>,
    ) -> (DVector<f64>, bool) {
        let tourne = self.rot0.len() == corps.len()
            && 6 * corps.len() <= dim
            && dim >= JAC_HERITE_MIN
            && self.tourne;
        let q: Vec<M3> = if tourne {
            corps
                .iter()
                .zip(&self.rot0)
                .map(|(c, r0)| c.rot * r0.transpose())
                .collect()
        } else {
            vec![]
        };
        let tourner = |v: &mut DVector<f64>, transpose: bool| {
            for (i, qi) in q.iter().enumerate() {
                for k in [6 * i, 6 * i + 3] {
                    let s = V3::new(v[k], v[k + 1], v[k + 2]);
                    let t = if transpose {
                        qi.transpose() * s
                    } else {
                        qi * s
                    };
                    v[k] = t.x;
                    v[k + 1] = t.y;
                    v[k + 2] = t.z;
                }
            }
        };
        if tourne && self.chrono {
            // MESURE : quels éléments ont vu leurs corps tourner l'un par rapport
            // à l'autre depuis l'assemblage (bloc approché), et de combien
            let mut pire = 0.0f64;
            let mut liste = vec![];
            for e in elems.iter() {
                let cs = e.corps();
                let qa = cs.first().copied().flatten().map(|i| q[i]);
                let qb = cs.get(1).copied().flatten().map(|i| q[i]);
                if let (Some(qa), Some(qb)) = (qa, qb) {
                    let rel = qa * qb.transpose();
                    let ang = ((rel.trace() - 1.0) / 2.0)
                        .clamp(-1.0, 1.0)
                        .acos()
                        .to_degrees();
                    if ang > pire {
                        pire = ang;
                    }
                    if ang > 0.05 {
                        liste.push(format!("{}:{:.2}°", e.nom(), ang));
                    }
                }
            }
            if pire > 3.0 {
                static UNE: std::sync::Once = std::sync::Once::new();
                UNE.call_once(|| {
                    eprintln!("éléments mixtes (dérive relative) : {}", liste.join(" · "))
                });
            }
        }
        let mut bt = b.clone();
        if tourne {
            tourner(&mut bt, true);
        }
        // Tourner en indices GLOBAUX avant restriction : les corps mobiles
        // peuvent être non contigus. Les lignes éliminées de Newton rendent
        // x_e = b_e, d'où le report A_le b_e dans le second membre local.
        let (mut y, svd) = if let Some(r) = &self.reduction {
            let rhs = r.rhs(&bt, &bt);
            let (y, svd) = if let Some(f) = &self.f {
                Modele::resout(f, &self.trip, rhs.len(), &rhs, !self.verifie)
            } else {
                (DVector::zeros(0), false)
            };
            (r.etend(&y, bt), svd)
        } else if let Some(f) = &self.f {
            Modele::resout(f, &self.trip, dim, &bt, !self.verifie)
        } else {
            (DVector::zeros(0), false)
        };
        if !svd {
            self.verifie = true;
        }
        if tourne {
            tourner(&mut y, false);
        }
        (y, svd)
    }
}

impl Modele {
    /// DIAGNOSTIC (VINKULUM_DIAG_J) : écart entre le jacobien tourné T·J₀·Tᵀ et le
    /// jacobien frais, agrégé par bloc (élément × corps), pour savoir quel bloc
    /// dérive vraiment quand un jacobien gardé se rafraîchit.
    fn diag_derive_j(&self, g: &JacGarde, frais: &Triplets, dim: usize, n: usize, t1: f64) {
        let dense = |tr: &Triplets| {
            let mut a = DMatrix::zeros(dim, dim);
            for t in tr {
                a[(t.row, t.col)] += t.val;
            }
            a
        };
        let trip0 = g.reduction.as_ref().map(|r| r.triplets_globaux(&g.trip));
        let j0 = dense(trip0.as_ref().unwrap_or(&g.trip));
        let jf = dense(frais);
        let q: Vec<M3> = self
            .corps
            .iter()
            .zip(&g.rot0)
            .map(|(c, r0)| c.rot * r0.transpose())
            .collect();
        let mut tr = DMatrix::identity(dim, dim);
        for (i, qi) in q.iter().enumerate() {
            for k in [6 * i, 6 * i + 3] {
                tr.view_mut((k, k), (3, 3)).copy_from(qi);
            }
        }
        // ESSAYÉ ET REFUSÉ (7 sept.) : tourner aussi les lignes de Φ des pivots au
        // bâti comme leur corps b — dérive totale 6,2 → 4,1 % et jac/pas 0,29 → 0,26,
        // mais Newton 4,0 → 4,4 et 0,141 → 0,153 ms/pas. Les lignes de Φ restent fixes.
        let jt = &tr * &j0 * tr.transpose();
        let d = &jf - &jt;
        // nom d'une ligne/colonne
        let mut lignes: Vec<String> = (0..n).map(|i| format!("corps{}", i / 6)).collect();
        let mut r = n;
        for e in &self.elems {
            for _ in 0..e.n() {
                lignes.push(format!("Φ:{}", e.nom()));
            }
            r += e.n();
        }
        while lignes.len() < dim {
            lignes.push(format!("ligne{}", lignes.len()));
        }
        let _ = r;
        let mut agg: std::collections::HashMap<(String, String), (f64, f64)> = Default::default();
        for i in 0..dim {
            for j in 0..dim {
                let k = (lignes[i].clone(), lignes[j].clone());
                let e = agg.entry(k).or_insert((0.0, 0.0));
                e.0 += d[(i, j)] * d[(i, j)];
                e.1 += jf[(i, j)] * jf[(i, j)];
            }
        }
        let mut v: Vec<_> = agg.into_iter().filter(|(_, (a, _))| *a > 0.0).collect();
        v.sort_by(|a, b| b.1 .0.partial_cmp(&a.1 .0).unwrap());
        eprintln!("dérive du jacobien tourné à t={t1:.4} (âge : depuis l'assemblage) — ‖Δ‖ / ‖J‖ par bloc :");
        for ((li, co), (a, b)) in v.iter().take(12) {
            eprintln!(
                "   {li:>22} × {co:<22} ‖Δ‖ {:.2e}  ‖J‖ {:.2e}  rel {:.1e}",
                a.sqrt(),
                b.sqrt(),
                a.sqrt() / b.sqrt().max(1e-300)
            );
        }
        // les huit plus grandes entrées : tournée contre fraîche
        let mut ent: Vec<(f64, usize, usize)> = vec![];
        for i in 0..dim {
            for j in 0..dim {
                if d[(i, j)].abs() > 1e-12 {
                    ent.push((d[(i, j)].abs(), i, j));
                }
            }
        }
        ent.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        for (dv, i, j) in ent.iter().take(10) {
            eprintln!(
                "   [{i:3},{j:3}] {:>22} × {:<22} tournée {:+.4e}  fraîche {:+.4e}  Δ {dv:.2e}",
                lignes[*i],
                lignes[*j],
                jt[(*i, *j)],
                jf[(*i, *j)]
            );
        }
        let (da, ja) = (d.norm(), jf.norm());
        eprintln!("   total ‖Δ‖ {da:.2e} / ‖J‖ {ja:.2e} = {:.1e}", da / ja);
    }
}

pub enum Facto {
    Dense(faer::linalg::solvers::PartialPivLu<f64>),
    Creux(Box<LuCreux>),
}

/// LU CREUX à parallélisme EXPLICITE.
///
/// Le solveur « haut niveau » de faer (`sparse::linalg::solvers::Lu`) lit le
/// réglage GLOBAL à chaque factorisation et à chaque descente, et le noyau le
/// basculait Seq ↔ rayon autour de chaque grande factorisation. Ce réglage est
/// un `static` du processus, partagé par tous les `Noyau` — et `simule` libère
/// le GIL précisément pour qu'un `ThreadPoolExecutor` en fasse tourner
/// plusieurs : course sur un état global (relu le 6 sept.). Ici `par` est un
/// champ, passé à la factorisation et à la descente par l'API bas niveau ; le
/// global reste Seq, posé une fois, et plus jamais touché.
pub struct LuCreux {
    symb: std::sync::Arc<faer::sparse::linalg::lu::SymbolicLu<usize>>,
    num: faer::sparse::linalg::lu::NumericLu<usize, f64>,
    par: faer::Par,
    mem_solve: std::sync::Mutex<faer::dyn_stack::MemBuffer>,
}

fn faer_seq() {
    // faer parallélise ses factorisations sur rayon par défaut : sur une
    // matrice de 143, le coût de synchronisation dépasse le calcul (mesuré :
    // 1,5 ms par pas au lieu de 0,05). Séquentiel — le parallélisme du
    // noyau est ailleurs (éléments), pas dans un LU. Posé UNE fois : le LU
    // creux choisit son parallélisme par appel (`LuCreux`), sans y toucher.
    static UNE_FOIS: std::sync::Once = std::sync::Once::new();
    UNE_FOIS.call_once(|| faer::set_global_parallelism(faer::Par::Seq));
}

impl Facto {
    #[cfg(test)]
    fn dense(trip: &Triplets, dim: usize) -> (DMatrix<f64>, Self) {
        let mut a = DMatrix::zeros(dim, dim);
        for t in trip {
            a[(t.row, t.col)] += t.val;
        }
        (a, Self::dense_seule(trip, dim))
    }
    /// la factorisation seule, assemblée directement dans la matrice faer
    fn dense_seule(trip: &Triplets, dim: usize) -> Self {
        #[cfg(test)]
        partition::tests::FACTOS.with(|f| f.borrow_mut().push(dim));
        faer_seq();
        let mut m = faer::Mat::<f64>::zeros(dim, dim);
        for t in trip {
            m[(t.row, t.col)] += t.val;
        }
        Facto::Dense(faer::linalg::solvers::PartialPivLu::new(m.as_ref()))
    }
    fn creux(
        trip: &Triplets,
        dim: usize,
        symb: &mut Option<std::sync::Arc<faer::sparse::linalg::lu::SymbolicLu<usize>>>,
    ) -> Result<Self, String> {
        #[cfg(test)]
        partition::tests::FACTOS.with(|f| f.borrow_mut().push(dim));
        use faer::dyn_stack::{MemBuffer, MemStack};
        use faer::sparse::linalg::lu::{factorize_symbolic_lu, NumericLu};
        faer_seq();
        let mat = faer::sparse::SparseColMat::<usize, f64>::try_new_from_triplets(dim, dim, trip)
            .map_err(|e| format!("triplets : {e:?}"))?;
        if symb.is_none() {
            *symb = Some(std::sync::Arc::new(
                factorize_symbolic_lu(mat.symbolic(), Default::default())
                    .map_err(|e| format!("LU symbolique : {e:?}"))?,
            ));
        }
        let symb = symb.clone().expect("symbolique posée ci-dessus");
        // le LU creux de faer sait paralléliser ses supernœuds : on l'y autorise
        // au-delà de CREUX_PAR inconnues (en dessous la synchronisation coûte plus
        // que le calcul — mesuré sur le dense à 143). Par APPEL, pas par réglage
        // global : cf. `LuCreux`.
        let par = if dim >= CREUX_PAR {
            faer::Par::rayon(0)
        } else {
            faer::Par::Seq
        };
        let mut num = NumericLu::new();
        let mut mem =
            MemBuffer::try_new(symb.factorize_numeric_lu_scratch::<f64>(par, Default::default()))
                .map_err(|e| format!("LU creux : mémoire ({e:?})"))?;
        symb.factorize_numeric_lu::<f64>(
            &mut num,
            mat.as_ref(),
            par,
            MemStack::new(&mut mem),
            Default::default(),
        )
        .map_err(|e| format!("LU creux : {e:?}"))?;
        let mem_solve = std::sync::Mutex::new(
            MemBuffer::try_new(symb.solve_in_place_scratch::<f64>(1, par))
                .map_err(|e| format!("LU creux : mémoire résolution ({e:?})"))?,
        );
        Ok(Facto::Creux(Box::new(LuCreux {
            symb,
            num,
            par,
            mem_solve,
        })))
    }
    fn solve(&self, b: &DVector<f64>) -> DVector<f64> {
        match self {
            Facto::Dense(lu) => {
                use faer::linalg::solvers::Solve;
                // en place sur une copie de b : une allocation au lieu de trois
                // (mesuré 7,2 µs par descente à 143 sur la tête S2, copies comprises)
                let mut x = b.clone();
                let n = x.len();
                lu.solve_in_place(faer::MatMut::from_column_major_slice_mut(
                    x.as_mut_slice(),
                    n,
                    1,
                ));
                x
            }
            Facto::Creux(lu) => {
                use faer::dyn_stack::MemStack;
                let mut mem = lu.mem_solve.lock().unwrap_or_else(|e| e.into_inner());
                let mut x = b.clone();
                let dim = x.len();
                faer::sparse::linalg::lu::LuRef::<usize, f64>::new_unchecked(&lu.symb, &lu.num)
                    .solve_in_place_with_conj(
                        faer::Conj::No,
                        faer::MatMut::from_column_major_slice_mut(x.as_mut_slice(), dim, 1),
                        lu.par,
                        MemStack::new(&mut mem),
                    );
                x
            }
        }
    }
}

/// produit jacobien × vecteur depuis les triplets (contrôle de la solution)
struct ChronoFac(std::time::Instant);
impl Drop for ChronoFac {
    fn drop(&mut self) {
        T_FAC.fetch_add(
            self.0.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }
}

fn trip_mv(trip: &Triplets, dim: usize, x: &DVector<f64>) -> DVector<f64> {
    let mut y = DVector::zeros(dim);
    for t in trip {
        y[t.row] += t.val * x[t.col];
    }
    y
}

/// Les constantes heuristiques de Newton (borne du prédicteur, rafraîchissement
/// du jacobien, plancher de stagnation, critère d'incrément) sont posées à la
/// main et MESURÉES sur les bancs. `VINKULUM_SEUILS=k` les multiplie toutes par
/// k : `bancs.seuils()` rejoue le corpus à ×0,5 et ×2 — un banc qui casse là
/// mesurait le seuil, pas la physique. Ce n'est pas un réglage utilisateur.
fn seuil(k: f64) -> f64 {
    static M: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    k * M.get_or_init(|| {
        std::env::var("VINKULUM_SEUILS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0)
    })
}

/// Un pas d'α-généralisé figé : ce qu'il faut pour évaluer le résidu en un x
/// quelconque, sur n'importe quelle copie du modèle (d'où le parallélisme).
struct Pas {
    am: f64,
    af: f64,
    gam: f64,
    bet: f64,
    sigma: f64,
    hh: f64,
    t1: f64,
    u0: DVector<f64>,
    poses0: Vec<Pose>,
    udot: DVector<f64>,
    a: DVector<f64>,
    /// GGL : G au DÉBUT du pas (lignes actives), et ∂Φ/∂t des cibles imposées
    /// pris à q₀ — O(h) sur la ligne, exact pour une cible à taux constant.
    ggl: Option<(Vec<(usize, usize, f64)>, DVector<f64>)>,
}

impl Pas {
    fn nz(&self) -> usize {
        self.ggl.as_ref().map_or(0, |(_, p)| p.len())
    }

    /// Pose le modèle à l'état (q₁, u₁) correspondant à x = (u̇₁, λ₁) ; rend (u̇₁, λ₁, a₁).
    fn etat(
        &self,
        me: &mut Modele,
        x: &DVector<f64>,
    ) -> Result<(DVector<f64>, DVector<f64>, DVector<f64>), String> {
        let (n, m) = (me.n(), me.m());
        let ud = x.rows(0, n).into_owned();
        let la = x.rows(n, m).into_owned();
        let (a1, dq) = self.a1_dq(x)?;
        let u1 = &self.u0 + (1.0 - self.gam) * self.hh * &self.a + self.gam * self.hh * &a1;
        me.avance(&self.poses0, &dq, self.hh);
        me.set_u(&u1);
        Ok((ud, la, a1))
    }

    /// Accélération algorithmique et incrément historique, avant sigma.
    fn a1_dq_classique(&self, x: &DVector<f64>) -> (DVector<f64>, DVector<f64>) {
        let ud = x.rows(0, self.u0.len());
        let a1 =
            ((1.0 - self.af) * ud + self.af * &self.udot - self.am * &self.a) / (1.0 - self.am);
        let mut dq = &self.u0 + (0.5 - self.bet) * self.hh * &self.a + self.bet * self.hh * &a1;
        if let Some((g0, phi_t)) = &self.ggl {
            // GGL : q̇ = u + Gᵀμ ; ζ = μ/(βh) pour que ∂Φ/∂ζ soit O(1) comme ∂Φ/∂u̇
            let ze = x.rows(self.u0.len() + phi_t.len(), phi_t.len());
            for &(i, j, v) in g0 {
                dq[j] += self.bet * self.hh * v * ze[i];
            }
        }
        (a1, dq)
    }

    /// Relation de configuration à sigma fixé, avec élimination locale exacte.
    fn a1_dq(&self, x: &DVector<f64>) -> Result<(DVector<f64>, DVector<f64>), String> {
        let (a1, mut dq) = self.a1_dq_classique(x);
        if self.sigma != 0.0 {
            let u1 = &self.u0 + (1.0 - self.gam) * self.hh * &self.a + self.gam * self.hh * &a1;
            let s = self.sigma * self.hh * self.bet / self.gam;
            for i in 0..self.u0.len() / 6 {
                let k = 6 * i + 3;
                let b = self.hh * V3::new(dq[k], dq[k + 1], dq[k + 2]);
                let w = V3::new(u1[k], u1[k + 1], u1[k + 2]);
                let theta = schema_lie::rotation(b, w, s)?;
                for c in 0..3 {
                    dq[k + c] = theta[c] / self.hh;
                }
            }
        }
        Ok((a1, dq))
    }

    /// Incrément classique : borne du prédicteur et J_l lorsque sigma=0.
    /// Pour sigma non nul, jacobien_ad recalcule la tangente implicite.
    fn dtheta(&self, x: &DVector<f64>, n: usize) -> Vec<V3> {
        let dq = self.a1_dq_classique(x).1;
        (0..n / 6)
            .map(|i| self.hh * V3::new(dq[6 * i + 3], dq[6 * i + 4], dq[6 * i + 5]))
            .collect()
    }

    /// Résidu en O(n + nnz) : M·u̇ par blocs de corps, Gᵀλ et Φ par éléments —
    /// aucune matrice dense (mesuré : à 1 800 inconnues le résidu dense
    /// coûtait 16 ms par pas, plus que le jacobien).
    fn residu(&self, me: &mut Modele, x: &DVector<f64>) -> DVector<f64> {
        // Un essai hors de la carte doit être rejeté par la recherche de pas,
        // sans modifier le modèle ni accepter un résidu d'un autre état.
        let Ok((ud, la, _)) = self.etat(me, x) else {
            return DVector::from_element(me.n() + me.m() + self.nz(), f64::NAN);
        };
        let (n, m) = (me.n(), me.m());
        let ze = x.rows(n + m, self.nz()).into_owned();
        me.residu_pose(
            &ud,
            &la,
            self.t1,
            self.bet * self.hh * self.hh,
            self.gam * self.hh,
            self.ggl.as_ref().map(|(_, phi_t)| (&ze, phi_t)),
        )
    }
}

impl Modele {
    /// Résidu sur l'état DÉJÀ POSÉ : m·v̇ − m·g, J_s·ω̇ + ω × J_s·ω − M_ext,
    /// moins les efforts, plus Gᵀλ ; puis Φ/`sc` sur les lignes de contrainte.
    /// Séparé de `Pas::residu` parce que le contrôle de pas a besoin de
    /// l'ÉQUILIBRE en un point quelconque du pas, pas du résidu du schéma.
    /// En O(n + nnz) : aucune matrice dense (mesuré, à 1 800 inconnues le
    /// résidu dense coûtait 16 ms par pas, plus que le jacobien).
    /// `sc` = βh² met à l'échelle les contraintes de POSITION ; `scv` = γh
    /// celles de VITESSE (non holonomes). Les deux sont choisies pour que la
    /// ligne rende le MÊME jacobien `G·c` — sans quoi Newton voit un gradient
    /// 2 000 fois trop petit sur les lignes non holonomes et diverge (mesuré).
    fn residu_pose(
        &self,
        ud: &DVector<f64>,
        la: &DVector<f64>,
        t: f64,
        sc: f64,
        scv: f64,
        ggl: Option<(&DVector<f64>, &DVector<f64>)>,
    ) -> DVector<f64> {
        let me = self;
        let (n, m) = (me.n(), me.m());
        let mut r = DVector::zeros(n + m + if ggl.is_some() { m } else { 0 });
        // corps : m·v̇ − m·g ; J_s·ω̇ + ω × J_s·ω − M_ext
        for (i, c) in me.corps.iter().enumerate() {
            if me.gele(i) {
                continue;
            }
            let k = 6 * i;
            let js = c.rot * c.j * c.rot.transpose();
            let wd = V3::new(ud[k + 3], ud[k + 4], ud[k + 5]);
            let rot = js * wd + c.w.cross(&(js * c.w));
            for a in 0..3 {
                r[k + a] = c.m * ud[k + a] - c.m * me.g[a];
                r[k + 3 + a] = rot[a];
            }
        }
        for (i, fo, mo) in &me.efforts {
            for a in 0..3 {
                r[6 * i + a] -= fo[a];
                r[6 * i + 3 + a] -= mo[a];
            }
        }
        for c in &me.couples {
            if !me.gel.is_empty() && [c.a, c.b].iter().flatten().all(|&i| me.gele(i)) {
                continue;
            }
            let tau = c.valeur(&me.corps, t);
            if let Some(ib) = c.b {
                for k in 0..3 {
                    r[6 * ib + 3 + k] -= tau[k];
                }
            }
            if let Some(ia) = c.a {
                for k in 0..3 {
                    r[6 * ia + 3 + k] += tau[k];
                }
            }
        }
        for p in &me.pales {
            let (fo, mo, _) = p.valeur(&me.corps, &me.inflows, me.vent);
            for k in 0..3 {
                r[6 * p.corps + k] -= fo[k];
                r[6 * p.corps + 3 + k] -= mo[k];
            }
        }
        for ct in &me.contacts {
            let ct_nl;
            let ct = if ct.nonlisse {
                ct_nl = Contact {
                    fn_impose: Some(la[ct.row].max(0.0)),
                    ..ct.clone()
                };
                &ct_nl
            } else {
                ct
            };
            let (fo, mo, fb, mb, _) = ct.valeur(&me.corps);
            for k in 0..3 {
                r[6 * ct.corps + k] -= fo[k];
                r[6 * ct.corps + 3 + k] -= mo[k];
            }
            if let Some(ib) = ct.b {
                for k in 0..3 {
                    r[6 * ib + k] -= fb[k];
                    r[6 * ib + 3 + k] -= mb[k];
                }
            }
        }
        for b in &me.poutres {
            if me.gele(b.a) && me.gele(b.b) {
                continue;
            }
            let fl = b.forces(&me.corps);
            for k in 0..6 {
                r[6 * b.a + k] -= fl[k];
                r[6 * b.b + k] -= fl[6 + k];
            }
        }
        for se in &me.supers {
            if !me.gel.is_empty() && se.noeuds.iter().all(|&i| me.gele(i)) {
                continue;
            }
            let fl = se.efforts(&me.corps, se.beta);
            for (i, &ni) in se.noeuds.iter().enumerate() {
                for k in 0..6 {
                    r[6 * ni + k] -= fl[6 * i + k];
                }
            }
        }
        // éléments : Φ et Gᵀλ, blocs locaux
        let mut row = 0;
        for (ie, e) in me.elems.iter().enumerate() {
            if me.gel_elem.get(ie).copied().unwrap_or(false) {
                let ne = e.n();
                r.rows_mut(n + row, ne).copy_from(&la.rows(row, ne));
                row += ne;
                continue;
            }
            let (phi, g) = match e {
                Elem::L(l) => l.phi_g(&me.corps, t),
                Elem::D(d) => d.phi_g(&me.corps),
                Elem::V(vi) => match vi.phi_g(&me.corps) {
                    Ok(v) => v,
                    Err(_) => return DVector::from_element(r.len(), 1e30),
                },
                Elem::C(c) => c.phi_g(&me.corps),
                Elem::K(k) => me.contacts[k.ct].phi_g(&me.corps),
                Elem::E(en) => match en.phi_g(&me.corps) {
                    Ok(v) => v,
                    // un angle d'engrenage qui saute d'un demi-tour dans un essai
                    // de Newton : résidu géant, Newton recule
                    Err(_) => return DVector::from_element(r.len(), 1e30),
                },
            };
            // Φ_t : demandé par la ligne GGL, ET par une liaison non holonome
            // PILOTÉE (cible_r/cible_t) — sans lui la cible y était ignorée et
            // la liaison imposait G·u = 0 (trouvé le 6 sept. : une section en
            // tangage imposé en vitesse ne bougeait pas)
            let pt_l = match (e, ggl) {
                (Elem::L(l), Some(_)) => Some(l.phi_dt(&me.corps, t)),
                (Elem::L(l), None) if l.nh && (l.cible_r.is_some() || l.cible_t.is_some()) => {
                    Some(l.phi_dt(&me.corps, t))
                }
                _ => None,
            };
            let ne = e.n();
            let cs = e.corps();
            let nh = matches!(e, Elem::L(l) if l.nh);
            let fbr = matches!(e, Elem::K(_));
            for i in 0..ne {
                let actif = me.actif.get(row + i).copied().unwrap_or(true);
                if let Some((ze, _)) = ggl {
                    // ligne GGL : G·u + Φ_t = 0 ; sans objet (ζ = 0) si la ligne
                    // est inactive ou déjà en vitesse (non holonome)
                    r[n + m + row + i] = ze[row + i];
                }
                if !actif {
                    r[n + row + i] = la[row + i];
                    continue;
                }
                let gu = || {
                    let mut s = 0.0;
                    for (kb, cb) in cs.iter().enumerate() {
                        if let Some(ib) = cb {
                            for c in 0..6 {
                                s += g[(i, 6 * kb + c)] * me.corps[*ib].u6(c);
                            }
                        }
                    }
                    s
                };
                // NON HOLONOME : on impose G·u + Φ_t = 0, pas Φ. Le jacobien
                // est le même ; c'est le NIVEAU qui change. Φ_t porte la cible
                // (vitesse imposée EXACTE sur le pas) ; nul sans cible.
                r[n + row + i] = if nh {
                    (gu() + pt_l.as_ref().map_or(0.0, |v| v[i])) / scv
                } else if fbr {
                    // Moreau–Jean : contact inactif ⇒ λ = 0 ; actif ⇒ Fischer–Burmeister
                    // sur (ġap₁ + e·ġap₀ + gap₀/h, λ) — nulle ssi a ≥ 0, λ ≥ 0, a·λ = 0
                    let Elem::K(k) = e else { unreachable!() };
                    let (actif_k, rhs) = me.k_pas.get(k.ct).copied().unwrap_or((false, 0.0));
                    if !actif_k {
                        r[n + row + i] = la[row + i];
                        continue;
                    }
                    let (a, b) = ((-gu() + rhs) / scv, la[row + i] * me.s_fb / sc);
                    a + b - (a * a + b * b).sqrt()
                } else {
                    phi[i] / sc
                };
                if let (Some((_, phi_t)), false) = (ggl, nh || fbr) {
                    let pt = pt_l.as_ref().map_or(phi_t[row + i], |v| v[i]);
                    r[n + m + row + i] = (gu() + pt) / scv;
                }
                let li = la[row + i];
                for (kb, cb) in cs.iter().enumerate() {
                    if let Some(ib) = cb {
                        for c in 0..6 {
                            r[6 * ib + c] += g[(i, 6 * kb + c)] * li;
                        }
                    }
                }
            }
            row += ne;
        }
        // multi-rythme : un corps gelé n'est pas une inconnue — sa ligne rend u̇
        for (i, _) in me.corps.iter().enumerate() {
            if me.gele(i) {
                for c in 0..6 {
                    r[6 * i + c] = ud[6 * i + c];
                }
            }
        }
        r
    }
}

// ── binding Python ───────────────────────────────────────────────────────────
fn valide(v: &[f64], nom: &str) -> PyResult<()> {
    fini(v, nom).map_err(pyo3::exceptions::PyValueError::new_err)
}
fn valide_axe(v: [f64; 3]) -> PyResult<()> {
    valide(&v, "axe")?;
    let n = v3(v).norm();
    if n == 0.0 || !n.is_finite() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "axe de norme finie non nulle requis",
        ));
    }
    Ok(())
}
fn v3(v: [f64; 3]) -> V3 {
    V3::new(v[0], v[1], v[2])
}
fn m3(m: [f64; 9]) -> M3 {
    M3::from_row_slice(&m)
}
fn loi_de(l: (String, Vec<f64>)) -> PyResult<Loi> {
    let (genre, p) = l;
    let taille_valide = match genre.as_str() {
        "constante" => p.len() == 1,
        "lineaire" => p.len() == 2,
        "table" => !p.is_empty() && p.len().is_multiple_of(2),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "loi inconnue : {genre}"
            )));
        }
    };
    if !taille_valide || p.iter().any(|x| !x.is_finite()) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "loi {genre} : paramètres finis requis ; constante = 1 valeur, \
             lineaire = 2 valeurs, table = au moins une paire (t, valeur)"
        )));
    }
    match genre.as_str() {
        "constante" => Ok(Loi::Constante(p[0])),
        "lineaire" => Ok(Loi::Lineaire {
            a0: p[0],
            taux: p[1],
        }),
        "table" => Loi::table(p.chunks(2).map(|c| (c[0], c[1])).collect())
            .map_err(pyo3::exceptions::PyValueError::new_err),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "loi inconnue : {genre}"
        ))),
    }
}

/// Historique du dernier calcul statique, distinct de l'état physique restaurable.
struct BilanStatique {
    statut: &'static str,
    tol: f64,
    strict: bool,
    residu_libre: Option<f64>,
    contraintes: Option<f64>,
    echelle_force: Option<f64>,
    evaluations: usize,
    tentatives: usize,
    paliers: usize,
    palier: usize,
    paliers_stagnation: usize,
    etat_restaure: bool,
    message: Option<String>,
}

/// Le modèle vu de Python : on déclare, on simule, on lit.
#[pyclass]
pub struct Noyau {
    mo: Modele,
    bilan_statique: Option<BilanStatique>,
}

#[pymethods]
impl Noyau {
    #[new]
    #[pyo3(signature = (g = [0.0, 0.0, -9.80665]))]
    fn new(g: [f64; 3]) -> PyResult<Self> {
        if g.iter().any(|x| !x.is_finite()) {
            // sondage du 7 sept. : une gravité NaN BLOQUAIT `simule` (SVD)
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "gravité non finie ({g:?})"
            )));
        }
        Ok(Noyau {
            mo: Modele::new(v3(g)),
            bilan_statique: None,
        })
    }

    /// Corps rigide ; `j` = tenseur 3×3 en ligne (axes corps, au CdM), `rot` = R en ligne.
    #[pyo3(signature = (nom, m, j, r, rot = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], v = [0.0; 3], w = [0.0; 3]))]
    fn corps(
        &mut self,
        nom: String,
        m: f64,
        j: [f64; 9],
        r: [f64; 3],
        rot: [f64; 9],
        v: [f64; 3],
        w: [f64; 3],
    ) -> PyResult<usize> {
        valide(&j, "inertie")?;
        let jm = m3(j);
        // UN TENSEUR D'INERTIE EST SYMÉTRIQUE. Ne pas le vérifier, c'était
        // lire la moitié de ce que l'utilisateur a écrit et jeter l'autre en
        // silence — trouvé par sondage le 3 sept. sur des modèles PLAUSIBLES
        // mal posés (et non absurdes), la classe de fautes qu'un utilisateur
        // nouveau commet le premier jour.
        let asym = (jm - jm.transpose()).amax();
        if asym > 1e-12 * jm.amax().max(1e-300) {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "{nom} : tenseur d'inertie NON SYMÉTRIQUE (écart {asym:.3e}) — \
                 un tenseur d'inertie l'est toujours ; vérifiez l'ordre des termes croisés"
            )));
        }
        let jm = jm * 0.5 + jm.transpose() * 0.5;
        // Une pose NON FINIE ne fait pas échouer le solveur : elle le fait
        // BOUCLER (mesuré le 3 sept. — un `sqrt` négatif dans un montage
        // d'inverseur de Peaucellier a produit un NaN, et `assemble` n'a
        // jamais rendu la main, même sur quatre corps). Les lecteurs URDF et
        // MJCF gardaient déjà leurs entrées ; l'API, non.
        for (v, quoi) in [(&r, "position"), (&v, "vitesse"), (&w, "rotation")] {
            if v.iter().any(|x| !x.is_finite()) {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "{nom} : {quoi} non finie ({v:?})"
                )));
            }
        }
        if rot.iter().any(|x| !x.is_finite()) || j.iter().any(|x| !x.is_finite()) {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "{nom} : orientation ou inertie non finie"
            )));
        }
        // ET UNE ROTATION EST ORTHOGONALE. Une matrice saisie ou arrondie à la
        // main ne l'est pas tout à fait, et le noyau la ré-orthonormalise en
        // douce à chaque pas : le modèle intégré n'est alors plus celui écrit.
        let rm = m3(rot);
        let ort = (rm.transpose() * rm - M3::identity()).amax();
        if ort > 1e-9 || (rm.determinant() - 1.0).abs() > 1e-9 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "{nom} : orientation hors SO(3) (‖RᵀR − I‖ = {ort:.3e}, det R = {:.3e}) — \
                 rotation orthogonale de déterminant +1 requise",
                rm.determinant()
            )));
        }
        if !(m > 0.0) || !m.is_finite() || jm.amax() == 0.0 || (jm / jm.amax()).cholesky().is_none()
        {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "{nom} : masse ou inertie non définie positive"
            )));
        }
        self.mo.corps.push(Corps {
            r_bas: crate::V3::zeros(),
            nom,
            m,
            j: jm,
            r: v3(r),
            rot: m3(rot),
            v: v3(v),
            w: v3(w),
        });
        Ok(self.mo.corps.len() - 1)
    }

    /// Valide un couple de références de corps AVANT de construire l'élément.
    ///
    /// Trouvé par sondage le 3 sept., et les deux fautes sont de nature
    /// différente : un indice hors bornes faisait PANIQUER Rust au travers de
    /// PyO3 (`index out of bounds`, pas d'exception Python) ; et `a == b` —
    /// bâti compris — était ACCEPTÉ en silence, en produisant une contrainte
    /// identiquement nulle. Le corps tombait librement sous une liaison qu'on
    /// croyait avoir posée : le pire des deux, parce que rien ne le dit.
    fn _refs(&self, nom: &str, a: Ref, b: Ref) -> PyResult<()> {
        let n = self.mo.corps.len();
        for (r, q) in [(a, "a"), (b, "b")] {
            if let Some(i) = r {
                if i >= n {
                    return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                        "« {nom} » : corps {q} = {i}, mais le modèle n'en a que {n}"
                    )));
                }
            }
        }
        if a == b {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "« {nom} » : a et b désignent le même corps ({}) — la contrainte serait \
                 identiquement nulle, donc muette",
                match a {
                    Some(i) => i.to_string(),
                    None => "le bâti".into(),
                }
            )));
        }
        Ok(())
    }

    /// Valide une référence de corps SIMPLE (pale, poutre, contact, effort).
    fn _ref1(&self, nom: &str, i: usize) -> PyResult<()> {
        let n = self.mo.corps.len();
        if i >= n {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "« {nom} » : corps {i}, mais le modèle n'en a que {n}"
            )));
        }
        Ok(())
    }

    /// LIAISON GÉNÉRIQUE — six degrés à la carte, entre `a` et `b`
    /// (`None` = le bâti).
    ///
    /// `pa` : le point d'attache, dans le repère LOCAL de `a` (ou en monde si
    /// `a` est le bâti). `ra` : le repère de la liaison, exprimé dans celui de
    /// `a` — c'est LUI qui définit les axes 0/1/2 que bloquent `bloque_t` et
    /// `bloque_r`, pas les axes du monde.
    ///
    /// `bloque_t` / `bloque_r` : les composantes de translation et de rotation
    /// relatives à annuler. Tout bloqué = encastrement ; `bloque_r = []` =
    /// rotule ; `bloque_r = [1, 2]` = pivot autour de l'axe 0 de `ra` ;
    /// `bloque_t = [1, 2]` = glissière. Choisir COMPOSANTE PAR COMPOSANTE est
    /// ce qui permet de décrire un mécanisme surcontraint sans le rendre
    /// surdéterminé — les benchmarks publiés en usent largement.
    ///
    /// `cible_t` / `cible_r` : (direction ou axe, loi) pour PILOTER la liaison
    /// au lieu de la bloquer à zéro. Loi = `("constante", [c])`,
    /// `("lineaire", [a0, taux])` ou `("table", [t0, v0, t1, v1, …])`. L'axe
    /// est exprimé dans le repère `ra`. Avec `nh=True` la cible est imposée en
    /// VITESSE (sa dérivée, pente de la table sur le pas) : c'est le mouvement
    /// imposé qu'il faut à une force qui dérive une vitesse (masse ajoutée) —
    /// imposée en position, la vitesse du corps oscille d'un pas à l'autre.
    ///
    /// ⚠ **LA LIAISON SE FERME SUR LA POSE COURANTE.** Le point homologue sur
    /// `b` n'est pas demandé : il est DÉDUIT pour que Φ soit nul à la
    /// création. Poser une liaison entre deux corps qui ne se rencontrent pas
    /// donne donc un lien rigide À DISTANCE, sans erreur — c'est voulu (on
    /// n'exige pas un assemblage parfait, et `assemble()` existe pour le
    /// reste), mais une faute de coordonnées passe en silence.
    ///
    /// Rend l'indice de l'élément ; `reactions()` donne ses multiplicateurs.
    #[pyo3(signature = (nom, a, b, pa = [0.0; 3], ra = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], bloque_t = vec![0, 1, 2], bloque_r = vec![0, 1, 2], cible_t = None, cible_r = None, nh = false))]
    #[allow(clippy::too_many_arguments)]
    fn liaison(
        &mut self,
        nom: String,
        a: Option<usize>,
        b: Option<usize>,
        pa: [f64; 3],
        ra: [f64; 9],
        bloque_t: Vec<usize>,
        bloque_r: Vec<usize>,
        cible_t: Option<([f64; 3], (String, Vec<f64>))>,
        cible_r: Option<([f64; 3], (String, Vec<f64>))>,
        nh: bool,
    ) -> PyResult<usize> {
        self._refs(&nom, a, b)?;
        if bloque_t.iter().chain(&bloque_r).any(|&i| i >= 3) {
            return Err(pyo3::exceptions::PyIndexError::new_err(
                "liaison : les axes de bloque_t et bloque_r sont 0, 1 ou 2",
            ));
        }
        if [&bloque_t, &bloque_r].iter().any(|axes| {
            axes.iter()
                .enumerate()
                .any(|(i, axe)| axes[..i].contains(axe))
        }) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "liaison : axes répétés dans bloque_t ou bloque_r",
            ));
        }
        valide(&pa, "point de liaison")?;
        valide(&ra, "repère de liaison")?;
        let rm = m3(ra);
        if (rm.transpose() * rm - M3::identity()).norm() > 1e-9
            || (rm.determinant() - 1.0).abs() > 1e-9
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "repère de liaison hors SO(3)",
            ));
        }
        for (d, _) in cible_t.iter().chain(cible_r.iter()) {
            valide_axe(*d)?;
        }
        let (ra_, rota) = pose_de(&self.mo.corps, a);
        let (rb_, rotb) = pose_de(&self.mo.corps, b);
        let pa_ = v3(pa);
        let ra_m = m3(ra);
        let pw = ra_ + rota * pa_;
        let pb = rotb.transpose()
            * ((pw - rb_) + partie_basse(&self.mo.corps, a) - partie_basse(&self.mo.corps, b));
        let rb = rotb.transpose() * (rota * ra_m);
        let ct = match cible_t {
            Some((d, l)) => Some((v3(d), loi_de(l)?)),
            None => None,
        };
        let cr = match cible_r {
            Some((ax, l)) => Some((v3(ax).normalize(), loi_de(l)?)),
            None => None,
        };
        self.mo.elems.push(Elem::L(Liaison {
            nom,
            a,
            b,
            pa: pa_,
            ra: ra_m,
            pb,
            rb,
            bt: bloque_t,
            br: bloque_r,
            cible_t: ct,
            cible_r: cr,
            nh,
        }));
        Ok(self.mo.elems.len() - 1)
    }

    /// Bielle à deux rotules ; `l` None = longueur initiale.
    #[pyo3(signature = (nom, a, b, pa, pb, l = None))]
    fn distance(
        &mut self,
        nom: String,
        a: Option<usize>,
        b: Option<usize>,
        pa: [f64; 3],
        pb: [f64; 3],
        l: Option<f64>,
    ) -> PyResult<usize> {
        self._refs(&nom, a, b)?;
        let (ra_, rota) = pose_de(&self.mo.corps, a);
        let (rb_, rotb) = pose_de(&self.mo.corps, b);
        valide(&pa, "point de bielle")?;
        valide(&pb, "point de bielle")?;
        let l0 = ((rb_ + rotb * v3(pb)) - (ra_ + rota * v3(pa)) + partie_basse(&self.mo.corps, b)
            - partie_basse(&self.mo.corps, a))
        .norm();
        if !l.unwrap_or(l0).is_finite() || l.unwrap_or(l0) <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "bielle : longueur finie > 0 requise",
            ));
        }
        self.mo.elems.push(Elem::D(Distance {
            nom,
            a,
            b,
            pa: v3(pa),
            pb: v3(pb),
            l: l.unwrap_or(l0),
        }));
        Ok(self.mo.elems.len() - 1)
    }

    /// CARDAN (croix de Hooke) : les axes `axe_a` (repère de a) et `axe_b`
    /// (repère de b) restent PERPENDICULAIRES — une contrainte scalaire. Le
    /// joint complet = cette contrainte PLUS une rotule au centre de la croix
    /// (`liaison` avec `bloque_t = [0,1,2]`, `bloque_r = []`).
    ///
    /// Refuse deux axes déjà colinéaires : la contrainte y est stationnaire
    /// (son gradient est le produit vectoriel, donc nul) et le jacobien perd
    /// son rang.
    #[pyo3(signature = (nom, a, b, axe_a, axe_b))]
    fn cardan(
        &mut self,
        nom: String,
        a: Option<usize>,
        b: Option<usize>,
        axe_a: [f64; 3],
        axe_b: [f64; 3],
    ) -> PyResult<usize> {
        self._refs(&nom, a, b)?;
        valide_axe(axe_a)?;
        valide_axe(axe_b)?;
        let (na, nb) = (v3(axe_a).normalize(), v3(axe_b).normalize());
        let (_, ra) = pose_de(&self.mo.corps, a);
        let (_, rb) = pose_de(&self.mo.corps, b);
        let c = (ra * na).cross(&(rb * nb)).norm();
        if c < 1e-6 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "cardan « {nom} » : axes colinéaires au montage (|n_a × n_b| = {c:.2e}) \
                 — le gradient de la contrainte est nul, le jacobien perd son rang"
            )));
        }
        self.mo.elems.push(Elem::C(Cardan { nom, a, b, na, nb }));
        Ok(self.mo.elems.len() - 1)
    }

    /// VIS–ÉCROU : la translation relative le long de `axe` est liée à la
    /// rotation relative autour du même axe. `pas` en MÈTRES PAR TOUR (la
    /// convention d'atelier ; le noyau travaille en m/rad et convertit).
    #[pyo3(signature = (nom, a, b, axe, pas))]
    fn vis(
        &mut self,
        nom: String,
        a: Option<usize>,
        b: Option<usize>,
        axe: [f64; 3],
        pas: f64,
    ) -> PyResult<usize> {
        self._refs(&nom, a, b)?;
        valide_axe(axe)?;
        valide(&[pas], "pas de vis")?;
        if pas == 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "vis : pas nul — la contrainte dégénérerait en simple glissière",
            ));
        }
        let (pa_, ra) = pose_de(&self.mo.corps, a);
        let (pb_, rb) = pose_de(&self.mo.corps, b);
        let n = v3(axe).normalize();
        let nw = ra * n;
        let pr = pas / (2.0 * std::f64::consts::PI);
        // c₀ pris sur la pose courante : Φ = 0 à la création, comme partout
        let c0 = nw.dot(&(pb_ - pa_))
            + nw.dot(&(partie_basse(&self.mo.corps, b) - partie_basse(&self.mo.corps, a)));
        self.mo.elems.push(Elem::V(Vis {
            nom,
            a,
            b,
            n,
            pas: pr,
            ref_rel: ra.transpose() * rb,
            c0,
            prev: 0.0,
        }));
        Ok(self.mo.elems.len() - 1)
    }

    /// Engrenage θ_a = rapport·θ_b, axes `axe_a`/`axe_b` dans le repère du
    /// porteur (`porteur` None = bâti). Rapport négatif pour deux roues
    /// extérieures. Rend l'indice de l'élément ; sa réaction (`reactions`) est
    /// le couple de denture sur a.
    #[pyo3(signature = (nom, a, b, axe_a, axe_b, rapport, porteur = None))]
    fn engrenage(
        &mut self,
        nom: String,
        a: Option<usize>,
        b: Option<usize>,
        axe_a: [f64; 3],
        axe_b: [f64; 3],
        rapport: f64,
        porteur: Option<usize>,
    ) -> PyResult<usize> {
        self._refs(&nom, a, b)?;
        valide_axe(axe_a)?;
        valide_axe(axe_b)?;
        valide(&[rapport], "rapport engrenage")?;
        if let Some(i) = porteur {
            self._ref1(&nom, i)?;
        }
        let (_, rc) = pose_de(&self.mo.corps, porteur);
        let (_, ra) = pose_de(&self.mo.corps, a);
        let (_, rb) = pose_de(&self.mo.corps, b);
        self.mo.elems.push(Elem::E(Engrenage {
            nom,
            a,
            b,
            c: porteur,
            na: v3(axe_a).normalize(),
            nb: v3(axe_b).normalize(),
            rapport,
            ref_a: rc.transpose() * ra,
            ref_b: rc.transpose() * rb,
            prev: (0.0, 0.0),
        }));
        Ok(self.mo.elems.len() - 1)
    }

    /// Couple à loi entre a et b autour de `axe` (repère de a ; monde si bâti).
    /// `loi` = ("constant", [τ]) | ("ressort", [k, c, θ0]) | ("pd", [kp, kd, q_max, w_nl]) |
    /// ("gouverneur", [kg, q_max]) ; `cible` = loi de consigne (angle pour pd, vitesse
    /// pour gouverneur) au format des cibles de liaison.
    #[pyo3(signature = (nom, a, b, axe, loi, cible = None))]
    fn couple(
        &mut self,
        nom: String,
        a: Option<usize>,
        b: Option<usize>,
        axe: [f64; 3],
        loi: (String, Vec<f64>),
        cible: Option<(String, Vec<f64>)>,
    ) -> PyResult<usize> {
        self._refs(&nom, a, b)?;
        valide_axe(axe)?;
        let (genre, p) = loi;
        valide(&p, "loi de couple")?;
        let attendu = match genre.as_str() {
            "constant" => 1,
            "ressort" => 3,
            "pd" | "butee" => 4,
            "gouverneur" => 2,
            _ => 0,
        };
        if attendu != 0 && p.len() != attendu {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "loi {genre} : {attendu} paramètres requis"
            )));
        }
        if (genre == "pd" && (p[2] <= 0.0 || p[3] <= 0.0)) || (genre == "gouverneur" && p[1] < 0.0)
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "couple : pd exige q_max et w_nl > 0 ; gouverneur exige q_max ≥ 0",
            ));
        }

        let cib = || -> PyResult<Loi> {
            loi_de(
                cible
                    .clone()
                    .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("cible requise"))?,
            )
        };
        let l = match genre.as_str() {
            "constant" => LoiCouple::Constant(p[0]),
            "ressort" => LoiCouple::Ressort {
                k: p[0],
                c: p[1],
                theta0: p[2],
            },
            "pd" => LoiCouple::PD {
                kp: p[0],
                kd: p[1],
                q_max: p[2],
                w_nl: p[3],
                cible: cib()?,
            },
            "gouverneur" => LoiCouple::Gouverneur {
                kg: p[0],
                q_max: p[1],
                cible: cib()?,
            },
            "butee" => {
                if p[2] >= p[3] {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "butée : θ_min doit être < θ_max (sinon la plage est vide et la liaison bloquée)"));
                }
                LoiCouple::Butee {
                    k: p[0],
                    c: p[1],
                    theta_min: p[2],
                    theta_max: p[3],
                }
            }
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "loi de couple inconnue : {genre}"
                )))
            }
        };
        let (_, ra) = pose_de(&self.mo.corps, a);
        let (_, rb) = pose_de(&self.mo.corps, b);
        self.mo.couples.push(Couple {
            nom,
            a,
            b,
            axe: v3(axe).normalize(),
            ref_rel: ra.transpose() * rb,
            prev: 0.0,
            loi: l,
        });
        Ok(self.mo.couples.len() - 1)
    }

    /// Modes propres du système contraint autour de l'état courant :
    /// [(fréquence Hz, forme modale sur les 6n degrés)], les `combien` premiers.
    #[pyo3(signature = (combien = 6, t = None))]
    fn modes(&self, combien: usize, t: Option<f64>) -> PyResult<Vec<(f64, Vec<f64>)>> {
        self.mo
            .modes(t.unwrap_or(self.mo.t), combien)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)
    }

    /// Modes locaux creux : valeurs propres signées proches du décalage,
    /// formes normalisées en masse, réactions et résidus. SciPy requis.
    /// Le décalage est en s⁻² ; None choisit un petit décalage négatif
    /// relatif à la raideur en métrique de masse, pour inclure les corps libres.
    /// Le problème utilise la partie symétrique de K, comme modes ; il ne
    /// remplace pas le spectre amorti ni une conclusion de stabilité globale.
    #[pyo3(signature = (combien = 6, t = None, decalage = None, tol = 1e-8, maxiter = 1000))]
    fn modes_creux(
        slf: PyRef<'_, Self>,
        py: Python<'_>,
        combien: usize,
        t: Option<f64>,
        decalage: Option<f64>,
        tol: f64,
        maxiter: usize,
    ) -> PyResult<Py<PyAny>> {
        let module = py.import("vinkulum.modes_creux")?;
        Ok(module
            .getattr("analyse")?
            .call1((slf, combien, t, decalage, tol, maxiter))?
            .unbind())
    }

    /// Domaine du pont adjoint mécanique : les états de fluide et les
    /// événements de contact ne figurent pas dans son état de schéma.
    fn _verifie_domaine_adjoint(&self) -> PyResult<()> {
        let mut raisons = self.mo.limites_domaine_spectre();
        if self
            .mo
            .elems
            .iter()
            .any(|e| matches!(e, Elem::L(l) if l.nh))
        {
            raisons.push("liaisons non holonomes : le pont dérive les contraintes de position");
        }
        if !raisons.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "PontNoyau hors domaine : {}",
                raisons.join(" ; ")
            )));
        }
        Ok(())
    }

    /// (K, C, M, G) en stockage CSC, autour de l'état courant.
    /// Chaque matrice est (n_lignes, n_colonnes, indptr, indices, donnees).
    /// Construire scipy.sparse.csc_matrix((donnees, indices, indptr),
    /// shape=(n_lignes, n_colonnes)) si SciPy est disponible.
    /// Même champ et mêmes conventions que k_c_m_z : précontrainte Gᵀλ,
    /// rotations perturbées à gauche, états internes figés. Aucun Z dense
    /// n'est construit ; les lignes de contact complémentaire restent nulles.
    /// Cette linéarisation locale ne dérive pas les changements de contact.
    #[pyo3(signature = (t = None))]
    fn k_c_m_g_creux(
        &self,
        t: Option<f64>,
    ) -> PyResult<(
        analyse::CscPython,
        analyse::CscPython,
        analyse::CscPython,
        analyse::CscPython,
    )> {
        self.mo
            .linearisation_creuse(t.unwrap_or(self.mo.t))
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)
    }

    /// Facteurs matériels des seules poutres, à l'état courant, sans K ni Z.
    /// Dictionnaire d, deformations, masse, contraintes, poutres, phi, t,
    /// infos_domaine, etat. Les matrices utilisent le même tuple CSC que
    /// k_c_m_g_creux. D a 6 lignes par poutre et 6 colonnes par corps ; aucune
    /// racine encastrée n'est retirée. poutres=[(nom,a,b)] suit l'ordre des
    /// lignes (gamma_x,y,z puis kappa_x,y,z, pondérés par sqrt(L C)).
    /// Colonnes : translations puis rotations monde perturbées à gauche.
    /// U_poutres=||deformations||²/2, f_poutres=-Dᵀ deformations.
    /// DᵀD est la contribution matérielle : même avec deformations non nul,
    /// l'extraction reste disponible, sans ajouter la tangente géométrique,
    /// les réactions, la gyroscopie ou les autres éléments. M est la masse
    /// spatiale complète, symétrisée coefficient par coefficient après R J Rᵀ ;
    /// G et phi gardent les lignes de contact non lisse
    /// nulles et n'appliquent pas le masque d'activité dynamique.
    /// infos_domaine signale les contributions et états hors périmètre ;
    /// cette liste ne certifie ni linéarité, ni équilibre, ni positivité.
    /// etat est une copie descriptive (corps, contacts, charges, contraintes,
    /// multiplicateurs), pas une sauvegarde de redémarrage. t ne change que
    /// l'évaluation des lois de contraintes ; les poses ne sont pas avancées.
    #[pyo3(signature = (t = None))]
    fn facteurs_materiels_poutres<'py>(
        &self,
        py: Python<'py>,
        t: Option<f64>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::PyDict;
        let t = t.unwrap_or(self.mo.t);
        let f = py
            .detach(|| self.mo.facteurs_materiels_poutres(t))
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let out = PyDict::new(py);
        out.set_item("d", f.d)?;
        out.set_item("deformations", f.deformations)?;
        out.set_item("masse", f.masse)?;
        out.set_item("contraintes", f.contraintes)?;
        out.set_item("poutres", f.poutres)?;
        out.set_item("phi", f.phi)?;
        out.set_item("t", t)?;
        out.set_item("infos_domaine", f.infos_domaine)?;

        let etat = PyDict::new(py);
        etat.set_item("temps_modele", self.mo.t)?;
        let mut corps = Vec::with_capacity(self.mo.corps.len());
        for (i, c) in self.mo.corps.iter().enumerate() {
            let b = PyDict::new(py);
            b.set_item("indice", i)?;
            b.set_item("nom", &c.nom)?;
            b.set_item("position", c.r.as_slice())?;
            b.set_item("position_basse", c.r_bas.as_slice())?;
            b.set_item(
                "rotation",
                (0..3)
                    .map(|j| (0..3).map(|k| c.rot[(j, k)]).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
            )?;
            b.set_item("vitesse", c.v.as_slice())?;
            b.set_item("vitesse_angulaire", c.w.as_slice())?;
            corps.push(b);
        }
        etat.set_item("corps", corps)?;
        etat.set_item("multiplicateurs", self.mo.lam.as_slice())?;
        etat.set_item("contraintes_actives", &self.mo.actif)?;
        etat.set_item("gravite", self.mo.g.as_slice())?;
        etat.set_item(
            "efforts",
            self.mo
                .efforts
                .iter()
                .map(|(i, f, m)| (*i, f.as_slice().to_vec(), m.as_slice().to_vec()))
                .collect::<Vec<_>>(),
        )?;
        let mut contraintes = Vec::with_capacity(self.mo.elems.len());
        let mut row = 0;
        for e in &self.mo.elems {
            let c = PyDict::new(py);
            let nature = match e {
                Elem::L(_) => "liaison",
                Elem::D(_) => "distance",
                Elem::E(_) => "engrenage",
                Elem::V(_) => "vis",
                Elem::C(_) => "cardan",
                Elem::K(_) => "contact_non_lisse",
            };
            c.set_item("nom", e.nom())?;
            c.set_item("nature", nature)?;
            c.set_item("premiere_ligne", row)?;
            c.set_item("nombre_lignes", e.n())?;
            c.set_item("corps", e.corps())?;
            c.set_item("non_holonome", matches!(e, Elem::L(l) if l.nh))?;
            c.set_item(
                "pilotee",
                matches!(e, Elem::L(l) if l.cible_t.is_some() || l.cible_r.is_some()),
            )?;
            contraintes.push(c);
            row += e.n();
        }
        etat.set_item("contraintes", contraintes)?;
        let mut contacts = Vec::with_capacity(self.mo.contacts.len());
        for ct in &self.mo.contacts {
            let c = PyDict::new(py);
            c.set_item("nom", &ct.nom)?;
            c.set_item("corps", (ct.corps, ct.b))?;
            c.set_item("non_lisse", ct.nonlisse)?;
            c.set_item("ligne", ct.nonlisse.then_some(ct.row))?;
            c.set_item("automatique", ct.auto)?;
            c.set_item("sortie", ct.sortie)?;
            c.set_item("force_normale_imposee", ct.fn_impose)?;
            c.set_item("point_maillage", ct.q_maille.as_slice())?;
            c.set_item("abscisse_axe", ct.s_axe)?;
            contacts.push(c);
        }
        etat.set_item("contacts", contacts)?;
        etat.set_item(
            "couples",
            self.mo
                .couples
                .iter()
                .map(|e| e.nom.as_str())
                .collect::<Vec<_>>(),
        )?;
        etat.set_item(
            "superelements",
            self.mo
                .supers
                .iter()
                .map(|e| e.nom.as_str())
                .collect::<Vec<_>>(),
        )?;
        etat.set_item(
            "pales",
            self.mo
                .pales
                .iter()
                .map(|e| e.nom.as_str())
                .collect::<Vec<_>>(),
        )?;
        etat.set_item("nombre_inflows", self.mo.inflows.len())?;
        etat.set_item("nombre_spheres_appariement", self.mo.spheres.len())?;
        etat.set_item("loi_contact_configuree", self.mo.loi_contact.is_some())?;
        out.set_item("etat", etat)?;
        Ok(out)
    }

    /// (K, M, base admissible Z) autour de l'état courant, en ligne — de quoi
    /// faire ses propres réductions et diagnostics en Python sans réécrire la
    /// raideur tangente. K est la raideur du champ de forces effectif
    /// (forces + précontrainte Gᵀλ + ω × J_s ω), M la masse spatiale, Z une
    /// base orthonormale du noyau des contraintes (6n × p).
    #[pyo3(signature = (t = None))]
    fn k_m_z(&self, t: Option<f64>) -> PyResult<(Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        let t = t.unwrap_or(self.mo.t);
        let k = self
            .mo
            .raideur(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let m = self.mo.masse_spatiale();
        let z = self
            .mo
            .base_admissible(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let lig = |a: &DMatrix<f64>| {
            (0..a.nrows())
                .map(|i| (0..a.ncols()).map(|j| a[(i, j)]).collect())
                .collect()
        };
        Ok((lig(&k), lig(&m), lig(&z)))
    }

    /// Résidu STATIQUE à l'état courant : −(forces + gravité), sur les 6n
    /// degrés physiques. C'est le `∂R/∂p` du calcul de sensibilité : on
    /// reconstruit le modèle au paramètre perturbé, on le POSE à l'état
    /// d'équilibre trouvé, et on lit de combien le résidu bouge — ce
    /// qu'ABAQUS appelle « holding the incremental displacement constant ».
    fn residu_statique(&self) -> PyResult<Vec<f64>> {
        let f = self.mo.forces();
        fini(f.as_slice(), "résidu statique").map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        Ok((0..f.len()).map(|i| -f[i]).collect())
    }

    /// Produits (∂R/∂p_i)ᵀ poids pour toutes les poutres, dans leur ordre.
    /// `poids` porte les 6n degrés physiques ; `quoi` est une direction de
    /// rigidité de d_residu_poutre. Chaque produit utilise les douze valeurs
    /// locales, sans matrice globale de dérivées ni vecteur global par poutre.
    /// Les directions ga et ei varient leurs deux composantes ; ei_e2 et
    /// ei_e3 les séparent. La condensation de l'intégrée est différenciée.
    #[pyo3(signature = (poids, quoi = "ei"))]
    fn d_residu_poutres_transpose(&self, poids: Vec<f64>, quoi: &str) -> PyResult<Vec<f64>> {
        if poids.len() != self.mo.n() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "poids : 6n valeurs requises",
            ));
        }
        fini(&poids, "poids du produit transposé")
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        let (cn, cm) =
            analyse::direction_poutre(quoi).map_err(pyo3::exceptions::PyValueError::new_err)?;
        let mut out = Vec::with_capacity(self.mo.poutres.len());
        for b in &self.mo.poutres {
            let local = b.forces_coefficients(&self.mo.corps, b.coefficients_derives(cn, cm), cm);
            let mut value = 0.0;
            for (j, node) in [b.a, b.b].iter().enumerate() {
                for k in 0..6 {
                    value -= poids[6 * node + k] * local[6 * j + k];
                }
            }
            out.push(value);
        }
        fini(&out, "produit transposé des sensibilités poutres")
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        Ok(out)
    }

    /// ∂(résidu statique)/∂(raideur de poutre), sans différence finie.
    /// La règle de chaîne dérive les coefficients de flexibilité condensée,
    /// puis les duaux calculent le gradient de cette énergie dérivée.
    /// `quoi` : "ea", "ga", "gj", "ei", "ei_e2" ou "ei_e3".
    fn d_residu_poutre(&self, i: usize, quoi: &str) -> PyResult<Vec<f64>> {
        if i >= self.mo.poutres.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("poutre inconnue"));
        }
        let (cn, cm) =
            analyse::direction_poutre(quoi).map_err(pyo3::exceptions::PyValueError::new_err)?;
        // Évaluer uniquement la poutre : les contacts et superéléments ne
        // dépendent pas de son paramètre, même lorsqu'ils sont déformés.
        let b = &self.mo.poutres[i];
        let local = b.forces_coefficients(&self.mo.corps, b.coefficients_derives(cn, cm), cm);
        fini(&local, "sensibilité poutre").map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let mut f = vec![0.0; self.mo.n()];
        for (j, node) in [b.a, b.b].iter().enumerate() {
            for k in 0..6 {
                f[6 * node + k] = -local[6 * j + k];
            }
        }
        Ok(f)
    }

    /// ∂K/∂(raideur de poutre), sans différence finie : coefficients
    /// condensés dérivés analytiquement, puis dérivée du champ de forces
    /// par duaux emboîtés. La masse ne dépend pas des raideurs.
    fn d_raideur_poutre(&self, i: usize, quoi: &str) -> PyResult<Vec<Vec<f64>>> {
        if i >= self.mo.poutres.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("poutre inconnue"));
        }
        let (cn, cm) =
            analyse::direction_poutre(quoi).map_err(pyo3::exceptions::PyValueError::new_err)?;
        let b = &self.mo.poutres[i];
        let local = b.tangente_coefficients(&self.mo.corps, b.coefficients_derives(cn, cm), cm);
        fini(local.as_slice(), "sensibilité de raideur poutre")
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let mut k = vec![vec![0.0; self.mo.n()]; self.mo.n()];
        for (i, ni) in [b.a, b.b].iter().enumerate() {
            for (j, nj) in [b.a, b.b].iter().enumerate() {
                for r in 0..6 {
                    for c in 0..6 {
                        k[6 * ni + r][6 * nj + c] = -local[(6 * i + r, 6 * j + c)];
                    }
                }
            }
        }
        Ok(k)
    }

    /// (K, C, M, Z, G) autour de l'état courant : raideur tangente,
    /// amortissement tangent (−∂f/∂u plus le terme gyroscopique), masse
    /// spatiale, base admissible et matrice de contrainte. De quoi assembler
    /// l'ÉQUATION VARIATIONNELLE du système et l'intégrer — ce qui donne la
    /// matrice de transition sans refaire 2m simulations non linéaires.
    #[pyo3(signature = (t = None))]
    fn k_c_m_z(
        &self,
        t: Option<f64>,
    ) -> PyResult<(
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
    )> {
        let t = t.unwrap_or(self.mo.t);
        let lig = |a: &DMatrix<f64>| -> Vec<Vec<f64>> {
            (0..a.nrows())
                .map(|i| (0..a.ncols()).map(|j| a[(i, j)]).collect())
                .collect()
        };
        let z = self
            .mo
            .base_admissible(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let g = self
            .mo
            .phi_g(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?
            .1;
        Ok((
            lig(&self
                .mo
                .raideur(t)
                .map_err(pyo3::exceptions::PyRuntimeError::new_err)?),
            lig(&self
                .mo
                .amortissement(t)
                .map_err(pyo3::exceptions::PyRuntimeError::new_err)?),
            lig(&self.mo.masse_spatiale()),
            lig(&z),
            lig(&g),
        ))
    }

    /// DIAGNOSTIC de l'affirmation « le transport sort de l'espace admissible » :
    /// rend (fuite relative de Ω·Z hors de l'image de Z, dim du noyau, ‖ΩZ‖).
    /// Une fuite ~0 réfuterait l'explication ; une fuite d'ordre 1 l'étaye.
    #[pyo3(signature = (t = None))]
    fn fuite_transport(&self, t: Option<f64>) -> PyResult<(f64, usize, f64)> {
        let t = t.unwrap_or(self.mo.t);
        let z = self
            .mo
            .base_admissible(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let n = self.mo.n();
        let mut omg = nalgebra::DMatrix::zeros(n, n);
        for (i, c) in self.mo.corps.iter().enumerate() {
            let sk = skew(&c.w);
            for a in 0..3 {
                for b in 0..3 {
                    omg[(6 * i + 3 + a, 6 * i + 3 + b)] = sk[(a, b)];
                }
            }
        }
        let oz = &omg * &z;
        let proj = &z * (z.transpose() * &oz); // Z orthonormale
        let hors = &oz - proj;
        let na = oz.norm();
        Ok((if na > 0.0 { hors.norm() / na } else { 0.0 }, z.ncols(), na))
    }

    /// DIAGNOSTIC : normes des blocs réduits (Kr, Cr, Mr, T) — pour comparer
    /// la raideur calculée à la raideur PHYSIQUE attendue.
    #[pyo3(signature = (t = None))]
    fn diag_modal(&self, t: Option<f64>) -> PyResult<(f64, f64, f64, f64)> {
        let t = t.unwrap_or(self.mo.t);
        let z = self
            .mo
            .base_admissible(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let k = self
            .mo
            .raideur(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let kr = z.transpose() * k * &z;
        let cr = z.transpose()
            * self
                .mo
                .amortissement(t)
                .map_err(pyo3::exceptions::PyRuntimeError::new_err)?
            * &z;
        let mr = z.transpose() * self.mo.masse_spatiale() * &z;
        let n = self.mo.n();
        let mut omg = nalgebra::DMatrix::zeros(n, n);
        for (i, c) in self.mo.corps.iter().enumerate() {
            let sk = skew(&c.w);
            for a in 0..3 {
                for b in 0..3 {
                    omg[(6 * i + 3 + a, 6 * i + 3 + b)] = sk[(a, b)];
                }
            }
        }
        let tr = z.transpose() * omg * &z;
        Ok((kr.norm(), cr.norm(), mr.norm(), tr.norm()))
    }

    /// Modes oscillants COMPLEXES : [(fréquence Hz, amortissement réduit ζ, σ)].
    /// σ > 0 indique une instabilité oscillante. Les racines purement réelles
    /// sont exclues ; cette liste ne suffit pas à conclure à la stabilité.
    #[pyo3(signature = (combien = 6, t = None))]
    fn modes_complexes(
        &self,
        py: Python<'_>,
        combien: usize,
        t: Option<f64>,
    ) -> PyResult<Vec<(f64, f64, f64)>> {
        let mo = &self.mo;
        py.detach(|| mo.modes_complexes(t.unwrap_or(mo.t), combien))
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)
    }

    /// Toutes les racines mécaniques locales : [(réel, imaginaire)] en s⁻¹.
    /// Multiplicités et conjugaisons conservées ; tri réel décroissant, puis
    /// imaginaire croissant. Date courante par défaut, sans modifier l'état.
    /// Contacts, appariement et aérodynamique refusés avec RuntimeError.
    /// Pour une trajectoire périodique, utiliser Floquet.
    #[pyo3(signature = (t = None))]
    fn spectre(&self, py: Python<'_>, t: Option<f64>) -> PyResult<Vec<(f64, f64)>> {
        py.detach(|| self.mo.spectre(t.unwrap_or(self.mo.t)))
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)
    }

    /// Bilan de la linéarisation mécanique locale, sans garantie de stabilité globale.
    /// Dictionnaire : t, racines [(réel, imaginaire)] en s⁻¹,
    /// abscisse_spectrale (max des parties réelles, None sans spectre), tol
    /// absolue en s⁻¹, statut et limites (liste de textes).
    /// Statuts : croissance_detectee (>tol), decroissance_spectrale (toutes
    /// <-tol), marginal_ou_indetermine, sans_ddl, hors_domaine.
    /// Hors domaine, racines=[] et abscisse_spectrale=None ; une erreur
    /// numérique reste un RuntimeError. tol doit être finie et >=0.
    #[pyo3(signature = (t = None, tol = 1e-8))]
    fn bilan_stabilite<'py>(
        &self,
        py: Python<'py>,
        t: Option<f64>,
        tol: f64,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        if !tol.is_finite() || tol < 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "tol doit être finie et >= 0",
            ));
        }
        let t = t.unwrap_or(self.mo.t);
        fini(&[t], "temps d'analyse").map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        self.mo
            .verifie_etat_fini()
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let mut limites = self.mo.limites_domaine_spectre();
        let hors_domaine = !limites.is_empty();
        limites.push("Linéarisation locale à l'état et à la date demandés ; aucune garantie de stabilité globale. Pour une trajectoire périodique, utiliser Floquet.");
        let racines = if hors_domaine {
            vec![]
        } else {
            self.spectre(py, Some(t))?
        };
        let abscisse = racines.iter().map(|v| v.0).reduce(f64::max);
        let statut = if hors_domaine {
            "hors_domaine"
        } else {
            match abscisse {
                None => "sans_ddl",
                Some(a) if a > tol => "croissance_detectee",
                Some(a) if a < -tol => "decroissance_spectrale",
                _ => "marginal_ou_indetermine",
            }
        };
        let out = pyo3::types::PyDict::new(py);
        out.set_item("t", t)?;
        out.set_item("racines", racines)?;
        out.set_item("abscisse_spectrale", abscisse)?;
        out.set_item("tol", tol)?;
        out.set_item("statut", statut)?;
        out.set_item("limites", limites)?;
        Ok(out)
    }

    /// (angle relatif déroulé, vitesse relative, couple) de chaque couple à loi, par nom
    fn couples(&self) -> Vec<(String, f64, f64, f64)> {
        self.mo
            .couples
            .iter()
            .map(|c| {
                let (th, w, tau) = c.etat(&self.mo.corps, self.mo.t);
                (c.nom.clone(), th, w, tau)
            })
            .collect()
    }

    /// Inflow uniforme : axe (monde), aire du disque, retard τ, masse volumique.
    /// Rend son indice.
    #[pyo3(signature = (axe, aire, tau, rho = 1.225))]
    fn inflow(&mut self, axe: [f64; 3], aire: f64, tau: f64, rho: f64) -> PyResult<usize> {
        valide_axe(axe)?;
        valide(&[aire, tau, rho], "inflow")?;
        if v3(axe).norm() == 0.0 || aire <= 0.0 || tau <= 0.0 || rho <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "inflow : axe non nul et paramètres > 0 requis",
            ));
        }
        self.mo.inflows.push(Inflow {
            axe: v3(axe).normalize(),
            aire,
            tau,
            rho,
            v_i: 0.0,
            poussee: 0.0,
            pitt: false,
            impose: false,
            profil: vec![],
            carte: (vec![], vec![], vec![]),
            centre: V3::zeros(),
            e1: V3::x(),
            vtip: 0.0,
            v1s: 0.0,
            v1c: 0.0,
            m_roul: 0.0,
            m_tang: 0.0,
        });
        Ok(self.mo.inflows.len() - 1)
    }

    /// Vitesse de l'air ambiant (repère monde). Le rotor de banc étant encastré,
    /// « avancer à V » s'écrit « vent de face à V » — c'est ce qui donne μ.
    fn vent(&mut self, v: [f64; 3]) -> PyResult<()> {
        valide(&v, "vent")?;
        self.mo.vent = v3(v);
        Ok(())
    }

    /// Bascule un inflow en PITT–PETERS à trois états. `centre` = centre du
    /// disque (monde), `e1` = référence d'azimut dans le plan, `vtip` = ΩR.
    /// Ces trois cotes ne se devinent pas depuis le multicorps : le solveur ne
    /// sait pas quel corps « est » le rotor, ni à quel régime il tourne. On les
    /// demande plutôt que de les inférer.
    #[pyo3(signature = (i, centre, e1, vtip, actif = true))]
    fn pitt_peters(
        &mut self,
        i: usize,
        centre: [f64; 3],
        e1: [f64; 3],
        vtip: f64,
        actif: bool,
    ) -> PyResult<()> {
        valide(&centre, "centre inflow")?;
        valide(&e1, "azimut inflow")?;
        valide(&[vtip], "vtip")?;
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "inflow {i}, mais le modèle n'en a que {}",
                self.mo.inflows.len()
            )));
        }
        if actif && !(vtip > 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "vtip = ΩR doit être > 0 : l'adimensionnement de Pitt–Peters en dépend",
            ));
        }
        let inf = &mut self.mo.inflows[i];
        let e = v3(e1);
        let e = e - inf.axe * e.dot(&inf.axe); // e1 doit être DANS le plan
        if e.norm() < 1e-9 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "e1 est colinéaire à l'axe : il ne définit aucun azimut",
            ));
        }
        inf.pitt = actif;
        inf.centre = v3(centre);
        inf.e1 = e.normalize();
        inf.vtip = vtip;
        Ok(())
    }

    /// UN PAS D'INFLOW SEUL, charges IMPOSÉES — le point d'entrée qui rend le
    /// modèle testable sans multicorps ni pale. C'est ce qui permet de le
    /// confronter à Froude, à sa constante de temps théorique et à Glauert :
    /// trois références analytiques qu'aucun couplage ne viendrait brouiller.
    fn pas_inflow(
        &mut self,
        i: usize,
        poussee: f64,
        m_roul: f64,
        m_tang: f64,
        h: f64,
    ) -> PyResult<()> {
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("inflow inconnu"));
        }
        let vent = self.mo.vent;
        let inf = &mut self.mo.inflows[i];
        inf.poussee = poussee;
        inf.m_roul = m_roul;
        inf.m_tang = m_tang;
        pitt_peters(inf, vent, h);
        Ok(())
    }

    /// (v_i, v1s, v1c) en m/s — les trois états d'inflow.
    fn etats_inflow(&self, i: usize) -> PyResult<(f64, f64, f64)> {
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("inflow inconnu"));
        }
        let f = &self.mo.inflows[i];
        Ok((f.v_i, f.v1s, f.v1c))
    }

    /// Force la vitesse induite de l'inflow `i` (m/s). Sert à partir d'un état
    /// convergé plutôt que de zéro — un rotor lâché à v_i = 0 traverse un
    /// transitoire qui n'a pas de sens physique.
    #[pyo3(signature = (i, v_i, impose = false))]
    fn pose_inflow(&mut self, i: usize, v_i: f64, impose: bool) -> PyResult<()> {
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "inflow {i}, mais le modèle n'en a que {}",
                self.mo.inflows.len()
            )));
        }
        self.mo.inflows[i].v_i = v_i;
        self.mo.inflows[i].impose = impose;
        Ok(())
    }

    /// Impose un PROFIL RADIAL d'inflow v(r̄) sur le disque `i` (m/s, > 0 vers
    /// le bas, ajouté à v_i), figé jusqu'au prochain appel — la sortie d'un
    /// sillage libre en stationnaire. `rbar` croissant dans [0, 1].
    /// Impose une CARTE d'inflow w(r̄, ψ) sur le disque `i` (m/s, > 0 vers le
    /// bas, ajoutée à v_i) : `rbar` croissant dans [0, 1], `psi` croissant dans
    /// [0, 2π) mesuré depuis `e1` dans le sens direct autour de l'axe, `w` en
    /// ligne (r̄ majeur : w[ir·n_psi + ip]). Bilinéaire, périodique en ψ. C'est
    /// l'entrée d'un Peters–He ou d'un sillage en avancement — toutes les
    /// harmoniques, là où `pose_inflow_harmoniques` n'en prend qu'une.
    #[pyo3(signature = (i, rbar, psi, w, centre, e1, vtip))]
    #[allow(clippy::too_many_arguments)]
    fn pose_inflow_carte(
        &mut self,
        i: usize,
        rbar: Vec<f64>,
        psi: Vec<f64>,
        w: Vec<f64>,
        centre: [f64; 3],
        e1: [f64; 3],
        vtip: f64,
    ) -> PyResult<()> {
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("inflow inconnu"));
        }
        if rbar.len() < 2
            || psi.len() < 2
            || w.len() != rbar.len() * psi.len()
            || rbar.windows(2).any(|a| a[1] <= a[0])
            || psi.windows(2).any(|a| a[1] <= a[0])
            || psi[0] < 0.0
            || *psi.last().unwrap() >= 2.0 * std::f64::consts::PI
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "carte : rbar et psi croissants (psi dans [0, 2π)), w de taille len(rbar)·len(psi)",
            ));
        }
        let inf = &mut self.mo.inflows[i];
        let e = v3(e1);
        let e = e - inf.axe * e.dot(&inf.axe);
        if e.norm() <= 0.0 || !(vtip > 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "e1 doit avoir une composante dans le plan du disque et vtip > 0",
            ));
        }
        inf.centre = v3(centre);
        inf.e1 = e.normalize();
        inf.vtip = vtip;
        inf.pitt = true;
        inf.impose = true;
        inf.carte = (rbar, psi, w);
        Ok(())
    }
    #[pyo3(signature = (i, rbar, v, centre, e1, vtip))]
    #[allow(clippy::too_many_arguments)]
    fn pose_inflow_profil(
        &mut self,
        i: usize,
        rbar: Vec<f64>,
        v: Vec<f64>,
        centre: [f64; 3],
        e1: [f64; 3],
        vtip: f64,
    ) -> PyResult<()> {
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("inflow inconnu"));
        }
        if rbar.len() != v.len() || rbar.len() < 2 || rbar.windows(2).any(|w| w[1] <= w[0]) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "profil : rbar croissant, même longueur que v, au moins deux points",
            ));
        }
        let inf = &mut self.mo.inflows[i];
        let e = v3(e1);
        let e = e - inf.axe * e.dot(&inf.axe);
        if e.norm() <= 0.0 || !(vtip > 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "e1 doit avoir une composante dans le plan du disque et vtip > 0",
            ));
        }
        inf.centre = v3(centre);
        inf.e1 = e.normalize();
        inf.vtip = vtip;
        inf.pitt = true;
        inf.impose = true;
        inf.profil = rbar.into_iter().zip(v).collect();
        Ok(())
    }

    /// Impose l'inflow d'un disque avec ses harmoniques 1/rev — ce qu'un
    /// sillage libre calculé DEHORS (`vinkulum.sillage`) rend au noyau :
    /// λ(r̄, ψ) = v_i + r̄·(v1c cos ψ + v1s sin ψ), figé jusqu'au prochain appel.
    /// `centre`, `e1`, `vtip` situent les stations comme pour Pitt–Peters.
    #[pyo3(signature = (i, v_i, v1c, v1s, centre, e1, vtip))]
    #[allow(clippy::too_many_arguments)]
    fn pose_inflow_harmoniques(
        &mut self,
        i: usize,
        v_i: f64,
        v1c: f64,
        v1s: f64,
        centre: [f64; 3],
        e1: [f64; 3],
        vtip: f64,
    ) -> PyResult<()> {
        if i >= self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("inflow inconnu"));
        }
        let inf = &mut self.mo.inflows[i];
        let e = v3(e1);
        let e = e - inf.axe * e.dot(&inf.axe);
        if e.norm() <= 0.0 || !(vtip > 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "e1 doit avoir une composante dans le plan du disque et vtip > 0",
            ));
        }
        inf.centre = v3(centre);
        inf.e1 = e.normalize();
        inf.vtip = vtip;
        inf.pitt = true;
        inf.impose = true;
        inf.v_i = v_i;
        inf.v1c = v1c;
        inf.v1s = v1s;
        Ok(())
    }

    /// Pale (théorie des tranches, polaire c81 mono-Mach). `p0` pied d'envergure
    /// et axes `es` (envergure), `ec` (avance), `en` (normale) en repère CORPS ;
    /// `polaire` = (alpha_deg, cl, cd, cm).
    #[pyo3(signature = (nom, corps, p0, es, ec, longueur, corde, polaire, inflow = None, rho = 1.225, gauss = 5, oye = false, lb = false, lb_const = None, masse_ajoutee = false, alpha_34 = false, retard_lineaire = false, a0_rad = None))]
    #[allow(clippy::too_many_arguments)]
    fn pale(
        &mut self,
        nom: String,
        corps: usize,
        p0: [f64; 3],
        es: [f64; 3],
        ec: [f64; 3],
        longueur: f64,
        corde: f64,
        polaire: (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>),
        inflow: Option<usize>,
        rho: f64,
        gauss: usize,
        oye: bool,
        lb: bool,
        lb_const: Option<(f64, f64, f64, f64)>,
        masse_ajoutee: bool,
        alpha_34: bool,
        retard_lineaire: bool,
        a0_rad: Option<f64>,
    ) -> PyResult<usize> {
        self._ref1(&nom, corps)?;
        if oye && lb {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "oye et lb : deux modèles de décrochage dynamique sur la même pale — en choisir un",
            ));
        }
        if let Some(i) = inflow {
            if i >= self.mo.inflows.len() {
                return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                    "« {nom} » : inflow {i}, mais le modèle n'en a que {}",
                    self.mo.inflows.len()
                )));
            }
        }
        let es_ = v3(es).normalize();
        let mut ec_ = v3(ec);
        ec_ -= es_ * ec_.dot(&es_);
        let ec_ = ec_.normalize();
        let en_ = es_.cross(&ec_);
        let (a, cl, cd, cm) = polaire;
        if !(a.len() == cl.len() && a.len() == cd.len() && a.len() == cm.len() && a.len() >= 2) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "polaire : tableaux de longueurs différentes",
            ));
        }
        if a.windows(2).any(|w| w[1] <= w[0]) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "polaire : incidences non croissantes",
            ));
        }
        // Øye : α₀ au premier passage par zéro de CL autour de 0°, a₀ par
        // moindres carrés sur |α − α₀| ≤ 5°, puis les deux tables dérivées
        let (mut alpha0, mut a0) = (0.0, 0.0);
        let (mut cl_fs, mut fst) = (vec![], vec![]);
        let mut lb_cn1 = 0.0;
        if oye || lb || retard_lineaire {
            let i0 = (0..a.len() - 1)
                .find(|&i| a[i] >= -15.0 && a[i + 1] <= 15.0 && cl[i] <= 0.0 && cl[i + 1] > 0.0)
                .ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(
                        "oye : pas de passage par CL = 0 entre −15° et 15°",
                    )
                })?;
            alpha0 = a[i0] - cl[i0] * (a[i0 + 1] - a[i0]) / (cl[i0 + 1] - cl[i0]);
            let (mut sxx, mut sxy) = (0.0, 0.0);
            for (x, y) in a.iter().zip(&cl) {
                if (x - alpha0).abs() <= 5.0 {
                    sxx += (x - alpha0) * (x - alpha0);
                    sxy += (x - alpha0) * y;
                }
            }
            a0 = sxy / sxx;
            if let Some(v) = a0_rad {
                a0 = v.to_radians(); // pente imposée (/rad → /°), déclarée par l'appelant
            }
            for (x, y) in a.iter().zip(&cl) {
                let att = a0 * (x - alpha0);
                let ratio = if att.abs() > 1e-9 {
                    (y / att).max(0.0)
                } else {
                    1.0
                };
                let f = (2.0 * ratio.sqrt() - 1.0).powi(2).clamp(0.0, 1.0);
                fst.push(f);
                cl_fs.push(if f < 1.0 - 1e-6 {
                    (y - f * att) / (1.0 - f)
                } else {
                    0.5 * y
                });
            }
            if lb {
                // C_N1 : le C_N attaché à la première incidence > α₀ où la polaire
                // statique dit f_st < 0,7 — le critère de Leishman–Beddoes (1989),
                // lu sur la table, pas posé. ESSAYÉ ET REFUSÉ le 6 sept. : le C_N
                // du CL max statique (usage pour les profils à décrochage de bord
                // de fuite) — sur le S809 mesuré (OSU 1995) il gagne 6 points à
                // 20° ± 5,5° (1,05 → 1,16 pour 1,63) et en PERD 17 à 8° ± 5,5°
                // (1,11 → 0,98 pour 1,22) : le seuil bas est le bon aux deux
                // premiers régimes, et aucun des deux ne donne le tourbillon
                // de 20° — c'est le modèle, pas le seuil.
                let i1 = a
                    .iter()
                    .zip(&fst)
                    .position(|(x, f)| *x > alpha0 && *f < 0.7)
                    .ok_or_else(|| {
                        pyo3::exceptions::PyValueError::new_err(
                            "lb : la polaire ne décroche pas (f_st ≥ 0,7 partout au-dessus de α₀)",
                        )
                    })?;
                lb_cn1 = a0 * (a[i1] - alpha0);
            }
        }
        let n_st = if gauss <= 3 { 3 } else { 5 };
        self.mo.pales.push(Pale {
            nom,
            corps,
            p0: v3(p0),
            es: es_,
            ec: ec_,
            en: en_,
            longueur,
            corde,
            polaire: std::sync::Arc::new(Polaire {
                alpha: a,
                cl,
                cd,
                cm,
            }),
            inflow,
            rho,
            gauss,
            sortie: (0.0, 0.0),
            instat: false,
            z: vec![(0.0, 0.0); n_st],
            oye,
            a0,
            alpha0,
            cl_fs,
            fst,
            f: vec![1.0; n_st],
            lb,
            lb_cn1,
            lb_const: lb_const.unwrap_or((LB_TP, LB_TF, LB_TV, LB_TVL)),
            nc: masse_ajoutee,
            nc_w: vec![0.0; n_st],
            nc_h: 0.0,
            a34: alpha_34,
            retard_lineaire,
            etats_lb: vec![
                EtatLb {
                    ff: 1.0,
                    ..Default::default()
                };
                n_st
            ],
        });
        Ok(self.mo.pales.len() - 1)
    }

    /// Bascule une pale en AÉRODYNAMIQUE INSTATIONNAIRE (Wagner/Jones).
    /// Deux états par station ; le modèle quasi-stationnaire reste le défaut.
    #[pyo3(signature = (i, actif = true))]
    fn instationnaire(&mut self, i: usize, actif: bool) -> PyResult<()> {
        if i >= self.mo.pales.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "pale {i}, mais le modèle n'en a que {}",
                self.mo.pales.len()
            )));
        }
        self.mo.pales[i].instat = actif;
        for z in self.mo.pales[i].z.iter_mut() {
            *z = (0.0, 0.0);
        }
        Ok(())
    }

    /// États de Leishman–Beddoes par station de la pale `i` :
    /// (f'' retardé, τ_v, C_N^v tourbillonnaire, C_N retardé C_N'), et le seuil C_N1.
    fn etats_lb(&self, i: usize) -> PyResult<(Vec<(f64, f64, f64, f64)>, f64)> {
        if i >= self.mo.pales.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("pale inconnue"));
        }
        let p = &self.mo.pales[i];
        Ok((
            p.etats_lb
                .iter()
                .map(|e| (e.ff, e.tau_v, e.cnv, e.cn_prev - e.dp))
                .collect(),
            p.lb_cn1,
        ))
    }

    /// États de Wagner (z₁, z₂) par station, pour la pale `i`.
    fn etats_wagner(&self, i: usize) -> PyResult<Vec<(f64, f64)>> {
        if i >= self.mo.pales.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("pale inconnue"));
        }
        Ok(self.mo.pales[i].z.clone())
    }

    /// (nom, poussée le long de l'axe d'inflow, couple autour de cet axe) par pale, et (v_i, poussée) par inflow
    fn aero(&self) -> (Vec<(String, f64, f64)>, Vec<(f64, f64)>) {
        (
            self.mo
                .pales
                .iter()
                .map(|p| (p.nom.clone(), p.sortie.0, p.sortie.1))
                .collect(),
            self.mo.inflows.iter().map(|i| (i.v_i, i.poussee)).collect(),
        )
    }

    /// Poutre géométriquement exacte entre les corps a et b (nœuds, à leur pose
    /// AU REPOS : l'axe de la poutre est la direction a → b, le repère matériel
    /// est celui de a) ; raideurs (EA, GA_y, GA_z) et (GJ, EI_y, EI_z).
    ///
    /// `ei3` et `ga3` valent `ei`/`ga` par défaut (section isotrope). Une
    /// section MINCE les sépare franchement — un ruban de 1 × 0,1 mm² a un
    /// facteur 100 sur EI et 4,7 sur GA — et c'est la séparation qui rend le
    /// flambement latéral possible.
    /// `formulation="milieu"` conserve l'élément historique. L'option
    /// `"integree"` intègre la rotation le long de l'élément et condense la
    /// flexibilité de flexion interne : meilleure précision sur les poutres
    /// régulières, mais encore expérimentale pour les câbles très souples.
    #[pyo3(signature = (nom, a, b, ea, ga, gj, ei, ei3 = None, ga3 = None, formulation = "milieu"))]
    #[allow(clippy::too_many_arguments)]
    fn poutre(
        &mut self,
        nom: String,
        a: usize,
        b: usize,
        ea: f64,
        ga: f64,
        gj: f64,
        ei: f64,
        ei3: Option<f64>,
        ga3: Option<f64>,
        formulation: &str,
    ) -> PyResult<usize> {
        let integree = match formulation {
            "milieu" => false,
            "integree" => true,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "formulation poutre : milieu ou integree",
                ))
            }
        };
        self._refs(&nom, Some(a), Some(b))?;
        let raideurs = [ea, ga, gj, ei, ei3.unwrap_or(ei), ga3.unwrap_or(ga)];
        valide(&raideurs, "raideur poutre")?;
        if raideurs.iter().any(|v| *v <= 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "poutre : raideurs > 0 requises",
            ));
        }
        let (ca, cb) = (&self.mo.corps[a], &self.mo.corps[b]);
        let (d, bas) = cb.difference_precise(ca);
        let d = d + bas;
        let l = d.norm();
        if !(l > 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "« {nom} » : les deux corps sont au même point — une poutre de longueur \
                 nulle n'a ni repère ni raideur"
            )));
        }
        // repère matériel au repos : e₁ = a → b, e₂ = le y du nœud a projeté, e₃ = e₁ × e₂
        let e1 = d / l;
        let ya = ca.rot.column(1).into_owned();
        let mut e2 = ya - e1 * ya.dot(&e1);
        if e2.norm() < 1e-9 {
            let za = ca.rot.column(2).into_owned();
            e2 = za - e1 * za.dot(&e1);
        }
        let e2 = e2.normalize();
        let e3 = e1.cross(&e2);
        let rmat = M3::from_columns(&[e1, e2, e3]);
        // DEUX inerties de flexion, et une pale en a BESOIN : sa flexion de
        // battement et sa flexion de traînée diffèrent d'un facteur 500 sur
        // une plaque cambrée. `ei` porte la flexion autour de e₂, `ei3` celle
        // autour de e₃ (défaut : la même, section isotrope).
        self.mo.poutres.push(Poutre {
            nom,
            integree,
            a,
            b,
            l,
            cn: [ea, ga, ga3.unwrap_or(ga)],
            cm: [gj, ei, ei3.unwrap_or(ei)],
            ra: ca.rot.transpose() * rmat,
            rb: cb.rot.transpose() * rmat,
        });
        Ok(self.mo.poutres.len() - 1)
    }

    /// (γ, κ) de chaque poutre — déformation et courbure matérielles courantes
    fn poutres(&self) -> Vec<(String, [f64; 3], [f64; 3])> {
        self.mo
            .poutres
            .iter()
            .map(|b| {
                let (g, k) = b.deformations(&self.mo.corps);
                (b.nom.clone(), g, k)
            })
            .collect()
    }

    /// Vue f64 : (t, [r], [R en ligne], [v], [ω], [v_i des inflows]).
    /// Les positions sont arrondies ; `etat_precis()` conserve aussi leur
    /// partie basse pour une restauration cinématique sans cette perte.
    #[allow(clippy::type_complexity)]
    fn etat(
        &self,
    ) -> (
        f64,
        Vec<[f64; 3]>,
        Vec<[f64; 9]>,
        Vec<[f64; 3]>,
        Vec<[f64; 3]>,
        Vec<f64>,
    ) {
        let r = self
            .mo
            .corps
            .iter()
            .map(|c| [c.r.x, c.r.y, c.r.z])
            .collect();
        let rot = self
            .mo
            .corps
            .iter()
            .map(|c| {
                let m = c.rot;
                [
                    m[(0, 0)],
                    m[(0, 1)],
                    m[(0, 2)],
                    m[(1, 0)],
                    m[(1, 1)],
                    m[(1, 2)],
                    m[(2, 0)],
                    m[(2, 1)],
                    m[(2, 2)],
                ]
            })
            .collect();
        let v = self
            .mo
            .corps
            .iter()
            .map(|c| [c.v.x, c.v.y, c.v.z])
            .collect();
        let w = self
            .mo
            .corps
            .iter()
            .map(|c| [c.w.x, c.w.y, c.w.z])
            .collect();
        (
            self.mo.t,
            r,
            rot,
            v,
            w,
            self.mo.inflows.iter().map(|i| i.v_i).collect(),
        )
    }

    /// État cinématique précis : les six champs de `etat()`, suivis de [r_bas].
    /// La position est r+r_bas, à conserver séparément. `pose_etat(*etat)`
    /// accepte ce tuple et préserve aussi les rotations telles quelles.
    /// Multiplicateurs, historiques de contact et états aérodynamiques autres
    /// que v_i ne font pas partie de ce tuple.
    fn etat_precis(
        &self,
    ) -> (
        f64,
        Vec<[f64; 3]>,
        Vec<[f64; 9]>,
        Vec<[f64; 3]>,
        Vec<[f64; 3]>,
        Vec<f64>,
        Vec<[f64; 3]>,
    ) {
        let (t, r, rot, v, w, vi) = self.etat();
        let bas = self.mo.corps.iter().map(|c| c.r_bas.into()).collect();
        (t, r, rot, v, w, vi, bas)
    }

    /// Restaure un état rendu par `etat`. Les rotations sont RÉ-ORTHONORMALISÉES
    /// (Gram-Schmidt) : une R perturbée composante par composante n'est plus dans
    /// SO(3), et l'y laisser fabrique une énergie qui n'existe pas.
    /// Avec r_bas (tuple `etat_precis()`), les positions sont normalisées
    /// en deux parties et les rotations validées sans être modifiées.
    /// Sans r_bas, la partie basse est remise à zéro : cette vue est arrondie.
    #[pyo3(signature = (t, r, rot, v, w, v_i = vec![], r_bas = None))]
    #[allow(clippy::too_many_arguments)]
    fn pose_etat(
        &mut self,
        t: f64,
        r: Vec<[f64; 3]>,
        rot: Vec<[f64; 9]>,
        v: Vec<[f64; 3]>,
        w: Vec<[f64; 3]>,
        v_i: Vec<f64>,
        r_bas: Option<Vec<[f64; 3]>>,
    ) -> PyResult<()> {
        let nb = self.mo.corps.len();
        if r.len() != nb || rot.len() != nb || v.len() != nb || w.len() != nb {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "pose_etat : {nb} corps, mais r/rot/v/w en donnent {}/{}/{}/{}",
                r.len(),
                rot.len(),
                v.len(),
                w.len()
            )));
        }
        if !t.is_finite()
            || r.iter()
                .flatten()
                .chain(v.iter().flatten())
                .chain(w.iter().flatten())
                .chain(rot.iter().flatten())
                .chain(&v_i)
                .any(|x| !x.is_finite())
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "pose_etat : temps, positions, orientations, vitesses et inflow doivent être finis",
            ));
        }
        if !v_i.is_empty() && v_i.len() != self.mo.inflows.len() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "pose_etat : v_i doit être vide ou contenir une valeur par inflow",
            ));
        }
        let precis = r_bas.is_some();
        let bas = r_bas.unwrap_or_else(|| vec![[0.0; 3]; nb]);
        if bas.len() != nb || bas.iter().flatten().any(|x| !x.is_finite()) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "pose_etat : r_bas doit contenir trois valeurs finies par corps",
            ));
        }
        let positions: Vec<_> = r
            .iter()
            .zip(&bas)
            .map(|(r, b)| numerique::deplace_compense(v3(*r), V3::zeros(), v3(*b)))
            .collect();
        if positions
            .iter()
            .any(|(r, b)| r.iter().chain(b.iter()).any(|x| !x.is_finite()))
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "pose_etat : position compensée non finie",
            ));
        }
        // Préparer TOUTES les rotations avant de modifier le modèle : une
        // erreur au dernier corps ne doit pas restaurer les premiers seuls.
        let rotations: Vec<M3> = rot
            .iter()
            .map(|r| {
                let m = m3(*r);
                if precis {
                    if !(m.determinant() > 0.0)
                        || !m.determinant().is_finite()
                        || (m.transpose() * m - M3::identity()).amax() > 1e-12
                    {
                        return Err(pyo3::exceptions::PyValueError::new_err(
                            "pose_etat précis : rotation non orthonormale ou réfléchie",
                        ));
                    }
                    return Ok(m);
                }
                let e1 = m.column(0).normalize();
                let e2 = (m.column(1) - e1 * e1.dot(&m.column(1))).normalize();
                let e3 = e1.cross(&e2);
                let rm = M3::from_columns(&[e1, e2, e3]);
                if !(m.determinant() > 0.0)
                    || !m.determinant().is_finite()
                    || rm.iter().any(|x| !x.is_finite())
                    || (rm.transpose() * rm - M3::identity()).amax() > 1e-9
                {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "pose_etat : orientation dégénérée ou réfléchie",
                    ));
                }
                Ok(rm)
            })
            .collect::<PyResult<_>>()?;
        for (i, c) in self.mo.corps.iter_mut().enumerate() {
            (c.r, c.r_bas) = positions[i];
            c.rot = rotations[i];
            c.v = v3(v[i]);
            c.w = v3(w[i]);
        }
        for (k, val) in v_i.iter().enumerate() {
            if k < self.mo.inflows.len() {
                self.mo.inflows[k].v_i = *val;
            }
        }
        self.mo.t = t;
        Ok(())
    }

    /// Arme la moyenne du torseur aérodynamique à partir de `t0` (remise à zéro).
    /// Ce qu'un rotor trime est une moyenne sur un tour : la lire ici évite de
    /// rapatrier la trajectoire entière à chaque évaluation du trim.
    fn moyenne_aero(&mut self, t0: f64) {
        self.mo.moy_aero = Some((t0, V3::zeros(), V3::zeros(), 0.0));
    }

    /// (F, M) aérodynamiques MOYENS depuis `t0`, réduits à l'origine du monde,
    /// et la durée intégrée.
    fn torseur_moyen(&self) -> PyResult<([f64; 3], [f64; 3], f64)> {
        match self.mo.moy_aero {
            Some((_, sf, sm, sh)) if sh > 0.0 => Ok(([sf.x / sh, sf.y / sh, sf.z / sh],
                                                    [sm.x / sh, sm.y / sh, sm.z / sh], sh)),
            _ => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "torseur_moyen : rien d'intégré — appeler moyenne_aero(t0) avant simule, avec t0 < t_end")),
        }
    }

    /// CONTACT par pénalité : une sphère portée par `corps` (centre `p0` en
    /// repère corps, rayon `r`) contre le demi-espace n·(p − origine) ≥ 0.
    ///
    /// `expo < 0` bascule sur la BARRIÈRE IPC au lieu de la loi de Hertz :
    /// B(d) = −(d − d̂)² ln(d/d̂), nulle au-delà de `d_hat`, et divergente en
    /// d → 0⁺ — l'interpénétration devient impossible au lieu d'être petite.
    ///
    /// Si `b` est donné, c'est un contact SPHÈRE/SPHÈRE entre deux corps
    /// mobiles (`pb`, `rayon_b`) et la force est appliquée aux deux, opposée.
    ///
    /// `k`, `expo` : raideur et exposant (1,5 = Hertz sphère/plan) ·
    /// `c` : amortissement de Hunt–Crossley, RELATIF (F = kδ^e(1 + c·δ̇)) —
    /// il part et revient à zéro avec l'enfoncement, ce qu'un amortisseur
    /// linéaire ne fait pas · `mu`, `v_eps` : Coulomb lissé par tanh.
    #[pyo3(signature = (nom, corps, p0, rayon, normale = [0.0, 0.0, 1.0], origine = [0.0; 3],
                        k = 1e6, expo = 1.5, c = 0.0, mu = 0.0, v_eps = 1e-3,
                        b = None, pb = [0.0; 3], rayon_b = 0.0, d_hat = 1e-3,
                        p1 = None, p1b = None, cable = false, demi = None, maille = None,
                        nonlisse = false, restitution = 0.0, cylindre = None))]
    #[allow(clippy::too_many_arguments)]
    fn contact(
        &mut self,
        nom: String,
        corps: usize,
        p0: [f64; 3],
        rayon: f64,
        normale: [f64; 3],
        origine: [f64; 3],
        k: f64,
        expo: f64,
        c: f64,
        mu: f64,
        v_eps: f64,
        b: Option<usize>,
        pb: [f64; 3],
        rayon_b: f64,
        d_hat: f64,
        p1: Option<[f64; 3]>,
        p1b: Option<[f64; 3]>,
        cable: bool,
        demi: Option<[f64; 3]>,
        maille: Option<usize>,
        nonlisse: bool,
        restitution: f64,
        cylindre: Option<([f64; 3], f64, f64)>,
    ) -> PyResult<usize> {
        if let Some((ax, r, l)) = cylindre {
            valide_axe(ax)?;
            valide(&[r, l], "dimensions cylindre")?;
            if !(r > 0.0) || !(l > 0.0) || v3(ax).norm() <= 0.0 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "cylindre : (axe non nul, rayon > 0, demi-longueur > 0)",
                ));
            }
            if b.is_none() {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "cylindre : il doit être porté par un second corps (b)",
                ));
            }
        }
        if let Some(mi) = maille {
            if mi >= self.mo.maillages.len() {
                return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                    "maillage {mi}, mais le modèle n'en a que {}",
                    self.mo.maillages.len()
                )));
            }
        }
        if let Some(d) = demi {
            if d.iter().any(|x| !(*x > 0.0) || !x.is_finite()) {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "boîte : les trois demi-dimensions doivent être > 0",
                ));
            }
            if b.is_none() {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "boîte : elle doit être portée par un second corps (b)",
                ));
            }
        }
        if corps >= self.mo.corps.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("corps inconnu"));
        }
        if let Some(ib) = b {
            if ib >= self.mo.corps.len() {
                return Err(pyo3::exceptions::PyIndexError::new_err(
                    "second corps inconnu",
                ));
            }
        }
        // SONDAGE DU 7 SEPT. : tout ce qui suit était ACCEPTÉ EN SILENCE — un
        // modèle faux qui simule, ou qui BLOQUE (normale nulle : `normalize()`
        // rend NaN, et la SVD de nalgebra ne converge jamais sur un NaN).
        let err = |m: String| pyo3::exceptions::PyValueError::new_err(m);
        let fini3 = |v: [f64; 3]| v.iter().all(|x| x.is_finite());
        if !fini3(p0)
            || !fini3(pb)
            || !fini3(origine)
            || !p1.is_none_or(fini3)
            || !p1b.is_none_or(fini3)
        {
            return Err(err(format!("contact « {nom} » : point non fini")));
        }
        if !(rayon >= 0.0) || !rayon.is_finite() || !(rayon_b >= 0.0) || !rayon_b.is_finite() {
            return Err(err(format!(
                "contact « {nom} » : rayon {rayon} / rayon_b {rayon_b} — un rayon est ≥ 0 et fini"
            )));
        }
        if b == Some(corps) {
            return Err(err(format!(
                "contact « {nom} » : les deux corps sont le même"
            )));
        }
        if b.is_none() && (!fini3(normale) || !(v3(normale).norm() > 0.0)) {
            return Err(err(format!(
                "contact « {nom} » : normale nulle ou non finie"
            )));
        }
        if b.is_none() {
            valide_axe(normale)?;
        }
        valide(&[v_eps, d_hat, restitution], "paramètres de contact")?;
        if !k.is_finite() || (!nonlisse && !(k > 0.0)) {
            return Err(err(format!("contact « {nom} » : k = {k} — la raideur doit être > 0 et finie (à k = 0 le contact n'agit pas)")));
        }
        if !expo.is_finite() || !(c >= 0.0) || !c.is_finite() || !(mu >= 0.0) || !mu.is_finite() {
            return Err(err(format!(
                "contact « {nom} » : expo, c ou mu non fini ou négatif"
            )));
        }
        if expo < 0.0 && !(d_hat > 0.0 && d_hat.is_finite()) {
            return Err(err(format!(
                "contact « {nom} » : barrière IPC à d_hat = {d_hat} (ln(d/d̂))"
            )));
        }
        if mu > 0.0 && !(v_eps > 0.0) {
            return Err(err(format!(
                "contact « {nom} » : v_eps doit être > 0 avec du frottement (tanh(v/0))"
            )));
        }
        if nonlisse && !(0.0..=1.0).contains(&restitution) {
            return Err(err(format!(
                "contact « {nom} » : restitution {restitution} hors de [0, 1]"
            )));
        }
        let n_prim = demi.is_some() as usize
            + cylindre.is_some() as usize
            + maille.is_some() as usize
            + (rayon_b > 0.0) as usize;
        if n_prim > 1 {
            return Err(err(format!("contact « {nom} » : une seule primitive sur le second corps (rayon_b, demi, cylindre ou maille)")));
        }
        // le maillage vit sur SON corps : b absent = ce corps, b autre = faux
        // (la bille tombait à travers le sol, sans un mot)
        let b = match maille {
            Some(mi) => {
                let cm = self.mo.maillages[mi].corps;
                match b {
                    None => Some(cm),
                    Some(ib) if ib != cm => {
                        return Err(err(format!("contact « {nom} » : le maillage {mi} est porté par le corps {cm}, pas par b = {ib}")));
                    }
                    _ => b,
                }
            }
            None => b,
        };
        if b == Some(corps) {
            return Err(err(format!(
                "contact « {nom} » : le maillage est porté par le corps en contact"
            )));
        }
        // p1 / p1b absents = sphère (segment dégénéré) : c'est le défaut, donc
        // tout modèle existant garde exactement son comportement.
        let m_avant = self.mo.m();
        // les contacts AUTOMATIQUES restent en queue : un contact déclaré après
        // une passe d'appariement s'insère devant eux, et les indices des
        // contacts manuels (`Elem::K::ct`) ne bougent pas
        let idx = self
            .mo
            .contacts
            .iter()
            .position(|c| c.auto)
            .unwrap_or(self.mo.contacts.len());
        if nonlisse {
            if expo < 0.0 || cable {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "nonlisse : incompatible avec la barrière IPC (expo < 0) et le câble",
                ));
            }
            self.mo.elems.push(Elem::K(Unilateral {
                nom: nom.clone(),
                ct: idx,
                a: corps,
                b,
            }));
            self.mo.actif = vec![];
        }
        self.mo.contacts.insert(
            idx,
            Contact {
                nom,
                corps,
                p0: v3(p0),
                p1: v3(p1.unwrap_or(p0)),
                rayon,
                b,
                pb: v3(pb),
                p1b: v3(p1b.unwrap_or(pb)),
                rayon_b,
                normale: v3(normale).normalize(),
                origine: v3(origine),
                k,
                expo,
                d_hat,
                c,
                mu,
                v_eps,
                demi: demi.map(v3),
                cylindre: cylindre.map(|(ax, r, l)| (v3(ax).normalize(), r, l)),
                maille,
                q_maille: V3::zeros(),
                s_axe: 0.0,
                cable,
                sortie: (0.0, 0.0),
                nonlisse,
                row: m_avant,
                fn_impose: None,
                restitution,
                auto: false,
            },
        );
        Ok(idx)
    }

    /// SUPERÉLÉMENT : la raideur d'un maillage EF quelconque, attachée à des
    /// nœuds qui sont des corps ordinaires. `k` est donnée à plat, ligne par
    /// ligne, en 6n × 6n — l'ordre des ddl est (tx,ty,tz,rx,ry,rz) par nœud.
    ///
    /// C'est le pont vers les coques et les solides 3D : on ne les maille pas
    /// ici, on importe leur K réduite (Craig–Bampton, que `vinkulum.reduction`
    /// fournit). Corotationnel : K vit dans le repère du PREMIER nœud.
    /// Efforts −BᵀKu, amortissement −βBᵀKBv, B = ∂u/∂q. `alpha` doit être
    /// nul ; son ancien paramètre était accepté sans effet. Les nœuds sont
    /// distincts, K finie symétrique, beta ≥ 0. Rotation relative à π :
    /// tangente non différentiable, refus explicite avant résolution.
    #[pyo3(signature = (nom, noeuds, k, alpha = 0.0, beta = 0.0))]
    fn superelement(
        &mut self,
        nom: String,
        noeuds: Vec<usize>,
        k: Vec<f64>,
        alpha: f64,
        beta: f64,
    ) -> PyResult<usize> {
        valide(&k, "raideur superélément")?;
        valide(&[alpha, beta], "amortissement superélément")?;
        if alpha != 0.0 || beta < 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "superélément : alpha non nul non pris en charge, beta ≥ 0 requis",
            ));
        }
        let uniques: std::collections::HashSet<_> = noeuds.iter().collect();
        if uniques.len() != noeuds.len() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "superélément : nœuds répétés",
            ));
        }
        let n = noeuds.len();
        if n < 2 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "superélément : au moins deux nœuds (avec un seul, K n'a aucun sens relatif)",
            ));
        }
        for &i in &noeuds {
            self._ref1(&nom, i)?;
        }
        let d = 6 * n;
        if k.len() != d * d {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "superélément « {nom} » : K devrait faire {d}×{d} = {} valeurs, en a {}",
                d * d,
                k.len()
            )));
        }
        let km = DMatrix::from_row_slice(d, d, &k);
        // K doit être SYMÉTRIQUE : une raideur non symétrique n'est pas une
        // énergie, et le noyau la propagerait sans le dire.
        let asym = (0..d)
            .map(|i| {
                (0..d)
                    .map(|j| (km[(i, j)] - km[(j, i)]).abs())
                    .fold(0.0f64, f64::max)
            })
            .fold(0.0f64, f64::max);
        let ech = km.amax().max(1.0);
        if asym > 1e-9 * ech {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "superélément « {nom} » : K n'est pas symétrique (écart {asym:.3e} pour une \
                 échelle {ech:.3e}) — une raideur non symétrique ne dérive d'aucune énergie"
            )));
        }
        let c0 = &self.mo.corps[noeuds[0]];
        let rt = c0.rot.transpose();
        let (ref_p, ref_r): (Vec<V3>, Vec<M3>) = noeuds
            .iter()
            .map(|&i| {
                let c = &self.mo.corps[i];
                let (d, bas) = c.difference_precise(c0);
                (rt * (d + bas), rt * c.rot)
            })
            .unzip();
        self.mo.supers.push(Superelement {
            nom,
            noeuds,
            k: 0.5 * (&km + km.transpose()),
            alpha,
            beta,
            ref_p,
            ref_r,
        });
        Ok(self.mo.supers.len() - 1)
    }

    /// MAILLAGE TRIANGULAIRE porté par un corps — la géométrie quelconque.
    /// `sommets` à plat (3 par point), `tris` à plat (3 index par triangle),
    /// dans le repère du corps. Le BVH se construit ICI, une fois.
    fn maillage(
        &mut self,
        nom: String,
        corps: usize,
        sommets: Vec<f64>,
        tris: Vec<usize>,
    ) -> PyResult<usize> {
        self._ref1(&nom, corps)?;
        if !sommets.len().is_multiple_of(3) || !tris.len().is_multiple_of(3) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "maillage : sommets et triangles se donnent à plat, par groupes de 3",
            ));
        }
        if sommets.iter().any(|x| !x.is_finite()) {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "maillage « {nom} » : sommet non fini"
            )));
        }
        let sv: Vec<V3> = sommets
            .chunks(3)
            .map(|c| V3::new(c[0], c[1], c[2]))
            .collect();
        let nt = sv.len();
        let mut tv = Vec::with_capacity(tris.len() / 3);
        for c in tris.chunks(3) {
            if c.iter().any(|&i| i >= nt) {
                return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                    "maillage « {nom} » : un triangle référence un sommet hors des {nt} donnés"
                )));
            }
            // un triangle dégénéré n'a pas de normale : `point_triangle` y divise
            // par zéro et le contact devient NaN, en silence
            let aire = (sv[c[1]] - sv[c[0]]).cross(&(sv[c[2]] - sv[c[0]])).norm();
            if c[0] == c[1] || c[1] == c[2] || c[0] == c[2] || !(aire > 0.0) {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "maillage « {nom} » : triangle dégénéré {c:?} (sommets répétés ou aire nulle)"
                )));
            }
            tv.push([c[0], c[1], c[2]]);
        }
        if tv.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "maillage : aucun triangle",
            ));
        }
        let (bvh, ordre) = Maillage::bati(&sv, &tv);
        self.mo.maillages.push(Maillage {
            nom,
            corps,
            sommets: sv,
            tris: tv,
            bvh,
            ordre,
        });
        Ok(self.mo.maillages.len() - 1)
    }

    /// Distance d'un point MONDE au maillage `i`, et le point le plus proche
    /// (monde). Sert aux contrôles et à toute requête géométrique.
    fn distance_maillage(&self, i: usize, p: [f64; 3]) -> PyResult<(f64, [f64; 3])> {
        if i >= self.mo.maillages.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("maillage inconnu"));
        }
        let m = &self.mo.maillages[i];
        let c = &self.mo.corps[m.corps];
        let pl = c.rot.transpose() * (v3(p) - c.r);
        let (d, q) = m.plus_proche(&pl);
        let qw = c.r + c.rot * q;
        Ok((d, [qw.x, qw.y, qw.z]))
    }

    /// Déclare une sphère candidate à l'APPARIEMENT automatique.
    fn sphere(&mut self, corps: usize, p0: [f64; 3], rayon: f64) -> PyResult<usize> {
        if corps >= self.mo.corps.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err("corps inconnu"));
        }
        if !(rayon > 0.0) || !rayon.is_finite() || p0.iter().any(|x| !x.is_finite()) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "sphère : centre fini et rayon fini strictement positif requis",
            ));
        }
        self.mo.spheres.push((corps, v3(p0), rayon));
        Ok(self.mo.spheres.len() - 1)
    }

    /// Active l'appariement automatique des sphères déclarées, avec une loi de
    /// contact commune. `marge` : distance sous laquelle une paire est gardée
    /// (elle doit couvrir ce qu'une sphère parcourt en un pas).
    #[pyo3(signature = (k = 1e6, expo = 1.5, c = 0.0, mu = 0.0, v_eps = 1e-3, marge = 0.0,
                        d_hat = 1e-3))]
    #[allow(clippy::too_many_arguments)]
    fn appariement(
        &mut self,
        k: f64,
        expo: f64,
        c: f64,
        mu: f64,
        v_eps: f64,
        marge: f64,
        d_hat: f64,
    ) -> PyResult<()> {
        valide(&[k, expo, c, mu, v_eps, marge, d_hat], "appariement")?;
        if k <= 0.0 || c < 0.0 || mu < 0.0 || v_eps <= 0.0 || marge < 0.0 || d_hat <= 0.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "appariement : k, v_eps, d_hat > 0 ; c, mu, marge ≥ 0 requis",
            ));
        }
        self.mo.loi_contact = Some((k, expo, c, mu, v_eps, marge, d_hat));
        for ct in self.mo.contacts.iter_mut().filter(|ct| ct.auto) {
            ct.k = k;
            ct.expo = expo;
            ct.c = c;
            ct.mu = mu;
            ct.v_eps = v_eps;
            ct.d_hat = d_hat;
        }
        self.mo.symb = None;
        Ok(())
    }

    /// (paires actives au dernier pas, secondes cumulées d'appariement)
    fn appariement_stats(&self) -> (usize, f64) {
        (self.mo.n_paires, self.mo.t_appar)
    }

    /// (nom, enfoncement, force normale) de chaque contact
    fn contacts(&self) -> Vec<(String, f64, f64)> {
        self.mo
            .contacts
            .iter()
            .map(|c| (c.nom.clone(), c.sortie.0, c.sortie.1))
            .collect()
    }

    /// Effort constant (F, M) au CdM du corps, repère monde.
    fn effort(&mut self, corps: usize, f: [f64; 3], m: [f64; 3]) -> PyResult<()> {
        self._ref1("effort", corps)?;
        valide(&f, "effort")?;
        valide(&m, "moment")?;
        self.mo.efforts.push((corps, v3(f), v3(m)));
        Ok(())
    }

    /// Intègre jusqu'à `t_end` au pas `h`. `adaptatif` = tolérance RELATIVE du
    /// résidu de demi-pas (None = pas FIXE, le défaut) ; `bornes` = (min, max)
    /// en multiples de `h`. Rend la trajectoire échantillonnée
    /// tous les `tous` pas : liste de (t, [r_i…], [R_i en ligne…], [v_i…], [w_i…], λ).
    /// `sigma_lie` dans [0,1] : correction géométrique expérimentale, 0 par
    /// défaut. Résout la cinématique implicite dans ||theta|| < pi ; voir
    /// docs/INTEGRATION_SIGMA.md pour les mesures et le domaine validé.
    #[pyo3(signature = (t_end, h, rho = 0.9, tol = 1e-12, newton_max = 25, tous = 1,
                       adaptatif = None, bornes = (1.0 / 64.0, 16.0),
                       pas_contact = None, ccd = false, moment = false, energie = false, ggl = false,
                       sigma_lie = 0.0))]
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    fn simule(
        &mut self,
        py: Python<'_>,
        t_end: f64,
        h: f64,
        rho: f64,
        tol: f64,
        newton_max: usize,
        tous: usize,
        adaptatif: Option<f64>,
        bornes: (f64, f64),
        pas_contact: Option<(f64, f64, f64)>,
        ccd: bool,
        moment: bool,
        energie: bool,
        ggl: bool,
        sigma_lie: f64,
    ) -> PyResult<
        Vec<(
            f64,
            Vec<[f64; 3]>,
            Vec<[f64; 9]>,
            Vec<[f64; 3]>,
            Vec<[f64; 3]>,
            Vec<f64>,
        )>,
    > {
        reglages(self.mo.t, t_end, h, tol, newton_max)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        if !(0.0..=1.0).contains(&sigma_lie) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "sigma_lie : valeur finie dans [0, 1] requise",
            ));
        }
        if sigma_lie != 0.0 && adaptatif.is_some() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "sigma_lie : estimateur adaptatif non validé, utiliser un pas fixe",
            ));
        }
        valide(&[bornes.0, bornes.1], "bornes")?;
        if bornes.0 <= 0.0 || bornes.0 > 1.0 || bornes.1 < 1.0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "0 < borne min ≤ 1 ≤ borne max requis",
            ));
        }
        if let Some(a) = adaptatif {
            valide(&[a], "adaptatif")?;
            if a <= 0.0 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "tolérance adaptative > 0 requise",
                ));
            }
        }
        if let Some((l, c, m)) = pas_contact {
            valide(&[l, c, m], "pas de contact")?;
            if c <= 0.0 || l < c || m < 0.0 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "0 < pas contact ≤ pas libre et marge ≥ 0 requis",
                ));
            }
        }
        if energie && moment {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "moment et energie ne se projettent pas ensemble",
            ));
        }
        if moment {
            if let Some(e) = self.mo.moment_invariant() {
                return Err(pyo3::exceptions::PyValueError::new_err(e));
            }
        }
        if energie {
            if let Some(e) = self.mo.energie_invariante() {
                return Err(pyo3::exceptions::PyValueError::new_err(e));
            }
        }
        if !(h > 0.0) || !h.is_finite() {
            // TROUVÉ PAR SONDAGE, 3 sept. : `h = 0` ne fait pas avancer t, donc
            // la boucle `while self.t < t_end` tourne À L'INFINI — pas une
            // erreur, un blocage, ce qui est pire (aucune trace, aucun
            // diagnostic, un processus à tuer). `h < 0` reculerait.
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "h = {h} : le pas doit être fini et strictement positif"
            )));
        }
        // sondage du 7 sept. : une fin NaN ou passée, ou `tous = 0`, rendaient
        // une trajectoire VIDE sans un mot
        if !t_end.is_finite() || t_end < self.mo.t - 1e-12 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "t_end = {t_end} : la fin doit être finie et ≥ t courant ({})",
                self.mo.t
            )));
        }
        if tous == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "tous = 0 : l'intervalle de sortie est ≥ 1 (tous = 10**9 pour ne rien garder)",
            ));
        }
        if !(0.0..1.0).contains(&rho) {
            // Domaine d'un ARGUMENT : ValueError, comme les gardes d'invariant.
            // Le solveur refuse aussi (String → RuntimeError) ; ici on refuse
            // AVANT, pour que tous les refus de domaine aient le même type.
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "rho = {rho} : ρ∞ doit être dans [0, 1[ — ρ∞ = 1 est instable en index 3"
            )));
        }
        self.mo.adapt = adaptatif;
        self.mo.adapt_bornes = bornes;
        self.mo.proj_moment = moment;
        if moment {
            // On refuse hors du domaine de Noether : projeter là où le moment
            // n'est PAS un invariant fabriquerait une conservation fausse.
            if let Some(cause) = self.mo.moment_invariant() {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "moment=True : le moment cinétique n'est pas un invariant de ce modèle ({cause})")));
            }
            self.mo.l_ref = Some(self.mo.moment());
        } else {
            self.mo.l_ref = None;
        }
        self.mo.proj_energie = energie;
        if energie && moment {
            // RÉSULTAT NÉGATIF, EXACT : pour un corps rigide, ∂E/∂ω = Jˢω et
            // ∂L/∂ω = Jˢ, donc la ligne d'énergie vaut ωᵀ fois les lignes du
            // moment — vérifié à 0.0, rang 3 sur 4. Les deux invariants ne
            // sont donc PAS imposables ensemble par une correction des seules
            // VITESSES : le système 3×3 devient singulier et la correction
            // part n'importe où (mesuré : |L| passe de 3e-7 à 3e-3, dix mille
            // fois pire). Les imposer ensemble demande de corriger aussi q,
            // c'est-à-dire un schéma énergie-moment dans la FORMULATION
            // (Simo–Tarnow), pas une projection.
            return Err(pyo3::exceptions::PyValueError::new_err(
                "moment et energie ne se projettent pas ensemble : pour un corps rigide \
                 ∂E/∂ω = ωᵀ·∂L/∂ω, la pile est de rang déficient et la correction diverge. \
                 Choisir l'invariant qui compte, ou un schéma énergie-moment dans la formulation.",
            ));
        }
        if energie {
            if let Some(cause) = self.mo.energie_invariante() {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "energie=True : l'énergie n'est pas un invariant de ce modèle ({cause})"
                )));
            }
            self.mo.e_ref = Some(self.mo.energie());
        } else {
            self.mo.e_ref = None;
        }
        self.mo.pas_contact = pas_contact;
        self.mo.ccd = ccd;
        self.mo.ggl = ggl;
        self.mo.sigma_lie = sigma_lie;
        let mut traj = vec![];
        let mut k = 0usize;
        // LE GIL EST LIBÉRÉ pendant l'intégration : un ThreadPoolExecutor
        // Python sur des Noyau indépendants tourne alors vraiment en
        // parallèle (mesuré 0,88× avant, GIL gardé). Rien de Python n'est
        // touché dans la boucle — la trajectoire est un Vec de f64.
        let mo = &mut self.mo;
        let res = py.detach(|| {
            mo.simule(t_end, h, rho, tol, newton_max, |t, mo| {
                k += 1;
                if k.is_multiple_of(tous) {
                    traj.push((
                        t,
                        mo.corps.iter().map(|c| [c.r.x, c.r.y, c.r.z]).collect(),
                        mo.corps
                            .iter()
                            .map(|c| {
                                let mut a = [0.0; 9];
                                for i in 0..3 {
                                    for j in 0..3 {
                                        a[3 * i + j] = c.rot[(i, j)];
                                    }
                                }
                                a
                            })
                            .collect(),
                        mo.corps.iter().map(|c| [c.v.x, c.v.y, c.v.z]).collect(),
                        mo.corps.iter().map(|c| [c.w.x, c.w.y, c.w.z]).collect(),
                        mo.lam.iter().copied().collect(),
                    ));
                }
            })
        });
        res.map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        Ok(traj)
    }

    /// AUDIT : jacobien AD de Newton contre différence finie du résidu, à
    /// l'état courant. Rend ([(bloc, écart relatif, max du bloc)], [(i, j, AD, DF)]).
    #[pyo3(signature = (h, rho = 0.9, ggl = false, sigma_lie = 0.0))]
    #[allow(clippy::type_complexity)]
    fn audit_jacobien(
        &mut self,
        h: f64,
        rho: f64,
        ggl: bool,
        sigma_lie: f64,
    ) -> PyResult<(Vec<(String, f64, f64)>, Vec<(usize, usize, f64, f64)>)> {
        if !(0.0..=1.0).contains(&sigma_lie) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "sigma_lie : valeur finie dans [0, 1] requise",
            ));
        }
        let avant = self.mo.sigma_lie;
        self.mo.sigma_lie = sigma_lie;
        let res = self.mo.audit_jacobien(h, rho, ggl);
        self.mo.sigma_lie = avant;
        res.map_err(pyo3::exceptions::PyRuntimeError::new_err)
    }

    /// SOUS-CYCLAGE MULTI-RYTHME : les corps `rapides` à h/k, les autres à h,
    /// couplage par les forces seulement (une liaison qui enjambe est refusée).
    /// Même sortie que `simule`, une image par pas LENT, avec les réactions
    /// des deux partitions. Systèmes résolus à la taille de chaque partition ;
    /// extrapolation des rapides, interpolation des lents, sans correcteur.
    /// Pas d'aérodynamique, contact non lisse, GGL ni adaptatif.
    #[pyo3(signature = (t_end, h, rapides, k, rho = 0.9, tol = 1e-12, newton_max = 25, tous = 1))]
    #[allow(clippy::too_many_arguments)]
    fn simule_multirythme(
        &mut self,
        py: Python<'_>,
        t_end: f64,
        h: f64,
        rapides: Vec<usize>,
        k: usize,
        rho: f64,
        tol: f64,
        newton_max: usize,
        tous: usize,
    ) -> PyResult<
        Vec<(
            f64,
            Vec<[f64; 3]>,
            Vec<[f64; 9]>,
            Vec<[f64; 3]>,
            Vec<[f64; 3]>,
            Vec<f64>,
        )>,
    > {
        reglages(self.mo.t, t_end, h, tol, newton_max)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        if tous == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err("tous ≥ 1 requis"));
        }
        if !(h > 0.0) || k == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "h doit être > 0 et k ≥ 1",
            ));
        }
        if !(0.0..1.0).contains(&rho) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "rho : ρ∞ doit être dans [0, 1[",
            ));
        }
        let nb = self.mo.corps.len();
        let mut fast = vec![false; nb];
        for i in rapides {
            if i >= nb {
                return Err(pyo3::exceptions::PyIndexError::new_err(
                    "corps rapide inconnu",
                ));
            }
            fast[i] = true;
        }
        if !self.mo.pales.is_empty()
            || !self.mo.inflows.is_empty()
            || self.mo.contacts.iter().any(|c| c.nonlisse)
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "multi-rythme : aérodynamique ou contact non lisse hors domaine",
            ));
        }
        for e in &self.mo.elems {
            let cs: Vec<_> = e.corps().into_iter().flatten().collect();
            let nr = cs.iter().filter(|&&i| fast[i]).count();
            if nr != 0 && nr != cs.len() {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "multi-rythme : la liaison « {} » enjambe les partitions",
                    e.nom()
                )));
            }
        }
        self.mo.adapt = None;
        self.mo.ggl = false;
        self.mo.sigma_lie = 0.0;
        self.mo.pas_contact = None;
        self.mo.ccd = false;
        self.mo.proj_moment = false;
        self.mo.proj_energie = false;
        self.mo.l_ref = None;
        self.mo.e_ref = None;
        let mut traj = vec![];
        let mut kk = 0usize;
        let mo = &mut self.mo;
        let res = py.detach(|| {
            mo.simule_multi(t_end, h, k, &fast, rho, tol, newton_max, |t, mo| {
                kk += 1;
                if kk.is_multiple_of(tous) {
                    traj.push((
                        t,
                        mo.corps.iter().map(|c| [c.r.x, c.r.y, c.r.z]).collect(),
                        mo.corps
                            .iter()
                            .map(|c| {
                                let mut a = [0.0; 9];
                                for i in 0..3 {
                                    for j in 0..3 {
                                        a[3 * i + j] = c.rot[(i, j)];
                                    }
                                }
                                a
                            })
                            .collect(),
                        mo.corps.iter().map(|c| [c.v.x, c.v.y, c.v.z]).collect(),
                        mo.corps.iter().map(|c| [c.w.x, c.w.y, c.w.z]).collect(),
                        mo.lam.iter().copied().collect(),
                    ));
                }
            })
        });
        res.map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        Ok(traj)
    }

    /// Intégrateur ÉNERGIE-MOMENT (Simo–Wong, Cayley au point milieu) — E et
    /// L conservés EXACTEMENT, là où `simule(moment=…)` ou `simule(energie=…)`
    /// n'en tient qu'un. Domaine : corps libres sous la seule gravité, refusé
    /// sinon. Même sortie que `simule` (λ vide). Voir `Modele::simule_em`.
    #[pyo3(signature = (t_end, h, tous = 1))]
    fn simule_em(
        &mut self,
        py: Python<'_>,
        t_end: f64,
        h: f64,
        tous: usize,
    ) -> PyResult<
        Vec<(
            f64,
            Vec<[f64; 3]>,
            Vec<[f64; 9]>,
            Vec<[f64; 3]>,
            Vec<[f64; 3]>,
            Vec<f64>,
        )>,
    > {
        reglages(self.mo.t, t_end, h, 1e-12, 30)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        if tous == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err("tous ≥ 1 requis"));
        }
        if !(h > 0.0) {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "h = {h} : le pas doit être strictement positif"
            )));
        }
        if let Some(e) = self.mo.energie_invariante() {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "énergie-moment : hors domaine ({e})"
            )));
        }
        let mut traj = vec![];
        let mut k = 0usize;
        let mo = &mut self.mo;
        let res = py.detach(|| {
            mo.simule_em(t_end, h, |t, mo| {
                k += 1;
                if k.is_multiple_of(tous) {
                    traj.push((
                        t,
                        mo.corps.iter().map(|c| [c.r.x, c.r.y, c.r.z]).collect(),
                        mo.corps
                            .iter()
                            .map(|c| {
                                let mut a = [0.0; 9];
                                for i in 0..3 {
                                    for j in 0..3 {
                                        a[3 * i + j] = c.rot[(i, j)];
                                    }
                                }
                                a
                            })
                            .collect(),
                        mo.corps.iter().map(|c| [c.v.x, c.v.y, c.v.z]).collect(),
                        mo.corps.iter().map(|c| [c.w.x, c.w.y, c.w.z]).collect(),
                        vec![],
                    ));
                }
            })
        });
        res.map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        Ok(traj)
    }

    /// Trace de diagnostic des moments effectivement portés par simule_em.
    /// Exécute une copie du modèle et conserve le point initial et chaque pas.
    /// Aucun certificat de trajectoire ni de propagation des arrondis.
    #[pyo3(signature = (t_end, h))]
    fn _trace_em(&self, py: Python<'_>, t_end: f64, h: f64) -> PyResult<Py<PyAny>> {
        reglages(self.mo.t, t_end, h, 1e-12, 30)
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        self.mo
            .verifie_etat_fini()
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        self.mo
            .verifie_options()
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        let mut mo = self.mo.clone();
        let matrice = |m: &M3| {
            let mut out = [0.0; 9];
            for i in 0..3 {
                for j in 0..3 {
                    out[3 * i + j] = m[(i, j)];
                }
            }
            out
        };
        let inerties: Vec<_> = mo.corps.iter().map(|c| matrice(&c.j)).collect();
        let mut trace = vec![];
        py.detach(|| {
            mo.simule_em_interne(t_end, h, |t, m, pis, _| {
                trace.push((
                    t,
                    pis.iter().map(|p| [p.x, p.y, p.z]).collect::<Vec<_>>(),
                    m.corps.iter().map(|c| matrice(&c.rot)).collect::<Vec<_>>(),
                    m.corps
                        .iter()
                        .map(|c| [c.w.x, c.w.y, c.w.z])
                        .collect::<Vec<_>>(),
                ));
            })
        })
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let out = pyo3::types::PyDict::new(py);
        out.set_item("schema", "vinkulum.trace.moments_em.1")?;
        out.set_item("inerties", inerties)?;
        out.set_item("echantillons", trace)?;
        out.set_item("propagation_certifiee", false)?;
        Ok(out.into_any().unbind())
    }

    /// Arme (ou désarme) l'HISTORIQUE DU SCHÉMA : (u̇, a, λ) à l'état initial
    /// puis après chaque pas accepté du prochain `simule`. Ce que `simule`
    /// rend par frame — (q, u, λ) — ne suffit pas à un adjoint en temps : le
    /// schéma traîne deux accélérations d'un pas à l'autre, et c'est d'elles
    /// que dépend le pas suivant.
    #[pyo3(signature = (actif = true))]
    fn enregistre_schema(&mut self, actif: bool) {
        self.mo.hist_schema = if actif { Some(vec![]) } else { None };
    }

    /// L'historique armé par `enregistre_schema` : liste de (u̇, a, λ), une
    /// entrée pour l'état initial puis une par pas accepté.
    fn schema(&self) -> Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        self.mo.hist_schema.clone().unwrap_or_default()
    }

    /// Pose les multiplicateurs λ de l'état courant (après `pose_etat`) : la
    /// raideur tangente `k_c_m_z` porte la précontrainte Gᵀλ, et sans λ posé
    /// elle le recalcule d'une accélération CONSISTANTE — qui n'est pas le λ
    /// du schéma à ce pas.
    fn pose_lam(&mut self, lam: Vec<f64>) -> PyResult<()> {
        if lam.len() != self.mo.m() {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "pose_lam : {} valeurs pour {} contraintes",
                lam.len(),
                self.mo.m()
            )));
        }
        valide(&lam, "multiplicateurs")?;
        self.mo.lam = DVector::from_vec(lam);
        Ok(())
    }

    /// Le temps courant du modèle, en secondes.
    fn t(&self) -> f64 {
        self.mo.t
    }
    /// chronos cumulés en secondes : (jacobien, résolution, résidus, fin de
    /// pas, TOTAL, tangente de exp). Le reste — total moins la somme des
    /// quatre premiers — est ce que le solveur fait sans l'avoir instrumenté ;
    /// il valait 82 % avant que ces deux-là ne soient posés (3 sept.).
    ///
    /// Le SIXIÈME est INCLUS dans le premier : c'est la part du jacobien
    /// passée à construire les J_l, posée le 5 sept. parce que le refus de
    /// Cayley reposait sur une supposition et non sur une mesure.
    /// Avec sigma_lie non nul, inclut aussi la résolution locale et les
    /// tangentes implicites de configuration nécessaires au jacobien.
    fn chronos(&self) -> (f64, f64, f64, f64, f64, f64) {
        if std::env::var("VINKULUM_CHRONO").is_ok() {
            let ns = |a: &std::sync::atomic::AtomicU64| {
                a.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9
            };
            eprintln!(
                "chronos fins : facto {:.2} s (dont motif {:.2}) · descente {:.2} s · contrôle Jx−b {:.2} s · nnz {}",
                ns(&T_FAC), ns(&T_HASH), ns(&T_SV), ns(&T_MV), NNZ.load(std::sync::atomic::Ordering::Relaxed)
            );
        }
        (
            self.mo.t_jac,
            self.mo.t_sol,
            self.mo.t_res,
            self.mo.t_fin,
            self.mo.t_pas,
            self.mo.t_exp,
        )
    }
    /// (itérations de Newton cumulées, replis SVD cumulés, jacobiens calculés, contraintes redondantes neutralisées)
    fn stats(&self) -> (usize, usize, usize, usize) {
        (
            self.mo.newton_total,
            self.mo.svd_total,
            self.mo.jacobien_total,
            self.mo.actif.iter().filter(|a| !**a).count(),
        )
    }
    /// (pas rejoués, dernier pas retenu, max de R_half, max du seuil, max de |Φ| à mi-pas relatif)
    fn adapt_stats(&self) -> (usize, f64, f64, f64, f64) {
        (
            self.mo.n_rejeu,
            self.mo.h_dernier,
            self.mo.r_half_max,
            self.mo.seuil_max,
            self.mo.phi_half_max,
        )
    }
    /// Moment cinétique total, au repère monde et autour de l'origine.
    ///
    /// Sans gravité ni effort extérieur, il est CONSERVÉ quelles que soient
    /// les liaisons internes : c'est un contrôle exact, applicable à
    /// n'importe quel modèle sans en connaître la solution.
    fn moment(&self) -> [f64; 3] {
        let l = self.mo.moment();
        [l.x, l.y, l.z]
    }

    /// Énergie mécanique totale : cinétique (translation + rotation, tenseur
    /// transporté en monde) plus potentielle de PESANTEUR. **N'inclut pas**
    /// l'énergie de déformation des poutres ni des contacts par pénalité :
    /// sur un système déformable, sa dérive ne mesure donc pas la dissipation
    /// du schéma. Sur un système de corps rigides liés, si.
    fn energie(&self) -> f64 {
        self.mo.energie()
    }
    /// (nom, multiplicateurs λ) par élément — les EFFORTS DE LIAISON, dans le
    /// repère de la liaison (`ra`) et dans l'ordre `bloque_t` puis `bloque_r` :
    /// forces en newtons, moments en N·m. C'est ce qui dimensionne un axe, un
    /// roulement ou une vis. Les lignes neutralisées comme redondantes portent 0.
    fn reactions(&self) -> PyResult<Vec<(String, Vec<f64>)>> {
        self.mo
            .reactions()
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)
    }
    /// Change un effort déjà déclaré. C'est ce qui rend la CONTINUATION
    /// possible : en grands déplacements, un Newton statique lancé d'un coup
    /// sur la charge finale diverge (mesuré sur la poutre de Princeton :
    /// convergent à θ = 0 où la flèche fait 10 mm, divergent à θ = 90 où elle
    /// en fait 146). On applique la charge par paliers, chacun partant de
    /// l'équilibre du précédent — la pratique de tout code non linéaire.
    fn pose_effort(&mut self, i: usize, f: [f64; 3], m: [f64; 3]) -> PyResult<()> {
        valide(&f, "effort")?;
        valide(&m, "moment")?;
        if i >= self.mo.efforts.len() {
            return Err(pyo3::exceptions::PyIndexError::new_err(format!(
                "effort {i}, mais le modèle n'en a que {}",
                self.mo.efforts.len()
            )));
        }
        let c = self.mo.efforts[i].0;
        self.mo.efforts[i] = (c, v3(f), v3(m));
        Ok(())
    }

    /// ANALYSE STATIQUE — l'équilibre, sans passer par le temps.
    ///
    /// Tout code multicorps en a une (ADAMS *equilibrium*, MBDyn `initial
    /// assembly`, Simscape *steady state*), et vinkulum ne l'avait pas : on
    /// atteignait l'équilibre en AMORTISSANT une trajectoire jusqu'à
    /// l'immobilité. C'est lent, et surtout le résultat dépend d'où l'on
    /// s'arrête — mesuré sur la poutre de Princeton : l'écart à la référence
    /// variait de façon NON MONOTONE avec le maillage (2,62 · 3,82 · 2,66 %),
    /// signature d'un bruit de relaxation et non d'une erreur de modèle.
    ///
    /// On résout donc le vrai système, par Newton :
    ///
    /// ```text
    ///     f(q) + Gᵀλ = 0        (équilibre des forces)
    ///     Φ(q)       = 0        (contraintes)
    /// ```
    ///
    /// Le jacobien est [[K, Gᵀ], [G, 0]] avec K = −∂f/∂q la raideur tangente
    /// que le noyau sait déjà former (précontrainte comprise). Newton amorti :
    /// un pas descend le mérite de forces ou réduit la correction de pose
    /// prédite tout en contrôlant séparément la fermeture des contraintes.
    ///
    /// `tol` est RELATIF au déséquilibre libre initial (norme infinie après
    /// projection des réactions, plancher 1). Un critère ABSOLU en newtons ne
    /// veut rien dire hors du modèle qui l'a vu naître : le même mécanisme
    /// coté en tonnes-millimètres le franchit ou pas selon l'unité.
    /// Les contraintes géométriques sont contrôlées séparément : ‖Φ‖∞ ≤ tol
    /// dans leurs unités (longueurs et radians), indépendamment de la charge.
    /// `t` est la date d'évaluation des contraintes et des couples à loi ;
    /// l'analyse statique ne fait pas avancer l'horloge du modèle.
    ///
    /// Rend (‖résidu‖ final, itérations).
    /// CONTINUATION AUTOMATIQUE. Newton d'une structure en grands
    /// déplacements ne converge pas depuis n'importe où : mesuré le 3 sept.
    /// sur une console de 1 m (EI = 100), `statique` tenait jusqu'à 0,5 N et
    /// STAGNAIT à 2 N — un cas parfaitement banal, et l'utilisateur devait
    /// deviner qu'il fallait monter la charge par paliers. Tous les codes de
    /// structure font cette continuation eux-mêmes (« automatic load
    /// stepping ») ; celui-ci la faisait faire à la main.
    ///
    /// On tente d'abord la charge pleine. En cas d'échec on subdivise —
    /// gravité et efforts multipliés par λ croissant — jusqu'à `paliers_max`
    /// subdivisions. Le nombre de paliers réellement utilisés est publié.
    /// `strict=True` refuse l'acceptation de secours sur stagnation.
    /// Sinon, cette acceptation émet un RuntimeWarning. `statique_info()`
    /// distingue la tolérance atteinte, la stagnation et l'échec, y compris
    /// les paliers intermédiaires. Le couple de retour historique est conservé.
    #[pyo3(signature = (tol = 1e-10, iters = 60, t = 0.0, amorti = true, paliers_max = 512, strict = false))]
    #[allow(clippy::too_many_arguments)]
    fn statique(
        &mut self,
        py: Python<'_>,
        tol: f64,
        iters: usize,
        t: f64,
        amorti: bool,
        paliers_max: usize,
        strict: bool,
    ) -> PyResult<(f64, usize)> {
        self.bilan_statique = None;
        reglages(t, t, 1.0, tol, iters).map_err(pyo3::exceptions::PyValueError::new_err)?;
        if paliers_max == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "paliers_max ≥ 1 requis",
            ));
        }
        self.mo
            .verifie_etat_fini()
            .map_err(pyo3::exceptions::PyValueError::new_err)?;
        self.bilan_statique = Some(BilanStatique {
            statut: "en_cours",
            tol,
            strict,
            residu_libre: None,
            contraintes: None,
            echelle_force: None,
            evaluations: 0,
            tentatives: 0,
            paliers: 1,
            palier: 0,
            paliers_stagnation: 0,
            etat_restaure: false,
            message: None,
        });
        let sauve = Sauvegarde::new(&self.mo);
        let g0 = self.mo.g;
        let ef0 = self.mo.efforts.clone();
        // ⚠ L'ÉTAT DE DÉPART SE RESTAURE À CHAQUE TENTATIVE. Sans ça, la
        // subdivision repart de la pose laissée par l'essai qui vient
        // d'échouer — donc d'un point pire — et il faut bien plus de paliers
        // (mesuré : 64 au lieu de 10 sur une console à 5 N).
        let poses0: Vec<Pose> = self.mo.corps.iter().map(Corps::pose_precise).collect();
        let mut np_ = 1usize;
        loop {
            let bilan = self.bilan_statique.as_mut().expect("bilan statique");
            bilan.paliers = np_;
            bilan.paliers_stagnation = 0;
            for (i, c) in self.mo.corps.iter_mut().enumerate() {
                c.r = poses0[i].0;
                c.r_bas = poses0[i].2;
                c.rot = poses0[i].1;
            }
            self.mo.lam = DVector::zeros(self.mo.m());
            let mut sortie = None;
            let mut ok = true;
            for k in 1..=np_ {
                self.bilan_statique.as_mut().expect("bilan statique").palier = k;
                let lam = k as f64 / np_ as f64;
                self.mo.g = g0 * lam;
                for (j, e) in self.mo.efforts.iter_mut().enumerate() {
                    e.1 = ef0[j].1 * lam;
                    e.2 = ef0[j].2 * lam;
                }
                // Priorité au mérite en forces, avec une acceptation
                // supplémentaire par corrections et filtre géométrique.
                // En cas d'échec, reprendre le même état avec le mérite
                // naturel historique, adapté aux raideurs étagées.
                let debut = amorti.then(|| Sauvegarde::new(&self.mo));
                let resultat = match self._statique_newton(tol, iters, t, amorti, false) {
                    Err(_) if amorti => {
                        debut.expect("sauvegarde du palier").restaure(&mut self.mo);
                        self._statique_newton(tol, iters, t, amorti, true)
                    }
                    resultat => resultat,
                };
                match resultat {
                    Ok(v) => {
                        let bilan = self.bilan_statique.as_mut().expect("bilan statique");
                        if bilan.statut == "stagnation" {
                            bilan.paliers_stagnation += 1;
                        }
                        sortie = Some(v);
                    }
                    Err(e) => {
                        ok = false;
                        if np_ >= paliers_max {
                            self.mo.g = g0;
                            self.mo.efforts = ef0;
                            sauve.restaure(&mut self.mo);
                            let bilan = self.bilan_statique.as_mut().expect("bilan statique");
                            bilan.statut = "echec";
                            bilan.etat_restaure = true;
                            bilan.message = Some(e.to_string());
                            return Err(pyo3::exceptions::PyRuntimeError::new_err(e.to_string()));
                        }
                        break;
                    }
                }
            }
            if ok {
                self.mo.g = g0;
                self.mo.efforts = ef0;
                let (r, it) = sortie.expect("au moins un palier");
                let bilan = self.bilan_statique.as_mut().expect("bilan statique");
                if bilan.paliers_stagnation > 0 {
                    bilan.statut = "stagnation";
                    let avertissement = PyErr::warn(
                        py, py.get_type::<pyo3::exceptions::PyRuntimeWarning>().as_any(),
                        c"statique : arret sur stagnation sans atteindre tol a tous les paliers ; consulter statique_info(). strict=True pour refuser.", 1);
                    if let Err(e) = avertissement {
                        sauve.restaure(&mut self.mo);
                        bilan.statut = "echec";
                        bilan.etat_restaure = true;
                        bilan.message = Some(e.to_string());
                        return Err(e);
                    }
                }
                return Ok((r, if np_ == 1 { it } else { it + 1000 * np_ }));
            }
            np_ = np_.saturating_mul(2).min(paliers_max);
        }
    }

    /// Rapport historique du dernier calcul statique ; None si aucun rapport
    /// n'est disponible. `statut` vaut "tolerance",
    /// "stagnation" (au moins un palier retenu) ou "echec". Les résidus et
    /// l'échelle sont ceux de la dernière évaluation Newton, avant restauration
    /// en cas d'échec ; ils ne diagnostiquent pas l'état courant du modèle.
    /// `evaluations` et `tentatives` incluent les reprises abandonnées.
    /// `paliers_stagnation` compte les paliers acceptés sur la dernière
    /// subdivision tentée, même si un échec ultérieur restaure tout le modèle.
    fn statique_info<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Option<Bound<'py, pyo3::types::PyDict>>> {
        let Some(b) = &self.bilan_statique else {
            return Ok(None);
        };
        let out = pyo3::types::PyDict::new(py);
        out.set_item("statut", b.statut)?;
        out.set_item("tol", b.tol)?;
        out.set_item("strict", b.strict)?;
        out.set_item("residu_libre", b.residu_libre)?;
        out.set_item("contraintes", b.contraintes)?;
        out.set_item("echelle_force", b.echelle_force)?;
        let relatif = b.residu_libre.zip(b.echelle_force).map(|(r, e)| r / e);
        out.set_item("residu_relatif", relatif)?;
        out.set_item(
            "tolerance_finale_atteinte",
            relatif
                .zip(b.contraintes)
                .map(|(r, p)| r <= b.tol && p <= b.tol),
        )?;
        out.set_item("evaluations", b.evaluations)?;
        out.set_item("tentatives", b.tentatives)?;
        out.set_item("paliers", b.paliers)?;
        out.set_item("palier", b.palier)?;
        out.set_item("paliers_stagnation", b.paliers_stagnation)?;
        out.set_item("etat_restaure", b.etat_restaure)?;
        out.set_item("message", &b.message)?;
        Ok(Some(out))
    }

    /// CONTRÔLE DU MODÈLE — « est-il bien posé ? », la question qu'un
    /// utilisateur nouveau ne sait pas encore poser.
    ///
    /// Rend la liste des anomalies, vide si tout va bien. Ce sont des choses
    /// qui ne peuvent PAS être des refus, parce qu'elles ont un usage
    /// légitime, mais qu'on ne doit pas non plus taire :
    ///
    /// · **inertie hors de l'inégalité triangulaire** (J₁ > J₂ + J₃) : aucun
    ///   solide réel ne l'a. Ce n'est pas un refus parce qu'un benchmark
    ///   publié le fait délibérément — `multibarmech` met un `1.` de
    ///   remplissage sur l'axe qui ne travaille pas ;
    /// · **contraintes non satisfaites** à l'état courant (Φ ≠ 0) : le
    ///   mécanisme n'est pas assemblé, `assemble()` est là pour ça ;
    /// · **vitesses incompatibles** (Φ̇ ≠ 0) : le premier pas les corrigera en
    ///   silence, en modifiant ce que l'utilisateur a écrit ;
    /// · **corps sans liaison, force ou interaction déclarée**.
    /// Date courante par défaut ; le contrôle ne modifie pas le modèle.
    #[pyo3(signature = (t = None))]
    fn controle(&self, t: Option<f64>) -> PyResult<Vec<String>> {
        let t = t.unwrap_or(self.mo.t);
        fini(&[t], "temps du diagnostic").map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        self.mo
            .verifie_etat_fini()
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let mut out = Vec::new();
        for c in &self.mo.corps {
            fini(c.j.as_slice(), "inertie du diagnostic")
                .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
            let scale = c.j.amax();
            if scale <= 0.0 {
                return Err(pyo3::exceptions::PyRuntimeError::new_err("inertie nulle"));
            }
            let j = c.j / scale;
            let e = numerique::propres_sym(DMatrix::from_column_slice(3, 3, j.as_slice()))
                .map_err(pyo3::exceptions::PyRuntimeError::new_err)?
                .eigenvalues;
            let (a, b, d) = (e[0], e[1], e[2]);
            for (x, y, z, q) in [(a, b, d, 1), (b, a, d, 2), (d, a, b, 3)] {
                if x > y + z + 1e-12 * x {
                    out.push(format!(
                        "{} : inertie hors de l'inégalité triangulaire (J{q}/max|J| = {x:.4e} > {:.4e}) — aucun solide réel n'a ce tenseur",
                        c.nom,
                        y + z
                    ));
                    break;
                }
            }
        }
        let phi = self
            .mo
            .phi_seul(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        fini(phi.as_slice(), "contraintes du diagnostic")
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        let e = self
            .mo
            .phi_holonome_max(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
        if e > 1e-9 * self.mo.echelle_l().max(1.0) {
            out.push(format!(
                        "contraintes non satisfaites : ‖Φ‖ = {e:.3e} — le mécanisme n'est pas assemblé (voir assemble())"
            ));
        }
        let d = self
            .mo
            .phi_dot_a(t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?
            .amax();
        if d > 1e-9 {
            out.push(format!(
                        "vitesses incompatibles avec les liaisons : ‖G·u + ∂Φ/∂t‖ = {d:.3e} (voir assemble(vitesses=True))"
            ));
        }
        let mut vus = vec![false; self.mo.corps.len()];
        for el in &self.mo.elems {
            for r in el.corps().into_iter().flatten() {
                vus[r] = true;
            }
        }
        for p in &self.mo.poutres {
            vus[p.a] = true;
            vus[p.b] = true;
        }
        for s in &self.mo.supers {
            for &i in &s.noeuds {
                vus[i] = true;
            }
        }
        for c in &self.mo.couples {
            for i in [c.a, c.b].into_iter().flatten() {
                vus[i] = true;
            }
        }
        for c in &self.mo.contacts {
            vus[c.corps] = true;
            if let Some(i) = c.b {
                vus[i] = true;
            }
        }
        for &(i, _, _) in &self.mo.spheres {
            vus[i] = true;
        }
        for &(i, _, _) in &self.mo.efforts {
            vus[i] = true;
        }
        for p in &self.mo.pales {
            vus[p.corps] = true;
        }
        for (i, v) in vus.iter().enumerate() {
            if !v {
                out.push(format!(
                    "{} : aucune liaison, force ou interaction déclarée sur ce corps",
                    self.mo.corps[i].nom
                ));
            }
        }
        Ok(out)
    }

    /// Résout les contraintes de position holonomes, puis G·u + ∂Φ/∂t = 0.
    ///
    /// Rend (résidu maximal de position, nombre de corrections). `iters`
    /// inclut la dernière correction contrôlée. `vitesses=False` limite
    /// l'opération aux positions ; les lignes `nh=True` agissent seulement
    /// sur les vitesses. La date `t` ne modifie pas l'horloge du modèle.
    ///
    /// La projection emploie les dérivées des lois : pente du segment gauche
    /// à un nœud intérieur d'une table, du segment adjacent aux extrémités,
    /// zéro hors table. Chaque résidu de vitesse est contrôlé à `tol` près,
    /// augmenté de 64 ε fois la somme absolue des termes de son équation.
    /// L'acceptation finale est vérifiée en arithmétique entière exacte sur
    /// les coefficients binary64 de G, u et ∂Φ/∂t. Ce contrôle ne borne pas
    /// l'erreur de ces coefficients ni la distance à la solution physique.
    /// Des commandes incompatibles ou une précision insuffisante provoquent
    /// une ValueError. Tout échec restitue l'état mécanique antérieur complet,
    /// même après une correction de position réussie.
    #[pyo3(signature = (t = 0.0, tol = 1e-12, iters = 60, vitesses = true))]
    fn assemble(
        &mut self,
        t: f64,
        tol: f64,
        iters: usize,
        vitesses: bool,
    ) -> PyResult<(f64, usize)> {
        let sauve = Sauvegarde::new(&self.mo);
        let resultat: Result<(f64, usize), String> = (|| {
            let r = self.mo.assemble_interne(t, tol, iters)?;
            if vitesses {
                self.mo.projette_vitesse(t, tol)?;
            }
            Ok(r)
        })();
        if resultat.is_err() {
            sauve.restaure(&mut self.mo);
        }
        resultat.map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Copie privée de géométrie pour qualifier le prototype de certificat.
    /// Toutes les contraintes originales sont conservées, y compris inactives.
    /// Aucun assemblage, projection ni changement d'état n'est effectué.
    fn _geometrie_pour_certificat<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::PyDict;
        let refus = |raison: &str| {
            pyo3::exceptions::PyValueError::new_err(format!("certification géométrique : {raison}"))
        };
        if !(1..=4).contains(&self.mo.corps.len()) {
            return Err(refus("nombre de corps hors domaine (1 à 4 requis)"));
        }
        if !self.mo.contacts.is_empty()
            || !self.mo.spheres.is_empty()
            || !self.mo.maillages.is_empty()
        {
            return Err(refus("contacts et maillages hors domaine"));
        }
        if !self.mo.poutres.is_empty() || !self.mo.supers.is_empty() {
            return Err(refus("poutres et superéléments hors domaine"));
        }
        if self.mo.gel.iter().any(|v| *v) {
            return Err(refus("corps gelés hors domaine"));
        }
        let matrice = |a: &M3| {
            (0..3)
                .map(|i| (0..3).map(|j| a[(i, j)]).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        };
        let mut poses = Vec::new();
        for c in &self.mo.corps {
            if c.r
                .iter()
                .chain(c.r_bas.iter())
                .chain(c.rot.iter())
                .any(|x| !x.is_finite())
            {
                return Err(refus("pose non finie"));
            }
            let p = PyDict::new(py);
            p.set_item("position", c.r.as_slice())?;
            p.set_item("position_basse", c.r_bas.as_slice())?;
            p.set_item("rotation", matrice(&c.rot))?;
            poses.push(p);
        }
        let mut contraintes = Vec::new();
        for elem in &self.mo.elems {
            let d = PyDict::new(py);
            match elem {
                Elem::L(l) => {
                    if l.nh {
                        return Err(refus("liaison non holonome hors domaine"));
                    }
                    if l.cible_t.is_some() || l.cible_r.is_some() {
                        return Err(refus("loi imposée hors domaine"));
                    }
                    if l.pa
                        .iter()
                        .chain(l.pb.iter())
                        .chain(l.ra.iter())
                        .chain(l.rb.iter())
                        .any(|x| !x.is_finite())
                    {
                        return Err(refus("point ou repère de liaison non fini"));
                    }
                    d.set_item("type", "liaison")?;
                    d.set_item("a", l.a)?;
                    d.set_item("b", l.b)?;
                    d.set_item("pa", l.pa.as_slice())?;
                    d.set_item("pb", l.pb.as_slice())?;
                    d.set_item("ra", matrice(&l.ra))?;
                    d.set_item("rb", matrice(&l.rb))?;
                    d.set_item("bt", &l.bt)?;
                    d.set_item("br", &l.br)?;
                }
                Elem::D(l) => {
                    if !l.l.is_finite()
                        || l.l <= 0.0
                        || l.pa.iter().chain(l.pb.iter()).any(|x| !x.is_finite())
                    {
                        return Err(refus("distance positive et points finis requis"));
                    }
                    d.set_item("type", "distance")?;
                    d.set_item("a", l.a)?;
                    d.set_item("b", l.b)?;
                    d.set_item("pa", l.pa.as_slice())?;
                    d.set_item("pb", l.pb.as_slice())?;
                    d.set_item("longueur", l.l)?;
                }
                _ => return Err(refus("type de contrainte hors domaine")),
            }
            contraintes.push(d);
        }
        let out = PyDict::new(py);
        out.set_item("corps", self.mo.corps.len())?;
        out.set_item("contraintes", contraintes)?;
        out.set_item("poses", poses)?;
        out.set_item("lignes_originales", self.mo.m())?;
        Ok(out)
    }

    /// Export privé vers le certificat Python : système exact de contributions
    /// et solution native. Toutes les évaluations portent sur un clone.
    #[pyo3(signature = (t = None, dimension_max = 64, redondances = false))]
    #[allow(clippy::type_complexity)]
    fn _initialisation_pour_certificat(
        &self,
        t: Option<f64>,
        dimension_max: usize,
        redondances: bool,
    ) -> PyResult<(
        Vec<(usize, usize, f64)>,
        Vec<f64>,
        Vec<f64>,
        usize,
        Vec<bool>,
        Vec<(usize, usize, f64)>,
        Vec<f64>,
    )> {
        let extrait = || -> Result<_, String> {
            if !(1..=128).contains(&dimension_max) {
                return Err(
                    "certificat initialisation : dimension_max doit être dans [1, 128]".into(),
                );
            }
            let dim = self
                .mo
                .n()
                .checked_add(self.mo.m())
                .ok_or("dimension excessive")?;
            if dim == 0 || dim > dimension_max {
                return Err(
                    "certificat initialisation : dimension vide ou supérieure au budget".into(),
                );
            }
            let date = t.unwrap_or(self.mo.t);
            fini(&[date], "temps du certificat")?;
            self.mo.verifie_options()?;
            self.mo.verifie_etat_fini()?;
            if self.mo.contacts.iter().any(|c| c.nonlisse)
                || self.mo.gel.iter().any(|&v| v)
                || self.mo.gel_elem.iter().any(|&v| v)
                || self.mo.schema_init.is_some()
            {
                return Err("certificat initialisation : contact non lisse, gel ou schéma imposé hors périmètre".into());
            }
            let mut mo = self.mo.clone();
            mo.t = date;
            if mo.actif.len() != mo.m() {
                mo.detecte_redondance(date)?;
            }
            if !redondances && mo.actif.iter().any(|&v| !v) {
                return Err("certificat initialisation : contraintes désactivées ou redondantes hors périmètre".into());
            }
            let mut temoin = None;
            mo.acc_init_temoin(date, Some(&mut temoin))?;
            let s = temoin.ok_or("certificat initialisation : témoin absent")?;
            Ok((
                s.a.iter().map(|t| (t.row, t.col, t.val)).collect(),
                s.b.as_slice().to_vec(),
                s.x.as_slice().to_vec(),
                mo.n(),
                s.actives,
                s.g_complet,
                s.c_complet.as_slice().to_vec(),
            ))
        };
        extrait().map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Φ̇ = G·u + ∂Φ/∂t à état fixé ; date courante par défaut, sans mutation.
    #[pyo3(signature = (t = None))]
    fn phi_dot(&self, t: Option<f64>) -> PyResult<Vec<f64>> {
        Ok(self
            .mo
            .phi_dot_a(t.unwrap_or(self.mo.t))
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?
            .iter()
            .copied()
            .collect())
    }
    /// Φ à l'état courant : la VIOLATION de chaque ligne de contrainte, dans
    /// l'ordre des éléments (longueurs en mètres, rotations en radians). Un
    /// modèle assemblé y est à la précision machine ; `assemble()` l'y ramène.
    fn phi(&self) -> PyResult<Vec<f64>> {
        Ok(self
            .mo
            .phi_seul(self.mo.t)
            .map_err(pyo3::exceptions::PyRuntimeError::new_err)?
            .iter()
            .copied()
            .collect())
    }
    /// (position, rotation en ligne) d'UN corps — le raccourci quand `etat()`,
    /// qui rend tout le modèle, est trop bavard.
    fn pose(&self, corps: usize) -> PyResult<([f64; 3], [f64; 9])> {
        self._ref1("pose", corps)?;
        let c = &self.mo.corps[corps];
        let mut a = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                a[3 * i + j] = c.rot[(i, j)];
            }
        }
        Ok(([c.r.x, c.r.y, c.r.z], a))
    }
}

#[pymodule(gil_used = false)]
fn _vinkulum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Noyau>()?;
    certificat::enregistrer(m)?;
    certificat_vitesse::enregistrer(m)?;
    // constantes de Leishman–Beddoes : lues par les bancs, jamais recopiées
    m.add("LB_TP", LB_TP)?;
    m.add("LB_TF", LB_TF)?;
    m.add("LB_TV", LB_TV)?;
    m.add("LB_TVL", LB_TVL)?;
    Ok(())
}

// ── contrôle en Rust : le pendule, pour que `cargo test` garde la brique seul ──
impl Noyau {
    fn _statique_newton(
        &mut self,
        tol: f64,
        iters: usize,
        t: f64,
        amorti: bool,
        merite_naturel: bool,
    ) -> PyResult<(f64, usize)> {
        let err = |e: String| pyo3::exceptions::PyValueError::new_err(e);
        let bilan = self.bilan_statique.as_mut().expect("bilan statique");
        bilan.tentatives += 1;
        bilan.statut = "en_cours";
        // En cas de refus avant la première évaluation, ne pas publier les
        // résidus de la tentative précédente comme ceux de cette tentative.
        bilan.residu_libre = None;
        bilan.contraintes = None;
        bilan.echelle_force = None;
        let strict = bilan.strict;
        let n = self.mo.n();
        let mut norme = f64::INFINITY;
        // Échelle figée du déséquilibre LIBRE initial, avec plancher 1.
        // Les efforts repris par les contraintes ne doivent pas agrandir
        // la tolérance : 1e12 N sur un axe bloqué faisait accepter sans
        // itération un déséquilibre de 1 N sur l'axe libre d'un ressort.
        // Le critère reste relatif au déséquilibre initial admissible ;
        // ce n'est ni une tolérance absolue ni une échelle de charge imposée.
        let mut ech = 0.0_f64;
        // Newton s'arrête aussi quand il ne PROGRESSE plus : la précision
        // machine est un plancher, pas un échec.
        //
        // ⚠ SUR UNE FENÊTRE, PAS D'UNE ITÉRATION À L'AUTRE. Première version
        // (3 sept.) : « moins de 1 % de progrès cinq fois de suite ». Elle
        // coupait un Newton qui convergeait très bien — un Newton AMORTI
        // descend par paliers avant d'entrer en régime quadratique, et sur ce
        // même cas il tombait de 3,3e4 à 3,6e-5 en trente itérations. Mesuré
        // par le banc d'invariance d'unités, qui l'a pris pour une limite
        // d'échelle alors que c'était le garde-fou lui-même.
        const FEN: usize = 10;
        let mut fen = [f64::INFINITY; FEN];
        // Une diminution de la correction peut s'accompagner d'une hausse
        // des forces sur les directions très raides. La fenêtre des forces
        // ne doit pas condamner ces pas acceptés dans la métrique de Newton.
        let mut fen_naturelle = [false; FEN];
        let autorise_contraction = !merite_naturel;
        let trace = std::env::var("VINKULUM_TRACE").is_ok();
        let mut symbolique = None;
        let mut motif_precedent = Vec::new();
        let mut dimension_precedente = 0;
        let mut projection = contraintes::Projection::default();
        let mut raffinement = raffinement::Cache::default();
        let mut equilibrage = equilibrage::Cache::default();
        for it in 0..iters {
            let chrono = trace.then(std::time::Instant::now);
            let temps = || chrono.map_or(0.0, |c| c.elapsed().as_secs_f64());
            let (phi, g) = self.mo.phi_g_locale(t).map_err(err)?;
            let m = phi.len();
            let f = self.mo.forces_a(t);
            fini(phi.as_slice(), "contraintes statiques").map_err(err)?;
            // Les lignes identiques et opposées, avec le même écart signé,
            // partagent une coordonnée orthonormale. On reconstruit toutes
            // les réactions avant les résidus et la raideur de précontrainte.
            let equivalentes = contraintes::Equivalentes::nouvelles(&g, &phi);
            let (gr, pr) = equivalentes
                .as_ref()
                .map_or((&g, &phi), |e| (&e.gradient, &e.phi));
            let mr = pr.len();
            // λ des moindres carrés sur Gᵀ, composante par composante.
            // La projection ne forme plus G dense ni la matrice normale GGᵀ.
            let mu = gr.reactions(&f, &mut projection).map_err(err)?;
            let lam = match &equivalentes {
                Some(e) => e.releve(&mu),
                None => mu,
            };
            fini(f.as_slice(), "forces statiques").map_err(err)?;
            fini(lam.as_slice(), "réactions statiques").map_err(err)?;
            let r_f = g.residu(&f, &lam); // équilibre hors contraintes
            fini(r_f.as_slice(), "résidu statique").map_err(err)?;
            if it == 0 {
                ech = r_f.amax().max(1.0);
            }
            norme = r_f.amax().max(phi.amax());
            let bilan = self.bilan_statique.as_mut().expect("bilan statique");
            bilan.evaluations += 1;
            bilan.residu_libre = Some(r_f.amax());
            bilan.contraintes = Some(phi.amax());
            bilan.echelle_force = Some(ech);
            // « pourquoi ma statique a-t-elle échoué » est LA question d'un
            // utilisateur, et le noyau n'y répondait que par un nombre final.
            // VINKULUM_TRACE=1 donne la descente terme par terme : on voit si
            // Newton stagne, remonte, ou se fait couper par le garde-fou.
            if trace {
                eprintln!(
                    "  statique it {it:3}  ‖r‖ {norme:.4e}  (forces {:.3e} · Φ {:.3e} · échelle libre {ech:.3e})",
                    r_f.amax(),
                    phi.amax()
                );
            }
            // ⚠ STOCKER λ AVANT DE RENDRE. Il l'était plus bas, après le
            // test de convergence : le λ de l'itération qui CONVERGE n'était
            // donc jamais gardé, et `reactions()` rendait celui d'AVANT — ou
            // rien du tout sur un mécanisme à zéro degré libre, où Newton
            // converge dès it = 0. Or la statique est précisément là où on
            // veut les efforts de liaison. (Mesuré le 3 sept. sur un hexapode
            // de Gough–Stewart : six vérins rendus à 0 N pour 251 N réels.)
            self.mo.lam = lam.clone();
            let merite = (r_f.amax() / ech).max(phi.amax());
            if merite <= tol {
                self.bilan_statique.as_mut().expect("bilan statique").statut = "tolerance";
                return Ok((norme, it));
            }
            // Plancher numérique : on accepte une stagnation seulement si
            // l'équilibre est déjà précis. Une remontée des forces pendant
            // la recherche ne prouve pas une stagnation des corrections de
            // pose (Princeton anisotrope) ; le budget borne les autres cas.
            let vieux = fen[it % FEN];
            fen[it % FEN] = merite;
            if !strict
                && it >= FEN
                && !fen_naturelle.iter().any(|&v| v)
                && merite > 0.5 * vieux
                && r_f.amax() <= 1e-6 * ech
                && phi.amax() <= tol
            {
                self.bilan_statique.as_mut().expect("bilan statique").statut = "stagnation";
                return Ok((norme, it));
            }
            if !merite_naturel
                && !fen_naturelle.iter().any(|&v| v)
                && it >= FEN
                && merite > 0.5 * vieux
            {
                return Err(err(format!(
                    "statique : Newton stagne à ‖r‖ = {norme:.3e} pour une échelle de force de {ech:.3e}"
                )));
            }
            // système augmenté [[K, Gᵀ], [G, 0]] · [δq ; δλ] = [r_f ; −Φ]
            //
            // ⚠ K CONTIENT LA PRÉCONTRAINTE ∂(Gᵀλ)/∂q, et elle se lit sur le λ
            // STOCKÉ. Sans cette ligne, la raideur est celle d'un λ périmé (ou
            // nul) : Newton part dans une direction qui n'est pas la sienne et
            // n'avance pas — mesuré, ‖r‖ bloqué à 6,4 sur la poutre de
            // Princeton, l'ordre de grandeur de la charge elle-même.
            let t_forces = temps();
            let mut trip = self.mo.raideur_locale(t).map_err(err)?;
            let t_raideur = temps();
            for v in &gr.trip {
                trip.push(faer::sparse::Triplet {
                    row: n + v.row,
                    col: v.col,
                    val: v.val,
                });
                trip.push(faer::sparse::Triplet {
                    row: v.col,
                    col: n + v.row,
                    val: v.val,
                });
            }
            analyse::regroupe(&mut trip);
            let mut b = DVector::zeros(n + mr);
            b.rows_mut(0, n).copy_from(&r_f);
            if mr > 0 {
                b.rows_mut(n, mr).copy_from(&(-pr));
            }
            // D⁻¹ A D⁻¹, à partir des lignes et colonnes : une diagonale
            // presque nulle ne doit pas amplifier les termes de liaison.
            // Les trois composantes de chaque vecteur spatial partagent
            // leur échelle ; les coefficients gardent leur valeur par
            // multiplication exacte par des puissances de deux.
            let d = equilibrage::statique(&mut trip, n, n + mr, &mut equilibrage).map_err(err)?;
            let b_s = equilibrage::vecteur(b, &d).map_err(err)?;
            let t_assemblage = temps();
            if trip.iter().any(|v| !v.val.is_finite()) {
                return Err(err("matrice statique équilibrée : valeur non finie".into()));
            }
            let dim = n + mr;
            let optimisee = if dim > DENSE_MAX {
                // Ne supprimer aucun petit coefficient : seul zéro est omis.
                // Le motif peut changer avec la configuration (rotation,
                // activation d'un contact), donc vérifier avant réutilisation.
                let motif: Vec<_> = trip.iter().map(|v| (v.row, v.col)).collect();
                if dim != dimension_precedente || motif != motif_precedent {
                    symbolique = None;
                    motif_precedent = motif;
                    dimension_precedente = dim;
                }
                if trip.len() * 4 < dim * dim {
                    Facto::creux(&trip, dim, &mut symbolique).ok()
                } else {
                    Some(Facto::dense_seule(&trip, dim))
                }
            } else {
                None
            };
            // Le refus du LU creux ne déclenche pas un second LU dense
            // avant d'essayer le raffinement. Seul le repli orthogonal
            // ci-dessous densifie un grand système non résolu.
            let a_s = (dim <= DENSE_MAX).then(|| analyse::dense(&trip, dim));
            let lu = a_s.as_ref().map(|a| a.clone().lu());
            let resout_direct = |rhs: &DVector<f64>| {
                if let Some(f) = &optimisee {
                    let x = f.solve(rhs);
                    if x.iter().all(|v| v.is_finite())
                        && (trip_mv(&trip, dim, &x) - rhs).norm() <= 1e-9 * (1.0 + rhs.norm())
                    {
                        Some(x)
                    } else {
                        None
                    }
                } else {
                    lu.as_ref().and_then(|lu| lu.solve(rhs)).filter(|x| {
                        x.iter().all(|v| v.is_finite())
                            && (trip_mv(&trip, dim, x) - rhs).norm() <= 1e-9 * (1.0 + rhs.norm())
                    })
                }
            };
            let mut solution = resout_direct(&b_s);
            // Régulariser uniquement la factorisation auxiliaire, puis
            // résoudre le système original par résidus compensés. Si le
            // raffinement ne converge pas, conserver le repli COD/SVD.
            let auxiliaire = if solution.is_none() && dim > DENSE_MAX && mr > 0 {
                raffinement::Facteur::nouveau(&trip, n, dim, &mut raffinement).ok()
            } else {
                None
            };
            if let Some(f) = &auxiliaire {
                solution = f.resout(&trip, &b_s);
            }
            let resout = |rhs: &DVector<f64>| {
                resout_direct(rhs).or_else(|| auxiliaire.as_ref()?.resout(&trip, rhs))
            };
            let naturel = solution.is_some();
            let dx = match solution {
                Some(v) => v,
                None => match moindres_carres::resout(
                    a_s.unwrap_or_else(|| analyse::dense(&trip, dim)),
                    &b_s,
                    1e-12,
                ) {
                    Ok(v) => v,
                    Err(_) => {
                        return Err(err(format!(
                        "statique : système singulier à l'itération {it} (‖r‖ = {norme:.3e}) — \
                         mécanisme sans raideur dans une direction libre ?")))
                    }
                },
            };
            let dx = equilibrage::vecteur(dx, &d).map_err(err)?;
            let dx = if let Some(e) = &equivalentes {
                let mut complet = DVector::zeros(n + m);
                complet.rows_mut(0, n).copy_from(&dx.rows(0, n));
                complet
                    .rows_mut(n, m)
                    .copy_from(&e.releve(&dx.rows(n, mr).into_owned()));
                complet
            } else {
                dx
            };
            fini(dx.as_slice(), "incrément statique").map_err(err)?;
            let t_resolution = temps();
            // NEWTON AMORTI : une structure en grands déplacements n'est pas
            // linéaire, et un pas plein peut faire empirer le résidu.
            let etat0 = (self.mo.poses(), self.mo.lam.clone());
            let mut sc = 1.0;
            // LE REBROUSSEMENT GARDAIT LE DERNIER ESSAI, PAS LE MEILLEUR :
            // `sc` se divise douze fois par deux, atteint 4,9e-4 sans passer
            // sous son seuil de 1e-4, donc la boucle s'épuisait en LAISSANT
            // appliqué un pas qui avait EMPIRÉ le résidu.
            // Mérite naturel : correction de Newton prédite par la même
            // factorisation. Une norme de force seule pénalise les directions
            // raides même lorsque leur correction de pose est minuscule.
            // Le critère final reste le résidu physique ci-dessus.
            let naturel0 = dx.rows(0, n).component_mul(&d.rows(0, n)).norm();
            let reference_merite = if merite_naturel && naturel {
                naturel0
            } else {
                merite
            };
            let mut meilleur: Option<(f64, f64)> = None;
            let mut contraction_naturelle = false;
            for _ in 0..12 {
                for (i, c) in self.mo.corps.iter_mut().enumerate() {
                    let (p0, r0, bas0) = &etat0.0[i];
                    c.deplace_depuis(
                        *p0,
                        *bas0,
                        sc * V3::new(dx[6 * i], dx[6 * i + 1], dx[6 * i + 2]),
                    );
                    c.rot = expm(&(sc * V3::new(dx[6 * i + 3], dx[6 * i + 4], dx[6 * i + 5]))) * r0;
                    c.rot = rotation::projette(&c.rot).map_err(err)?;
                }
                if !amorti {
                    break;
                }
                let (p2, g2) = self.mo.phi_g_locale(t).map_err(err)?;
                let f2 = self.mo.forces_a(t);
                if !g2.est_fini() || !f2.iter().chain(p2.iter()).all(|v| v.is_finite()) {
                    sc *= 0.5;
                    continue;
                }
                let contrainte2 = p2.amax();
                // ⚠ LE MÉRITE DOIT JUGER CE QUE LE PAS A RÉSOLU.
                //
                // Le système augmenté rend [δq ; δλ], et δλ était JETÉ : le
                // mérite recalculait λ par moindres carrés à chaque essai. Le
                // pas est donc dérivé à λ FIGÉ et jugé à λ RECALCULÉ — deux
                // fonctions différentes. Sur un problème bien conditionné
                // elles coïncident au premier ordre ; sur une raideur tangente
                // quasi singulière (câble, chaîne), non — et le rebroussement
                // rejette alors une direction qui était bonne. Mesuré le
                // 3 sept. : à n = 54 le résidu plafonne à 3,168 AVEC UN PAS
                // NUL, alors que la résolution linéaire est exacte à 1e-13 et
                // que ‖dx‖ vaut 56.
                let n2 = if p2.is_empty() {
                    f2.amax() / ech
                } else {
                    let l2 = &etat0.1 + sc * dx.rows(n, m);
                    (g2.residu(&f2, &l2).amax() / ech).max(p2.amax())
                };
                // Pour un pas refusé par les forces, mesurer la correction
                // prédite avec la factorisation COURANTE, déjà disponible.
                // Une contraction forte peut justifier ce pas même quand la
                // courbure injecte un petit déplacement dans un axe très raide.
                let faisable = contrainte2 <= tol.max((1.0 - 0.25 * sc) * phi.amax());
                let correction = if naturel
                    && (merite_naturel || (autorise_contraction && faisable && n2 >= merite))
                {
                    let pr2 = if let Some(e) = &equivalentes {
                        match e.reduit_phi(&p2) {
                            Some(p) => p,
                            None => {
                                sc *= 0.5;
                                continue;
                            }
                        }
                    } else {
                        p2
                    };
                    let l2 = &etat0.1 + sc * dx.rows(n, m);
                    let rf2 = g2.residu(&f2, &l2);
                    let mut rhs2 = DVector::zeros(n + mr);
                    rhs2.rows_mut(0, n).copy_from(&rf2);
                    if mr > 0 {
                        rhs2.rows_mut(n, mr).copy_from(&(-pr2));
                    }
                    Some(
                        equilibrage::vecteur(rhs2, &d)
                            .ok()
                            .and_then(|rhs| resout(&rhs))
                            .map_or(f64::INFINITY, |v| v.rows(0, n).norm()),
                    )
                } else {
                    None
                };
                // Le modèle linéaire prédit (1-alpha)*dx. Exiger au moins
                // un quart de cette diminution dans la norme des corrections,
                // y compris pour alpha<1, sans accepter un simple petit pas.
                contraction_naturelle = autorise_contraction
                    && naturel0 > 0.0
                    && naturel0.is_finite()
                    // Une petite correction préconditionnée peut masquer
                    // l'ouverture d'une liaison non linéaire. Le filtre
                    // contrôle séparément la faisabilité, dans les unités
                    // et à la tolérance des contraintes originales.
                    && faisable
                    && correction.is_some_and(|c| {
                        c.is_finite() && c <= (1.0 - 0.25 * sc) * naturel0
                    });
                if trace {
                    let corr =
                        correction.map_or_else(|| "non_evaluee".into(), |c| format!("{c:.6e}"));
                    eprintln!("  statique recherche alpha={sc:.6e} forces={n2:.6e} contraintes={contrainte2:.6e} correction={corr} correction0={naturel0:.6e} contraction={contraction_naturelle}");
                }
                let n2 = if merite_naturel {
                    correction.unwrap_or(n2)
                } else {
                    n2
                };
                if meilleur.is_none_or(|(_, b)| n2 < b) {
                    meilleur = Some((sc, n2));
                }
                if contraction_naturelle || n2 < reference_merite || sc < 1e-4 {
                    break;
                }
                sc *= 0.5;
            }
            if amorti {
                let (sb, nb) =
                    meilleur.ok_or_else(|| err("statique : aucun essai admissible".into()))?;
                let s_ret = if contraction_naturelle {
                    // La contraction valide cet essai précis, pas l'essai
                    // précédent qui minimisait éventuellement les forces.
                    sc
                } else if nb < reference_merite {
                    sb
                } else {
                    0.0
                };
                fen_naturelle[it % FEN] = contraction_naturelle;
                if (s_ret - sc).abs() > 0.0 {
                    for (i, c) in self.mo.corps.iter_mut().enumerate() {
                        let (p0, r0, bas0) = &etat0.0[i];
                        c.deplace_depuis(
                            *p0,
                            *bas0,
                            s_ret * V3::new(dx[6 * i], dx[6 * i + 1], dx[6 * i + 2]),
                        );
                        c.rot =
                            expm(&(s_ret * V3::new(dx[6 * i + 3], dx[6 * i + 4], dx[6 * i + 5])))
                                * r0;
                        c.rot = rotation::projette(&c.rot).map_err(err)?;
                    }
                }
            }
            if trace {
                eprintln!("  statique temps forces={:.6} raideur={:.6} assemblage={:.6} resolution={:.6} recherche={:.6}",
                    t_forces, t_raideur-t_forces, t_assemblage-t_raideur,
                    t_resolution-t_assemblage, temps()-t_resolution);
            }
        }
        Err(err(format!(
            "statique : pas convergé en {iters} itérations, ‖r‖ = {norme:.3e} pour une \
             échelle de force de {ech:.3e} ({:.1e} en relatif)",
            norme / ech
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pendule_periode() {
        let (l, th0, g) = (1.0_f64, 5.0_f64.to_radians(), 9.80665_f64);
        let mut mo = Modele::new(V3::new(0.0, 0.0, -g));
        mo.corps.push(Corps {
            r_bas: crate::V3::zeros(),
            nom: "m".into(),
            m: 0.2,
            j: M3::identity() * 1e-6,
            r: V3::new(l * th0.sin(), 0.0, -l * th0.cos()),
            rot: M3::identity(),
            v: V3::zeros(),
            w: V3::zeros(),
        });
        let pb = mo.corps[0].rot.transpose() * (V3::zeros() - mo.corps[0].r);
        mo.elems.push(Elem::L(Liaison {
            nom: "rotule".into(),
            a: None,
            b: Some(0),
            pa: V3::zeros(),
            ra: M3::identity(),
            pb,
            rb: M3::identity(),
            bt: vec![0, 1, 2],
            br: vec![],
            cible_t: None,
            cible_r: None,
            nh: false,
        }));
        let t_th = 2.0 * std::f64::consts::PI * (l / g).sqrt() * (1.0 + th0 * th0 / 16.0);
        let mut xs = vec![];
        mo.simule(10.0 * t_th, t_th / 400.0, 0.9, 1e-12, 25, |t, m| {
            xs.push((t, m.corps[0].r.x))
        })
        .unwrap();
        let mut z = vec![];
        for w in xs.windows(2) {
            if w[0].1 > 0.0 && w[1].1 <= 0.0 {
                z.push(w[0].0 - w[0].1 * (w[1].0 - w[0].0) / (w[1].1 - w[0].1));
            }
        }
        let per = (z[z.len() - 1] - z[0]) / (z.len() - 1) as f64;
        assert!((per / t_th - 1.0).abs() < 2e-4, "période {per} vs {t_th}");
        assert!((mo.corps[0].r.norm() - l).abs() < 1e-10, "contrainte");
    }

    fn corps_test(nom: &str, r: V3, w: V3) -> Corps {
        Corps {
            r_bas: crate::V3::zeros(),
            nom: nom.into(),
            m: 1.0,
            j: M3::identity() * 1e-3,
            r,
            rot: M3::identity(),
            v: V3::zeros(),
            w,
        }
    }

    /// L'appariement automatique change la liste des contacts EN COURS de
    /// `simule` ; le masque `actif` des lignes de contrainte doit y survivre.
    /// Jusqu'au 6 sept. `apparie` le vidait et `jacobien_ad` l'indexait à
    /// vide au pas suivant : panique dès qu'une liaison accompagnait les sphères.
    #[test]
    fn ccd_capsule_balayee() {
        // triangle dans le plan x = 0 ; un segment qui le traverse est à 0,
        // un segment parallèle à x = 0,3 est à 0,3 — minimum plat, section dorée
        let (a, b, c) = (
            V3::new(0.0, -1.0, -1.0),
            V3::new(0.0, 1.0, -1.0),
            V3::new(0.0, 0.0, 1.0),
        );
        let f = |p: &V3| (point_triangle(p, &a, &b, &c) - p).norm();
        assert!(
            min_convexe_segment(f, &V3::new(-1.0, 0.0, 0.0), &V3::new(1.0, 0.0, 0.0)).0 < 1e-12
        );
        let d = min_convexe_segment(f, &V3::new(0.3, -2.0, 0.0), &V3::new(0.3, 2.0, 0.0)).0;
        assert!((d - 0.3).abs() < 1e-12, "{d}");
        // boîte unité : segment vertical au-dessus du coin (2, 2) → √2
        let d = min_convexe_segment(
            |p| dist_boite(p, &V3::new(1.0, 1.0, 1.0)),
            &V3::new(2.0, 2.0, -3.0),
            &V3::new(2.0, 2.0, 3.0),
        )
        .0;
        assert!((d - 2f64.sqrt()).abs() < 1e-12, "{d}");
        // cylindre d'axe z, R 1, L 1 : segment tangent à x = 1,5 → 0,5 ;
        // et au-dessus du fond, à z = 1,25 sur l'axe → 0,25
        let cyl = |p: &V3| dist_cylindre(p, &V3::z(), 1.0, 1.0);
        let d = min_convexe_segment(cyl, &V3::new(1.5, -3.0, 0.0), &V3::new(1.5, 3.0, 0.0)).0;
        assert!((d - 0.5).abs() < 1e-12, "{d}");
        let d = min_convexe_segment(cyl, &V3::new(-3.0, 0.0, 1.25), &V3::new(3.0, 0.0, 1.25)).0;
        assert!((d - 0.25).abs() < 1e-12, "{d}");
        // maillage : la capsule de rayon 0,1 à x = 0,05 touche, celle à x = 0,3 non,
        // et celle qui passe HORS du triangle (y = 5) non plus malgré x = 0
        let (bvh, ordre) = Maillage::bati(&[a, b, c], &[[0, 1, 2]]);
        let m = Maillage {
            nom: "t".into(),
            corps: 0,
            sommets: vec![a, b, c],
            tris: vec![[0, 1, 2]],
            bvh,
            ordre,
        };
        assert!(m.traverse(&V3::new(0.05, -2.0, 0.0), &V3::new(0.05, 2.0, 0.0), 0.1));
        assert!(!m.traverse(&V3::new(0.3, -2.0, 0.0), &V3::new(0.3, 2.0, 0.0), 0.1));
        assert!(!m.traverse(&V3::new(-1.0, 5.0, 0.0), &V3::new(1.0, 5.0, 0.0), 0.1));
        // seg_tri en forme fermée : traversée (0, s = 0,5), parallèle à x = 0,3
        // (0,3), au-delà d'une arête (distance à l'arête), segment dégénéré
        let (dd, sx, _) = seg_tri(
            &V3::new(-1.0, 0.0, 0.0),
            &V3::new(1.0, 0.0, 0.0),
            &a,
            &b,
            &c,
        );
        assert!(dd < 1e-15 && (sx - 0.5).abs() < 1e-15, "{dd} {sx}");
        let (dd, _, _) = seg_tri(
            &V3::new(0.3, -2.0, 0.0),
            &V3::new(0.3, 2.0, 0.0),
            &a,
            &b,
            &c,
        );
        assert!((dd - 0.3).abs() < 1e-15, "{dd}");
        // arête [a, b] à z = −1 : segment à z = −1,5, x = 0 → 0,5, s = 0,5, point (0, 0, −1)
        let (dd, sx, q) = seg_tri(
            &V3::new(0.0, -3.0, -1.5),
            &V3::new(0.0, 3.0, -1.5),
            &a,
            &b,
            &c,
        );
        assert!(
            (dd - 0.5).abs() < 1e-15
                && (1.0 / 3.0 - 1e-15..=2.0 / 3.0 + 1e-15).contains(&sx)
                && q.x.abs() < 1e-15
                && (q.z + 1.0).abs() < 1e-15
                && q.y.abs() <= 1.0,
            "{dd} {sx} {q}"
        );
        let (dd, _, _) = seg_tri(&V3::new(0.2, 0.0, 0.0), &V3::new(0.2, 0.0, 0.0), &a, &b, &c);
        assert!((dd - 0.2).abs() < 1e-15, "{dd}");
        // plus_proche_segment sur le triangle seul = seg_tri
        let (dd, sx, _) =
            m.plus_proche_segment(&V3::new(0.0, -3.0, -1.5), &V3::new(0.0, 3.0, -1.5));
        assert!((dd - 0.5).abs() < 1e-15 && (1.0 / 3.0 - 1e-15..=2.0 / 3.0 + 1e-15).contains(&sx));
    }

    #[test]
    fn appariement_avec_liaison_ne_panique_pas() {
        let mut mo = Modele::new(V3::new(0.0, 0.0, -9.81));
        mo.corps.push(corps_test("a", V3::zeros(), V3::zeros()));
        mo.corps
            .push(corps_test("b", V3::new(0.0, 0.0, 0.5), V3::zeros()));
        mo.elems.push(Elem::L(Liaison {
            nom: "enc".into(),
            a: None,
            b: Some(0),
            pa: V3::zeros(),
            ra: M3::identity(),
            pb: V3::zeros(),
            rb: M3::identity(),
            bt: vec![0, 1, 2],
            br: vec![0, 1, 2],
            cible_t: None,
            cible_r: None,
            nh: false,
        }));
        mo.spheres.push((0, V3::zeros(), 0.1));
        mo.spheres.push((1, V3::zeros(), 0.1));
        mo.loi_contact = Some((1e5, 1.5, 0.0, 0.0, 1e-3, 0.0, 1e-3));
        let n = mo.simule(0.4, 1e-3, 0.9, 1e-12, 25, |_, _| {}).unwrap();
        assert_eq!(n, 400);
        assert_eq!(mo.actif.len(), mo.m());
        // la paire est apparue vers t = 0,25 s et la bille est repartie vers le haut
        assert!(mo.corps[1].v.z > 0.0, "v_b = {}", mo.corps[1].v.z);
    }

    /// Un engrenage dont l'angle saute d'un quart de tour sur l'itéré prédit :
    /// le résidu vaut 1e30, le jacobien n'existe pas. `simule` doit rendre
    /// l'erreur de l'élément, pas paniquer (6 sept.).
    #[test]
    fn engrenage_angle_saute_rend_err() {
        let mut mo = Modele::new(V3::zeros());
        mo.corps
            .push(corps_test("a", V3::zeros(), V3::new(0.0, 0.0, 100.0)));
        mo.corps.push(corps_test(
            "b",
            V3::new(0.1, 0.0, 0.0),
            V3::new(0.0, 0.0, -100.0),
        ));
        mo.elems.push(Elem::E(Engrenage {
            nom: "e".into(),
            a: Some(0),
            b: Some(1),
            c: None,
            na: V3::z(),
            nb: V3::z(),
            rapport: -1.0,
            ref_a: M3::identity(),
            ref_b: M3::identity(),
            prev: (0.0, 0.0),
        }));
        // h·ω = 2 rad > π/2 : l'angle ne se déroule plus sur ce pas
        let r = mo.simule(0.1, 0.02, 0.9, 1e-12, 25, |_, _| {});
        let msg = r.expect_err("Err attendu");
        assert!(msg.contains("quart de tour"), "{msg}");
        // en adaptatif, le même modèle est rejoué à pas plus court et passe
        let mut mo2 = Modele::new(V3::zeros());
        mo2.corps
            .push(corps_test("a", V3::zeros(), V3::new(0.0, 0.0, 100.0)));
        mo2.corps.push(corps_test(
            "b",
            V3::new(0.1, 0.0, 0.0),
            V3::new(0.0, 0.0, -100.0),
        ));
        mo2.elems = mo.elems.clone();
        mo2.adapt = Some(1e-3);
        let n = mo2.simule(0.1, 0.02, 0.9, 1e-12, 25, |_, _| {}).unwrap();
        assert!(mo2.n_rejeu > 0 && n > 5, "rejoués {} pas {n}", mo2.n_rejeu);
    }

    /// Un état de DÉPART où l'angle d'engrenage est à plus d'un quart de tour
    /// de `prev` : `detecte_redondance` et `acc_init` paniquaient sur leur
    /// `expect("Φ initial")` ; `simule` doit rendre l'erreur de l'élément.
    #[test]
    fn engrenage_prev_incoherent_rend_err() {
        let mut mo = Modele::new(V3::zeros());
        mo.corps.push(corps_test("a", V3::zeros(), V3::zeros()));
        mo.corps
            .push(corps_test("b", V3::new(0.1, 0.0, 0.0), V3::zeros()));
        mo.elems.push(Elem::E(Engrenage {
            nom: "e".into(),
            a: Some(0),
            b: Some(1),
            c: None,
            na: V3::z(),
            nb: V3::z(),
            rapport: -1.0,
            ref_a: M3::identity(),
            ref_b: M3::identity(),
            prev: (3.0, 0.0),
        }));
        let msg = mo
            .simule(0.01, 1e-3, 0.9, 1e-12, 25, |_, _| {})
            .expect_err("Err attendu");
        assert!(msg.contains("quart de tour"), "{msg}");
        // et les tangentes exposées (raideur → modes) rendent Err aussi
        assert!(mo.modes(0.0, 3).is_err());
    }

    /// Le LU creux à parallélisme EXPLICITE rend la même solution que le
    /// dense, en séquentiel (dim < CREUX_PAR) comme en parallèle (dim ≥
    /// CREUX_PAR), la symbolique gardée resert, et le réglage global de faer
    /// n'a pas bougé — c'est lui que l'ancien `creux` basculait (6 sept.).
    #[test]
    fn facto_creux_par_explicite() {
        let tri = |dim: usize| -> Triplets {
            let t = |r, c, v| faer::sparse::Triplet {
                row: r,
                col: c,
                val: v,
            };
            let mut trip = Vec::new();
            for i in 0..dim {
                trip.push(t(i, i, 4.0 + (i % 7) as f64));
                if i + 1 < dim {
                    trip.push(t(i, i + 1, -1.0));
                    trip.push(t(i + 1, i, -1.3));
                }
                if i + 17 < dim {
                    trip.push(t(i, i + 17, 0.5));
                }
            }
            trip
        };
        for dim in [200usize, CREUX_PAR + 50] {
            let trip = tri(dim);
            let x_ref = DVector::from_fn(dim, |i, _| 1.0 + (i as f64) * 1e-3);
            let b = trip_mv(&trip, dim, &x_ref);
            let mut symb = None;
            let f = Facto::creux(&trip, dim, &mut symb).unwrap();
            let x = f.solve(&b);
            let ecart = (&x - &x_ref).amax();
            assert!(ecart < 1e-9, "dim {dim} : écart {ecart}");
            let f2 = Facto::creux(&trip, dim, &mut symb).unwrap();
            assert!((f2.solve(&b) - &x_ref).amax() < 1e-9);
            if dim <= 400 {
                let (_, fd) = Facto::dense(&trip, dim);
                assert!((fd.solve(&b) - &x).amax() < 1e-9);
            }
            assert!(matches!(faer::get_global_parallelism(), faer::Par::Seq));
        }
    }

    /// L'appariement ne doit pas écraser les contacts déclarés à la main :
    /// jusqu'au 6 sept. `apparie` remplaçait la liste ENTIÈRE, et le sol
    /// disparaissait au premier pas.
    #[test]
    fn appariement_garde_les_contacts_manuels() {
        let mut mo = Modele::new(V3::new(0.0, 0.0, -9.81));
        mo.corps.push(corps_test("a", V3::zeros(), V3::zeros()));
        mo.corps
            .push(corps_test("b", V3::new(0.0, 0.0, 0.5), V3::zeros()));
        // sol manuel sous a : plan z = −0,1, sphère de a de rayon 0,1 posée dessus
        mo.contacts.push(Contact {
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
            origine: V3::new(0.0, 0.0, -0.1),
            k: 1e5,
            expo: 1.5,
            d_hat: 1e-3,
            c: 0.0,
            mu: 0.0,
            v_eps: 1e-3,
            demi: None,
            cylindre: None,
            maille: None,
            q_maille: V3::zeros(),
            s_axe: 0.0,
            cable: false,
            sortie: (0.0, 0.0),
            nonlisse: false,
            row: 0,
            fn_impose: None,
            restitution: 0.0,
            auto: false,
        });
        mo.spheres.push((0, V3::zeros(), 0.1));
        mo.spheres.push((1, V3::zeros(), 0.1));
        mo.loi_contact = Some((1e5, 1.5, 0.0, 0.0, 1e-3, 0.0, 1e-3));
        mo.simule(0.4, 1e-3, 0.9, 1e-12, 25, |_, _| {}).unwrap();
        assert_eq!(mo.contacts[0].nom, "sol");
        assert!(!mo.contacts[0].auto);
        assert!(mo.contacts.iter().skip(1).all(|c| c.auto));
        // a est restée sur le sol : le contact manuel a agi tout le long
        assert!(mo.corps[0].r.z > -0.05, "z_a = {}", mo.corps[0].r.z);
    }
}
