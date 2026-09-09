use crate::{ad, M3, V3};
use ad::{Dual, Scalar, M, V};
use nalgebra::{DMatrix, DVector};

#[derive(Clone, Debug)]
pub struct Pose {
    pub r: V3,
    pub rot: M3,
}

#[derive(Clone, Debug)]
pub struct Interne {
    pub rot: M3,
    pub pente: V3,
    /// r_milieu = (r_A+r_B)/2 + L*bulle ; coordonnées monde.
    pub bulle: V3,
}

#[derive(Clone, Debug)]
pub struct Element {
    pub ordre: usize,
    pub l: f64,
    pub cn: [f64; 3],
    pub cm: [f64; 3],
    pub repere: M3,
}

pub struct Jet {
    pub u: f64,
    pub g: DVector<f64>,
    pub h: DMatrix<f64>,
}

pub struct Condense {
    pub u: f64,
    pub g: DVector<f64>,
    pub k: DMatrix<f64>,
    pub residu_interne: f64,
    pub iterations: usize,
}

pub fn exp(v: V3) -> M3 {
    let m = ad::expm::<f64>(v.into());
    M3::from_fn(|i, j| m[i][j])
}

/// Branche principale, refus explicite au voisinage de la coupure de pi.
/// Série asin(s)/s en zéro, sans racine singulière pour les jets emboîtés.
fn log<T: Scalar>(r: &M<T>) -> V<T> {
    let w = [
        (r[2][1] - r[1][2]).scale(0.5),
        (r[0][2] - r[2][0]).scale(0.5),
        (r[1][0] - r[0][1]).scale(0.5),
    ];
    let c = (ad::m_trace(r) - T::one()).scale(0.5);
    if c.val() < -0.999999 {
        return [T::from_f64(f64::NAN); 3];
    }
    let s2 = ad::v_dot(w, w);
    let factor = if c.val() > 0.0 && s2.val() < 1e-6 {
        T::one()
            + s2.scale(1.0 / 6.0)
            + (s2 * s2).scale(3.0 / 40.0)
            + (s2 * s2 * s2).scale(5.0 / 112.0)
    } else {
        let s = s2.sqrt();
        s.atan2(c) / s
    };
    ad::v_scale(factor, w)
}

fn rotation<T: Scalar>(rot: &M3, x: [T; 3]) -> M<T> {
    ad::m_mul_cst(&ad::expm(x), rot)
}

impl Interne {
    pub fn incremente(&mut self, x: &DVector<f64>, alpha: f64) {
        self.rot = exp(alpha * V3::new(x[0], x[1], x[2])) * self.rot;
        if x.len() == 9 {
            self.pente += alpha * V3::new(x[3], x[4], x[5]);
            self.bulle += alpha * V3::new(x[6], x[7], x[8]);
        }
    }
}

impl Element {
    pub fn initial(&self) -> Interne {
        Interne {
            rot: self.repere,
            pente: V3::zeros(),
            bulle: V3::zeros(),
        }
    }

    /// Les moments sont éliminés analytiquement par leur matrice de masse
    /// polynomiale H. Energie de flexion = (1/2) b^T H^-1 C_m b.
    /// b comprend la courbure volumique et les sauts aux deux extrémités.
    fn energie<T: Scalar>(
        &self,
        poses: &[Pose; 2],
        interne: &Interne,
        seed: impl Fn(usize) -> T,
    ) -> T {
        let ra = rotation(&(poses[0].rot * self.repere), [seed(3), seed(4), seed(5)]);
        let rb = rotation(&(poses[1].rot * self.repere), [seed(9), seed(10), seed(11)]);
        let rm = rotation(&interne.rot, [seed(12), seed(13), seed(14)]);
        let a: V<T> = if self.ordre == 2 {
            std::array::from_fn(|i| T::from_f64(interne.pente[i]) + seed(15 + i))
        } else {
            [T::zero(); 3]
        };
        let bulle: V<T> = if self.ordre == 2 {
            std::array::from_fn(|i| T::from_f64(interne.bulle[i]) + seed(18 + i))
        } else {
            [T::zero(); 3]
        };
        let corde: V<T> = std::array::from_fn(|i| {
            (T::from_f64(poses[1].r[i] - poses[0].r[i]) + seed(6 + i) - seed(i)).scale(1.0 / self.l)
        });
        let mut u = T::zero();
        // Intégration réduite à k points de Gauss de la partie extension/cisaillement.
        let points: &[f64] = if self.ordre == 1 {
            &[0.0]
        } else {
            &[-0.5773502691896258, 0.5773502691896258]
        };
        for &xi in points {
            let re = ad::m_mul(&rm, &ad::expm(ad::v_scale(T::from_f64(xi), a)));
            let dr = ad::v_sub(corde, ad::v_scale(T::from_f64(4.0 * xi), bulle));
            let mut gamma = ad::m_tv(&re, dr);
            gamma[0] = gamma[0] - T::one();
            for (i, gi) in gamma.into_iter().enumerate() {
                u = u + (gi * gi).scale(0.5 * self.l * self.cn[i] / points.len() as f64);
            }
        }
        let re_a = ad::m_mul(&rm, &ad::expm(ad::v_scale(T::from_f64(-1.0), a)));
        let re_b = ad::m_mul(&rm, &ad::expm(a));
        let ta = log(&ad::m_mul(&ad::m_t(&re_a), &ra));
        let tb = log(&ad::m_mul(&ad::m_t(&re_b), &rb));
        let (b, hinv): (Vec<V<T>>, Vec<Vec<f64>>) = if self.ordre == 1 {
            (
                vec![ad::v_scale(T::from_f64(-1.0), ta), tb],
                vec![vec![4., -2.], vec![-2., 4.]],
            )
        } else {
            (
                vec![
                    ad::v_sub(ad::v_scale(T::from_f64(1. / 3.), a), ta),
                    ad::v_scale(T::from_f64(4. / 3.), a),
                    ad::v_add(ad::v_scale(T::from_f64(1. / 3.), a), tb),
                ],
                vec![
                    vec![9., -1.5, 3.],
                    vec![-1.5, 2.25, -1.5],
                    vec![3., -1.5, 9.],
                ],
            )
        };
        for i in 0..b.len() {
            for j in 0..b.len() {
                for (k, coefficient) in self.cm.iter().enumerate() {
                    u = u + (b[i][k] * b[j][k]).scale(0.5 * hinv[i][j] * coefficient / self.l);
                }
            }
        }
        u
    }

    fn jet<const N: usize>(&self, p: &[Pose; 2], y: &Interne, local: bool) -> Jet {
        type D<const N: usize> = Dual<Dual<f64, N>, N>;
        let u = self.energie(p, y, |i| {
            if local && i < 12 {
                D::<N>::from_f64(0.)
            } else {
                let j = if local { i - 12 } else { i };
                D::<N>::var(Dual::var(0., j), j)
            }
        });
        Jet {
            u: u.v.v,
            g: DVector::from_fn(N, |i, _| u.v.d[i]),
            h: DMatrix::from_fn(N, N, |i, j| u.d[i].d[j]),
        }
    }

    pub fn complet(&self, p: &[Pose; 2], y: &Interne) -> Jet {
        if self.ordre == 1 {
            self.jet::<15>(p, y, false)
        } else {
            self.jet::<21>(p, y, false)
        }
    }

    fn gradient<const N: usize>(&self, p: &[Pose; 2], y: &Interne) -> DVector<f64> {
        let u = self.energie(p, y, |i| {
            if i < 12 {
                Dual::<f64, N>::cst(0.)
            } else {
                Dual::<f64, N>::var(0., i - 12)
            }
        });
        DVector::from_column_slice(&u.d)
    }

    fn residu(&self, g: &DVector<f64>) -> f64 {
        g.iter()
            .enumerate()
            .map(|(i, v)| v.abs() / if i >= 6 { self.l } else { 1. })
            .fold(0., f64::max)
    }

    pub fn condense(&self, p: &[Pose; 2], y: &mut Interne) -> Result<Condense, String> {
        let mut iterations = 0;
        let mut converged = false;
        for it in 0..40 {
            iterations = it + 1;
            let j = if self.ordre == 1 {
                self.jet::<3>(p, y, true)
            } else {
                self.jet::<9>(p, y, true)
            };
            if !j.u.is_finite() || !j.g.iter().chain(j.h.iter()).all(|x| x.is_finite()) {
                return Err("poutre mixte : branche de rotation ou état non fini".into());
            }
            let scale = j.h.diagonal().map(|x| x.abs().max(1e-30).sqrt().recip());
            let hs = DMatrix::from_fn(j.g.len(), j.g.len(), |a, b| {
                scale[a] * j.h[(a, b)] * scale[b]
            });
            let gs = j.g.component_mul(&scale);
            // Résolution locale dans des coordonnées équilibrées, avec contrôle
            // sur le gradient original et sur l'incrément prédit.
            let eig = hs.clone().symmetric_eigen().eigenvalues.min();
            let mut regularized = hs;
            if eig < 1e-8 {
                regularized += DMatrix::identity(j.g.len(), j.g.len()) * (1e-8 - eig);
            }
            let lu = regularized.lu();
            let dx = lu
                .solve(&(-gs))
                .ok_or("poutre mixte : résolution locale singulière")?
                .component_mul(&scale);
            if dx.amax() < 1e-13 && self.residu(&j.g) < 1e-9 && eig > 0. {
                converged = true;
                break;
            }
            let slope = j.g.dot(&dx);
            let mut accepted = None;
            let mut alpha = 1.0;
            for _ in 0..20 {
                let mut next = y.clone();
                next.incremente(&dx, alpha);
                let u = self.energie::<f64>(p, &next, |_| 0.);
                if u.is_finite() {
                    let gn = if self.ordre == 1 {
                        self.gradient::<3>(p, &next)
                    } else {
                        self.gradient::<9>(p, &next)
                    };
                    let correction = lu
                        .solve(&gn.component_mul(&scale))
                        .ok_or("correction locale singulière")?
                        .component_mul(&scale);
                    // La différence d'énergies devient mal résolue près de
                    // l'équilibre d'une poutre raide. Contrôle complémentaire
                    // dans la métrique de la factorisation locale courante.
                    if self.residu(&gn) < 1e-9
                        || correction.norm() < (1. - 1e-4 * alpha) * dx.norm()
                        || u <= j.u + 1e-4 * alpha * slope
                            && slope.abs() > 1e-12 * j.u.abs().max(1.)
                    {
                        accepted = Some(next);
                        break;
                    }
                }
                alpha *= 0.5;
            }
            *y = accepted.ok_or_else(|| {
                format!(
                    "poutre mixte : recherche locale refusée à {it}, gradient {}",
                    j.g.amax()
                )
            })?;
        }
        let j = self.complet(p, y);
        let ni = j.g.len() - 12;
        let gy = j.g.rows(12, ni).into_owned();
        if !converged {
            return Err(format!(
                "poutre mixte : budget local épuisé, gradient {}",
                gy.amax()
            ));
        }
        let hyy = j.h.view((12, 12), (ni, ni)).into_owned();
        let hyx = j.h.view((12, 0), (ni, 12)).into_owned();
        let solved = hyy
            .lu()
            .solve(&hyx)
            .ok_or("poutre mixte : condensation singulière")?;
        let g = j.g.rows(0, 12).into_owned();
        let mut k = j.h.view((0, 0), (12, 12)).into_owned() - hyx.transpose() * solved;
        // Hessienne de la rétraction -> dérivée du gradient spatial.
        // Pour les rotations gauches : Dg = Hess(E)-[g]x/2.
        for off in [3, 9] {
            let skew = V3::new(g[off], g[off + 1], g[off + 2]).cross_matrix();
            for a in 0..3 {
                for b in 0..3 {
                    k[(off + a, off + b)] -= 0.5 * skew[(a, b)];
                }
            }
        }
        Ok(Condense {
            u: j.u,
            g,
            k,
            residu_interne: self.residu(&gy),
            iterations,
        })
    }
}
