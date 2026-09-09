//! Efforts corotationnels obtenus par travail virtuel : f = -Bᵀ K u,
//! B = ∂u/∂q. L'amortissement emploie la même cinématique : -β Bᵀ K B v.
use crate::ad::{self, Dual, Scalar, M, V};
use crate::{Corps, Superelement};
use nalgebra::{DMatrix, DVector};

impl Superelement {
    fn eval<T: Scalar>(
        &self,
        poses: &[(V<T>, M<T>)],
        bas: &[V<T>],
        vit: &[(V<T>, V<T>)],
        beta: f64,
    ) -> Vec<T> {
        let nb = poses.len();
        let (p0, r0) = poses[0];
        let (v0, w0) = vit[0];
        let mut u = vec![T::zero(); 6 * nb];
        let mut invs = vec![ad::m_ident(); nb];
        for i in 1..nb {
            let (p, r) = poses[i];
            let d = ad::v_sub(p, p0);
            let db = ad::v_sub(bas[i], bas[0]);
            let dp = ad::v_add(
                ad::v_sub(ad::m_tv(&r0, d), ad::v_from(self.ref_p[i].into())),
                ad::m_tv(&r0, db),
            );
            let rr = ad::m_mul(
                &ad::m_mul(&ad::m_t(&r0), &r),
                &ad::m_from(&self.ref_r[i].transpose()),
            );
            let th = crate::tangent::logm(&rr);
            let a2 = ad::v_dot(th, th);
            let coeff = if a2.val() < 1e-6 {
                T::from_f64(1.0 / 12.0) + a2.scale(1.0 / 720.0) + (a2 * a2).scale(1.0 / 30240.0)
            } else {
                let half = a2.sqrt().scale(0.5);
                (T::one() - half * half.cos() / half.sin()) / a2
            };
            let s = ad::skew(th);
            let s2 = ad::m_mul(&s, &s);
            for a in 0..3 {
                for b in 0..3 {
                    invs[i][a][b] = invs[i][a][b] - s[a][b].scale(0.5) + coeff * s2[a][b];
                }
            }
            let dv = ad::m_tv(
                &r0,
                ad::v_sub(
                    ad::v_sub(ad::v_sub(vit[i].0, v0), ad::v_cross(w0, d)),
                    ad::v_cross(w0, db),
                ),
            );
            let dw = ad::m_v(&invs[i], ad::m_tv(&r0, ad::v_sub(vit[i].1, w0)));
            for a in 0..3 {
                u[6 * i + a] = dp[a] + dv[a].scale(beta);
                u[6 * i + 3 + a] = th[a] + dw[a].scale(beta);
            }
        }
        let mut stress = vec![T::zero(); u.len()];
        for (i, si) in stress.iter_mut().enumerate() {
            for (j, uj) in u.iter().enumerate() {
                *si = *si + uj.scale(self.k[(i, j)]);
            }
        }
        let mut f = vec![T::zero(); u.len()];
        for i in 1..nb {
            let ft = ad::m_v(
                &r0,
                [-stress[6 * i], -stress[6 * i + 1], -stress[6 * i + 2]],
            );
            let mt = ad::m_v(
                &r0,
                ad::m_tv(
                    &invs[i],
                    [-stress[6 * i + 3], -stress[6 * i + 4], -stress[6 * i + 5]],
                ),
            );
            let reaction = ad::v_add(
                ad::v_add(mt, ad::v_cross(ad::v_sub(poses[i].0, p0), ft)),
                ad::v_cross(ad::v_sub(bas[i], bas[0]), ft),
            );
            for a in 0..3 {
                f[6 * i + a] = ft[a];
                f[6 * i + 3 + a] = mt[a];
                f[a] = f[a] - ft[a];
                f[3 + a] = f[3 + a] - reaction[a];
            }
        }
        f
    }

    pub fn forces(&self, corps: &[Corps]) -> DVector<f64> {
        self.efforts(corps, 0.0)
    }

    pub fn amortissement(&self, corps: &[Corps]) -> DVector<f64> {
        if self.beta == 0.0 {
            return DVector::zeros(6 * self.noeuds.len());
        }
        self.efforts(corps, self.beta) - self.efforts(corps, 0.0)
    }

    pub(crate) fn efforts(&self, corps: &[Corps], beta: f64) -> DVector<f64> {
        let (poses, bas) = self.poses_t::<f64>(corps);
        let vit: Vec<_> = self
            .noeuds
            .iter()
            .map(|&i| (corps[i].v.into(), corps[i].w.into()))
            .collect();
        DVector::from_vec(self.eval(&poses, &bas, &vit, beta))
    }

    fn poses_t<T: Scalar>(&self, corps: &[Corps]) -> (Vec<(V<T>, M<T>)>, Vec<V<T>>) {
        self.noeuds
            .iter()
            .map(|&i| {
                let (r, bas) = corps[i].difference_precise(&corps[self.noeuds[0]]);
                (
                    (ad::v_from(r.into()), ad::m_from(&corps[i].rot)),
                    ad::v_from(bas.into()),
                )
            })
            .unzip()
    }

    pub(crate) fn verifie_branche(&self, corps: &[Corps]) -> Result<(), String> {
        let r0 = corps[self.noeuds[0]].rot;
        for (i, &ni) in self.noeuds.iter().enumerate().skip(1) {
            let r = r0.transpose() * corps[ni].rot * self.ref_r[i].transpose();
            if !r.iter().all(|v| v.is_finite()) || (r.trace() + 1.0).abs() < 1e-12 {
                return Err(format!(
                    "superélément « {} » : rotation relative à π, tangente non différentiable",
                    self.nom
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn tangentes(
        &self,
        corps: &[Corps],
    ) -> Result<(DMatrix<f64>, DMatrix<f64>), String> {
        self.verifie_branche(corps)?;
        type D = Dual<f64, 12>;
        let dim = 6 * self.noeuds.len();
        let mut k = DMatrix::zeros(dim, dim);
        let mut c = k.clone();
        for (j, &nj) in self.noeuds.iter().enumerate() {
            let (mut poses, bas) = self.poses_t::<D>(corps);
            let mut vit: Vec<_> = self
                .noeuds
                .iter()
                .map(|&i| {
                    (
                        ad::v_from::<D>(corps[i].v.into()),
                        ad::v_from::<D>(corps[i].w.into()),
                    )
                })
                .collect();
            poses[j].0 = std::array::from_fn(|a| D::var(poses[j].0[a].val(), a));
            let dr = ad::expm(std::array::from_fn(|a| D::var(0.0, 3 + a)));
            poses[j].1 = ad::m_mul(&dr, &poses[j].1);
            vit[j].0 = std::array::from_fn(|a| D::var(corps[nj].v[a], 6 + a));
            vit[j].1 = std::array::from_fn(|a| D::var(corps[nj].w[a], 9 + a));
            for (i, fi) in self.eval(&poses, &bas, &vit, self.beta).iter().enumerate() {
                for a in 0..6 {
                    k[(i, 6 * j + a)] = fi.d[a];
                    c[(i, 6 * j + a)] = fi.d[6 + a];
                }
            }
        }
        Ok((k, c))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{expm, log_so3, M3, V3};

    fn modele() -> (Superelement, Vec<Corps>) {
        let cs: Vec<_> = (0..3)
            .map(|i| Corps {
                r_bas: crate::V3::zeros(),
                nom: format!("{i}"),
                m: 1.0,
                j: M3::identity(),
                r: V3::new(i as f64, 0.2 * i as f64, 0.0),
                rot: expm(&V3::new(0.1 * i as f64, 0.0, 0.2)),
                v: V3::new(0.1, 0.4 * i as f64, -0.2),
                w: V3::new(0.1 * i as f64, 0.2, 0.3),
            })
            .collect();
        let rt = cs[0].rot.transpose();
        let a = DMatrix::from_fn(18, 18, |i, j| ((i * 7 + j * 3 + 1) as f64).sin());
        let se = Superelement {
            nom: "test".into(),
            noeuds: vec![0, 1, 2],
            k: a.transpose() * a + DMatrix::identity(18, 18),
            alpha: 0.0,
            beta: 0.03,
            ref_p: cs.iter().map(|c| rt * (c.r - cs[0].r)).collect(),
            ref_r: cs.iter().map(|c| rt * c.rot).collect(),
        };
        (se, cs)
    }

    fn energie(se: &Superelement, cs: &[Corps]) -> f64 {
        let mut u = DVector::zeros(18);
        let rt = cs[0].rot.transpose();
        for i in 1..3 {
            u.rows_mut(6 * i, 3)
                .copy_from(&(rt * (cs[i].r - cs[0].r) - se.ref_p[i]));
            u.rows_mut(6 * i + 3, 3)
                .copy_from(&log_so3(&(rt * cs[i].rot * se.ref_r[i].transpose())));
        }
        0.5 * u.dot(&(&se.k * &u))
    }

    fn perturbe(cs: &mut [Corps], j: usize, eps: f64) {
        let i = j / 6;
        let a = j % 6;
        if a < 3 {
            cs[i].r[a] += eps;
        } else {
            let mut v = V3::zeros();
            v[a - 3] = eps;
            cs[i].rot = expm(&v) * cs[i].rot;
        }
    }

    #[test]
    fn travail_virtuel_resultantes_et_tangentes() {
        let (se, mut cs) = modele();
        cs[1].r += V3::new(0.1, 0.02, -0.04);
        cs[2].rot = expm(&V3::new(-0.1, 0.2, 0.3)) * cs[2].rot;
        let f = se.forces(&cs);
        let (k, c) = se.tangentes(&cs).unwrap();
        let mut force = V3::zeros();
        let mut moment = V3::zeros();
        for i in 0..3 {
            let fi = V3::from_column_slice(f.rows(6 * i, 3).as_slice());
            force += fi;
            moment += cs[i].r.cross(&fi) + V3::from_column_slice(f.rows(6 * i + 3, 3).as_slice());
        }
        assert!(force.norm() < 1e-12 && moment.norm() < 1e-12);
        for eps in [1e-5, 1e-6] {
            for j in 0..18 {
                let mut p = cs.clone();
                let mut m = cs.clone();
                perturbe(&mut p, j, eps);
                perturbe(&mut m, j, -eps);
                let df = (se.efforts(&p, se.beta) - se.efforts(&m, se.beta)) / (2.0 * eps);
                assert!((df - k.column(j)).norm() < 1e-6 * (1.0 + k.column(j).norm()));
                assert!(
                    ((energie(&se, &p) - energie(&se, &m)) / (2.0 * eps) + f[j]).abs()
                        < 1e-7 * (1.0 + f[j].abs())
                );
                let mut p = cs.clone();
                let mut m = cs.clone();
                if j % 6 < 3 {
                    p[j / 6].v[j % 6] += eps;
                    m[j / 6].v[j % 6] -= eps;
                } else {
                    p[j / 6].w[j % 6 - 3] += eps;
                    m[j / 6].w[j % 6 - 3] -= eps;
                }
                let df = (se.efforts(&p, se.beta) - se.efforts(&m, se.beta)) / (2.0 * eps);
                assert!((df - c.column(j)).norm() < 1e-6 * (1.0 + c.column(j).norm()));
            }
        }
        let damp = se.amortissement(&cs);
        let power: f64 = cs
            .iter()
            .enumerate()
            .map(|(i, c)| {
                (0..3)
                    .map(|a| damp[6 * i + a] * c.v[a] + damp[6 * i + 3 + a] * c.w[a])
                    .sum::<f64>()
            })
            .sum();
        assert!(power < 0.0);
    }
}
