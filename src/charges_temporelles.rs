//! External time-dependent world-frame wrench at a body-fixed application point.
//! No user callback runs during Newton. Laws are immutable native values.

use crate::{fini, skew, Corps, Loi, M3, V3};

#[derive(Clone, Debug)]
pub struct EffortTemporel {
    pub nom: String,
    pub corps: usize,
    pub point: V3,
    pub force: [Loi; 3],
    pub moment: [Loi; 3],
}

impl EffortTemporel {
    pub fn verifie(&self, corps: &[Corps]) -> Result<(), String> {
        if self.corps >= corps.len() {
            return Err(format!("{} : corps de charge inexistant", self.nom));
        }
        fini(self.point.as_slice(), "point d'application")?;
        for loi in self.force.iter().chain(self.moment.iter()) {
            match loi {
                Loi::Constante(c) => fini(&[*c], "charge constante")?,
                Loi::Lineaire { a0, taux } => fini(&[*a0, *taux], "charge linéaire")?,
                Loi::Table(pts) => {
                    Loi::table(pts.clone())?;
                }
            }
        }
        Ok(())
    }

    pub fn valeur(&self, corps: &[Corps], t: f64) -> (V3, V3) {
        let f = V3::from_fn(|i, _| self.force[i].valeur(t));
        let m = V3::from_fn(|i, _| self.moment[i].valeur(t));
        let bras = corps[self.corps].rot * self.point;
        (f, m + bras.cross(&f))
    }

    /// ∂M/∂δθ for the spatial perturbation R' = exp(δθ) R.
    pub fn tangente_moment(&self, corps: &[Corps], t: f64) -> M3 {
        let (f, _) = self.valeur(corps, t);
        let bras = corps[self.corps].rot * self.point;
        skew(&f) * skew(&bras)
    }
}

pub fn loi_info(loi: &Loi) -> (String, Vec<f64>) {
    match loi {
        Loi::Constante(c) => ("constante".into(), vec![*c]),
        Loi::Lineaire { a0, taux } => ("lineaire".into(), vec![*a0, *taux]),
        Loi::Table(pts) => (
            "table".into(),
            pts.iter().flat_map(|&(t, y)| [t, y]).collect(),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expm;

    fn body() -> Corps {
        Corps {
            nom: "body".into(),
            m: 1.,
            j: M3::identity(),
            r: V3::zeros(),
            r_bas: V3::zeros(),
            rot: expm(&V3::new(0.2, -0.4, 0.1)),
            v: V3::zeros(),
            w: V3::zeros(),
        }
    }

    #[test]
    fn excentrage_tangente_et_puissance() {
        let bodies = vec![body()];
        let load = EffortTemporel {
            nom: "force".into(),
            corps: 0,
            point: V3::new(1., 2., -0.5),
            force: [
                Loi::Lineaire { a0: 1., taux: 2. },
                Loi::Constante(-3.),
                Loi::Constante(2.),
            ],
            moment: std::array::from_fn(|_| Loi::Constante(0.)),
        };
        load.verifie(&bodies).unwrap();
        let t = 0.7;
        let tangent = load.tangente_moment(&bodies, t);
        for axis in 0..3 {
            let mut plus = bodies.clone();
            let mut minus = bodies.clone();
            let delta = V3::from_fn(|i, _| if i == axis { 1e-6 } else { 0. });
            plus[0].rot = expm(&delta) * bodies[0].rot;
            minus[0].rot = expm(&(-delta)) * bodies[0].rot;
            let derivative = (load.valeur(&plus, t).1 - load.valeur(&minus, t).1) / 2e-6;
            assert!((derivative - tangent.column(axis)).norm() < 1e-8);
        }
        let v = V3::new(2., -1., 3.);
        let w = V3::new(0.4, -0.5, 0.1);
        let (f, m) = load.valeur(&bodies, t);
        let point_velocity = v + w.cross(&(bodies[0].rot * load.point));
        assert!((f.dot(&point_velocity) - f.dot(&v) - m.dot(&w)).abs() < 1e-12);
    }
}
