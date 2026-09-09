use crate::{
    element::{exp, Element, Interne, Pose},
    M3, V3,
};
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

pub struct Systeme {
    pub poses: Vec<Pose>,
    pub elements: Vec<(usize, usize, Element)>,
    pub internes: Vec<Interne>,
    pub evaluations: usize,
    pub iterations_locales: usize,
}

impl Systeme {
    pub fn console(ne: usize, ordre: usize, l: f64, cn: [f64; 3], cm: [f64; 3]) -> Self {
        let poses = (0..=ne)
            .map(|i| Pose {
                r: V3::new(l * i as f64 / ne as f64, 0., 0.),
                rot: M3::identity(),
            })
            .collect();
        let elements: Vec<_> = (0..ne)
            .map(|i| {
                (
                    i,
                    i + 1,
                    Element {
                        ordre,
                        l: l / ne as f64,
                        cn,
                        cm,
                        repere: M3::identity(),
                    },
                )
            })
            .collect();
        let internes = elements.iter().map(|(_, _, e)| e.initial()).collect();
        Self {
            poses,
            elements,
            internes,
            evaluations: 0,
            iterations_locales: 0,
        }
    }

    pub fn assemble(&mut self) -> Result<(f64, DVector<f64>, DMatrix<f64>, f64), String> {
        self.evaluations += 1;
        let n = 6 * self.poses.len();
        let mut g = DVector::zeros(n);
        let mut k = DMatrix::zeros(n, n);
        let mut u = 0.;
        let mut ri: f64 = 0.;
        for (idx, (a, b, e)) in self.elements.iter().enumerate() {
            let c = e.condense(
                &[self.poses[*a].clone(), self.poses[*b].clone()],
                &mut self.internes[idx],
            )?;
            self.iterations_locales += c.iterations;
            u += c.u;
            ri = ri.max(c.residu_interne);
            let map: Vec<_> = (0..6)
                .map(|j| 6 * a + j)
                .chain((0..6).map(|j| 6 * b + j))
                .collect();
            for i in 0..12 {
                g[map[i]] += c.g[i];
                for j in 0..12 {
                    k[(map[i], map[j])] += c.k[(i, j)];
                }
            }
        }
        Ok((u, g, k, ri))
    }

    /// Prototype : racine de l'équilibre avec bâti encastré et charge au bout.
    /// Le mérite emploie la correction prédite ; les résidus physiques restent
    /// le critère final. Le code de production Vinkulum n'est pas modifié.
    pub fn equilibre(&mut self, force: V3, moment: V3, tol: f64) -> Result<f64, String> {
        let n = 6 * (self.poses.len() - 1);
        let mut load = DVector::zeros(n);
        load.rows_mut(n - 6, 3).copy_from(&force);
        load.rows_mut(n - 3, 3).copy_from(&moment);
        self.equilibre_charges(&load, tol)
    }

    pub fn equilibre_charges(&mut self, load: &DVector<f64>, tol: f64) -> Result<f64, String> {
        let n = 6 * (self.poses.len() - 1);
        if load.len() != n {
            return Err("dimensions des charges incompatibles".into());
        }
        for it in 0..60 {
            let (_, g, k, _) = self.assemble()?;
            let residual = g.rows(6, n).into_owned() - load;
            if residual.amax() < tol {
                return Ok(residual.amax());
            }
            let matrix = k.view((6, 6), (n, n)).into_owned();
            let scale = matrix.diagonal().map(|v| v.abs().max(1e-30).sqrt().recip());
            let ks = DMatrix::from_fn(n, n, |i, j| scale[i] * matrix[(i, j)] * scale[j]);
            let lu = ks.lu();
            let dx = lu
                .solve(&(-residual.component_mul(&scale)))
                .ok_or("équilibre global singulier")?
                .component_mul(&scale);
            let base_poses = self.poses.clone();
            let base_y = self.internes.clone();
            let mut alpha = 1.;
            let mut accepted = false;
            for _ in 0..18 {
                self.poses = base_poses.clone();
                self.internes = base_y.clone();
                for i in 1..self.poses.len() {
                    let off = 6 * (i - 1);
                    self.poses[i].r += alpha * V3::new(dx[off], dx[off + 1], dx[off + 2]);
                    self.poses[i].rot = exp(alpha * V3::new(dx[off + 3], dx[off + 4], dx[off + 5]))
                        * self.poses[i].rot;
                }
                if let Ok((_, g2, _, _)) = self.assemble() {
                    let res2 = g2.rows(6, n).into_owned() - load;
                    let corr = lu
                        .solve(&res2.component_mul(&scale))
                        .ok_or("correction globale singulière")?
                        .component_mul(&scale);
                    if corr.norm() < (1. - 1e-4 * alpha) * dx.norm() || res2.amax() < tol {
                        accepted = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if !accepted {
                self.poses = base_poses;
                self.internes = base_y;
                return Err(format!(
                    "recherche globale refusée à {it}, résidu {}",
                    residual.amax()
                ));
            }
        }
        Err("budget global épuisé".into())
    }
}

pub fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    let ordre = args
        .get(1)
        .map_or(Ok(2), |s| s.parse::<usize>())
        .map_err(|e| e.to_string())?;
    let ne = args
        .get(2)
        .map_or(Ok(4), |s| s.parse::<usize>())
        .map_err(|e| e.to_string())?;
    let cas = args.get(3).map_or("princeton", String::as_str);
    let paliers = args
        .get(4)
        .map_or(Ok(50), |s| s.parse::<usize>())
        .map_err(|e| e.to_string())?;
    let profil = args.get(5).is_some_and(|s| s == "profil");
    if ![1, 2].contains(&ordre) || ne == 0 || paliers == 0 {
        return Err("ordre 1 ou 2, éléments et paliers positifs requis".into());
    }
    let princeton = cas == "princeton";
    let mut sys = if princeton {
        Systeme::console(
            ne,
            ordre,
            0.508,
            [2.84191e6, 6.40131e5, 9.03881e5],
            [3.10338, 36.2794, 2.42873],
        )
    } else {
        Systeme::console(ne, ordre, 1., [1e6, 1e5, 1e5], [50., 100., 100.])
    };
    let start = Instant::now();
    let mut residual = 0.;
    let mut max_residual: f64 = 0.;
    for step in 1..=paliers {
        let load = 0.5 * (1. - (std::f64::consts::PI * step as f64 / paliers as f64).cos());
        let (f, m) = match cas {
            "princeton" => (
                V3::new(0., 8.896 * load / 2_f64.sqrt(), 8.896 * load / 2_f64.sqrt()),
                V3::zeros(),
            ),
            "moment" => (V3::zeros(), V3::new(0., 0., 100. * load)),
            "force" => (V3::new(0., load, 0.), V3::zeros()),
            _ => return Err("cas princeton, moment ou force requis".into()),
        };
        residual = sys.equilibre(f, m, 1e-8)?;
        max_residual = max_residual.max(residual);
    }
    let elapsed = start.elapsed().as_secs_f64();
    let (_, gradient, _, internal_residual) = sys.assemble()?;
    let p = &sys.poses[ne];
    print!("{{\"ordre\":{ordre},\"elements\":{ne},\"cas\":\"{cas}\",\"paliers\":{paliers},\"temps_s\":{elapsed},\"residu\":{residual},\"residu_max_paliers\":{max_residual},\"residu_interne\":{internal_residual},\"evaluations\":{},\"iterations_locales\":{},\"position\":[{},{},{}],\"angle_z\":{},\"reaction_racine\":[{},{},{},{},{},{}]",
        sys.evaluations,sys.iterations_locales,p.r.x,p.r.y,p.r.z,p.rot[(1,0)].atan2(p.rot[(0,0)]),gradient[0],gradient[1],gradient[2],gradient[3],gradient[4],gradient[5]);
    if profil {
        print!(",\"positions\":[");
        for i in 0..=960 {
            let location = ne as f64 * i as f64 / 960.;
            let index = (location.floor() as usize).min(ne - 1);
            let xi = 2. * (location - index as f64) - 1.;
            let r = (1. - xi) * 0.5 * sys.poses[index].r
                + (1. + xi) * 0.5 * sys.poses[index + 1].r
                + (1. - xi * xi) * sys.elements[index].2.l * sys.internes[index].bulle;
            if i != 0 {
                print!(",");
            }
            print!("[{},{},{}]", r.x, r.y, r.z);
        }
        print!("]");
    }
    println!("}}");
    Ok(())
}
