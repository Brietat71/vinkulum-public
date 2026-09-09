use crate::{
    element::{exp, Element, Pose},
    systeme::Systeme,
    M3, V3,
};
use nalgebra::{DMatrix, DVector};

fn pile(f: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024)
        .spawn(f)
        .unwrap()
        .join()
        .unwrap();
}

fn element(ordre: usize) -> (Element, [Pose; 2]) {
    (
        Element {
            ordre,
            l: 1.3,
            cn: [1e6, 8e4, 2e5],
            cm: [70., 110., 35.],
            repere: M3::identity(),
        },
        [
            Pose {
                r: V3::zeros(),
                rot: M3::identity(),
            },
            Pose {
                r: V3::new(1.3, 0., 0.),
                rot: M3::identity(),
            },
        ],
    )
}

#[test]
fn flexibilite_timoshenko_anisotrope_et_modes_rigides() {
    pile(|| {
        for ordre in [1, 2] {
            let (e, p) = element(ordre);
            let c = e.condense(&p, &mut e.initial()).unwrap();
            let kb = c.k.view((6, 6), (6, 6)).into_owned();
            let got = kb.try_inverse().unwrap();
            let mut exact = DMatrix::zeros(6, 6);
            let l = e.l;
            exact[(0, 0)] = l / e.cn[0];
            exact[(3, 3)] = l / e.cm[0];
            for (translation, rotation, rigidite, signe) in
                [(1, 5, e.cm[2], 1.), (2, 4, e.cm[1], -1.)]
            {
                exact[(translation, translation)] =
                    l / e.cn[translation] + l.powi(3) / (3. * rigidite);
                exact[(rotation, rotation)] = l / rigidite;
                exact[(translation, rotation)] = signe * l * l / (2. * rigidite);
                exact[(rotation, translation)] = exact[(translation, rotation)];
            }
            assert!(
                (&got - &exact).amax() < 2e-11,
                "ordre {ordre}: {}",
                (&got - &exact).amax()
            );
            for axis in 0..3 {
                let mut translation = DVector::zeros(12);
                translation[axis] = 1.;
                translation[6 + axis] = 1.;
                assert!((&c.k * translation).amax() < 1e-8);
                let mut rotation = DVector::zeros(12);
                rotation[3 + axis] = 1.;
                rotation[9 + axis] = 1.;
                let mut w = V3::zeros();
                w[axis] = 1.;
                rotation.rows_mut(6, 3).copy_from(&w.cross(&p[1].r));
                assert!((&c.k * rotation).amax() < 1e-8);
            }
        }
    });
}

#[test]
fn objectivite_apres_condensation_non_lineaire() {
    pile(|| {
        for ordre in [1, 2] {
            let (e, mut p) = element(ordre);
            p[1].r += V3::new(0.002, 0.04, -0.01);
            p[1].rot = exp(V3::new(0.03, 0.02, 0.05));
            let mut y = e.initial();
            let c = e.condense(&p, &mut y).unwrap();
            let q = exp(V3::new(0.5, -0.7, 0.2));
            let t = V3::new(2., -3., 1.);
            let pp = p.map(|p| Pose {
                r: q * p.r + t,
                rot: q * p.rot,
            });
            y.rot = q * y.rot;
            y.bulle = q * y.bulle;
            let other = e.condense(&pp, &mut y).unwrap();
            assert!((other.u - c.u).abs() < 2e-8);
            for off in [0, 3, 6, 9] {
                let expected = q * c.g.rows(off, 3).into_owned();
                assert!((other.g.rows(off, 3) - expected).amax() < 2e-7);
            }
        }
    });
}

#[test]
fn tangente_condensee_contre_difference_du_gradient_spatial() {
    pile(|| {
        for ordre in [1, 2] {
            let (e, mut p) = element(ordre);
            p[1].r += V3::new(0.001, 0.02, -0.01);
            p[1].rot = exp(V3::new(0.03, 0.01, 0.04));
            let mut y = e.initial();
            let c = e.condense(&p, &mut y).unwrap();
            let mut errors = Vec::new();
            for h in [1e-4, 2e-5, 4e-6, 8e-7] {
                let mut fd = DMatrix::zeros(12, 12);
                for col in 0..12 {
                    let mut gradients = Vec::new();
                    for sign in [-1., 1.] {
                        let mut pp = p.clone();
                        let node = col / 6;
                        let component = col % 6;
                        if component < 3 {
                            pp[node].r[component] += sign * h;
                        } else {
                            let mut w = V3::zeros();
                            w[component - 3] = sign * h;
                            pp[node].rot = exp(w) * pp[node].rot;
                        }
                        gradients.push(e.condense(&pp, &mut y.clone()).unwrap().g);
                    }
                    fd.set_column(col, &((&gradients[1] - &gradients[0]) / (2. * h)));
                }
                errors.push((&fd - &c.k).norm() / c.k.norm());
            }
            println!("tangente ordre {ordre}, h=[1e-4,2e-5,4e-6,8e-7] : {errors:?}");
            assert!(
                errors.iter().copied().fold(f64::INFINITY, f64::min) < 2e-7,
                "ordre {ordre}: {errors:?}"
            );
            assert!(
                errors[1] < errors[0] / 10.,
                "absence de convergence des différences finies : {errors:?}"
            );
        }
    });
}

#[test]
fn arc_sous_moment_pur_converge_en_ordre_quatre() {
    pile(|| {
        let exact = V3::new(1_f64.sin(), 1. - 1_f64.cos(), 0.);
        let mut errors = Vec::new();
        for ne in [1, 2, 4] {
            let mut s = Systeme::console(ne, 2, 1., [1e6, 1e5, 1e5], [50., 100., 100.]);
            for step in 1..=10 {
                s.equilibre(V3::zeros(), V3::new(0., 0., 10. * step as f64), 1e-8)
                    .unwrap();
            }
            let p = &s.poses[ne];
            let error = (p.r - exact).norm();
            errors.push(error);
            assert!((p.rot[(1, 0)].atan2(p.rot[(0, 0)]) - 1.).abs() < 2e-10);
        }
        assert!(
            errors[0] / errors[1] > 15. && errors[1] / errors[2] > 15.,
            "{errors:?}"
        );
    });
}

#[test]
fn jonction_ramifiee_equilibre_les_moments() {
    pile(|| {
        let positions = [
            V3::zeros(),
            V3::new(1., 0., 0.),
            V3::new(2., 0.5, 0.),
            V3::new(2., -0.5, 0.),
        ];
        let poses = positions
            .into_iter()
            .map(|r| Pose {
                r,
                rot: M3::identity(),
            })
            .collect();
        let elements: Vec<_> = [(0, 1), (1, 2), (1, 3)]
            .into_iter()
            .map(|(a, b)| {
                let tangent = positions[b] - positions[a];
                let l = tangent.norm();
                let x = tangent / l;
                let y = V3::new(0., 0., 1.);
                (
                    a,
                    b,
                    Element {
                        ordre: 2,
                        l,
                        cn: [1e6, 1e5, 1e5],
                        cm: [50., 100., 100.],
                        repere: M3::from_columns(&[x, y, x.cross(&y)]),
                    },
                )
            })
            .collect();
        let internes = elements.iter().map(|(_, _, e)| e.initial()).collect();
        let mut s = Systeme {
            poses,
            elements,
            internes,
            evaluations: 0,
            iterations_locales: 0,
        };
        let mut load = DVector::zeros(18);
        for step in 1..=10 {
            load[8] = step as f64 / 10.;
            load[14] = -step as f64 / 10.;
            s.equilibre_charges(&load, 1e-8).unwrap();
        }
        let (_, g, _, internal) = s.assemble().unwrap();
        assert!(internal < 1e-9);
        assert!(g.rows(6, 6).amax() < 1e-8, "résidu au nœud de jonction");
        assert!(
            g.rows(0, 3).amax() < 1e-8,
            "la paire de forces a une résultante nulle"
        );
        let torque = (s.poses[2].r - s.poses[3].r).cross(&V3::new(0., 0., 1.));
        assert!(
            (g.rows(3, 3) + torque).amax() < 1e-8,
            "bilan spatial des moments"
        );
    });
}
