//! Diagnostics de conditionnement du prototype, sans verdict de performance.
//! Compare la flexibilité condensée à la solution linéaire de Timoshenko.

#[allow(dead_code)]
#[path = "../src/ad.rs"]
mod ad;
#[allow(dead_code)]
#[path = "poutre_mixte/element.rs"]
mod element;

use element::{Element, Pose};
use nalgebra::{DMatrix, Matrix3, Vector3};
type V3 = Vector3<f64>;
type M3 = Matrix3<f64>;

fn sonde(ordre: usize, ratio: f64, rigidite_flexion: f64) {
    let e = Element {
        ordre,
        l: 1.,
        cn: [10. * ratio, ratio, ratio],
        cm: [rigidite_flexion; 3],
        repere: M3::identity(),
    };
    let p = [
        Pose {
            r: V3::zeros(),
            rot: M3::identity(),
        },
        Pose {
            r: V3::new(1., 0., 0.),
            rot: M3::identity(),
        },
    ];
    print!("{{\"ordre\":{ordre},\"rigidite_cisaillement\":{ratio},\"rigidite_flexion\":{rigidite_flexion}");
    match e.condense(&p, &mut e.initial()) {
        Err(error) => print!(",\"condensation\":false,\"diagnostic\":{error:?}"),
        Ok(c) => {
            let kb = c.k.view((6, 6), (6, 6)).into_owned();
            let eigen_min = kb.clone().symmetric_eigen().eigenvalues.min();
            print!(",\"condensation\":true,\"valeur_propre_min\":{eigen_min}");
            if let Some(flexibility) = kb.try_inverse() {
                if rigidite_flexion > 0. {
                    let mut exact = DMatrix::<f64>::zeros(6, 6);
                    exact[(0, 0)] = 1. / e.cn[0];
                    exact[(3, 3)] = 1. / rigidite_flexion;
                    for (i, j, sign) in [(1, 5, 1.), (2, 4, -1.)] {
                        exact[(i, i)] = 1. / ratio + 1. / (3. * rigidite_flexion);
                        exact[(j, j)] = 1. / rigidite_flexion;
                        exact[(i, j)] = sign / (2. * rigidite_flexion);
                        exact[(j, i)] = exact[(i, j)];
                    }
                    let mut error: f64 = 0.;
                    for i in 0..6 {
                        for j in 0..6 {
                            error = error.max(
                                (flexibility[(i, j)] - exact[(i, j)]).abs()
                                    / (exact[(i, i)] * exact[(j, j)]).sqrt(),
                            );
                        }
                    }
                    print!(",\"erreur_flexibilite_equilibree\":{error}");
                }
            } else {
                print!(",\"diagnostic\":\"raideur condensée singulière\"");
            }
        }
    }
    println!("}}");
}

fn main() {
    std::thread::Builder::new()
        .stack_size(16 * 1024 * 1024)
        .spawn(|| {
            for ordre in [1, 2] {
                for ratio in [1., 1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14, 1e16] {
                    sonde(ordre, ratio, 1.);
                }
                // Au repos sans tension, le câble a des mécanismes physiques.
                // Le rejet ne démontre rien sur un câble tendu non linéaire.
                sonde(ordre, 1e5, 0.);
            }
        })
        .expect("thread de diagnostic")
        .join()
        .expect("diagnostic sans panique");
}
