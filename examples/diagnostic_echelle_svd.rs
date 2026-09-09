//! Spectre connu sous rotations denses et changements d'échelle.
//! Compare le moteur faer brut et l'adaptateur normalisé du noyau.
use nalgebra::{DMatrix, Matrix3, Vector3};

fn main() {
    let q = vinkulum::expm(&Vector3::new(0.6, -0.7, 1.8));
    let p = vinkulum::expm(&Vector3::new(-1.2, 0.4, 0.8));
    print!("[");
    for (i, scale) in [1e-300, 1e-200, 1., 1e200, 1e300].into_iter().enumerate() {
        let a = q
            * Matrix3::from_diagonal(&Vector3::new(2. * scale, scale, -0.5 * scale))
            * p.transpose();
        let brute = match faer::Mat::from_fn(3, 3, |i, j| a[(i, j)]).thin_svd() {
            Ok(s) => {
                let sigma: Vec<_> = (0..3).map(|i| s.S()[i]).collect();
                let u: Vec<_> = (0..9).map(|i| s.U()[(i % 3, i / 3)]).collect();
                let vt: Vec<_> = (0..9).map(|i| s.V()[(i / 3, i % 3)]).collect();
                if sigma.iter().chain(&u).chain(&vt).all(|x| x.is_finite()) {
                    format!(
                        "{{\"erreur\":null,\"sigma\":{sigma:?},\"u_col\":{u:?},\"vt_col\":{vt:?}}}"
                    )
                } else {
                    "{\"erreur\":\"facteurs non finis\"}".into()
                }
            }
            Err(e) => format!("{{\"erreur\":{:?}}}", format!("{e:?}")),
        };
        let courante = vinkulum::svd_sure(DMatrix::from_column_slice(3, 3, a.as_slice())).unwrap();
        println!("{}{{\"echelle\":{scale:?},\"a_col\":{:?},\"brute\":{brute},\"sigma\":{:?},\"u_col\":{:?},\"vt_col\":{:?}}}",
            if i == 0 { "" } else { "," }, a.as_slice(), courante.singular_values.as_slice(),
            courante.u.unwrap().as_slice(), courante.v_t.unwrap().as_slice());
    }
    println!("]");
}
