//! Frontière sûre des certificats C++ : données possédées et CSR vérifié.
//! Les seules zones unsafe sont les appels à l'ABI bornée, sans pointeur Python.
use pyo3::exceptions::{PyMemoryError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

type Sparse = (Vec<i64>, Vec<i64>, Vec<f64>);
type Retour = (
    i32,
    Py<PyBytes>,
    Py<PyBytes>,
    Py<PyBytes>,
    Vec<i64>,
    Vec<f64>,
    Vec<f64>,
);
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Pivot {
    indices: [i64; 4],
    valeurs: [f64; 8],
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Summary {
    valeurs: [i64; 8],
}
const _: () = assert!(std::mem::size_of::<Pivot>() == 96 && std::mem::size_of::<Summary>() == 64);

extern "C" {
    fn vinkulum_inertia_binary64(
        n: i64,
        s: i64,
        nr: i64,
        dp: *const i64,
        dj: *const i64,
        dv: *const f64,
        mp: *const i64,
        mj: *const i64,
        mv: *const f64,
        b: *const f64,
        gamma: f64,
        perm: *const i64,
        ops: i64,
        coeff: i64,
        diagonal: bool,
        pm: *mut Pivot,
        pk: *mut Pivot,
        summaries: *mut Summary,
        diagnostic: *mut i64,
        bad: *mut f64,
        phases: *mut f64,
    ) -> i32;
    fn vinkulum_compare_mass(
        n: i64,
        mp: *const i64,
        mj: *const i64,
        mv: *const f64,
        ell: *const f64,
        alpha: f64,
        beta: f64,
        perm: *const i64,
        ops: i64,
        coeff: i64,
        lower: *mut Pivot,
        upper: *mut Pivot,
        summaries: *mut Summary,
        diagnostic: *mut i64,
        bad: *mut f64,
    ) -> i32;
    fn vinkulum_interval_binary64(
        op: i32,
        alo: f64,
        ahi: f64,
        blo: f64,
        bhi: f64,
        out: *mut f64,
    ) -> i32;
}
fn erreur(message: &str) -> PyErr {
    PyValueError::new_err(message.to_owned())
}
fn valider_csr(a: &Sparse, rows: usize, cols: usize) -> PyResult<()> {
    let (p, j, v) = a;
    if p.len()
        != rows
            .checked_add(1)
            .ok_or_else(|| erreur("dimension CSR excessive"))?
        || p.first() != Some(&0)
        || p.last().copied() != i64::try_from(j.len()).ok()
        || j.len() != v.len()
        || p.windows(2).any(|w| w[0] < 0 || w[0] > w[1])
        || j.iter().any(|&x| x < 0 || x as usize >= cols)
        || v.iter().any(|x| !x.is_finite())
    {
        return Err(erreur("stockage CSR ou coefficients invalides"));
    }
    for row in 0..rows {
        if j[p[row] as usize..p[row + 1] as usize]
            .windows(2)
            .any(|w| w[0] >= w[1])
        {
            return Err(erreur("CSR canonique trié sans doublon requis"));
        }
    }
    Ok(())
}
fn valider_masse(m: &Sparse, n: usize) -> PyResult<bool> {
    valider_csr(m, n, n)?;
    let (p, j, v) = m;
    for row in 0..n {
        for k in p[row] as usize..p[row + 1] as usize {
            let col = j[k] as usize;
            let start = p[col] as usize;
            let end = p[col + 1] as usize;
            let other = j[start..end]
                .binary_search(&(row as i64))
                .ok()
                .map_or(0., |offset| v[start + offset]);
            if v[k] != other {
                return Err(erreur("masse exactement symétrique requise"));
            }
        }
    }
    let diagonal = j.len() == n && (0..n).all(|i| p[i] == i as i64 && j[i] == i as i64);
    if diagonal && v.iter().any(|&x| x <= 0.) {
        return Err(erreur("masse diagonale positive requise"));
    }
    Ok(diagonal)
}
fn valider_communs(n: usize, perm: &[i64], ops: i64, coeff: i64) -> PyResult<()> {
    if n == 0 || i64::try_from(n).is_err() || ops < 1 || coeff < 1 || n as u64 > coeff as u64 {
        return Err(erreur("dimensions ou budgets invalides"));
    }
    if perm.len() != n {
        return Err(erreur("permutation complète requise"));
    }
    let mut sorted = perm.to_vec();
    sorted.sort_unstable();
    if sorted.iter().enumerate().any(|(i, &v)| v != i as i64) {
        return Err(erreur("permutation invalide"));
    }
    Ok(())
}
fn allouer<T: Clone + Default>(n: usize) -> PyResult<Vec<T>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(n)
        .map_err(|_| PyMemoryError::new_err("allocation des preuves indisponible"))?;
    result.resize(n, T::default());
    Ok(result)
}
struct Brut {
    status: i32,
    first: Vec<Pivot>,
    second: Vec<Pivot>,
    summaries: [Summary; 2],
    diag: [i64; 3],
    bad: [f64; 2],
    phases: [f64; 3],
}
impl Brut {
    fn new(n: usize, total: usize) -> PyResult<Self> {
        Ok(Self {
            status: 0,
            first: allouer(n)?,
            second: allouer(total)?,
            summaries: [Summary::default(); 2],
            diag: [0; 3],
            bad: [0.; 2],
            phases: [0.; 3],
        })
    }
    fn python(self, py: Python<'_>) -> Retour {
        fn bytes(p: &[Pivot]) -> Vec<u8> {
            p.iter()
                .flat_map(|x| {
                    x.indices
                        .iter()
                        .map(|v| v.to_ne_bytes())
                        .chain(x.valeurs.iter().map(|v| v.to_ne_bytes()))
                })
                .flatten()
                .collect()
        }
        let summaries: Vec<u8> = self
            .summaries
            .iter()
            .flat_map(|s| s.valeurs.iter().flat_map(|v| v.to_ne_bytes()))
            .collect();
        (
            self.status,
            PyBytes::new(py, &bytes(&self.first)).unbind(),
            PyBytes::new(py, &bytes(&self.second)).unbind(),
            PyBytes::new(py, &summaries).unbind(),
            self.diag.to_vec(),
            self.bad.to_vec(),
            self.phases.to_vec(),
        )
    }
}

#[pyfunction]
fn _certificat_inertie_binaire(
    py: Python<'_>,
    n: usize,
    rang: usize,
    d: Sparse,
    m: Sparse,
    b: Vec<f64>,
    gamma: f64,
    permutation: Vec<i64>,
    budget_operations: i64,
    budget_coefficients: i64,
    budget_rectangulaire: usize,
) -> PyResult<Retour> {
    valider_communs(n, &permutation, budget_operations, budget_coefficients)?;
    let total = n
        .checked_add(rang)
        .ok_or_else(|| erreur("dimension conservée excessive"))?;
    let count = n
        .checked_mul(rang)
        .ok_or_else(|| erreur("contraintes trop grandes"))?;
    if rang == 0
        || rang >= n
        || total as u64 > budget_coefficients as u64
        || count > budget_rectangulaire
        || count > i64::MAX as usize
        || b.len() != count
        || b.iter().any(|v| !v.is_finite())
        || !gamma.is_finite()
        || gamma <= 0.
    {
        return Err(erreur("rang, seuil ou budget de contraintes invalide"));
    }
    let nr =
        d.0.len()
            .checked_sub(1)
            .ok_or_else(|| erreur("pointeurs D absents"))?;
    valider_csr(&d, nr, n)?;
    let diagonal = valider_masse(&m, n)?;
    let mut result = Brut::new(n, total)?;
    let result = py.detach(move || {
        // Sécurité : CSR, permutations, dimensions et B sont vérifiés ci-dessus.
        // Les vecteurs sont possédés, immuables, et vivants durant tout l'appel.
        // Le C++ écrit au plus n et n+rang pivots et les tableaux fixes de l'ABI.
        result.status = unsafe {
            vinkulum_inertia_binary64(
                n as i64,
                rang as i64,
                nr as i64,
                d.0.as_ptr(),
                d.1.as_ptr(),
                d.2.as_ptr(),
                m.0.as_ptr(),
                m.1.as_ptr(),
                m.2.as_ptr(),
                b.as_ptr(),
                gamma,
                permutation.as_ptr(),
                budget_operations,
                budget_coefficients,
                diagonal,
                result.first.as_mut_ptr(),
                result.second.as_mut_ptr(),
                result.summaries.as_mut_ptr(),
                result.diag.as_mut_ptr(),
                result.bad.as_mut_ptr(),
                result.phases.as_mut_ptr(),
            )
        };
        result
    });
    Ok(result.python(py))
}
#[pyfunction]
fn _certificat_comparaison_masse(
    py: Python<'_>,
    n: usize,
    m: Sparse,
    ell: Vec<f64>,
    alpha: f64,
    beta: f64,
    permutation: Vec<i64>,
    budget_operations: i64,
    budget_coefficients: i64,
) -> PyResult<Retour> {
    valider_communs(n, &permutation, budget_operations, budget_coefficients)?;
    valider_masse(&m, n)?;
    if ell.len() != n
        || ell.iter().any(|&x| !x.is_finite() || x <= 0.)
        || !alpha.is_finite()
        || !beta.is_finite()
        || alpha <= 0.
        || beta <= alpha
    {
        return Err(erreur("comparaison et diagonale invalides"));
    }
    let mut result = Brut::new(n, n)?;
    let result = py.detach(move || {
        // Sécurité : toutes les longueurs et tous les indices ont été vérifiés.
        // Les sorties possèdent n pivots par preuve et deux résumés fixes.
        result.status = unsafe {
            vinkulum_compare_mass(
                n as i64,
                m.0.as_ptr(),
                m.1.as_ptr(),
                m.2.as_ptr(),
                ell.as_ptr(),
                alpha,
                beta,
                permutation.as_ptr(),
                budget_operations,
                budget_coefficients,
                result.first.as_mut_ptr(),
                result.second.as_mut_ptr(),
                result.summaries.as_mut_ptr(),
                result.diag.as_mut_ptr(),
                result.bad.as_mut_ptr(),
            )
        };
        result
    });
    Ok(result.python(py))
}
#[pyfunction]
fn _intervalle_binaire(op: i32, alo: f64, ahi: f64, blo: f64, bhi: f64) -> (i32, Vec<f64>) {
    let mut out = [0.; 2];
    // Sécurité : l'ABI scalaire valide ses opérandes et écrit deux flottants.
    let status = unsafe { vinkulum_interval_binary64(op, alo, ahi, blo, bhi, out.as_mut_ptr()) };
    (status, out.to_vec())
}
#[pyfunction]
fn _construction_certificat() -> (&'static str, &'static str, &'static str, &'static str) {
    (
        env!("VINKULUM_CERTIFICAT_COMPILATEUR"),
        env!("VINKULUM_CERTIFICAT_OPTIONS"),
        include_str!("certificat/inertie_binaire_native.cpp"),
        include_str!("certificat/comparaison_masse_native.cpp"),
    )
}
pub fn enregistrer(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_certificat_inertie_binaire, m)?)?;
    m.add_function(wrap_pyfunction!(_certificat_comparaison_masse, m)?)?;
    m.add_function(wrap_pyfunction!(_intervalle_binaire, m)?)?;
    m.add_function(wrap_pyfunction!(_construction_certificat, m)?)?;
    Ok(())
}
