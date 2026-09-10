//! Reusable native CPU pools. A pool can be shared by independent models.
//!
//! `install` scopes nested Rayon assembly and faer work to the same pool. It
//! never changes the process-wide Rayon or faer configuration. Sharing a pool
//! bounds the workers even when several calling threads submit calculations.

use pyo3::prelude::*;
use pyo3::types::{PyBool, PyInt};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

/// Run in an explicit pool, retaining the historical executor when absent.
pub(crate) fn run<F, R>(pool: Option<&ExecutionPool>, operation: F) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    match pool {
        Some(pool) => pool.install(operation),
        None => operation(),
    }
}

/// A fixed-size native Rayon executor, shareable between independent `Noyau`.
/// Create `ExecutionPool(threads)` and pass it as `Noyau(executor=pool)`.
/// The pool limits native workers, not the amount of parallel work in a model.
#[pyclass(frozen, skip_from_py_object, module = "vinkulum._vinkulum")]
#[derive(Clone)]
pub struct ExecutionPool {
    pool: Arc<rayon::ThreadPool>,
    id: u64,
}

impl ExecutionPool {
    pub fn with_threads(threads: usize) -> Result<Self, String> {
        if threads == 0 || threads > rayon::max_num_threads() {
            return Err(format!(
                "ExecutionPool requires 1–{} threads; received {threads}",
                rayon::max_num_threads()
            ));
        }
        static NEXT: AtomicU64 = AtomicU64::new(1);
        let id = NEXT.fetch_add(1, Ordering::Relaxed);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .thread_name(move |index| format!("vinkulum-cpu-{id}-{index}"))
            .build()
            .map_err(|error| format!("Could not create native CPU pool: {error}"))?;
        Ok(Self {
            pool: Arc::new(pool),
            id,
        })
    }

    /// Nested native parallel operations use this pool's workers.
    pub fn install<F, R>(&self, operation: F) -> R
    where
        F: FnOnce() -> R + Send,
        R: Send,
    {
        self.pool.install(operation)
    }
}

#[pymethods]
impl ExecutionPool {
    #[new]
    fn new(threads: &Bound<'_, PyAny>) -> PyResult<Self> {
        if !threads.is_instance_of::<PyInt>() || threads.is_instance_of::<PyBool>() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "ExecutionPool requires an integer thread count, excluding bool",
            ));
        }
        Self::with_threads(threads.extract()?).map_err(pyo3::exceptions::PyValueError::new_err)
    }

    /// Number of native worker threads in this pool.
    #[getter]
    pub fn threads(&self) -> usize {
        self.pool.current_num_threads()
    }

    /// Process-local identity; clones and models sharing a pool retain it.
    #[getter]
    pub fn pool_id(&self) -> u64 {
        self.id
    }

    fn __repr__(&self) -> String {
        format!(
            "ExecutionPool(threads={}, pool_id={})",
            self.threads(),
            self.id
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::prelude::*;
    use std::collections::HashSet;

    #[test]
    fn native_workers_are_distinct_and_nested_work_stays_in_its_pool() {
        let pool = ExecutionPool::with_threads(2).unwrap();
        let workers = pool.pool.broadcast(|_| {
            assert_eq!(rayon::current_num_threads(), 2);
            std::thread::current().id()
        });
        assert_eq!(workers.into_iter().collect::<HashSet<_>>().len(), 2);
        let result = pool.install(|| {
            (0..128)
                .into_par_iter()
                .map(|i| {
                    assert_eq!(rayon::current_num_threads(), 2);
                    assert!(pool.pool.current_thread_index().is_some());
                    i * i
                })
                .collect::<Vec<_>>()
        });
        assert_eq!(result, (0..128).map(|i| i * i).collect::<Vec<_>>());
    }

    #[test]
    fn sharing_and_concurrent_pools_do_not_reconfigure_each_other() {
        let a = ExecutionPool::with_threads(1).unwrap();
        let b = ExecutionPool::with_threads(2).unwrap();
        let shared = b.clone();
        assert!(Arc::ptr_eq(&b.pool, &shared.pool));
        assert_eq!(b.pool_id(), shared.pool_id());
        assert_ne!(a.pool_id(), b.pool_id());
        std::thread::scope(|scope| {
            let first = scope.spawn(|| a.install(rayon::current_num_threads));
            let second = scope.spawn(|| shared.install(rayon::current_num_threads));
            assert_eq!(first.join().unwrap(), 1);
            assert_eq!(second.join().unwrap(), 2);
        });
        assert!(ExecutionPool::with_threads(0).is_err());
        assert!(ExecutionPool::with_threads(usize::MAX).is_err());
    }
}
