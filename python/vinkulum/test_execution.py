"""Native pool identity, numerical equivalence and concurrent model isolation."""

from concurrent.futures import ThreadPoolExecutor
import unittest

import numpy as np

from vinkulum import ExecutionPool, Noyau


def chain(count, pool):
    model = Noyau(executor=pool)
    length, mass = .05, .02
    inertia = np.diag([1e-8, mass * length**2 / 12, mass * length**2 / 12])
    for i in range(count):
        model.corps(str(i), mass, inertia.ravel().tolist(), [length * (i + .5), 0, 0])
        model.liaison(str(i), i - 1 if i else None, i,
                      pa=[length / 2 if i else 0, 0, 0], bloque_r=[])
    return model


def run_chain(pool, count=64):
    model = chain(count, pool)
    model.assemble()
    model.simule(.01, .001, tous=10**9)
    return model.etat(), model.stats(), model.phi()


class ExecutionContracts(unittest.TestCase):
    def test_strict_counts_and_model_ownership(self):
        for invalid in (True, False, 1.5, '2', None):
            with self.subTest(value=invalid), self.assertRaises(TypeError):
                ExecutionPool(invalid)
        for invalid in (0, -1, 10**100):
            with self.subTest(value=invalid), self.assertRaises((ValueError, OverflowError)):
                ExecutionPool(invalid)
        pool = ExecutionPool(2)
        first, second = Noyau(executor=pool), Noyau(executor=pool)
        self.assertEqual(first.execution_pool_id, pool.pool_id)
        self.assertEqual(second.execution_pool_id, pool.pool_id)
        del pool
        self.assertEqual(first.execution_threads, 2)
        self.assertEqual(second.execution_threads, 2)
        self.assertEqual(Noyau().execution_pool_id, 0)

    def test_parallel_assembly_keeps_same_state_iterations_and_constraints(self):
        # 64 bodies exceeds the assembly threshold (48); 256 also exercises
        # faer's explicit sparse parallel path (2304 saddle-point unknowns).
        for count in (64, 256):
            reference = run_chain(ExecutionPool(1), count)
            for threads in (2, 4):
                actual = run_chain(ExecutionPool(threads), count)
                self.assertEqual(actual[0], reference[0])
                self.assertEqual(actual[1], reference[1])
                np.testing.assert_array_equal(actual[2], reference[2])
                self.assertLess(np.linalg.norm(actual[2]), 1e-9)

    def test_concurrent_models_share_bounded_pool_and_failure_does_not_poison_it(self):
        pool = ExecutionPool(2)
        reference = run_chain(pool)
        with ThreadPoolExecutor(max_workers=3) as callers:
            results = list(callers.map(lambda _: run_chain(pool), range(3)))
        for result in results:
            self.assertEqual(result[0], reference[0])
            self.assertEqual(result[1], reference[1])
            np.testing.assert_array_equal(result[2], reference[2])
        failed = chain(64, pool)
        with self.assertRaises(RuntimeError):
            failed.simule(.01, .001, newton_max=1, tol=1e-30)
        self.assertEqual(run_chain(pool)[0], reference[0])

    def test_all_time_integrators_accept_shared_pool(self):
        for method in ('simule', 'simule_em', 'simule_multirythme'):
            outputs = []
            for threads in (1, 2):
                model = Noyau(g=[0, 0, 0], executor=ExecutionPool(threads))
                for i in range(2):
                    model.corps(str(i), 1, [1, 0, 0, 0, 2, 0, 0, 0, 3],
                                [i, 0, 0], w=[.1, .2, .3])
                kwargs = dict(rapides=[0], k=2) if method == 'simule_multirythme' else {}
                getattr(model, method)(.01, .001, **kwargs)
                outputs.append(model.etat())
            self.assertEqual(*outputs)


if __name__ == '__main__':
    unittest.main()
