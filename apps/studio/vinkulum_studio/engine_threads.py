"""Captured thread settings for optional native engines; no Qt or engine imports."""

import os
import re
from dataclasses import asdict

from .execution import ExecutionPlan


def thread_environment(threads, *, engine):
    ExecutionPlan(threads, threads)
    values = {
        "OMP_NUM_THREADS": str(threads),
        "OMP_THREAD_LIMIT": str(threads),
        "OMP_MAX_ACTIVE_LEVELS": "1",
        "OMP_DYNAMIC": "FALSE",
        "OPENBLAS_NUM_THREADS": "1",
        "GOTO_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "PYTHONUNBUFFERED": "1",
    }
    if engine == "calculix":
        values.update(
            {
                name: str(threads)
                for name in (
                    "NUMBER_OF_CPUS",
                    "CCX_NPROC_STIFFNESS",
                    "CCX_NPROC_RESULTS",
                    "CCX_NPROC_EQUATION_SOLVER",
                )
            }
        )
    elif engine == "occt":
        values["VINKULUM_OCCT_THREADS"] = str(threads)
    elif engine not in ("gmsh", "pinocchio"):
        raise ValueError("Unsupported native engine thread policy.")
    return values


def process_environment(threads, *, engine):
    return {**os.environ, **thread_environment(threads, engine=engine)}


def calculix_threads(log, plan):
    """Reported phase maxima, not a claim that every requested CPU was active."""
    phases = []
    for match in re.finditer(r"Using up to\s+(\d+)\s+cpu\(s\)\s+for\s+([^\r\n]+)", log):
        count, label = int(match[1]), match[2].strip().rstrip(".")
        if not 1 <= count <= plan.threads:
            raise ValueError("CalculiX reported a thread count outside its allocation.")
        phases.append({"phase": label, "threads": count})
    if not phases:
        raise ValueError("CalculiX did not report its native thread allocation.")
    return {
        "allocation": asdict(plan),
        "backend": "calculix",
        "reported_phases": phases,
    }


def configure_occt_threads(threads):
    """Call once in a fresh CAD worker, before importing build123d or making shapes."""
    ExecutionPlan(threads, threads)
    from OCP.OSD import OSD_Parallel, OSD_ThreadPool

    OSD_Parallel.SetUseOcctThreads_s(True)
    pool = OSD_ThreadPool.DefaultPool_s(threads)
    if pool.NbThreads() != threads:
        pool.Init(threads)
    pool.SetNbDefaultThreadsToLaunch(threads)
    if pool.NbThreads() != threads or not OSD_Parallel.ToUseOcctThreads_s():
        raise RuntimeError("OCCT did not accept its native CPU allocation.")
    return {
        "backend": "occt-native",
        "threads": pool.NbThreads(),
        "default_threads": pool.NbDefaultThreadsToLaunch(),
    }
