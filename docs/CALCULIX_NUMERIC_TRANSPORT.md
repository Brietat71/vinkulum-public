# CalculiX input precision and the geometry actually solved

This developing source change introduces result schema 3 / adapter 0.3.0.
It is not in the published Studio 0.5.0 Linux binary.

## The failure

CalculiX reads node coordinates, concentrated loads and elastic material values
through 20-character Fortran numeric fields. The previous `.17g` exporter could
exceed that width. Scientific notation could be truncated into an invalid
exponent, or a valid but different number. The issue became visible on a real
Gmsh mesh, with coordinates such as `-3.552713678800501e-18`.

This concerns the external CalculiX adapter's decimal serialization, not the
Rust multibody equations. The reader is visible in `nodes.f`, `cloads.f` and
`elastics.f` in the [official CalculiX source](https://www.dhondt.de/ccx_2.23.src.tar.bz2).
Local solver qualification uses the installed CalculiX 2.21 executable.

## Transport contract

Every exported floating-point field fits in 20 characters. The formatter first
tries to preserve the binary64 value exactly, removing redundant characters
and using Fortran's permitted implicit exponent notation. For example,
`-3552713678800501-33` fits and represents the coordinate above without changing
its binary64 value. The exact deck remains an ASCII file.

If an exact representation does not fit, the formatter considers nearby decimal
representations with at least 15 significant digits. It accepts a conversion
only when the exact rational difference is at most **64 ULPs of the requested
value**, including subnormal inputs. Nonfinite values are refused. The fixed
budget is an admission rule, not a claim about how that perturbation affects a
poorly conditioned mechanical problem.

The requested study remains unchanged in `study.json`. `StaticStudy.solver_study`
derives the numerical coordinates, nodal forces and material values actually
encoded by the new deck. If any values change, the derived study passes the
same element geometry, connectivity, material and rigid-mode admission checks.
This happens before a calculation directory is created. In particular, the
exact local C3D10 Jacobian check is repeated on the transported coordinates.
A thin tetrahedron whose distinct coordinates collapse during decimal
conversion is rejected before the solver starts.

Result validation computes reactions, moment balance, work and integrated energy
using the transported study. The result view and its load table use those same
inputs. The requested CAD mesh and its boundary-condition intent remain available
in the original study capture; rounding does not rewrite or reattach that intent.

## Captured evidence

Schema-3 results include `input_transport`, with:

- the `calculix-f20-v1` format and the enforced field/ULP budgets;
- changed-value counts and maximum absolute changes for coordinates in metres,
  nodal forces in newtons, Young modulus in pascals and Poisson ratio;
- a SHA-256 fingerprint of the derived solver study.

Reopening regenerates the transport, revalidates its study, compares the deck and
the captured transport contract, and rechecks the raw solver output. It does not
run an engine. As before, fingerprints establish archive consistency, not
authorship or the historical identity of a solver installation.

Schema-1 and schema-2 result folders retain their original deck and result
interpretation. Reading them does not silently upgrade their numerical contract.
New runs use schema 3 even when every input value survives transport unchanged.

## Independent checks and remaining limits

`apps/studio/tests/test_calculix_numbers.py` exercises randomly sampled finite
binary64 values, extreme exponents, signed zero and subnormals. An independent
Fortran program is compiled with `gfortran` and uses the upstream reader's
`F20.0` format; its decoded binary64 bits must match the transported values.
This test is skipped explicitly when the compiler is unavailable.

The suite also runs a translated tension specimen through real CalculiX,
compares its energy against the analytic reference, reopens the original request
and result, checks published legacy archives, and refuses geometry that would
collapse during conversion. These checks do not certify all Fortran runtimes,
global mesh injectivity, stiffness conditioning or discretisation error.

The full CAD meshing interface and its publication are separate work. These
transport checks are required before using arbitrary mesher coordinates in
that workflow.
