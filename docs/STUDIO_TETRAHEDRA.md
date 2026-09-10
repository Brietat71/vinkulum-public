# Tetrahedral CalculiX studies and curved-element checks

**Studio 0.6.0.dev3 source; experimental.** The statics workspace now reads,
calculates, displays and reopens studies using C3D4 or C3D10 tetrahedra, including
curved quadratic geometry. Affine C3D8 studies remain supported. This extends
the element contract now used by the [0.6.0.dev4 CAD meshing workspace](STUDIO_CAD_MESHING.md).
The public **Studio 0.5.0 Linux binary** predates
this extension.

The [installed-package qualification](bancs/studio-tetrahedra-060/README.md)
contains the actual calculation, test log and source fingerprints.

![Studio 0.6.0.dev3: six quadratic tetrahedra under pure bending, displacement field and integration-point values](assets/studio-tetra-bending.png)

## Try a complete calculation

Open [the quadratic pure-bending study](../examples/studio/fem/pure-bending-c3d10.ccx.json)
in **Run → Linear statics · CalculiX… → Open study…**, then run the installed
`ccx` executable. Select **Integration-point stress / energy** to inspect the
four material points per element. Geometry deformation is amplified only for
display; exported displacements and stresses retain SI values.

The example contains six quadratic tetrahedra in a 1 × 0.2 × 0.3 m specimen,
with E = 210 GPa and ν = 0.3. Opposite end tractions generate pure bending with
σxx = −10⁶(y − 0.1) Pa. Supports remove rigid motions while allowing the exact
Poisson contraction. This is not a fully clamped cantilever test.

For κ = 10⁶/E, the independent displacement field is:

```text
ux = −κ x(y − 0.1)
uy = κ/2 [x² + ν((y − 0.1)² − z² − 0.1²)]
uz = νκ(y − 0.1)z
```

Its strain energy is `10¹² I L / (2E)`, with
`I = 0.3 × 0.2³ / 12 m⁴`: approximately **0.00047619047619 J**.
Quadratic displacement interpolation represents this polynomial on straight
C3D10 geometry. Six elements reproduce the reference within the rounding
tolerance of the solver's ASCII tables. This exactness is specific to the
chosen polynomial solution.

The same load case exposes the stiffness of linear tetrahedra:

| C3D4 elements | Relative strain-energy error, CalculiX 2.21 |
|---:|---:|
| 30 | 76.0501% |
| 240 | 48.1062% |
| 1,920 | 20.4617% |

These errors decrease with refinement but remain large. The reference is useful
for diagnosing element behaviour, not for claiming a general solver speedup.
The [CalculiX manual, sections 6.2.6–6.2.7](https://www.dhondt.de/ccx_2.23.pdf)
documents one integration point for C3D4 and four for C3D10, and warns about
the stiffness of coarse linear-tetrahedron meshes. Local execution here uses
2.21; the manual's 2.23 version is not a tested executable claim.

## A bound over a curved element

A quadratic tetrahedral map has an affine Jacobian matrix. Therefore its
determinant is a polynomial of degree at most three in barycentric coordinates.
Positive values at the four integration points do not establish positivity
throughout the element. A regression moves one edge node so all four values
remain positive while the determinant at a vertex becomes negative.

The local checker follows the Bernstein bounding principle described by
[Johnen, Remacle and Geuzaine](https://gmsh.info/doc/preprints/gmsh_curved_preprint.pdf).
Bernstein basis functions are nonnegative and sum to one on the reference
simplex, so their coefficients bound the polynomial. Subdivision can tighten
an inconclusive bound. Vinkulum's implementation uses exact integer arithmetic
for binary64 geometry. The [0.6.0.dev4 numeric transport](CALCULIX_NUMERIC_TRANSPORT.md)
repeats admission on coordinates actually encoded in the solver input deck.

The implementation constructs its cubic coefficients without a floating-point
interpolation solve. Write the affine matrix as `J(λ) = Σ λi Ji`, where each
`Ji` is the matrix at a reference vertex. Expand the determinant by multilinearity
in its three columns. Group the 64 products by their barycentric multi-index;
divide each grouped sum by its positive multinomial count to obtain the 20
Bernstein coefficients. A shared power-of-two denominator makes every sign
decision integral. Longest-reference-edge bisection averages vertex matrices;
the same common-denominator representation remains exact on every child.

`certify_tetra10` returns one of three outcomes:

- **positive**: every examined leaf has strictly positive coefficients; the
  result includes an exact rational lower bound on the physical determinant.
- **nonpositive**: an exact vertex evaluation supplies a witness, including a
  vertex introduced by subdivision.
- **unresolved**: positivity could not be established within the depth/cell
  budget. The element is refused, without being labelled mathematically inverted.

Default limits are depth 10 and 255 examined cells per element. A separate
numerical quality rule requires the admitted determinant lower bound, scaled
by the corner-edge length cubed, to exceed 10⁻¹². The bound applies to the
**local determinant of the binary64 geometry**. It does not establish global
injectivity, absence of intersections between distant elements, fidelity to
the original CAD boundary, stiffness conditioning or solution accuracy.

## Quadrature, connectivity and archives

Each study declares one element family. C3D10 edge nodes follow CalculiX/VTK
order: `(1,2), (2,3), (3,1), (1,4), (2,4), (3,4)`. Shared edges must use the same
node identity, and an edge node cannot also be a corner node. Meshes must be
face connected and retain the existing node/element and file budgets.

Energy uses the physical weight `det(J(ξg)) wg` of each integration point.
Curved-element weights are generally unequal; averaging energy density and
multiplying by a corner-tetrahedron volume is incorrect. Tests include curved
traction and bending cases. Applied nodal forces come from a separate polynomial
boundary integral, independent of the implementation's shape-gradient and
volume-quadrature routines. Force, moment and work/energy balance checks remain
active on real CalculiX output.

The 0.6.0.dev3 studies use schema 2 with an explicit `element_type`; its captured
results use schema 2 / adapter 0.2.0 and identify their integration-point count. Legacy
schema-1 C3D8 studies and result folders remain readable. A tetrahedral study
cannot be downgraded to schema 1. Results remain **NotAssessed**.

Source version 0.6.0.dev4 uses result schema 3 for explicit numeric transport and
study schema 3 for CAD mesh/face provenance. Mixed-family meshes, warped C3D8
geometry, nonlinear materials, contact and discretisation-error bounds remain
separate work. The CAD meshing service requires Gmsh built with OCCT 8 or newer;
the local Gmsh 4.12.1 / OCC 7.6.3 binary is refused.
