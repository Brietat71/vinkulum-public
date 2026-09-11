> The maintained adapter package (development 0.6.1.dev2) no longer requires PySide6 or VTK.
> A new engine installation checks that neither is importable. FreeCAD keeps
> its own Qt in the host process. Existing engine environments are not modified;
> create a new environment to benefit from the smaller dependency set.
> See the [installed engine qualification](bancs/freecad-headless-engine-2026/README.md).

# Install the separate FreeCAD engine on Linux

The FreeCAD extension sends captured geometry to a separate Python process.
FreeCAD keeps its own OCCT and Python; the engine uses Vinkulum 0.20 and OCCT 8.
Run the installation below in a normal Linux terminal, outside FreeCAD's Python
console. It does not start the former Studio GUI.

## One installation command

Use a checkout containing [the installation script](../ci/freecad_engine.py),
with `uv`, a Rust/C build toolchain, `git` and `patch` available on `PATH`.
Internet access is required for Python, package downloads and source builds.
The qualified platform is Linux x86-64 / Ubuntu 24.04.

From the repository root:

```sh
python3 ci/freecad_engine.py "$HOME/vinkulum-freecad-engine"
```

The destination must be new. The command refuses an existing directory, including
an existing virtual environment; it never upgrades or removes it. Choose a new
destination for each installation. `--python /path/to/python3.14` selects an
existing interpreter instead of asking `uv` to provide Python 3.14.

The command creates `venv`, checks the upstream build123d and ocpsvg source
SHA-256 hashes, applies the repository's reviewed OCCT 8 adaptations, and builds
and installs the native kernel and backend adapters from this checkout.
Dependency constraints remain enforced. The engine installation omits Qt and
VTK and checks that neither is importable. The optional `legacy-desktop` extra
exists only for reproducing archived Studio work.

Before reporting success, it checks installed dependencies, imports OCCT 8 and
the kernel, runs an actual mechanics calculation against the independent pendulum
reference, and verifies rejection of an inconsistent captured mass. The tested
paths are not a general kernel or CAD certification.

## Connect FreeCAD

Install the [FreeCAD extension](../apps/freecad/README.md#install), then open
**Vinkulum → Motion analysis**. Set **Engine Python** to the path printed by the
installer, normally:

```text
/home/your-user/vinkulum-freecad-engine/venv/bin/python
```

Choose the output folder for captured calculations and click **Run motion**.
Keep the virtual-environment path as printed; resolving its Python symlink to a
system interpreter can select the wrong environment. Install no OCCT 8 bindings
into FreeCAD's embedded Python.

## Records and failure handling

The destination retains `installation.json`, `installation.log`, the installer
source and its hash, `runtime.json`, patched CAD sources and the virtual
environment. Records identify the source commit, whether the checkout was dirty,
executed commands and installed package versions. Primary CAD versions and source
hashes are pinned; the complete transitive package set is recorded after
resolution, not supplied as a globally locked dependency set.

On failure, the partial directory and logs remain available for diagnosis. Correct
the prerequisite or network problem and choose a new destination. An existing
working engine remains usable. Deleting an installation is a separate user action;
the installer does not perform cleanup of existing environments.

For actual FreeCAD verification, use a dedicated process and a new output directory:

```sh
python3 ci/freecad_extension.py \
  --freecad /absolute/path/to/FreeCAD/AppRun \
  --engine-python "$HOME/vinkulum-freecad-engine/venv/bin/python" \
  --archive /absolute/path/to/Vinkulum-FreeCAD-0.1.0a2.zip \
  --output /tmp/new-freecad-engine-check
```

The qualification opens and closes its own example documents and exits its
dedicated FreeCAD process. See the [document and lifecycle record](bancs/freecad-analysis-010a2/README.md)
for its scope. That record covers the released single-solid path. The
[development first-run guide](FREECAD_FIRST_RUN.md) describes the separately
qualified native Assembly and linear-static tasks; other operating systems remain
unqualified.

The [fresh-installation qualification record](bancs/freecad-engine-setup-2026/README.md)
contains the executed installer, runtime inventory and actual FreeCAD checks
against the newly created environment.

The development static workflow in extension 0.1.0a3.dev6 requires adapter
0.6.1.dev2 for captured CAD face witnesses. Create a fresh engine environment;
older installations do not provide this contract. See the explicit
[geometry persistence requirements and patched FreeCAD host recipe](FREECAD_GEOMETRY_IDENTITY.md).
