"""Verify installed scientific adapters without any Studio desktop dependency."""

import importlib
import importlib.metadata
import importlib.util
import json
import sys
from pathlib import Path

GUI_MODULES = ("PySide6", "shiboken6", "vtk", "vtkmodules")
ADAPTER_MODULES = (
    "cad",
    "calculix",
    "document",
    "engine_threads",
    "execution",
    "mechanism",
    "mesh_binding",
    "meshing",
    "model",
    "pinocchio_backend",
)


def verify():
    repository = Path(__file__).resolve().parents[1]
    for name in GUI_MODULES:
        if importlib.util.find_spec(name) is not None:
            raise RuntimeError(
                f"Desktop dependency present in engine environment: {name}"
            )
    imported = {}
    for name in ADAPTER_MODULES:
        module = importlib.import_module(f"vinkulum_studio.{name}")
        path = Path(module.__file__).resolve()
        if path.is_relative_to(repository):
            raise RuntimeError(
                f"Expected installed adapter outside the checkout: {path}"
            )
        imported[name] = str(path)
    import build123d
    import OCP
    import vinkulum

    if OCP.__version__.split(".")[0] != "8":
        raise RuntimeError("The engine requires OCCT 8")
    if any(name in sys.modules for name in GUI_MODULES):
        raise RuntimeError("An engine import loaded a desktop dependency")
    return {
        "status": "passed",
        "python": sys.version,
        "adapter_version": importlib.metadata.version("vinkulum-studio"),
        "ocp_binding": OCP.__version__,
        "absent_modules": list(GUI_MODULES),
        "installed_adapters": imported,
        "build123d": str(Path(build123d.__file__).resolve()),
        "vinkulum": str(Path(vinkulum.__file__).resolve()),
    }


if __name__ == "__main__":
    result = json.dumps(verify(), indent=2) + "\n"
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(result)
    print(result, end="")
