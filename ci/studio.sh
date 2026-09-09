#!/usr/bin/env bash
# Desktop tests use a real Qt platform plugin and VTK render window.
set -euo pipefail
VINKULUM_SOURCE_ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-python}
export VINKULUM_3D_TESTS=1
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-2}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
case "$(uname -s)" in
  Linux)
    export QT_QPA_PLATFORM=xcb
    xvfb-run -a -s '-screen 0 1600x1100x24' "$PY" -X faulthandler -m unittest discover -s "$VINKULUM_SOURCE_ROOT/apps/studio/tests" -v
    ;;
  Darwin)
    export QT_QPA_PLATFORM=cocoa
    "$PY" -X faulthandler -m unittest discover -s "$VINKULUM_SOURCE_ROOT/apps/studio/tests" -v
    ;;
  *) echo "Plateforme Studio non qualifiée." >&2; exit 1 ;;
esac
