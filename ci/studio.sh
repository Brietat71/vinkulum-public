#!/usr/bin/env bash
# Desktop tests use a real Qt platform plugin and VTK render window.
set -euo pipefail
VINKULUM_SOURCE_ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=${PY:-python}
export VINKULUM_3D_TESTS=1
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-2}
export VINKULUM_STUDIO_CPUS=${VINKULUM_STUDIO_CPUS:-2}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
case "$(uname -s)" in
  Linux)
    export QT_QPA_PLATFORM=xcb
    # Probe a fresh process: another library must not mask missing Qt/X11
    # dependencies by loading its own bundled copies before QApplication.
    xvfb-run -a "$PY" -c 'from PySide6.QtWidgets import QApplication; app = QApplication([]); print("Qt/X11 startup OK before importing CAD or VTK")'
    xvfb-run -a -s '-screen 0 1600x1100x24' "$PY" -X faulthandler "$VINKULUM_SOURCE_ROOT/ci/studio_test_runner.py" discover -s "$VINKULUM_SOURCE_ROOT/apps/studio/tests" -v
    ;;
  Darwin)
    export QT_QPA_PLATFORM=cocoa
    # Keyboard shortcuts require a foreground Cocoa application. Disabling
    # Qt's normal activation policy leaves a plain interpreter without focus.
    unset QT_MAC_DISABLE_FOREGROUND_APPLICATION_TRANSFORM
    "$PY" -X faulthandler "$VINKULUM_SOURCE_ROOT/ci/studio_test_runner.py" discover -s "$VINKULUM_SOURCE_ROOT/apps/studio/tests" -v
    ;;
  *) echo "Plateforme Studio non qualifiée." >&2; exit 1 ;;
esac
