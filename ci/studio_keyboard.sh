#!/usr/bin/env bash
# Issue #2: keyboard-only CAD qualification on an isolated Linux/X11 desktop.
set -euo pipefail
VINKULUM_SOURCE_ROOT=$(cd "$(dirname "$0")/.." && pwd)
export VINKULUM_SOURCE_ROOT
export PY=${PY:-python}
export VINKULUM_KEYBOARD_TESTS=1
export QT_QPA_PLATFORM=xcb
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 RAYON_NUM_THREADS=2

if [[ ${1:-} == --session ]]; then
  openbox --sm-disable --config-file /etc/xdg/openbox/rc.xml > "$VINKULUM_KEYBOARD_SCREENSHOTS/openbox.log" 2>&1 &
  VINKULUM_WM_PID=$!
  trap 'kill "$VINKULUM_WM_PID" 2>/dev/null || true; wait "$VINKULUM_WM_PID" 2>/dev/null || true' EXIT
  ready=0
  for attempt in {1..100}; do
    if [[ $(xprop -root _NET_SUPPORTING_WM_CHECK 2>/dev/null) == *"window id"* ]]; then
      ready=1
      break
    fi
    sleep 0.05
  done
  [[ $ready == 1 ]] || { echo "Openbox did not start." >&2; exit 1; }
  "$PY" -X faulthandler -m unittest discover -s "$VINKULUM_SOURCE_ROOT/apps/studio/tests" -p test_keyboard_cad.py -v
  exit
fi

[[ $(uname -s) == Linux ]] || { echo "This recipe qualifies Linux/X11 only." >&2; exit 1; }
for tool in xvfb-run openbox xprop; do
  command -v "$tool" >/dev/null || { echo "Missing $tool (install xvfb, xauth, openbox and x11-utils)." >&2; exit 1; }
done
output=${1:-$(mktemp -d /tmp/vinkulum-keyboard-XXXXXX)}
mkdir -p "$output"
output=$(cd "$output" && pwd)
"$PY" - <<'PY' > "$output/versions.txt"
import platform
from importlib.metadata import version
print(platform.platform())
print('Python', platform.python_version())
for name in ('vinkulum', 'PySide6', 'vtk', 'cadquery-ocp-novtk', 'build123d'):
    print(name, version(name))
import vinkulum_studio
print('Studio', vinkulum_studio.__version__, vinkulum_studio.__file__)
PY
openbox --version >> "$output/versions.txt"
for scale in 1 2; do
  export QT_SCALE_FACTOR=$scale
  export VINKULUM_KEYBOARD_SCREENSHOTS="$output/$((100 * scale))"
  mkdir -p "$VINKULUM_KEYBOARD_SCREENSHOTS"
  xvfb-run -a -s "-screen 0 $((1600 * scale))x$((1100 * scale))x24" \
    bash "$VINKULUM_SOURCE_ROOT/ci/studio_keyboard.sh" --session \
    2>&1 | tee "$VINKULUM_KEYBOARD_SCREENSHOTS/recipe.log"
done
printf 'Keyboard CAD qualification: %s\n' "$output"
