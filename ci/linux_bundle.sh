#!/usr/bin/env bash
# Build locally with an environment containing the native kernel, Studio and PyInstaller.
set -euo pipefail
[[ "$(uname -sm)" == "Linux x86_64" ]] || { echo 'Linux x86_64 required.' >&2; exit 1; }
VINKULUM_SOURCE_ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$VINKULUM_SOURCE_ROOT"
PY=${PY:-python}
PY=$("$PY" -c 'import sys; print(sys.executable)')
VERSION=$("$PY" -c 'import tomllib; print(tomllib.load(open("apps/studio/pyproject.toml", "rb"))["project"]["version"])')
"$PY" -c 'import importlib.metadata, sys; assert importlib.metadata.version("vinkulum-studio") == sys.argv[1], "Reinstall Studio before packaging"' "$VERSION"
VINKULUM_LINUX_OUT=${VINKULUM_LINUX_OUT:-"$VINKULUM_SOURCE_ROOT/dist/linux"}
mkdir -p "$VINKULUM_LINUX_OUT"
VINKULUM_LINUX_OUT=$(cd "$VINKULUM_LINUX_OUT" && pwd)
BUILD_DIR=$(mktemp -d /tmp/vinkulum-linux.XXXXXX)
trap 'rm -rf -- "$BUILD_DIR"' EXIT
"$PY" -m PyInstaller --noconfirm --clean --distpath "$BUILD_DIR/dist" \
  --workpath "$BUILD_DIR/build" apps/studio/packaging/studio.spec
APP="$BUILD_DIR/dist/Vinkulum Studio"
cp apps/studio/packaging/INSTALLATION-LINUX.txt "$APP/INSTALLATION.txt"
mkdir -p "$APP/Licences"
cp -R apps/studio/packaging/licenses/. "$APP/Licences/"
cp LICENSE NOTICE THIRD_PARTY_NOTICES.md "$APP/Licences/"
"$PY" - "$APP/build-info.json" <<'PY'
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys

report = {
    "platform": platform.platform(),
    "libc": platform.libc_ver(),
    "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "source_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
    "packages": {name: importlib.metadata.version(name) for name in
                 ("vinkulum", "vinkulum-studio", "PySide6", "vtk", "numpy", "pyinstaller", "cadquery-ocp-novtk", "build123d", "ocpsvg", "ocp_gordon")},
}
Path(sys.argv[1]).write_text(json.dumps(report, indent=2) + "\n")
PY
# Test the archived payload after extraction, outside the source/build trees.
ARCHIVE="Vinkulum-Studio-$VERSION-linux-x86_64.tar.gz"
tar -czf "$BUILD_DIR/$ARCHIVE" -C "$BUILD_DIR/dist" 'Vinkulum Studio'
mkdir "$BUILD_DIR/extracted"
tar -xzf "$BUILD_DIR/$ARCHIVE" -C "$BUILD_DIR/extracted"
cd "$BUILD_DIR/extracted"
COMMAND=(env -u PYTHONPATH -u PYTHONHOME -u VIRTUAL_ENV
  XDG_CONFIG_HOME="$BUILD_DIR/config" QT_QPA_PLATFORM=xcb
  "$BUILD_DIR/extracted/Vinkulum Studio/Vinkulum Studio")
STARTUP=("${COMMAND[@]}" --startup-check "$BUILD_DIR/check")
CHECK=("${COMMAND[@]}" --bundle-check "$BUILD_DIR/check")
if [[ -n "${DISPLAY:-}" ]]; then
  timeout 30s "${STARTUP[@]}"
  "${CHECK[@]}"
else
  xvfb-run -a -s '-screen 0 1600x1100x24' timeout 30s "${STARTUP[@]}"
  xvfb-run -a -s '-screen 0 1600x1100x24' "${CHECK[@]}"
fi
"$PY" - "$BUILD_DIR/check/bundle-check.json" "$VERSION" <<'PY'
import json
import sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
assert report["status"] == "passed" and report["frozen"]
assert report["machine"] == "x86_64"
assert report["manifest"]["app_version"] == sys.argv[2]
from pathlib import Path
startup = json.loads(Path(sys.argv[1]).with_name("startup-check.json").read_text())
assert startup["status"] == "passed" and startup["frozen"]
assert startup["app_version"] == sys.argv[2] and startup["locale"] == "en_GB"
PY
# Only deliver an archive whose extracted executable passed the check.
cp "$BUILD_DIR/$ARCHIVE" "$VINKULUM_LINUX_OUT/"
mkdir -p "$VINKULUM_LINUX_OUT/check"
cp "$BUILD_DIR/check/bundle-check.json" "$BUILD_DIR/check/scene.png" "$VINKULUM_LINUX_OUT/check/"
cp "$BUILD_DIR/check/startup-check.json" "$BUILD_DIR/check/startup-scene.png" "$VINKULUM_LINUX_OUT/check/"
cp "$APP/build-info.json" "$VINKULUM_LINUX_OUT/"
cd "$VINKULUM_LINUX_OUT"
sha256sum "$ARCHIVE" > SHA256SUMS
sha256sum --check SHA256SUMS
echo "Linux bundle tested and ready: $VINKULUM_LINUX_OUT/$ARCHIVE"
