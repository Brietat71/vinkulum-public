#!/usr/bin/env bash
# Build on native Apple Silicon, then qualify the application copied from the DMG.
set -euo pipefail
[[ "$(uname -sm)" == "Darwin arm64" ]] || { echo 'Native Apple Silicon macOS required.' >&2; exit 1; }
PY=${PY:-python}
PY=$("$PY" -c 'import sys; print(sys.executable)')
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
VERSION=$("$PY" -c 'import tomllib; print(tomllib.load(open("apps/studio/pyproject.toml", "rb"))["project"]["version"])')
"$PY" -c 'import importlib.metadata, sys; assert importlib.metadata.version("vinkulum-studio") == sys.argv[1], "Reinstall Studio before packaging"' "$VERSION"
OUT="$ROOT/dist/macos"
mkdir -p "$OUT"
BUILD_DIR=$(mktemp -d "${RUNNER_TEMP:-/tmp}/vinkulum-macos.XXXXXX")
MOUNT="$BUILD_DIR/mounted"
cleanup() {
  status=$?
  if mount | /usr/bin/grep -Fq " on $MOUNT "; then
    hdiutil detach "$MOUNT" || true
  fi
  if (( status == 0 )); then
    rm -rf -- "$BUILD_DIR"
  else
    echo "Build or qualification failed; diagnostics preserved in $BUILD_DIR" >&2
  fi
}
trap cleanup EXIT
export VINKULUM_MACOS_PAYLOAD="$BUILD_DIR/payload"
mkdir -p "$VINKULUM_MACOS_PAYLOAD"
"$PY" ci/studio_bundle_examples.py "$VINKULUM_MACOS_PAYLOAD/Examples"
"$PY" - "$VINKULUM_MACOS_PAYLOAD/build-info.json" <<'PY'
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import subprocess
import sys

report = {
    "platform": platform.platform(),
    "machine": platform.machine(),
    "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "source_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
    "workflow_run": os.environ.get("GITHUB_RUN_ID"),
    "signing": "ad-hoc; not notarised",
    "external_engines_bundled": [],
    "packages": {name: importlib.metadata.version(name) for name in
                 ("vinkulum", "vinkulum-studio", "PySide6", "vtk", "numpy", "pyinstaller", "cadquery-ocp-novtk", "build123d", "ocpsvg", "ocp_gordon")},
}
# Packaging must use exactly the reviewed Studio modules, not an older installed wheel.
import vinkulum_studio
installed = Path(vinkulum_studio.__file__).parent
for source in Path("apps/studio/vinkulum_studio").rglob("*.py"):
    assert source.read_bytes() == (installed / source.relative_to("apps/studio/vinkulum_studio")).read_bytes(), source
Path(sys.argv[1]).write_text(json.dumps(report, indent=2) + "\n")
PY
"$PY" -m PyInstaller --noconfirm --clean --distpath "$BUILD_DIR/dist" \
  --workpath "$BUILD_DIR/build" apps/studio/packaging/studio.spec
APP="$BUILD_DIR/dist/Vinkulum Studio.app"
codesign --verify --deep --strict "$APP"
lipo "$APP/Contents/MacOS/Vinkulum Studio" -verify_arch arm64
mkdir "$BUILD_DIR/dmg-root"
/usr/bin/ditto "$APP" "$BUILD_DIR/dmg-root/Vinkulum Studio.app"
ln -s /Applications "$BUILD_DIR/dmg-root/Applications"
cp apps/studio/packaging/INSTALLATION.txt "$BUILD_DIR/dmg-root/"
cp -R "$VINKULUM_MACOS_PAYLOAD/Examples" "$BUILD_DIR/dmg-root/Examples"
cp -R apps/studio/packaging/licenses "$BUILD_DIR/dmg-root/Licences"
cp LICENSE NOTICE THIRD_PARTY_NOTICES.md "$BUILD_DIR/dmg-root/Licences/"
DMG="Vinkulum-Studio-$VERSION-apple-silicon.dmg"
hdiutil create -volname 'Vinkulum Studio' -srcfolder "$BUILD_DIR/dmg-root" -ov -format UDZO "$BUILD_DIR/$DMG"
hdiutil verify "$BUILD_DIR/$DMG"
mkdir "$MOUNT" "$BUILD_DIR/installed"
hdiutil attach "$BUILD_DIR/$DMG" -readonly -nobrowse -mountpoint "$MOUNT"
/usr/bin/ditto "$MOUNT/Vinkulum Studio.app" "$BUILD_DIR/installed/Vinkulum Studio.app"
hdiutil detach "$MOUNT"
APP="$BUILD_DIR/installed/Vinkulum Studio.app"
codesign --verify --deep --strict "$APP"
lipo "$APP/Contents/MacOS/Vinkulum Studio" -verify_arch arm64
# The installed payload runs outside the source/build trees without Python paths.
cd "$BUILD_DIR/installed"
COMMAND=(env -u PYTHONPATH -u PYTHONHOME -u VIRTUAL_ENV
  -u VINKULUM_BUNDLE_EXAMPLES -u VINKULUM_BUNDLE_GMSH -u VINKULUM_BUNDLE_CALCULIX
  XDG_CONFIG_HOME="$BUILD_DIR/config" QT_QPA_PLATFORM=cocoa
  "$APP/Contents/MacOS/Vinkulum Studio")
"$PY" - "${COMMAND[@]}" --startup-check "$BUILD_DIR/check" <<'PY'
import subprocess
import sys
subprocess.run(sys.argv[1:], check=True, timeout=45)
PY
"${COMMAND[@]}" --bundle-check "$BUILD_DIR/check"
"$PY" - "$BUILD_DIR/check" "$VERSION" "$APP/Contents/Resources/build-info.json" <<'PY'
import json
from pathlib import Path
import sys

checks = Path(sys.argv[1])
report = json.loads((checks / "bundle-check.json").read_text())
assert report["status"] == "passed" and report["frozen"]
assert report["machine"] == "arm64"
assert report["manifest"]["app_version"] == sys.argv[2]
assert report["calculix"]["status"] == "not_requested"
previews = report["source_previews"]
for name in ("pinocchio_capture", "cad_history", "cad_mesh_capture", "quadratic_static_capture"):
    assert previews[name]["status"] == "passed", name
assert previews["example_files_unchanged"]
assert previews["gmsh"]["status"] == "not_requested"
startup = json.loads((checks / "startup-check.json").read_text())
assert startup["status"] == "passed" and startup["frozen"]
assert startup["app_version"] == sys.argv[2] and startup["locale"] == "en_GB"
build = json.loads(Path(sys.argv[3]).read_text())
assert build["machine"] == "arm64"
assert build["packages"]["vinkulum-studio"] == sys.argv[2]
assert not build["source_dirty"]
PY
# Only expose a DMG after its copied application passed every check.
cp "$BUILD_DIR/$DMG" "$OUT/"
cp "$VINKULUM_MACOS_PAYLOAD/build-info.json" "$OUT/"
cp -R "$BUILD_DIR/check" "$OUT/"
cd "$OUT"
shasum -a 256 "$DMG" > SHA256SUMS
shasum -a 256 --check SHA256SUMS
echo "Apple Silicon DMG tested and ready: $OUT/$DMG"
