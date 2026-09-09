#!/usr/bin/env bash
set -euo pipefail
[[ "$(uname -sm)" == "Darwin arm64" ]] || { echo 'Native Apple Silicon macOS required.' >&2; exit 1; }
PY=${PY:-python}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
"$PY" -m PyInstaller --noconfirm --clean --distpath dist/macos --workpath build/studio-macos apps/studio/packaging/studio.spec
APP="$ROOT/dist/macos/Vinkulum Studio.app"
codesign --verify --deep --strict "$APP"
lipo -verify_arch arm64 "$APP/Contents/MacOS/Vinkulum Studio"
# Launch the executable as installed, from outside the source tree.
cd "${RUNNER_TEMP:-/tmp}"
"$APP/Contents/MacOS/Vinkulum Studio" --bundle-check "$ROOT/dist/macos/check"
cd "$ROOT"
mkdir -p dist/dmg-root
/usr/bin/ditto "$APP" 'dist/dmg-root/Vinkulum Studio.app'
ln -s /Applications dist/dmg-root/Applications
cp apps/studio/packaging/INSTALLATION.txt dist/dmg-root/
hdiutil create -volname 'Vinkulum Studio' -srcfolder dist/dmg-root -ov -format UDZO dist/macos/Vinkulum-Studio-0.2.0-apple-silicon.dmg
hdiutil verify dist/macos/Vinkulum-Studio-0.2.0-apple-silicon.dmg
shasum -a 256 dist/macos/*.dmg > dist/macos/SHA256SUMS
