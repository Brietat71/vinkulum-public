#!/usr/bin/env bash
set -euo pipefail
[[ "$(uname -sm)" == "Darwin arm64" ]] || { echo 'Native Apple Silicon macOS required.' >&2; exit 1; }
PY=${PY:-python}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
"$PY" -m PyInstaller --noconfirm --clean --distpath dist/macos --workpath build/studio-macos apps/studio/packaging/studio.spec
APP="$ROOT/dist/macos/Vinkulum Studio.app"
codesign --verify --deep --strict "$APP"
lipo "$APP/Contents/MacOS/Vinkulum Studio" -verify_arch arm64
# Launch the executable as installed, from outside the source tree.
cd "${RUNNER_TEMP:-/tmp}"
"$APP/Contents/MacOS/Vinkulum Studio" --bundle-check "$ROOT/dist/macos/check"
cd "$ROOT"
mkdir -p dist/dmg-root
/usr/bin/ditto "$APP" 'dist/dmg-root/Vinkulum Studio.app'
ln -s /Applications dist/dmg-root/Applications
cp apps/studio/packaging/INSTALLATION.txt dist/dmg-root/
cp -R apps/studio/packaging/licenses dist/dmg-root/Licences
cp LICENSE NOTICE THIRD_PARTY_NOTICES.md dist/dmg-root/Licences/
hdiutil create -volname 'Vinkulum Studio' -srcfolder dist/dmg-root -ov -format UDZO dist/macos/Vinkulum-Studio-0.3.0-apple-silicon.dmg
hdiutil verify dist/macos/Vinkulum-Studio-0.3.0-apple-silicon.dmg
shasum -a 256 dist/macos/*.dmg > dist/macos/SHA256SUMS
